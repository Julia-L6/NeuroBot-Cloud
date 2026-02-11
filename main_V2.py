import os
import csv
import datetime
import time
import smtplib
import re
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.header import Header
import google.generativeai as genai
from Bio import Entrez

# 核心作用：本代码实现了一个全自动化的神经科学科研情报工具，每日定时从 PubMed 检索并筛选 Alzheimer's Disease 与 Microglia 领域的最新高分文献，利用 Gemini AI 生成包含“重点深度解读”与“快速浏览列表”的专业中文简报，最终自动以邮件形式（附带数据 CSV）发送给用户。

# --- 1. 全局配置 ---
Entrez.email = "julia_light@msn.cn"
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
SMTP_SERVER = os.getenv("SMTP_SERVER")
SMTP_PORT = os.getenv("SMTP_PORT")
EMAIL_SENDER = os.getenv("EMAIL_SENDER")
raw_password = os.getenv("EMAIL_PASSWORD")
EMAIL_RECEIVER = "julia_light@msn.cn"

if raw_password:
    EMAIL_PASSWORD = raw_password.replace(' ', '').replace('\xa0', '').strip()
else:
    EMAIL_PASSWORD = None

# 检索配置
SEARCH_WINDOW_DAYS = 3  # Day 0 - Day 3
TARGET_TOTAL_COUNT = 40 # 最终保留的 Top 文献数

# --- 2. 期刊数据库加载 ---
def load_journal_db():
    db = {}
    csv_path = 'journals.csv'
    if not os.path.exists(csv_path):
        return {}
    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if len(row) >= 2:
                    name = row[0].lower().strip()
                    try: if_val = float(row[1].strip())
                    except: if_val = 0.0
                    zone = row[2].strip() if len(row) > 2 else "?"
                    db[name] = {"if": if_val, "q": zone}
    except: pass
    return db

JOURNAL_DB = load_journal_db()

def get_journal_metrics(journal_name):
    if not journal_name: return 0.0, "?", False
    clean_name = journal_name.lower().split('(')[0].strip().replace('.', '')
    data = {"if": 0.0, "q": "?"}
    
    if clean_name in JOURNAL_DB:
        data = JOURNAL_DB[clean_name]
    elif clean_name.startswith("the ") and clean_name[4:] in JOURNAL_DB:
        data = JOURNAL_DB[clean_name[4:]]
        
    is_high = data["if"] > 10.0 or "Q1" in data["q"]
    return data["if"], data["q"], is_high

# --- 3. 辅助提取函数 ---
def extract_country(text):
    if not text: return "Unknown"
    # 优先匹配常见科研大国
    countries = ["USA", "China", "UK", "Germany", "Japan", "Canada", "Australia", "France", "Italy", "Spain", "Netherlands", "Switzerland", "Korea", "Singapore"]
    text_u = text.upper()
    if "USA" in text_u or "UNITED STATES" in text_u: return "USA"
    if "CHINA" in text_u: return "China"
    if "UK" in text_u or "UNITED KINGDOM" in text_u: return "UK"
    
    for c in countries:
        if c.upper() in text_u: return c
    return "Global"

def format_author_details(author_entry):
    """提取格式: 姓名 [机构, 国家]"""
    if not author_entry: return "Unknown"
    
    last = author_entry.get('LastName', '')
    initials = author_entry.get('Initials', '')
    name = f"{last} {initials}"
    
    aff_info = author_entry.get('AffiliationInfo', [])
    if aff_info:
        raw_aff = aff_info[0].get('Affiliation', '')
        country = extract_country(raw_aff)
        # 简化机构名: 取逗号前的内容，或者 University of X
        short_aff = raw_aff.split(',')[0]
        # 简单的清洗逻辑
        short_aff = re.sub(r'Department of [A-Za-z\s]+', '', short_aff).strip()
        if len(short_aff) > 40: short_aff = short_aff[:37] + "..."
        
        return f"{name} [{short_aff}, {country}]"
    return name

# --- 4. PubMed 核心逻辑 ---
def build_search_queries():
    TERMS_AD = '("Alzheimer Disease"[Mesh] OR "Alzheimer\'s"[Title/Abstract] OR "Alzheimers"[Title/Abstract] OR "AD"[Title/Abstract] OR "Amyloid beta"[Title/Abstract] OR "Tau"[Title/Abstract])'
    TERMS_MICROGLIA = '("Microglia"[Mesh] OR "microglia"[Title/Abstract] OR "microglial"[Title/Abstract] OR "Neuroinflammation"[Title/Abstract] OR "DAM"[Title/Abstract] OR "TREM2"[Title/Abstract])'
    TERMS_NEURO_BROAD = '("Neuroscience"[Title/Abstract] OR "Neurology"[Title/Abstract] OR "Brain"[Title/Abstract] OR "CNS"[Title/Abstract])'

    return [
        f"({TERMS_AD} AND {TERMS_MICROGLIA})", # 核心 Q1-Q3
        f"{TERMS_AD}",                          # AD 高分
        f"({TERMS_MICROGLIA} AND {TERMS_NEURO_BROAD})" # Microglia 高分
    ]

def search_and_fetch_all():
    queries = build_search_queries()
    all_pmids = set()
    
    print(f"🔍 执行高级检索 (Window: {SEARCH_WINDOW_DAYS} days)...")
    for i, q in enumerate(queries):
        try:
            handle = Entrez.esearch(db="pubmed", term=q, retmax=100, sort="date", reldate=SEARCH_WINDOW_DAYS, datetype="pdat")
            ids = Entrez.read(handle)["IdList"]
            all_pmids.update(ids)
            time.sleep(1)
        except Exception as e:
            print(f"  Query {i} Error: {e}")

    unique_ids = list(all_pmids)
    if not unique_ids: return []

    print(f"📥 下载 {len(unique_ids)} 篇文献详情...")
    batch_size = 50
    final_papers = []
    
    for i in range(0, len(unique_ids), batch_size):
        try:
            handle = Entrez.efetch(db="pubmed", id=unique_ids[i:i+batch_size], rettype="xml", retmode="xml")
            data = Entrez.read(handle)
            handle.close()
            
            for article in data.get('PubmedArticle', []):
                try:
                    med = article['MedlineCitation']
                    art = med['Article']
                    
                    title = art.get('ArticleTitle', 'No Title')
                    journal = art.get('Journal', {}).get('Title', 'Unknown')
                    if_val, zone, is_high_impact = get_journal_metrics(journal)
                    
                    # 日期处理
                    d = art.get('ArticleDate', [])
                    date_str = f"{d[0]['Year']}-{d[0]['Month']}-{d[0]['Day']}" if d else "Recent"
                    pub_status = "Online" if d else "Print"
                    
                    # 作者处理
                    authors = art.get('AuthorList', [])
                    first_auth = format_author_details(authors[0]) if authors else "Unknown"
                    senior_auth = format_author_details(authors[-1]) if len(authors) > 1 else first_auth
                    
                    abstract = " ".join(art.get('Abstract', {}).get('AbstractText', []))
                    if not abstract: abstract = "No Abstract"
                    
                    # DOI & Types
                    doi = next((str(x) for x in med.get('Article', {}).get('ELocationID', []) if x.attributes.get('EIdType')=='doi'), "")
                    types = [t for t in art.get('PublicationTypeList', []) if "Journal" not in t and "Support" not in t]
                    type_str = ", ".join(types) if types else "Research"

                    # 评分筛选逻辑
                    score = if_val
                    txt = (title + abstract).lower()
                    is_core = ("alzheimer" in txt or "amyloid" in txt) and ("microglia" in txt or "neuroinflammation" in txt)
                    
                    if is_core: score += 50
                    elif is_high_impact: score += 30
                    
                    # 过滤: 非高分 且 非核心 且 IF<3
                    if not is_high_impact and not is_core and if_val < 3.0:
                        continue

                    final_papers.append({
                        "PMID": str(med['PMID']),
                        "Title": title,
                        "Journal": journal,
                        "IF": if_val,
                        "Date": date_str,
                        "Status": pub_status,
                        "Type": type_str,
                        "First_Author": first_auth,
                        "Senior_Author": senior_auth,
                        "Abstract": abstract,
                        "DOI": doi,
                        "Score": score
                    })
                except: continue
        except: pass
            
    final_papers.sort(key=lambda x: x["Score"], reverse=True)
    return final_papers[:TARGET_TOTAL_COUNT]

# --- 5. Gemini AI 分析 ---
def analyze_with_gemini(papers):
    if not GOOGLE_API_KEY: return "❌ 无 API Key"
    
    # 强制使用 Flash 模型
    model_name = 'gemini-flash-latest'
    
    # 配置生成参数：降低 Temperature 防止幻觉
    generation_config = {
        "temperature": 0.2,
        "top_p": 0.95,
        "top_k": 64,
        "max_output_tokens": 8192,
    }
    
    genai.configure(api_key=GOOGLE_API_KEY)
    model = genai.GenerativeModel(model_name, generation_config=generation_config)
    
    # 构建输入数据
    csv_block = "PMID | Title | Journal (IF) | First Author [Affiliation] | Senior Author [Affiliation] | Abstract\n"
    for p in papers:
        # 摘要截断以防超长，但Flash通常没问题
        abs_clean = p['Abstract'][:1200]
        csv_block += f"{p['PMID']} | {p['Title']} | {p['Journal']} (IF:{p['IF']}) | {p['First_Author']} | {p['Senior_Author']} | {abs_clean}\n"

    prompt = f"""
# Role Assignment
你是一位拥有 20 年经验的资深神经科学领域研究员，擅长快速阅读英文学术文献，并将其核心价值转化为逻辑严密、通俗易懂的中文技术简报。你对 Alzheimer's disease、Microglia 以及其他神经退行性疾病领域前沿科研知识与技术有深刻的理解。

# Task
我将提供一批最新的文献信息（CSV格式数据）。请你阅读并分析这些文献，生成一份中文“每日科研简报”。

# Constraints
1. **必须使用中文进行输出**，但必须保留必要的英文专业术语，或在中文后括号内保留英文专业术语（例如：靶点蛋白如 TREM2，模型名称如 5xFAD，关键表型如 phagocytosis 等）。
2. **严禁直接翻译原文摘要**，必须基于理解进行重述和概括。
3. **语气保持客观、专业**，避免使用营销式夸张词汇。
4. **"创新点"部分必须具体**，指出该论文解决了什么具体痛点，不仅是罗列功能。
5. **必须一次性完整输出**：请合理分配篇幅，务必在一条回复中同时包含【🌟 重点关注】和【📂 快速浏览】两个版块。

# Analysis Requirements

## Section 1: 🌟 重点关注 (High Priority)
筛选标准：高分期刊 (IF>10) 或 机制创新极强（如发现新靶点/新通路）的研究。
请对每一篇重点文献按以下格式进行深读：

- **标题**：(中文翻译)
- **期刊/IF**：(保留原名)
- **👥 关键作者**：
  - 一作：(提取自 CSV First Author 列，保留姓名、机构、国家)
  - 通讯：(提取自 CSV Senior Author 列，保留姓名、机构、国家)
- **🏷️ 类型**：[in vivo机制 / in vitro机制 / 生信 / 临床 / 综述] (可多选，如涵盖多种以+号连接)
- **🧐 核心发现**：(One-Liner 亮点)
- **✨ 创新点**：(具体指出该论文解决了什么具体痛点，或填补了什么空白)
- **🧬 关键机制/靶点**：(指明研究聚焦的具体通路、靶点、表型、细胞功能等，及研究使用的主要方法)
- **💡 简评**：(专业评价其科学价值)

---

## Section 2: 📂 快速浏览 (Quick Browse)
筛选标准：验证性研究、低分期刊或纯临床统计文章。
**请务必以表格形式展示，不要分段描述：**

| 序号 | 标题 (中文) | 期刊 | 类型 | 核心发现 (一句话) |
| :--- | :--- | :--- | :--- | :--- |
| 1 | ... | ... | ... | ... |

# Input Data
{csv_block}
"""
    
    print(f"🤖 请求 Gemini ({model_name}) 分析 {len(papers)} 篇文献...")
    try:
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        return f"❌ AI 分析生成失败: {e}"

# --- 6. 邮件发送 ---
def send_email_report(content, csv_file):
    if not EMAIL_PASSWORD: return

    msg = MIMEMultipart()
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER
    today = datetime.date.today()
    msg['Subject'] = Header(f"🧠 NeuroBot 深度简报 - {today}", 'utf-8')

    msg.attach(MIMEText(f"NeuroBot 自动分析完成。\n\n{content}", 'plain', 'utf-8'))

    try:
        with open(csv_file, "rb") as f:
            att = MIMEApplication(f.read(), Name=os.path.basename(csv_file))
            att['Content-Disposition'] = f'attachment; filename="{os.path.basename(csv_file)}"'
            msg.attach(att)
    except: pass

    try:
        server = smtplib.SMTP(SMTP_SERVER, int(SMTP_PORT))
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        server.quit()
        print("✅ 邮件发送成功！")
    except Exception as e:
        print(f"❌ 邮件发送失败: {e}")

# --- 7. 主程序 ---
def main():
    print(f"🚀 NeuroBot V2 启动 - {datetime.date.today()}")
    
    # 1. 检索与清洗
    papers = search_and_fetch_all()
    if not papers:
        print("📭 今日无符合条件数据")
        return

    # 2. 生成 CSV
    csv_name = f"NeuroBot_Data_{datetime.date.today()}.csv"
    headers = ["PMID", "Title", "Journal", "IF", "Date", "Status", "Type", "First_Author", "Senior_Author", "DOI", "Abstract", "Score"]
    with open(csv_name, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(papers)

    # 3. AI 分析
    report = analyze_with_gemini(papers)
    
    # 4. 发送
    send_email_report(report, csv_name)

if __name__ == "__main__":
    main()
