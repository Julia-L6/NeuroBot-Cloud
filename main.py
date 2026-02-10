import os
import csv
import datetime
import time
import smtplib
from email.mime.text import MIMEText
from email.header import Header
import google.generativeai as genai
from Bio import Entrez

# --- 1. 配置区域 ---
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

# --- 2. 动态期刊数据库 ---
def load_journal_db():
    """读取仓库中的 journals.csv 文件"""
    db = {}
    csv_path = 'journals.csv'
    
    if not os.path.exists(csv_path):
        print("⚠️ 未找到 journals.csv，IF 功能将受限")
        return {}

    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader, None) # 跳过标题行
            for row in reader:
                if len(row) >= 2:
                    # key: 归一化后的期刊名 (小写, 去空格)
                    name = row[0].lower().strip()
                    if_val = row[1].strip()
                    zone = row[2].strip() if len(row) > 2 else "?"
                    db[name] = {"if": if_val, "q": zone}
        print(f"📚 已加载 {len(db)} 条期刊数据")
    except Exception as e:
        print(f"❌ 读取 CSV 失败: {e}")
    return db

# 全局加载一次
JOURNAL_DB = load_journal_db()

def get_journal_metrics(journal_name):
    """
    匹配期刊 IF。
    策略：
    1. 完全匹配 (忽略大小写)
    2. 如果匹配不到，返回 N/A，不乱猜
    """
    if not journal_name: return "N/A", "N/A"
    
    # 清洗：Journal of X (London) -> journal of x
    clean_name = journal_name.lower().split('(')[0].strip()
    
    # 1. 直接查表
    if clean_name in JOURNAL_DB:
        return JOURNAL_DB[clean_name]['if'], JOURNAL_DB[clean_name]['q']
    
    # 2. 尝试移除 "The" (例如 "The Lancet" vs "Lancet")
    if clean_name.startswith("the "):
        alt_name = clean_name[4:].strip()
        if alt_name in JOURNAL_DB:
            return JOURNAL_DB[alt_name]['if'], JOURNAL_DB[alt_name]['q']

    return "N/A", ""

def setup_gemini():
    if not GOOGLE_API_KEY:
        print("❌ 错误: 缺少 GOOGLE_API_KEY")
        return None
    genai.configure(api_key=GOOGLE_API_KEY)
    return genai.GenerativeModel('gemini-flash-latest')

def search_pubmed_ids(query, max_results):
    print(f"🔍 检索: {query[:50]}... (Max: {max_results})")
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="date", reldate=7, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"⚠️ 检索失败: {e}")
        return []

def fetch_and_parse(id_list, tag_label):
    if not id_list: return []
    print(f"📥 [{tag_label}] 下载 {len(id_list)} 篇详情...")
    
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="xml")
        xml_data = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"❌ 下载失败: {e}")
        return []

    parsed_list = []
    articles = xml_data.get('PubmedArticle', [])
    
    for article in articles:
        try:
            medline = article.get('MedlineCitation', {})
            art = medline.get('Article', {})
            
            title = art.get('ArticleTitle')
            if not title: continue 

            abstract_list = art.get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if abstract_list else "【无摘要】"
            
            journal_info = art.get('Journal', {})
            journal_name = journal_info.get('Title', 'Unknown')
            
            # 获取 IF 和分区
            if_val, zone = get_journal_metrics(journal_name)
            
            # 只有当 IF 是数字时才用于排序，否则设为 -1
            try:
                sort_score = float(if_val)
            except:
                sort_score = -1.0

            pub_date = art.get('ArticleDate', [])
            date_str = "Recent"
            if pub_date:
                date_str = f"{pub_date[0].get('Year')}-{pub_date[0].get('Month')}-{pub_date[0].get('Day')}"

            authors = art.get('AuthorList', [])
            first_auth = "Unknown"
            if authors:
                f = authors[0]
                aff = ""
                if f.get('AffiliationInfo'):
                    aff = f['AffiliationInfo'][0].get('Affiliation', '')
                flag = ""
                if "China" in aff: flag = "[🇨🇳China]"
                elif "USA" in aff: flag = "[🇺🇸USA]"
                first_auth = f"{f.get('LastName')} {f.get('ForeName')} {flag}"

            parsed_list.append({
                "title": title,
                "abstract": abstract,
                "journal": journal_name,
                "if": if_val,
                "zone": zone,
                "sort_score": sort_score, # 专门用于排序的数字
                "date": date_str,
                "author": first_auth,
                "tag": tag_label
            })
            
        except Exception as e:
            continue

    return parsed_list

def analyze_with_ai(model, paper):
    print(f"🤖 AI 阅读: {paper['title'][:30]}...")
    
    # 构造 IF 显示字符串
    if_str = paper['if']
    if paper['zone']:
        if_str += f" ({paper['zone']})"
    
    prompt = f"""
    你是一位资深神经科学研究员。请简要分析这篇文献。
    
    标题: {paper['title']}
    摘要: {paper['abstract']}
    
    请严格按Markdown格式输出:
    ### {paper['tag']} | {paper['title']}
    > 📅 {paper['date']} | 📖 {paper['journal']} | 📊 IF: {if_str} | 👤 {paper['author']}
    
    - **🏷️ 类型**: [综述/动物/细胞/临床/数据]
    - **🧐 核心**: (简要总结，关注定量数据)
    - **🔬 机制**: (关键分子/通路/药味)
    ---
    """
    
    for _ in range(3):
        try:
            res = model.generate_content(prompt)
            return res.text.replace('\xa0', ' ')
        except:
            time.sleep(5)
    return ""

def send_email(subject, content):
    if not EMAIL_PASSWORD:
        return
    
    msg = MIMEText(content, 'plain', 'utf-8')
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER
    msg['Subject'] = Header(subject, 'utf-8')
    
    try:
        s = smtplib.SMTP(SMTP_SERVER, int(SMTP_PORT))
        s.starttls()
        s.login(EMAIL_SENDER, EMAIL_PASSWORD)
        s.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        s.quit()
        print("✅ 邮件发送成功")
    except Exception as e:
        print(f"❌ 发送失败: {e}")

def main():
    model = setup_gemini()
    if not model: return

    # 1. Neuro 通道
    q_neuro = '(Alzheimer\'s disease[Title/Abstract] AND (microglia[Title/Abstract] OR neuroinflammation[Title/Abstract]))'
    # 抓取 40 篇，增加筛选池
    neuro_raw = fetch_and_parse(search_pubmed_ids(q_neuro, 40), "🧠Neuro")
    
    # 排序：优先按 IF 分数高低，Unknown 的排后面
    neuro_raw.sort(key=lambda x: x['sort_score'], reverse=True)
    final_neuro = neuro_raw[:20]

    # 2. TCM 通道
    q_tcm = '((Traditional Chinese Medicine[Title/Abstract] OR Herbal[Title/Abstract] OR Acupuncture[Title/Abstract]) AND (Alzheimer[Title/Abstract] OR Brain[Title/Abstract]))'
    tcm_raw = fetch_and_parse(search_pubmed_ids(q_tcm, 15), "🌿TCM")
    
    tcm_raw.sort(key=lambda x: x['sort_score'], reverse=True)
    final_tcm = tcm_raw[:5]

    total = len(final_neuro) + len(final_tcm)
    if total == 0:
        print("📭 无数据")
        return

    content = f"🧠 NeuroBot 日报 ({datetime.date.today()})\n"
    content += f"📚 IF数据源: 本地数据库 (匹配失败显示 N/A)\n\n"

    if final_tcm:
        content += "--- 🌿 TCM 特别关注 ---\n\n"
        for p in final_tcm:
            content += analyze_with_ai(model, p)
            time.sleep(3)

    content += "\n--- 🧠 神经科学核心推荐 ---\n\n"
    for p in final_neuro:
        content += analyze_with_ai(model, p)
        time.sleep(3)

    send_email(f"NeuroBot日报 ({total}篇) - {datetime.date.today()}", content)

if __name__ == "__main__":
    main()
