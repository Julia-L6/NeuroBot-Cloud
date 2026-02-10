import os
import datetime
import time
import smtplib
import re
from email.mime.text import MIMEText
from email.header import Header
import google.generativeai as genai
from Bio import Entrez

# --- 1. 配置区域 ---
# 必须配置 email，PubMed 要求用于追踪
Entrez.email = "your_email@example.com" 
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
# 邮箱配置
SMTP_SERVER = os.getenv("SMTP_SERVER")
SMTP_PORT = os.getenv("SMTP_PORT")
EMAIL_SENDER = os.getenv("EMAIL_SENDER")
raw_password = os.getenv("EMAIL_PASSWORD")
EMAIL_RECEIVER = "julia_light@msn.cn"

# 清洗密码
if raw_password:
    EMAIL_PASSWORD = raw_password.replace(' ', '').replace('\xa0', '').strip()
else:
    EMAIL_PASSWORD = None

# --- 2. 期刊 IF 简易数据库 (可自行补充) ---
# 由于 IF 是商业数据，无法免费通过 API 实时获取。
# 这里采用“本地字典”方式，覆盖常见 Neuro 和 TCM 核心期刊。
JOURNAL_DB = {
    "nature": {"if": 64.8, "q": "Q1"},
    "science": {"if": 56.9, "q": "Q1"},
    "cell": {"if": 45.5, "q": "Q1"},
    "nature neuroscience": {"if": 21.2, "q": "Q1"},
    "neuron": {"if": 16.2, "q": "Q1"},
    "molecular neurodegeneration": {"if": 15.1, "q": "Q1"},
    "alzheimer's & dementia": {"if": 14.0, "q": "Q1"},
    "acta neuropathologica": {"if": 12.7, "q": "Q1"},
    "brain": {"if": 10.6, "q": "Q1"},
    "journal of neuroinflammation": {"if": 9.3, "q": "Q1"},
    "glia": {"if": 8.0, "q": "Q1"},
    "phytomedicine": {"if": 6.7, "q": "Q1"}, # TCM 强刊
    "journal of ethnopharmacology": {"if": 5.4, "q": "Q1"}, # TCM
    "frontiers in immunology": {"if": 5.7, "q": "Q1"},
    "aging cell": {"if": 7.8, "q": "Q1"},
}

def get_journal_metrics(journal_name):
    """根据期刊名查找 IF 和分区"""
    if not journal_name: return "N/A", "N/A"
    clean_name = journal_name.lower().strip()
    # 模糊匹配
    for k, v in JOURNAL_DB.items():
        if k in clean_name:
            return v["if"], v["q"]
    return "Unknown", "?"

def setup_gemini():
    if not GOOGLE_API_KEY:
        print("❌ 错误: 缺少 GOOGLE_API_KEY")
        return None
    genai.configure(api_key=GOOGLE_API_KEY)
    return genai.GenerativeModel('gemini-flash-latest')

def search_pubmed(query, max_retries=3):
    """在 PubMed 中搜索并获取详细信息"""
    print(f"🔍 正在 PubMed 检索: {query}")
    
    # 1. 搜索 ID (过去 5 天，确保包含 Online First)
    # reldate=5 表示最近5天，datetype="pdat" 表示出版日期(含电子出版)
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="date", reldate=7, datetype="pdat")
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            break
        except Exception as e:
            print(f"⚠️ 连接 PubMed 失败 ({attempt+1}/{max_retries}): {e}")
            time.sleep(5)
    else:
        return []

    if not id_list:
        return []

    print(f"📥 获取到 {len(id_list)} 篇文献 ID，正在下载详情...")

    # 2. 获取详情
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="xml")
        papers = Entrez.read(handle)
        handle.close()
        return papers['PubmedArticle']
    except Exception as e:
        print(f"❌ 下载详情失败: {e}")
        return []

def extract_paper_info(article):
    """从 XML 中解析复杂的论文元数据"""
    try:
        medline = article['MedlineCitation']
        article_data = medline['Article']
        journal_info = article_data['Journal']
        
        # 1. 基础信息
        title = article_data.get('ArticleTitle', 'No Title')
        abstract_list = article_data.get('Abstract', {}).get('AbstractText', [])
        abstract = " ".join(abstract_list) if abstract_list else "No Abstract"
        
        # 2. 期刊与 IF
        journal_name = journal_info.get('Title', 'Unknown Journal')
        if_score, quartile = get_journal_metrics(journal_name)
        
        # 3. 日期 (优先用 PubMed PubDate)
        pub_date = article_data.get('ArticleDate', [])
        if pub_date:
            date_str = f"{pub_date[0].get('Year')}-{pub_date[0].get('Month')}-{pub_date[0].get('Day')}"
        else:
            date_str = "Recent"

        # 4. 作者与机构 (难点)
        authors = article_data.get('AuthorList', [])
        first_author = "Unknown"
        first_affil = "Unknown"
        corresp_author = "Unknown"
        corresp_affil = "Unknown"

        if authors:
            # 一作
            f_auth = authors[0]
            first_author = f"{f_auth.get('LastName', '')} {f_auth.get('ForeName', '')}"
            if 'AffiliationInfo' in f_auth and f_auth['AffiliationInfo']:
                first_affil = f_auth['AffiliationInfo'][0].get('Affiliation', 'Unknown')

            # 通讯作者 (通常是列表最后一个，或带有 Affiliation 的)
            l_auth = authors[-1]
            corresp_author = f"{l_auth.get('LastName', '')} {l_auth.get('ForeName', '')}"
            if 'AffiliationInfo' in l_auth and l_auth['AffiliationInfo']:
                corresp_affil = l_auth['AffiliationInfo'][0].get('Affiliation', 'Unknown')
        
        # 简单的国家提取正则
        def extract_country(text):
            common_countries = ["China", "USA", "United States", "UK", "Germany", "Japan", "Canada", "Australia"]
            for c in common_countries:
                if c.lower() in text.lower(): return c
            return "Global"

        first_country = extract_country(first_affil)
        
        return {
            "title": title,
            "abstract": abstract,
            "journal": journal_name,
            "if": if_score,
            "q": quartile,
            "date": date_str,
            "first_author": f"{first_author} ({first_country})",
            "corresp_author": corresp_author,
            "first_affil": first_affil[:50] + "..." if len(first_affil)>50 else first_affil, # 截断一下防止太长
            "doi": f"10.1038/..." # 简化，实际需遍历ELocationID
        }
    except Exception as e:
        print(f"⚠️ 解析单篇出错: {e}")
        return None

def analyze_paper_gemini(model, paper_info):
    print(f"🤖 正在阅读: {paper_info['title'][:30]}...")
    prompt = f"""
    你是一位神经科学与中医药领域的审稿专家。请分析这篇文献。
    
    标题: {paper_info['title']}
    摘要: {paper_info['abstract']}
    期刊: {paper_info['journal']} (IF: {paper_info['if']})
    
    请用以下中文 Markdown 格式输出 (不要使用 \xa0 特殊空格):
    
    ### 📄 {paper_info['title']}
    - **📚 期刊**: {paper_info['journal']} | **IF: {paper_info['if']} ({paper_info['q']})** | 📅 {paper_info['date']}
    - **👤 作者**: 一作 {paper_info['first_author']} | 通讯 {paper_info['corresp_author']}
    - **🏫 机构**: {paper_info['first_affil']}
    - **🧐 核心发现**: (一句话概括，如果是TCM相关请重点突出中药成分)
    - **🔬 机制/靶点**: (关键通路、分子或技术)
    - **💡 创新评级**: (根据IF和内容给出 ⭐~⭐⭐⭐⭐⭐)
    ---
    """
    for _ in range(3):
        try:
            response = model.generate_content(prompt)
            return response.text.replace('\xa0', ' ')
        except Exception as e:
            time.sleep(10)
    return "❌ AI 分析超时\n\n"

def send_email(subject, content):
    if not EMAIL_PASSWORD:
        print("⚠️ 邮箱配置缺失")
        return
    
    content = content.replace('\xa0', ' ')
    msg = MIMEText(content, 'plain', 'utf-8')
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER
    msg['Subject'] = Header(subject, 'utf-8')
    
    try:
        server = smtplib.SMTP(SMTP_SERVER, int(SMTP_PORT))
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        server.quit()
        print(f"✅ 邮件已发送")
    except Exception as e:
        print(f"❌ 发送失败: {e}")

def main():
    model = setup_gemini()
    if not model: return
    
    # --- 构造 PubMed 查询 ---
    # 逻辑: (AD + 小胶质) OR (AD + TCM/中药)
    # 这样既能包含纯神经机制，也能包含中药干预
    query_neuro = '(Alzheimer[Title/Abstract] AND microglia[Title/Abstract])'
    query_tcm = '(Alzheimer[Title/Abstract] AND (Traditional Chinese Medicine[Title/Abstract] OR herbal[Title/Abstract] OR acupuncture[Title/Abstract]))'
    
    full_query = f"({query_neuro} OR {query_tcm})"
    
    # 1. 获取原始数据
    raw_articles = search_pubmed(full_query)
    
    parsed_papers = []
    for art in raw_articles:
        info = extract_paper_info(art)
        if info:
            parsed_papers.append(info)
            
    if not parsed_papers:
        print("📭 今日无新文献")
        send_email("NeuroBot: 今日无新文献", "PubMed 检索结果为空。")
        return

    # 2. 排序与筛选 (核心算法)
    # 优先看 IF (如果是 Unknown 设为 0)，其次看日期
    def sort_key(p):
        try:
            score = float(p['if'])
        except:
            score = 0.0
        return score

    # 按 IF 降序排列
    parsed_papers.sort(key=sort_key, reverse=True)
    
    # 取前 20 篇
    top_papers = parsed_papers[:20]
    
    print(f"📊 已筛选 Top {len(top_papers)} 篇 (按 IF 排序)")

    # 3. AI 分析与生成
    content = f"🧠 NeuroBot 每日精选 ({datetime.date.today()})\n"
    content += f"🔍 来源: PubMed (Online First & Published)\n"
    content += f"📊 策略: IF 优先排序 | 覆盖 AD, Microglia, TCM\n\n"
    
    for paper in top_papers:
        summary = analyze_paper_gemini(model, paper)
        content += summary
        time.sleep(3) # 避免 Gemini 速率限制
        
    send_email(f"NeuroBot日报 (Top {len(top_papers)}) - {datetime.date.today()}", content)

if __name__ == "__main__":
    main()
