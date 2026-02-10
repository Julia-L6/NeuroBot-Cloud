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

# 🔴 策略调整区
# 每天只深度分析多少篇？(建议不超过15篇以避免429错误)
MAX_AI_ANALYSIS_NEURO = 10 
MAX_AI_ANALYSIS_TCM = 3
# 检索最近几天？(建议2-3天，避免0结果)
SEARCH_WINDOW_DAYS = 2 

if raw_password:
    EMAIL_PASSWORD = raw_password.replace(' ', '').replace('\xa0', '').strip()
else:
    EMAIL_PASSWORD = None

# --- 2. 动态期刊数据库 ---
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
                    db[name] = {"if": row[1].strip(), "q": row[2].strip() if len(row) > 2 else "?"}
    except:
        pass
    return db

JOURNAL_DB = load_journal_db()

def get_journal_metrics(journal_name):
    if not journal_name: return "N/A", "N/A"
    clean_name = journal_name.lower().split('(')[0].strip()
    if clean_name in JOURNAL_DB:
        return JOURNAL_DB[clean_name]['if'], JOURNAL_DB[clean_name]['q']
    if clean_name.startswith("the "):
        alt = clean_name[4:].strip()
        if alt in JOURNAL_DB: return JOURNAL_DB[alt]['if'], JOURNAL_DB[alt]['q']
    return "N/A", ""

def setup_gemini():
    if not GOOGLE_API_KEY:
        print("❌ 无 API KEY")
        return None
    genai.configure(api_key=GOOGLE_API_KEY)
    return genai.GenerativeModel('gemini-flash-latest')

def search_pubmed_ids(query, max_results):
    print(f"🔍 检索(近{SEARCH_WINDOW_DAYS}天): {query[:30]}...")
    try:
        # ✅ 关键修改：reldate 使用配置变量
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="date", reldate=SEARCH_WINDOW_DAYS, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"⚠️ 检索失败: {e}")
        return []

def fetch_and_parse(id_list, tag_label):
    if not id_list: return []
    print(f"📥 [{tag_label}] 下载 {len(id_list)} 篇...")
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="xml")
        xml_data = Entrez.read(handle)
        handle.close()
    except:
        return []

    parsed_list = []
    for article in xml_data.get('PubmedArticle', []):
        try:
            art = article.get('MedlineCitation', {}).get('Article', {})
            title = art.get('ArticleTitle')
            if not title: continue
            
            abstract = " ".join(art.get('Abstract', {}).get('AbstractText', [])) or "【无摘要】"
            j_name = art.get('Journal', {}).get('Title', 'Unknown')
            if_val, zone = get_journal_metrics(j_name)
            
            try: sort_score = float(if_val)
            except: sort_score = -1.0
            
            d = art.get('ArticleDate', [])
            date_str = f"{d[0]['Year']}-{d[0]['Month']}-{d[0]['Day']}" if d else "Recent"
            
            authors = art.get('AuthorList', [])
            auth = "Unknown"
            if authors:
                f = authors[0]
                aff = f['AffiliationInfo'][0]['Affiliation'] if f.get('AffiliationInfo') else ""
                flag = "[🇨🇳CN]" if "China" in aff else ("[🇺🇸US]" if "USA" in aff else "")
                auth = f"{f.get('LastName')} {flag}"

            parsed_list.append({
                "title": title, "abstract": abstract, "journal": j_name,
                "if": if_val, "zone": zone, "sort_score": sort_score,
                "date": date_str, "author": auth, "tag": tag_label
            })
        except: continue
    return parsed_list

def analyze_with_ai(model, paper):
    print(f"🤖 AI 阅读: {paper['title'][:20]}...")
    
    prompt = f"""
    你是神经科学专家。请用中文简述这篇文献。
    标题: {paper['title']}
    摘要: {paper['abstract']}
    
    输出格式(Markdown):
    ### {paper['tag']} | {paper['title']}
    > 📅 {paper['date']} | 📖 {paper['journal']} (IF:{paper['if']}) | 👤 {paper['author']}
    - **🏷️ 类型**: [综述/实验/临床]
    - **🧐 核心**: (一句话结论,含数据)
    - **🔬 机制**: (通路/靶点)
    """
    
    # ✅ 关键修改：更稳健的重试逻辑
    for attempt in range(3):
        try:
            res = model.generate_content(prompt)
            # 成功后，强制休息 15 秒 (避免RPM超限)
            time.sleep(15) 
            return res.text.replace('\xa0', ' ')
        except Exception as e:
            err_str = str(e)
            if "429" in err_str:
                print(f"⚠️ 触发限流 (429)，冷却 60秒...")
                time.sleep(60) # 罚站 60s
            else:
                print(f"⚠️ 其他错误: {e}")
                time.sleep(5)
    
    return f"❌ {paper['title']} (分析失败)\n\n"

def format_simple_list(papers):
    """不经过AI，只列出标题，节省额度"""
    if not papers: return ""
    txt = "\n#### 📋 其他新收录文献 (仅列表)\n"
    for p in papers:
        txt += f"- **{p['date']}** | {p['title']} | *{p['journal']}*\n"
    return txt + "\n"

def main():
    model = setup_gemini()
    
    # Neuro
    neuro_ids = search_pubmed_ids('(Alzheimer\'s disease[Title/Abstract] AND (microglia[Title/Abstract] OR neuroinflammation[Title/Abstract]))', 30)
    neuro_papers = fetch_and_parse(neuro_ids, "🧠")
    neuro_papers.sort(key=lambda x: x['sort_score'], reverse=True)
    
    # TCM
    tcm_ids = search_pubmed_ids('((Traditional Chinese Medicine[Title/Abstract] OR Herbal[Title/Abstract] OR Acupuncture[Title/Abstract]) AND (Alzheimer[Title/Abstract] OR Brain[Title/Abstract]))', 10)
    tcm_papers = fetch_and_parse(tcm_ids, "🌿")
    tcm_papers.sort(key=lambda x: x['sort_score'], reverse=True)

    if not neuro_papers and not tcm_papers:
        print("📭 今日无数据")
        return

    # === 分级处理 ===
    # 1. 精选 (AI分析)
    top_neuro = neuro_papers[:MAX_AI_ANALYSIS_NEURO]
    top_tcm = tcm_papers[:MAX_AI_ANALYSIS_TCM]
    
    # 2. 列表 (仅标题)
    rest_neuro = neuro_papers[MAX_AI_ANALYSIS_NEURO:]
    rest_tcm = tcm_papers[MAX_AI_ANALYSIS_TCM:]

    content = f"🧠 NeuroBot 日报 ({datetime.date.today()})\n"
    content += f"⏱️ 检索范围: 过去 {SEARCH_WINDOW_DAYS} 天 | 📊 收录: {len(neuro_papers)+len(tcm_papers)} 篇\n\n"

    # TCM 板块
    if top_tcm:
        content += "--- 🌿 TCM 精选 ---\n\n"
        for p in top_tcm: content += analyze_with_ai(model, p)
        content += format_simple_list(rest_tcm)

    # Neuro 板块
    if top_neuro:
        content += "\n--- 🧠 Neuro 精选 ---\n\n"
        for p in top_neuro: content += analyze_with_ai(model, p)
        content += format_simple_list(rest_neuro)

    send_email(f"NeuroBot日报 - {datetime.date.today()}", content)

if __name__ == "__main__":
    main()
