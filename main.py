import os
import datetime
import time
import smtplib
from email.mime.text import MIMEText
from email.header import Header
import google.generativeai as genai
from Bio import Entrez

# --- 1. 配置区域 ---
# ✅ 这里填写你的真实邮箱即可，符合 NCBI 规范
Entrez.email = "lifuyaojulia@gmail.com"
 
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
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

# --- 2. 核心期刊数据库 (包含 Neuro & TCM 高分刊) ---
# 用于粗略判断 IF，未收录的期刊将标记为 Unknown
JOURNAL_DB = {
    "nature": 64.8, "science": 56.9, "cell": 45.5,
    "nature neuroscience": 21.2, "neuron": 16.2, 
    "nature reviews neuroscience": 38.7, "molecular neurodegeneration": 15.1,
    "alzheimer's & dementia": 14.0, "brain": 10.6, "acta neuropathologica": 12.7,
    "journal of neuroinflammation": 9.3, "glia": 8.0, 
    "molecular psychiatry": 11.0, "biological psychiatry": 10.6,
    "autophagy": 13.3, "redox biology": 10.7,
    # TCM & Pharma
    "phytomedicine": 6.7, "journal of ethnopharmacology": 5.4,
    "pharmacological research": 9.3, "british journal of pharmacology": 7.3,
    "chinese medicine": 4.9, "journal of advanced research": 10.7
}

def get_journal_if(journal_name):
    """模糊匹配获取 IF"""
    if not journal_name: return 0.0
    name_clean = journal_name.lower()
    for k, v in JOURNAL_DB.items():
        if k in name_clean:
            return v
    return 0.0 # 未知期刊默认 0 分，排在后面

def setup_gemini():
    if not GOOGLE_API_KEY:
        print("❌ 错误: 缺少 GOOGLE_API_KEY")
        return None
    genai.configure(api_key=GOOGLE_API_KEY)
    return genai.GenerativeModel('gemini-flash-latest')

def search_pubmed_ids(query, max_results):
    """基础检索函数，返回 ID 列表"""
    print(f"🔍 检索策略: {query} (Max: {max_results})")
    try:
        # reldate=7: 过去7天; datetype='pdat': 包含Online First和正式出版
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="date", reldate=7, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"⚠️ 检索失败: {e}")
        return []

def fetch_details(id_list):
    """批量下载文献详情"""
    if not id_list: return []
    print(f"📥 下载 {len(id_list)} 篇文献元数据...")
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="xml")
        papers = Entrez.read(handle)
        handle.close()
        return papers['PubmedArticle']
    except Exception as e:
        print(f"❌ 下载详情失败: {e}")
        return []

def parse_paper(article, section_tag):
    """解析 XML，提取关键元数据"""
    try:
        medline = article['MedlineCitation']
        art = medline['Article']
        
        title = art.get('ArticleTitle', 'No Title')
        abstract_list = art.get('Abstract', {}).get('AbstractText', [])
        abstract = " ".join(abstract_list) if abstract_list else "No Abstract"
        journal = art['Journal'].get('Title', 'Unknown')
        
        # 提取日期
        pub_date = art.get('ArticleDate', [])
        date_str = "Recent"
        if pub_date:
            date_str = f"{pub_date[0].get('Year')}-{pub_date[0].get('Month')}-{pub_date[0].get('Day')}"
            
        # 提取通讯/一作国家
        authors = art.get('AuthorList', [])
        first_auth_str = "Unknown"
        corresp_auth_str = "Unknown"
        
        if authors:
            f = authors[0]
            first_auth_str = f"{f.get('LastName')} {f.get('ForeName')} "
            if f.get('AffiliationInfo'):
                # 提取国家 (简单正则逻辑)
                aff = f['AffiliationInfo'][0]['Affiliation']
                if "China" in aff: first_auth_str += "[🇨🇳China]"
                elif "USA" in aff: first_auth_str += "[🇺🇸USA]"
                else: first_auth_str += "[🌍Global]"

            l = authors[-1]
            corresp_auth_str = f"{l.get('LastName')} {l.get('ForeName')}"

        if_score = get_journal_if(journal)
        
        return {
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "if": if_score,
            "date": date_str,
            "authors": f"1st: {first_auth_str} | Rep: {corresp_auth_str}",
            "tag": section_tag # 标记是 Neuro 还是 TCM
        }
    except:
        return None

def analyze_with_ai(model, paper):
    print(f"🤖 AI 阅读 ({paper['tag']}): {paper['title'][:30]}...")
    
    # --- 资深学者专用 Prompt ---
    prompt = f"""
    你是一位神经科学领域的资深审稿人。请根据摘要快速提炼关键信息。
    
    文章标题: {paper['title']}
    摘要: {paper['abstract']}
    
    请严格按照以下 Markdown 格式输出（不要使用不换行空格）：
    
    ### {paper['tag']} | {paper['title']}
    > 📅 {paper['date']} | 📖 {paper['journal']} (IF: {paper['if']})
    > 👥 {paper['authors']}
    
    - **🏷️ 研究类型**: [请判断: 综述 / 临床试验 / 动物实验(In Vivo) / 细胞实验(In Vitro) / 数据挖掘]
    - **🧐 核心发现**: (一句话总结，尽量包含定量数据，如 p<0.05 或 变化幅度)
    - **🔬 关键机制**: (如果是TCM请指出具体化合物和靶点；如果是Neuro请指出通路)
    - **💡 简评**: (一句话评价其实际意义)
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
    msg = MIMEText(content.replace('\xa0', ' '), 'plain', 'utf-8')
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
        print(f"❌ 邮件发送失败: {e}")

def main():
    model = setup_gemini()
    if not model: return

    # --- 策略 A: 神经科学核心 (Target: AD, Microglia) ---
    # 逻辑: (AD AND Microglia) 优先，其次泛 AD
    query_neuro = '(Alzheimer\'s disease[Title/Abstract] AND (microglia[Title/Abstract] OR neuroinflammation[Title/Abstract]))'
    ids_neuro = search_pubmed_ids(query_neuro, 30) # 多抓一点用来按 IF 排序

    # --- 策略 B: TCM 专区 (Target: TCM + Neuro) ---
    # 逻辑: (中药/针灸) AND (神经/脑/AD)
    query_tcm = '((Traditional Chinese Medicine[Title/Abstract] OR Herbal[Title/Abstract] OR Acupuncture[Title/Abstract] OR formula[Title/Abstract]) AND (Neuroscience[Title/Abstract] OR Alzheimer[Title/Abstract] OR Brain[Title/Abstract]))'
    ids_tcm = search_pubmed_ids(query_tcm, 10) # 抓 10 篇选 5 篇

    # --- 下载详情 ---
    # 去重
    all_ids = list(set(ids_neuro + ids_tcm))
    if not all_ids:
        print("📭 今日无匹配文献")
        return

    raw_xmls = fetch_details(all_ids)
    
    neuro_papers = []
    tcm_papers = []

    # --- 分类与解析 ---
    for xml in raw_xmls:
        # 判断这篇文章属于哪个 ID 列表 (粗略反推)
        # PubMed返回的顺序可能乱，这里根据 Title 简单二次归类，或者直接解析
        # 为了简单，我们先把所有解析出来，再根据 ID 归类，但 XML 里很难直接对 ID
        # 简单策略：先解析，然后根据内容关键词打标
        p = parse_paper(xml, "待定")
        if not p: continue
        
        # 简单的关键词打标逻辑
        text_for_tag = (p['title'] + p['abstract']).lower()
        is_tcm = any(k in text_for_tag for k in ['chinese medicine', 'herbal', 'acupuncture', 'decoction', 'ginsenoside'])
        
        if is_tcm:
            p['tag'] = "🌿 [TCM]"
            tcm_papers.append(p)
        else:
            p['tag'] = "🧠 [Neuro]"
            neuro_papers.append(p)

    # --- 排序与截取 ---
    # 1. Neuro: 按 IF 降序，取前 20
    neuro_papers.sort(key=lambda x: x['if'], reverse=True)
    final_neuro = neuro_papers[:20]
    
    # 2. TCM: 按 IF 降序 (或日期)，取前 5
    tcm_papers.sort(key=lambda x: x['if'], reverse=True)
    final_tcm = tcm_papers[:5]
    
    total_papers = final_neuro + final_tcm
    print(f"📊 最终入选: Neuro {len(final_neuro)} 篇, TCM {len(final_tcm)} 篇")

    # --- 生成报告 ---
    content = f"🧠 NeuroBot 深度简报 ({datetime.date.today()})\n"
    content += f"🎯 策略: Neuro (IF Sort) + TCM (Special)\n"
    content += f"📌 总计: {len(total_papers)} 篇 | 来源: PubMed Online First\n\n"
    
    # 先展示 TCM (如果有)，让你眼前一亮
    if final_tcm:
        content += "--- 🌿 TCM 特别关注 ---\n\n"
        for p in final_tcm:
            content += analyze_with_ai(model, p)
            time.sleep(3)
            
    content += "\n--- 🧠 神经科学/AD 核心文献 ---\n\n"
    for p in final_neuro:
        content += analyze_with_ai(model, p)
        time.sleep(3)

    send_email(f"NeuroBot 深度日报 - {datetime.date.today()}", content)

if __name__ == "__main__":
    main()
