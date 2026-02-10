import os
import csv
import datetime
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from Bio import Entrez

# --- 配置区 ---
Entrez.email = "julia_light@msn.cn"
SMTP_SERVER = os.getenv("SMTP_SERVER")
SMTP_PORT = os.getenv("SMTP_PORT")
EMAIL_SENDER = os.getenv("EMAIL_SENDER")
raw_password = os.getenv("EMAIL_PASSWORD")
EMAIL_PASSWORD = raw_password.replace(' ', '').strip() if raw_password else None
EMAIL_RECEIVER = "julia_light@msn.cn"

# 筛选目标数量
TARGET_COUNT_NEURO = 40
TARGET_COUNT_TCM = 10
SEARCH_DAYS = 2  # 检索最近2天

# --- 1. 加载期刊库 (用于识别高分文章) ---
def load_journal_db():
    db = {}
    if os.path.exists('journals.csv'):
        try:
            with open('journals.csv', 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                next(reader, None)
                for row in reader:
                    if len(row) >= 2:
                        # 格式: Name, IF, Zone
                        db[row[0].lower().strip()] = float(row[1]) if row[1].replace('.','').isdigit() else 0.0
        except:
            pass
    return db

JOURNAL_DB = load_journal_db()

def get_if(journal_name):
    """获取期刊IF，未匹配则返回0"""
    if not journal_name: return 0.0
    name = journal_name.lower().strip()
    return JOURNAL_DB.get(name, JOURNAL_DB.get(name.split('(')[0].strip(), 0.0))

# --- 2. PubMed 核心功能 ---
def search_pubmed(query, max_ret):
    """基础检索函数"""
    print(f"🔍 Searching: {query[:50]}...")
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_ret, sort="date", reldate=SEARCH_DAYS, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error searching: {e}")
        return []

def fetch_details(id_list):
    """批量下载详情"""
    if not id_list: return []
    print(f"📥 Fetching details for {len(id_list)} papers...")
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="xml")
        data = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error fetching: {e}")
        return []

    papers = []
    for article in data.get('PubmedArticle', []):
        try:
            med = article.get('MedlineCitation', {})
            art = med.get('Article', {})
            
            # 提取基础信息
            pmid = med.get('PMID', '')
            title = art.get('ArticleTitle', '')
            journal = art.get('Journal', {}).get('Title', 'Unknown')
            pub_date_list = art.get('ArticleDate', [])
            pub_date = f"{pub_date_list[0]['Year']}-{pub_date_list[0]['Month']}-{pub_date_list[0]['Day']}" if pub_date_list else "Recent"
            
            # 提取摘要
            abstract = " ".join(art.get('Abstract', {}).get('AbstractText', []))
            
            # 提取DOI
            doi = ""
            for el in med.get('Article', {}).get('ELocationID', []):
                if el.attributes.get('EIdType') == 'doi':
                    doi = str(el)
            
            # 提取作者
            authors_list = art.get('AuthorList', [])
            first_author = ""
            if authors_list:
                fa = authors_list[0]
                aff = fa.get('AffiliationInfo', [{}])[0].get('Affiliation', '')
                first_author = f"{fa.get('LastName')} {fa.get('ForeName')} ({aff[:30]}...)"

            # 计算IF
            impact_factor = get_if(journal)

            papers.append({
                'PMID': str(pmid),
                'Title': title,
                'Journal': journal,
                'IF': impact_factor,
                'Date': pub_date,
                'FirstAuthor': first_author,
                'DOI': doi,
                'Abstract': abstract,
                'Category': 'Unsorted'
            })
        except:
            continue
    return papers

# --- 3. 核心筛选逻辑 (您的定制需求) ---
def logic_neuro_process():
    """
    逻辑:
    1. AD + Microglia (最高优先级)
    2. AD (高分期刊 IF>10)
    3. Neurodegen + Microglia (高分期刊 IF>10)
    """
    pool = []
    
    # 检索策略 A: AD + Microglia (核心关注)
    ids_a = search_pubmed('(Alzheimer\'s disease[Title/Abstract] AND (microglia[Title/Abstract] OR neuroinflammation[Title/Abstract]))', 50)
    papers_a = fetch_details(ids_a)
    for p in papers_a: 
        p['Category'] = 'Neuro_Core(AD+Microglia)'
        p['Priority'] = 100 + p['IF'] # 基础分100，叠加IF
    
    # 检索策略 B: AD Broad (仅关注高分)
    # 为了避免检索量过大，我们只取前50篇最新的，然后通过代码滤掉低分的
    ids_b = search_pubmed('(Alzheimer\'s disease[Title/Abstract])', 60)
    papers_b = fetch_details(ids_b)
    for p in papers_b:
        if p['IF'] >= 7.0: # 设定子刊级别门槛，例如 IF>7
            p['Category'] = 'Neuro_HighImpact(AD)'
            p['Priority'] = 80 + p['IF']
        else:
            p['Priority'] = 0 # 丢弃

    # 检索策略 C: Other Neuro + Microglia (仅关注高分)
    ids_c = search_pubmed('((Parkinson[Title/Abstract] OR ALS[Title/Abstract] OR Neurodegenerative[Title/Abstract]) AND microglia[Title/Abstract])', 40)
    papers_c = fetch_details(ids_c)
    for p in papers_c:
        if p['IF'] >= 7.0:
            p['Category'] = 'Neuro_HighImpact(Other+Microglia)'
            p['Priority'] = 70 + p['IF']
        else:
            p['Priority'] = 0

    # 合并去重
    unique_pool = {p['PMID']: p for p in (papers_a + papers_b + papers_c) if p.get('Priority', 0) > 0}
    final_list = list(unique_pool.values())
    
    # 排序：按优先级（即核心程度+IF）倒序
    final_list.sort(key=lambda x: x['Priority'], reverse=True)
    return final_list[:TARGET_COUNT_NEURO]

def logic_tcm_process():
    """
    逻辑: TCM 领域，优先高分
    """
    ids = search_pubmed('((Traditional Chinese Medicine[Title/Abstract] OR Herbal[Title/Abstract] OR Acupuncture[Title/Abstract]) AND (Brain[Title/Abstract] OR Alzheimer[Title/Abstract]))', 40)
    papers = fetch_details(ids)
    
    for p in papers:
        p['Category'] = 'TCM'
        p['Priority'] = p['IF'] # 纯按IF排序
        
    papers.sort(key=lambda x: x['Priority'], reverse=True)
    return papers[:TARGET_COUNT_TCM]

# --- 4. 生成 CSV 与 发送邮件 ---
def save_csv(papers, filename):
    headers = ['PMID', 'Category', 'IF', 'Date', 'Journal', 'Title', 'DOI', 'FirstAuthor', 'Abstract']
    with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for p in papers:
            # 仅写入需要的字段
            clean_p = {k: p.get(k, '') for k in headers}
            writer.writerow(clean_p)
    return filename

def send_email_with_attachment(subject, body, filepath):
    msg = MIMEMultipart()
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER
    msg['Subject'] = Header(subject, 'utf-8')
    
    msg.attach(MIMEText(body, 'plain', 'utf-8'))
    
    # 添加附件
    with open(filepath, 'rb') as f:
        part = MIMEApplication(f.read(), Name=os.path.basename(filepath))
    part['Content-Disposition'] = f'attachment; filename="{os.path.basename(filepath)}"'
    msg.attach(part)
    
    try:
        s = smtplib.SMTP(SMTP_SERVER, int(SMTP_PORT))
        s.starttls()
        s.login(EMAIL_SENDER, EMAIL_PASSWORD)
        s.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        s.quit()
        print("✅ Email sent successfully.")
    except Exception as e:
        print(f"❌ Email failed: {e}")

def main():
    print("🚀 Starting Raw Data Collection...")
    
    # 1. 获取数据
    neuro_list = logic_neuro_process()
    tcm_list = logic_tcm_process()
    all_papers = neuro_list + tcm_list
    
    print(f"📊 Collected: {len(neuro_list)} Neuro + {len(tcm_list)} TCM = {len(all_papers)} Total")
    
    if not all_papers:
        print("No papers found today.")
        return

    # 2. 保存 CSV
    today_str = datetime.date.today().strftime('%Y-%m-%d')
    csv_filename = f"NeuroBot_Raw_{today_str}.csv"
    save_csv(all_papers, csv_filename)
    
    # 3. 发送邮件
    email_body = f"""
    NeuroBot 原始数据采集报告
    日期: {today_str}
    
    包含数据:
    - Neuro/AD 核心: {len(neuro_list)} 篇
    - TCM 精选: {len(tcm_list)} 篇
    
    请查收附件 CSV 文件。
    """
    send_email_with_attachment(f"NeuroBot数据源 ({today_str})", email_body, csv_filename)

if __name__ == "__main__":
    main()
