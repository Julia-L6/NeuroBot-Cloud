import os
import csv
import datetime
import time
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.header import Header
from Bio import Entrez

# --- 配置区域 ---
Entrez.email = "julia_light@msn.cn"
SMTP_SERVER = os.getenv("SMTP_SERVER")
SMTP_PORT = os.getenv("SMTP_PORT")
EMAIL_SENDER = os.getenv("EMAIL_SENDER")
raw_password = os.getenv("EMAIL_PASSWORD")
EMAIL_PASSWORD = raw_password.replace(' ', '').replace('\xa0', '').strip() if raw_password else None
EMAIL_RECEIVER = "julia_light@msn.cn"

TARGET_COUNT_NEURO = 40
TARGET_COUNT_TCM = 10
SEARCH_WINDOW_DAYS = 3 

# --- 核心功能函数 ---
def search_pubmed(query, max_ret):
    print(f"🔍 检索关键词: {query[:60]}...")
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_ret, sort="date", reldate=SEARCH_WINDOW_DAYS, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
        time.sleep(2)
        return record["IdList"]
    except Exception as e:
        print(f"⚠️ 检索失败: {e}")
        return []

def extract_author_info(article_dict):
    """提取第一作者和最后一位作者（假设为通讯）的信息"""
    authors = article_dict.get('AuthorList', [])
    if not authors:
        return "Unknown", "Unknown"
    
    def format_author(auth_entry):
        if not auth_entry: return ""
        name = f"{auth_entry.get('LastName', '')} {auth_entry.get('Initials', '')}"
        aff_info = auth_entry.get('AffiliationInfo', [])
        aff = aff_info[0]['Affiliation'] if aff_info else "No Affiliation"
        # 简化机构信息，只保留国家或前50个字符
        return f"{name} ({aff[:60]}...)"
    
    first_auth = format_author(authors[0])
    last_auth = format_author(authors[-1]) if len(authors) > 1 else ""
    
    return first_auth, last_auth

def fetch_details(id_list, category_label):
    if not id_list: return []
    print(f"📥 正在获取 {len(id_list)} 篇 [{category_label}] 的详情...")
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="xml")
        data = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"❌ 下载失败: {e}")
        return []

    results = []
    for article in data.get('PubmedArticle', []):
        try:
            med = article['MedlineCitation']
            art = med['Article']
            journal = art['Journal']['Title']
            title = art.get('ArticleTitle', 'No Title')
            
            d = art.get('ArticleDate', [])
            date_str = f"{d[0]['Year']}-{d[0]['Month']}-{d[0]['Day']}" if d else "Recent"
            
            abs_text = " ".join(art.get('Abstract', {}).get('AbstractText', []))
            
            # 提取 DOI
            doi = ""
            for id_obj in med.get('Article', {}).get('ELocationID', []):
                if id_obj.attributes.get('EIdType') == 'doi':
                    doi = str(id_obj)

            # --- 新增：提取作者信息 ---
            first_author, last_author = extract_author_info(art)

            results.append({
                "PMID": str(med['PMID']),
                "Date": date_str,
                "Journal": journal,
                "Title": title,
                "Category": category_label,
                "First_Author": first_author,  # 新增字段
                "Senior_Author": last_author,  # 新增字段
                "DOI": f"https://doi.org/{doi}" if doi else "",
                "Abstract": abs_text
            })
        except Exception as e: 
            print(f"解析错误: {e}")
            continue
    return results

def send_raw_data_email(filepath, summary_text):
    if not EMAIL_PASSWORD:
        print("❌ 邮箱密码缺失")
        return

    today = datetime.date.today()
    msg = MIMEMultipart()
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER
    msg['Subject'] = Header(f"NeuroBot 原始文献数据 - {today}", 'utf-8')
    msg.attach(MIMEText(summary_text, 'plain', 'utf-8'))

    try:
        with open(filepath, "rb") as f:
            attachment = MIMEApplication(f.read(), _subtype="csv")
            attachment.add_header('Content-Disposition', 'attachment', filename=os.path.basename(filepath))
            msg.attach(attachment)
    except Exception as e:
        print(f"❌ 附件打包失败: {e}")

    try:
        server = smtplib.SMTP(SMTP_SERVER, int(SMTP_PORT))
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        server.quit()
        print("✅ 包含原始数据的邮件已发送！")
    except Exception as e:
        print(f"❌ 邮件发送失败: {e}")

def main():
    print(f"🚀 开始采集当日原始文献数据 ({datetime.date.today()})...")
    
    ids_neuro_core = search_pubmed('(Alzheimer[Title/Abstract] AND microglia[Title/Abstract])', 30)
    ids_neuro_high = search_pubmed('(Alzheimer[Title/Abstract] AND (Nature[Journal] OR Science[Journal] OR Cell[Journal] OR Lancet[Journal]))', 20)
    ids_tcm = search_pubmed('((Traditional Chinese Medicine[Title/Abstract] OR herbal[Title/Abstract]) AND Brain[Title/Abstract])', 20)

    neuro_papers = fetch_details(list(set(ids_neuro_core + ids_neuro_high)), "Neuroscience")
    tcm_papers = fetch_details(ids_tcm, "TCM")
    all_data = neuro_papers[:TARGET_COUNT_NEURO] + tcm_papers[:TARGET_COUNT_TCM]

    if not all_data:
        print("📭 今日未检索到符合条件的文献")
        return

    filename = f"NeuroBot_Raw_{datetime.date.today()}.csv"
    if all_data:
        keys = all_data[0].keys()
        with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_data)

    summary = f"今日数据采集完成：\n- Neuro: {len(neuro_papers[:TARGET_COUNT_NEURO])} 篇\n- TCM: {len(tcm_papers[:TARGET_COUNT_TCM])} 篇\n\n含作者信息 CSV 已发送。"
    send_raw_data_email(filename, summary)

if __name__ == "__main__":
    main()
