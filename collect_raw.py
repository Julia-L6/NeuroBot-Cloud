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

# --- 1. 配置区域 ---
Entrez.email = "julia_light@msn.cn"
SMTP_SERVER = os.getenv("SMTP_SERVER")
SMTP_PORT = os.getenv("SMTP_PORT")
EMAIL_SENDER = os.getenv("EMAIL_SENDER")
raw_password = os.getenv("EMAIL_PASSWORD")
EMAIL_PASSWORD = raw_password.replace(' ', '').replace('\xa0', '').strip() if raw_password else None
EMAIL_RECEIVER = "julia_light@msn.cn"

# 检索配置
TARGET_COUNT_NEURO = 40
TARGET_COUNT_TCM = 10
SEARCH_WINDOW_DAYS = 3 # 检索近3天，确保不漏掉周一的数据

# --- 2. 核心功能函数 ---
def search_pubmed(query, max_ret):
    print(f"🔍 检索关键词: {query[:60]}...")
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_ret, sort="date", reldate=SEARCH_WINDOW_DAYS)
        record = Entrez.read(handle)
        handle.close()
        time.sleep(2) # 强制延迟，规避限流
        return record["IdList"]
    except Exception as e:
        print(f"⚠️ 检索失败: {e}")
        return []

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
            
            # 提取日期
            d = art.get('ArticleDate', [])
            date_str = f"{d[0]['Year']}-{d[0]['Month']}-{d[0]['Day']}" if d else "Recent"
            
            # 提取摘要
            abs_text = " ".join(art.get('Abstract', {}).get('AbstractText', []))
            
            # 提取DOI
            doi = ""
            for id_obj in med.get('Article', {}).get('ELocationID', []):
                if id_obj.attributes.get('EIdType') == 'doi':
                    doi = str(id_obj)

            results.append({
                "PMID": str(med['PMID']),
                "Date": date_str,
                "Journal": journal,
                "Title": title,
                "Category": category_label,
                "DOI": f"https://doi.org/{doi}" if doi else "",
                "Abstract": abs_text
            })
        except: continue
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

    # 邮件正文
    msg.attach(MIMEText(summary_text, 'plain', 'utf-8'))

    # 添加附件
    try:
        with open(filepath, "rb") as f:
            attachment = MIMEApplication(f.read(), _subtype="csv")
            attachment.add_header('Content-Disposition', 'attachment', filename=os.path.basename(filepath))
            msg.attach(attachment)
    except Exception as e:
        print(f"❌ 附件打包失败: {e}")

    # 发送
    try:
        server = smtplib.SMTP(SMTP_SERVER, int(SMTP_PORT))
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        server.quit()
        print("✅ 包含原始数据的邮件已发送！")
    except Exception as e:
        print(f"❌ 邮件发送失败: {e}")

# --- 3. 主逻辑 ---
def main():
    print(f"🚀 开始采集当日原始文献数据 ({datetime.date.today()})...")

    # 逻辑 1: Neuro 领域 (AD + Microglia)
    ids_neuro_core = search_pubmed('(Alzheimer[Title/Abstract] AND microglia[Title/Abstract])', 30)
    # 逻辑 2: Neuro 领域 (其他 AD 高分 - 简化为关键词检索)
    ids_neuro_high = search_pubmed('(Alzheimer[Title/Abstract] AND (Nature[Journal] OR Science[Journal] OR Cell[Journal] OR Lancet[Journal]))', 20)
    
    # 逻辑 3: TCM 领域
    ids_tcm = search_pubmed('((Traditional Chinese Medicine[Title/Abstract] OR herbal[Title/Abstract]) AND Brain[Title/Abstract])', 20)

    # 汇总并去重
    neuro_papers = fetch_details(list(set(ids_neuro_core + ids_neuro_high)), "Neuroscience")
    tcm_papers = fetch_details(ids_tcm, "TCM")
    all_data = neuro_papers[:TARGET_COUNT_NEURO] + tcm_papers[:TARGET_COUNT_TCM]

    if not all_data:
        print("📭 今日未检索到符合条件的文献")
        return

    # 生成 CSV
    filename = f"NeuroBot_Raw_{datetime.date.today()}.csv"
    keys = all_data[0].keys()
    with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(all_data)

    # 发送
    summary = f"今日文献数据采集完成：\n- Neuro 领域: {len(neuro_papers[:TARGET_COUNT_NEURO])} 篇\n- TCM 领域: {len(tcm_papers[:TARGET_COUNT_TCM])} 篇\n\n数据已作为 CSV 附件发送，请下载后使用本地 AI 或 Zotero 处理。"
    send_raw_data_email(filename, summary)

if __name__ == "__main__":
    main()
