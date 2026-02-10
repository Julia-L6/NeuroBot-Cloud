import os
import datetime
import time
import smtplib
from email.mime.text import MIMEText
from email.header import Header
import google.generativeai as genai
import arxiv # 导入整个 arxiv 库
from arxiv import Search, SortCriterion

# --- 1. 获取密钥 ---
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
SMTP_SERVER = os.getenv("SMTP_SERVER")
SMTP_PORT = os.getenv("SMTP_PORT")
EMAIL_SENDER = os.getenv("EMAIL_SENDER")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")

# --- 2. 接收邮箱 ---
EMAIL_RECEIVER = "julia_light@msn.cn"

def setup_gemini():
    if not GOOGLE_API_KEY:
        print("❌ 错误: 缺少 GOOGLE_API_KEY")
        return None
    genai.configure(api_key=GOOGLE_API_KEY)
    return genai.GenerativeModel('gemini-flash-latest')

def get_latest_papers(topics):
    print(f"🔍 正在检索: {topics}")
    
    # ✅ 改进点1：配置 ArXiv 客户端，设置延迟和重试
    client = arxiv.Client(
        page_size=3,
        delay_seconds=10.0, # 每次请求强制间隔10秒，对服务器更友好
        num_retries=5       # 库内部自动重试5次
    )

    search = Search(
        query=topics,
        max_results=3, 
        sort_by=SortCriterion.SubmittedDate
    )
    
    # ✅ 改进点2：外层手动重试循环，专门对抗 429 错误
    for attempt in range(3): # 给它3次“死里复活”的机会
        try:
            # 使用新版写法 client.results
            return list(client.results(search))
        except Exception as e:
            print(f"⚠️ 检索遭遇拥堵 (尝试 {attempt+1}/3): {e}")
            if "429" in str(e):
                print("⏳ 触发 ArXiv 限流，强制休眠 30 秒...")
                time.sleep(30) # 休息30秒再试
            else:
                time.sleep(5)
    
    print("❌ 三次尝试均失败，今日跳过检索。")
    return []

def analyze_paper(model, paper):
    print(f"🤖 正在阅读: {paper.title}")
    prompt = f"""
    你是一位神经科学专家。请阅读摘要并写出中文日报。
    
    标题: {paper.title}
    摘要: {paper.summary}
    
    请严格使用以下 Markdown 格式:
    ### 📄 {paper.title}
    > *{paper.published.strftime('%Y-%m-%d')}*
    - **🧐 核心发现**: (一句话概括)
    - **🔬 机制/方法**: (关键技术或靶点)
    - **💡 创新点**: (1-2个亮点)
    - **🔗 原文**: {paper.entry_id}
    ---
    """
    for _ in range(3):
        try:
            response = model.generate_content(prompt)
            # 清洗特殊字符
            safe_text = response.text.replace('\xa0', ' ')
            return safe_text
        except Exception as e:
            if "429" in str(e): time.sleep(20)
            else: return f"❌ Error: {str(e)}\n\n"
    return "❌ Error: Rate Limit\n\n"

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
        print(f"✅ 邮件已成功发送至 {EMAIL_RECEIVER}")
    except Exception as e:
        print(f"❌ 发送失败: {repr(e)}")

def main():
    model = setup_gemini()
    if not model: return
    
    keywords = '(ti:"Alzheimer" OR abs:"Alzheimer") AND (ti:"microglia" OR abs:"microglia")'
    papers = get_latest_papers(keywords)
    
    if not papers:
        # 如果是因为网络拥堵没拿到论文，发个邮件通知一下，而不是报错
        print("📭 本次运行未获取到论文 (可能是无新文，也可能是网络拥堵)")
        return

    content = f"🧠 NeuroBot 日报 ({datetime.date.today()})\n\n"
    for paper in papers:
        summary = analyze_paper(model, paper)
        content += summary
        time.sleep(5)
        
    send_email(f"NeuroBot日报 - {datetime.date.today()}", content)

if __name__ == "__main__":
    main()
