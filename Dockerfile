FROM python:3.9.2
WORKDIR /

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

CMD ["streamlit", "app.py"]
