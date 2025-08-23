FROM python:3.9.2

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --branch v1.3.6 https://github.com/karstenkunneman/P2F-HydrateCalc.git /app
WORKDIR /app

RUN pip install -r requirements.txt
RUN pip install streamlit && streamlit --1.40.1

EXPOSE 8000

HEALTHCHECK CMD curl --fail http://localhost:8000/_stcore/health || exit 1

CMD ["streamlit", "run", "app.py", "--server.port=8000", "--server.address=0.0.0.0"]