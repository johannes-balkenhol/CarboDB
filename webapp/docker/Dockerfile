FROM python:3.11-slim
RUN apt-get update && apt-get install -y hmmer wget && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN python -c "import esm; esm.pretrained.esm2_t33_650M_UR50D()" || echo "ESM-2 will download at startup"
COPY app/ ./app/
COPY static/ ./static/
COPY models/ ./models/
RUN mkdir -p /app/jobs
EXPOSE 8080
ENV MODELS_DIR=/app/models JOBS_DIR=/app/jobs ESM2_DEVICE=cpu
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080", "--workers", "2"]
