FROM python:3.10

WORKDIR /app


# Install system-level dependencies
RUN apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*
# Copy the requirements.txt file into the current working directory (/app)
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy all content from the current directory into /app in the container
COPY . .

EXPOSE 8001

CMD ["python", "run.py"]
