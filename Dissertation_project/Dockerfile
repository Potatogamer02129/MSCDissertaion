# 1. Use an official, lightweight Python 3.12 image
FROM python:3.12

# 2. Set the working directory inside the container
WORKDIR /app

#RUN apt-get update && apt-get install -y \
#    libxrender1 \
#    libxext6 \
#    libsm6 \
#    && rm -rf /var/lib/apt/lists/*

# 3. Copy the requirements file and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copy all your actual code (main.py, pages/, etc.) into the container
COPY . .

# 5. Tell Docker that Streamlit runs on port 8501
EXPOSE 8501

# 6. The command to run your app when the container starts
CMD ["streamlit", "run", "main.py", "--server.port=8501", "--server.address=0.0.0.0"]