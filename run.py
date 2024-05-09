import subprocess
import sys
import os

def run_fastapi():
    fastapi_process = subprocess.Popen(["uvicorn", "app:app", "--reload", "--port", "8001"])
    fastapi_process.wait()

def run_streamlit():
    streamlit_process = subprocess.Popen(["streamlit", "run", "streamlit_app.py"])
    streamlit_process.wait()

if __name__ == "__main__":
    fastapi_process = subprocess.Popen(["uvicorn", "app:app", "--reload", "--port", "8001"])
    streamlit_process = subprocess.Popen(["streamlit", "run", "streamlit_app.py"])

    try:
        fastapi_process.wait()
        streamlit_process.wait()
    except KeyboardInterrupt:
        fastapi_process.terminate()
        streamlit_process.terminate()
        sys.exit(0)

