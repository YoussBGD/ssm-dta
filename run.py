import subprocess
import sys
import os


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

