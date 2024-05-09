from fastapi import FastAPI, UploadFile, File, HTTPException
from pydantic import BaseModel
from typing import Optional  
import main as script

app = FastAPI()

class InputData(BaseModel):
    molecules: str
    proteins: str
    mode: str
    labels: Optional[str] = None 
    model: str
    output_fn: str
    batch_size: int = 1
    num_threads: int = 4

@app.post("/process")
async def process_files(input_data: InputData):
    args = input_data.dict()
    
    if args["mode"] == "Evaluation" and args["labels"] is None:
        raise HTTPException(status_code=400, detail="Labels file is required for evaluation mode.")
    
    results = script.main(args)
    return {"results": results}

@app.post("/upload_molecules")
async def upload_molecules(file: UploadFile = File(...)):
    with open(file.filename, "wb") as f:
        f.write(await file.read())
    return {"message": f"Molecules file {file.filename} uploaded successfully"}

@app.post("/upload_proteins")
async def upload_proteins(file: UploadFile = File(...)):
    with open(file.filename, "wb") as f:
        f.write(await file.read())
    return {"message": f"Proteins file {file.filename} uploaded successfully"}

@app.post("/upload_labels")
async def upload_labels(file: UploadFile = File(...)):
    with open(file.filename, "wb") as f:
        f.write(await file.read())
    return {"message": f"Labels file {file.filename} uploaded successfully"}

