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



