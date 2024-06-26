# SSM DRUG-TARGET Affinity Prediction APP<br><br><br>

## Introduction

This repository contains the code and resources for using **SSM-DTA**, a Deep learning model developed by QIZHI PEI et al. in this research paper: https://academic.oup.com/bib/article/24/6/bbad386/7333673?login=false.

The provided code launches a user-friendly graphical interface built using **Streamlit**, an open-source app framework for Machine Learning and Data Science. This interface allows for easy **PREDICTION** of drug-target activity scores **(Ki, Kd, IC50, KIBA score)** and also enables the **EVALUATION** and benchmarking of the SSM-DTA model on labeled data. 
The backend of the application is powered by **FastAPI**.<br><br><br><br><br>


## A- Input Files Format ( **** IMPORTANT ****) read this before using the application

To use the SSM-DTA model, the input files must follow the required format as the examples files in **Examples** repository.

**Examples/ :**  This repository contains example files demonstrating the required format for the drug-target paired input data.

- example.mol: Example file containing molecule SMILES strings, with each SMILES on a separate line.
- example.pro: Example file containing protein sequences, with each sequence on a separate line.
- example.label: Example file containing binding affinity scores (labels) of the protein-target interaction,  with each score on a separate line. <br>

    
**Labels file :**  The labels file is mandatory for evaluation mode. In prediction mode, only the molecule and protein files are required.

The scores in the label file can be one of the following:

- **KI values**: Inhibition constant, should be in **nM** and transformed using **$-\log_{10}(\frac{Ki}{10^9})$**.
- **IC50 values**: Half maximal inhibitory concentration, should be in **nM** and transformed using **$-\log_{10}(\frac{IC50}{10^9})$**.
- **Kd values**: Dissociation constant, should be in **nM** and transformed using **$-\log_{10}(\frac{Kd}{10^9})$**.
- **KIBA scores**: An aggregation of IC50, Ki, and Kd measurements. KIBA scores should remain as they are and do not require any transformation.

If you possess SDF files or a CSV file containing CIDs, protein sequences, and affinity values, you can utilize the **formatfiles.py** script. This script will interactively guide you through the process and automatically convert your data into the appropriate format using the command provided below.

    python formatfiles.py
To run formatfiles.py you need the following libraries 
- rdkit
- pandas
- numpy
- pubchempy

<br><br><br><br>


## B- Getting Started <br><br>

**1- Clone this repository:**

    git clone https://YoussBGD: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@github.com/YoussBGD/ssm-dta.git

#### 2- Download fairseq and ckpt repertories  **** IMPORTANT ****
from the links below, donwload the repertories **fairseq.zip** and **ckpt.zip**, paste them in the **cloned directory ssm-dta** and unzip them **without changing their names**.

   - ckpt : https://drive.google.com/file/d/1bg2zBfAWx-fct8BRyEqznxqWexlI53-c/view?usp=drive_link<br>
   - fairseq :  https://drive.google.com/file/d/1ImPyGGYzTBIS-iq7Lt2vZCx-PbZeu9nU/view?usp=drive_link<br><br><br>

### To run the code with docker<br><br>

### Prerequisites

- Docker installed on your system<br><br>


**1- Navigate to the project directory:**

    cd ssm-dta

**2- Build the Docker image:**

    docker build -t ssm-dta .

**3- Run the Docker container:**

    docker run ssm-dta 

**- Or run the Docker container with more CPUs to increase the power.**

    docker run --cpus {CPUs number} ssm-dta

**4- To access the application, Open a web browser and navigate to the "Network URL" showed in the terminal after running the Docker container**.<br><br><br><br>  


### To run the code with conda<br><br>

### Prerequisites

- Conda installed on your system<br><br>


**1- Navigate to the project directory:**

    cd ssm-dta

**2- Create a new Conda environment:**

    conda create --name ssm-dta python=3.10 

**3- Activate the Conda environment:**

    conda activate ssm-dta

**4- Install the required dependencies:**

    pip install build-essential
    pip install -r requirements.txt

    

**5- Run the application:**

    python run.py


<br><br><br>


## C- How to Use the Interface

The SSM-DTA application provides a user-friendly web interface for predicting drug-target binding affinity. Here's a step-by-step guide on how to use the interface:

#### 1- Upload Input Files:
- Click on the "Browse files" button next to "Upload Molecules File" to select and upload a file containing molecule SMILES strings. (see example.mol)
- Click on the "Browse files" button next to "Upload Proteins File" to select and upload a file containing protein sequences. (see example.pro).
- If you want to perform evaluation, click on the "Browse files" button next to "Upload Labels File" to select and upload a file containing binding affinity scores (labels). (see example.label).<br><br>

#### 2- Select Prediction Mode:
-Choose the desired prediction mode from the "Select Mode" dropdown menu. There are two options:
- "Prediction": Use this mode if you only want to predict the binding affinity scores for the given molecule-protein pairs.
- "Evaluation": Use this mode if you have a labels file and want to evaluate the performance of the SSM-DTA model by comparing the predicted scores with the true scores.<br><br>
                you can download the **Benchmark/** directory from the link below, it contains test sets from DAVIS, KIBA, BindingDB_IC50 and BIndingDB_KI databases, you can use them directly for benchmarking the SSM-DTA model as they are alredy in the good format file.

- Benchmark databases: https://drive.google.com/file/d/1Ih2_UATxAkvFskAgCtO90P0Kj_b6S84P/view?usp=drive_link<br><br>

#### 3- Select Model: 
-Choose the appropriate model from the "Model" dropdown menu based on the type of binding affinity scores you want to predict or evaluate. The available options are:
- "DAVIS": Select this model if your want to predict dissociation constant (Kd) values.
- "KIBA ": Select this model if your want to predict KIBA scores (an aggregation of IC50, Ki, and Kd measurements).
- "BindingDB_IC50": Select this model if your want to predict half maximal inhibitory concentration (IC50) values.
- "BindingDB_Ki": Select this model if your want to predict inhibition constant (Ki) values.<br><br>

#### 4- Specify Output File Name: 
-Enter a name for the output file in the "Output File Name" text input field. This file will contain the predicted binding affinity scores.<br><br>

#### 5- Adjust Number of Threads (Optional):
-If you want to speed up the processing, you can increase the number of threads used by the application. The default value is 4, but you can adjust it based on the number of CPU cores available on your machine. Note that using too many threads may cause performance issues.<br><br>

#### 6- Start Processing: 
-Click the "**Process**" button to start the prediction or evaluation process. The application will display a progress message indicating that the processing has begun.<br><br>

#### 7- View Results: 
Once the processing is complete, the application will display the results.
- If you performed prediction, you will see a download link for the output file containing the predicted binding affinity scores and a table showing the predicted scores and the Ki, Kd or IC50 values. Click on the links to download the files. 
- If you performed evaluation, you will see the evaluation metrics such as MSE, RMSE, Pearson correlation, and C-index. Additionally, you will find a download link for the output file containing the predicted scores and a table comparing the true and predicted scores. <br><br>


That's it! You can repeat the process with different input files or adjust the settings as needed. <br><br><br>

## D- Overview of the code files 

-The SSM-DTA project consists of several Python files that work together to provide the functionality for drug-target binding affinity prediction. Here's an overview of each file:

- app.py: This file contains the FastAPI application that defines the API endpoints for processing files and handling file uploads.
- main.py: This file contains the main script that performs the preprocessing, feature extraction, and affinity prediction using the SSM-DTA model.
- streamlit_app.py: This file contains the Streamlit application that provides the user interface for interacting with the SSM-DTA model.
- run.py: This script is responsible for running both the FastAPI and Streamlit applications concurrently.

- requirements.txt: This file lists all the Python dependencies required for the SSM-DTA project
- Dockerfile: This file contains the instructions for building a Docker image of the SSM-DTA application.
