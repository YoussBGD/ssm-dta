import streamlit as st
import requests
import asyncio
import nest_asyncio
import numpy as np
import base64
import pandas as pd
import os


#nest_asyncio.apply()

st.title("SSM DRUG-TARGET Affinity Prediction")

current_dir = os.path.dirname(os.path.abspath(__file__))

molecules_file = st.file_uploader("Upload Molecules File")
proteins_file = st.file_uploader("Upload Proteins File")
labels_file = st.file_uploader("Upload Labels File (required for Evaluation)")



if molecules_file:
    molecules_path = os.path.join(current_dir, molecules_file.name)
    with open(molecules_path, "wb") as f:
        f.write(molecules_file.getvalue())
    st.success(f"Molecules file {molecules_file.name} uploaded successfully")

if proteins_file:
    proteins_path = os.path.join(current_dir, proteins_file.name)
    with open(proteins_path, "wb") as f:
        f.write(proteins_file.getvalue())
    st.success(f"Proteins file {proteins_file.name} uploaded successfully")

if labels_file:
    labels_path = os.path.join(current_dir, labels_file.name)
    with open(labels_path, "wb") as f:
        f.write(labels_file.getvalue())
    st.success(f"Labels file {labels_file.name} uploaded successfully")

mode = st.selectbox("Select Mode", ["Prediction", "Evaluation"])

model_descriptions = {
    "DAVIS": " To predict Kd",
    "KIBA": "To predict KIBA score (Aggregation of IC50, Ki and Kd measurements)",
    "BindingDB_IC50": " To predict IC50",
    "BindingDB_Ki": "To predict Ki"
}

model_options = [f"{model} - {description}" for model, description in model_descriptions.items()]

selected_model_option = st.selectbox("Model", model_options)

model = selected_model_option.split(" - ")[0]

output_fn = st.text_input("Output File Name (required)")
batch_size = 1
num_threads = st.number_input("Number of Threads (increasing this may speed up processing)", value=4, min_value=1, max_value=16, step=1,
                              help="Recommended: Use a number of threads equal to or slightly higher than the number of CPU cores available on your machine. Using too many threads may cause performance issues.")


def download_link(object_to_download, download_filename, download_link_text):
    if isinstance(object_to_download, pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)

    b64 = base64.b64encode(object_to_download.encode()).decode()

    return f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>'

if st.button("Process", disabled=(not output_fn)):
    if not molecules_file or not proteins_file:
        st.error("Please upload both molecules and proteins files.")
    elif mode.lower() == "evaluation" and not labels_file:
        st.error("Please upload the labels file for evaluation mode.")
    else:
        input_data = {
            "molecules": molecules_file.name if molecules_file else None,
            "proteins": proteins_file.name if proteins_file else None,
            "mode": mode.lower(),
            "labels": labels_file.name if labels_file else None,
            "model": model,
            "output_fn": output_fn,
            "batch_size": batch_size,
            "num_threads": num_threads
        }
        st.info("Processing in progress. Please check the terminal where you executed the code to follow the progress.")

        response = requests.post("http://localhost:8001/process", json=input_data)
        
        if response.status_code == 200:
            results = response.json()["results"]

            if "error" in results:
                st.error(f"Error: {results['error']}")
            else:
                st.header(f"Prediction with model trained on {model} database")

                st.subheader("Download Output File")
                st.markdown(f'<a href="data:file/txt;base64,{results["output_file"]}" download="{output_fn}">Click here to download {output_fn} (a text file containing predicted affinity scores of your molecules-target pairs )</a>', unsafe_allow_html=True)

                if mode.lower() == "evaluation":
                    if "mse" in results:
                        st.subheader("\nModel Evaluation Results")
                        st.write(f"MSE: {results['mse']}")
                        st.write(f"RMSE: {results['rmse']}")
                        st.write(f"Pearson: {results['pearson']}")
                        st.write(f"C-index: {results['c_index']}")

                        st.subheader("\nEvaluation Scores")
                        if model == "DAVIS":
                            st.latex(r"Score = -\log_{10}\left(\frac{K_d}{10^9}\right)")
                        elif model == "BindingDB_IC50":
                            st.latex(r"Score = -\log_{10}\left(\frac{IC_{50}}{10^9}\right)")
                        elif model == "BindingDB_Ki":
                            st.latex(r"Score = -\log_{10}\left(\frac{K_i}{10^9}\right)")
                        elif model == "KIBA":
                            st.write("Score = KIBA score")

                        true_scores = results["true_scores"]
                        pred_scores = results["pred_scores"]
                        
                        transformed_true_scores = []
                        transformed_pred_scores = []
                        
                        if model in ["DAVIS", "BindingDB_IC50", "BindingDB_Ki"]:
                            transformed_true_scores = 10 ** (9 - np.array(true_scores))
                            transformed_pred_scores = 10 ** (9 - np.array(pred_scores))
                        
                        if model == "DAVIS":
                            score_name = "Kd"
                        elif model == "BindingDB_IC50":
                            score_name = "IC50"
                        elif model == "BindingDB_Ki":
                            score_name = "Ki"
                        else:
                            score_name = "KIBA"
                        
                        eval_scores = {
                            "True Scores": true_scores,
                            "Predicted Scores": pred_scores,
                        }
                        
                        if len(transformed_true_scores) > 0 and len(transformed_pred_scores) > 0:
                            eval_scores[f"True {score_name} (nM)"] = transformed_true_scores
                            eval_scores[f"Predicted {score_name} (nM)"] = transformed_pred_scores

                        eval_scores_df = pd.DataFrame(eval_scores)
                        st.markdown(download_link(eval_scores_df, f"{output_fn}_eval_scores.csv", "Download this table"), unsafe_allow_html=True)
                        st.table(eval_scores_df)
                    else:
                        st.error("Evaluation results not found. Please make sure you have provided the correct labels file.")
                else:
                    st.subheader("\nPrediction Scores")
                    if model == "DAVIS":
                        st.latex(r"Score = -\log_{10}\left(\frac{K_d}{10^9}\right)")
                    elif model == "BindingDB_IC50":
                        st.latex(r"Score = -\log_{10}\left(\frac{IC_{50}}{10^9}\right)")
                    elif model == "BindingDB_Ki":
                        st.latex(r"Score = -\log_{10}\left(\frac{K_i}{10^9}\right)")
                    elif model == "KIBA":
                        st.write("Score = KIBA score (An aggregation of IC50, Ki and Kd measurements)")
                    
                    pred_scores = results["pred_scores"]
                    transformed_scores = []
                        
                    if model in ["DAVIS", "BindingDB_IC50", "BindingDB_Ki"]:
                        transformed_scores = 10 ** (9 - np.array(pred_scores))
                    
                    if model == "DAVIS":
                        score_name = "Kd"
                    elif model == "BindingDB_IC50":
                        score_name = "IC50"
                    elif model == "BindingDB_Ki":
                        score_name = "Ki"
                    else:
                        score_name = "KIBA"
                    
                    scores_table = {
                        "Predicted Scores": pred_scores,
                    }
                    
                    if len(transformed_scores) > 0:
                        scores_table[f"Predicted {score_name} (nM)"] = transformed_scores
                    
                    scores_table_df = pd.DataFrame(scores_table)
                    st.markdown(download_link(scores_table_df, f"{output_fn}_pred_scores.csv", "Download this table"), unsafe_allow_html=True)
                    st.table(scores_table_df)
                    
        elif response.status_code == 400:
            st.error(f"Error: {response.json()['detail']}")
        else:
            st.error(f"Error processing files: {response.status_code} - {response.text}")
