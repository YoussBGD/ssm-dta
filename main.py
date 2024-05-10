import argparse
import os
import sys
import logging
from tqdm import tqdm
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from lifelines.utils import concordance_index
from fairseq.models.roberta import RobertaModel
from torch.nn.utils.rnn import pad_sequence
import torch
import base64

def main(args):
    num_threads = args.get("num_threads", 4)
    torch.set_num_threads(num_threads)

    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=os.environ.get("LOGLEVEL", "INFO").upper(),
        stream=sys.stdout,
    )

    canon_command = f"python preprocess/canonicalize.py {args['molecules']} --output-fn {args['molecules']}.can --workers 1"
    tokenize_command = f"python preprocess/tokenize_re.py {args['molecules']}.can --output-fn {args['molecules']}.can.re --workers 1"
    add_space_command = f"python preprocess/add_space.py {args['proteins']} --output-fn {args['proteins']}.addspace --workers 1"

    for command in [canon_command, tokenize_command, add_space_command]:
        print("Executing command:", command)
        if os.system(command) != 0:
            print("Command preprocessing failed, check if your input files are in the required format.")
            return

    checkpoint_path = f"./ckpt/{args['model']}/checkpoint_best_20221021.pt"

    roberta = RobertaModel.from_pretrained(
        os.path.dirname(checkpoint_path),
        checkpoint_file=os.path.basename(checkpoint_path),
        data_name_or_path="dict",
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    roberta.to(device)
    roberta.eval()

    total = len(open(args["molecules"] + ".can.re", 'r').readlines())
    pbar = tqdm(total=total, desc='Predicting')

    with open(f'{args["molecules"]}.can.re', 'r') as mol_in, \
         open(f'{args["proteins"]}.addspace', 'r') as pro_in, \
         open(args["output_fn"], 'w') as out_f:
        batch_mol_buf = []
        batch_pro_buf = []
        for i, mol_pro in enumerate(zip(mol_in, pro_in)):
            mol, pro = mol_pro
            tmp_mol, tmp_pro = roberta.myencode_separate(mol.strip(), pro.strip())
            batch_mol_buf.append(tmp_mol[:512])
            batch_pro_buf.append(tmp_pro[:1024])
            if ((i+1) % args["batch_size"] == 0) or ((i+1) == total):
                tokens_0 = pad_sequence(batch_mol_buf, batch_first=True, padding_value=1).to(device)
                tokens_1 = pad_sequence(batch_pro_buf, batch_first=True, padding_value=1).to(device)
                predictions = roberta.myextract_features_separate(tokens_0, tokens_1)
                for result in predictions:
                    out_f.write(f'{str(result.item())}\n')
                batch_mol_buf.clear()
                batch_pro_buf.clear()
            pbar.update(1)
        pbar.close()
        
    results = {}
    if args["mode"] == 'evaluation':
        if args["labels"] is not None:
            try:
                pred_scores = [float(line.strip()) for line in open(args["output_fn"], 'r').readlines()]
                true_scores = [float(line.strip()) for line in open(args["labels"], 'r').readlines()]
                
                if len(true_scores) != len(pred_scores):
                    results = {
                        "error": "Mismatch in the number of true scores and predicted scores."
                    }
                else:
                    mse = mean_squared_error(true_scores, pred_scores)
                    rmse = np.sqrt(mse)
                    pearson = pearsonr(true_scores, pred_scores)[0]
                    c_index = concordance_index(true_scores, pred_scores)
                    
                    results = {
                        "mse": mse,
                        "rmse": rmse,
                        "pearson": pearson,
                        "c_index": c_index,
                        "true_scores": true_scores,
                        "pred_scores": pred_scores
                    }

            except FileNotFoundError:
                results = {
                    "error": "Labels file not found."
                }
            except Exception as e:
                results = {
                    "error": f"An error occurred while processing the labels file: {str(e)}"
                }
        else:
            results = {
                "error": "Labels file not provided for evaluation."
            }
    else:
        pred_scores = [float(line.strip()) for line in open(args["output_fn"], 'r').readlines()]
        results = {
            "pred_scores": pred_scores
        }
    
    with open(args["output_fn"], "rb") as f:
        output_file = base64.b64encode(f.read()).decode("utf-8")
    
    results["output_file"] = output_file

    return results

