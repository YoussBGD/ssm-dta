import os
from rdkit import Chem
import pandas as pd
import numpy as np
from pubchempy import get_compounds

# Ask the user for the file type
file_type = input("What is the file type (sdf, csv, tsv)? ")

# Ask the user for the file path
file_path = input("Enter the path to the file: ")

# Get the file extension
file_extension = os.path.splitext(file_path)[1].lower()

# Check if the declared file type matches the file extension
if file_type.lower() == 'sdf' and file_extension!= '.sdf':
    print("\nError: The file extension does not match the declared file type (SDF). Please check the file path and try again.\n")
elif file_type.lower() in ['csv', 'tsv'] and file_extension not in ['.csv', '.tsv']:
    print("\nError: The file extension does not match the declared file type (CSV or TSV). Please check the file path and try again.\n")
else:
    # Handling SDF files
    if file_type.lower() == 'sdf':
        output_file = input("Enter the output file name for SMILES: ")
        
        mols = Chem.SDMolSupplier(file_path)
        smiles_list = [Chem.MolToSmiles(mol) for mol in mols if mol is not None]
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(smiles_list))
        
        print(f"\nSMILES have been written to the file {output_file}.\n")

    # Handling CSV or TSV files
    elif file_type.lower() in ['csv', 'tsv']:
        separator = ',' if file_type.lower() == 'csv' else '\t'
        
        # Ask the user about each column separately for clarity
        has_molecule = input("Does the file contain a column for molecules (yes/no)? ")
        has_protein = input("Does the file contain a column for proteins (yes/no)? ")
        has_label = input("Does the file contain a column for labels (yes/no)? ")
        
        # Read the CSV/TSV file
        df = pd.read_csv(file_path, sep=separator)
        
        # Handling the molecule column
        if has_molecule.lower() == 'yes':
            mol_col = input("Enter the column name for molecules: ")
            mol_type = input("Are the molecules in CID (cid) or SMILES (smiles) format? ")
            
            if mol_type.lower() == 'cid':
                # Fetch molecules from PubChem and convert to SMILES
                smiles_list = []
                print("\nProcessing...\n")
                for cid in df[mol_col]:
                    compounds = get_compounds(cid, 'cid')
                    if compounds:
                        mol = compounds[0]
                        smiles = mol.isomeric_smiles
                        smiles_list.append(smiles)
            
            elif mol_type.lower() == 'smiles':
                smiles_list = df[mol_col].tolist()
            
            output_file = input("Enter the output file name for SMILES: ")
            with open(output_file, 'w') as f:
                f.write('\n'.join(smiles_list))
            
            print(f"\nSMILES have been written to the file {output_file}.\n")
        
        # Handling the protein column
        if has_protein.lower() == 'yes':
            protein_col = input("Enter the column name for proteins: ")
            output_file = input("Enter the output file name for proteins: ")
            
            with open(output_file, 'w') as f:
                f.write('\n'.join(df[protein_col].tolist()))
            
            print(f"\nProteins have been written to the file {output_file}.\n")
        
        # Handling the label column
        if has_label.lower() == 'yes':
            label_col = input("Enter the column name for labels: ")
            label_unit = input("What is the type of your labels (Ki, Kd, IC50, score KIBA)? ")
            
            if label_unit.lower() == 'score kiba':
                # If the labels are KIBA scores, write them directly to the file
                output_file = input("Enter the output file name for labels: ")
                df[label_col].to_csv(output_file, index=False, header=False)
                print(f"\nLabels have been written to the file {output_file}.\n")
            else:
                # For Ki, Kd, IC50, check if they are in -log format
                is_log = input(f"Are the labels already in -log10({label_unit}/10^9) format? (yes/no): ")
                
                # Convert to nM if necessary
                if label_unit.lower() in ['ki', 'kd', 'ic50'] and is_log.lower() == "no":
                    current_unit = input("What is the current unit of the labels (nM, µM, mM, M)? ")
                    if current_unit.lower() == 'µm':
                        df[label_col] = df[label_col] * 1000  # Convert µM to nM
                    elif current_unit.lower() == 'mm':
                        df[label_col] = df[label_col] * 1e6  # Convert mM to nM
                    elif current_unit.lower() == 'm':
                        df[label_col] = df[label_col] * 1e9  # Convert M to nM
                
                # Convert to -log(....) if not already in -log format
                if is_log.lower() == 'no':
                    df[label_col] = -np.log10(df[label_col] / 1e9)
                
                output_file = input("Enter the output file name for labels: ")
                df[label_col].to_csv(output_file, index=False, header=False)
                
                print(f"\nLabels have been written to the file {output_file}.\n")

    else:
        print("\nUnrecognized file type. Please enter 'sdf', 'csv', or 'tsv'.\n")
