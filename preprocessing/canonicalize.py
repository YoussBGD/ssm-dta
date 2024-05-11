import re
from rdkit import Chem

def rm_map_number(smiles):
    return re.sub(':\d*', '', smiles)

def canonicalize(smiles, keep_atommap=False):
    if not keep_atommap:
        smiles = rm_map_number(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    else:
        return Chem.MolToSmiles(mol)

