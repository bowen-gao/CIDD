import argparse
import os

import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from tqdm.auto import tqdm
from glob import glob
from collections import Counter
import subprocess



# add the path of the evaluation folder to the sys path

import sys

current_dir = os.path.dirname(os.path.abspath(__file__))

parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from evaluation.docking_qvina import QVinaDockingTask
from evaluation.docking_vina import VinaDockingTask



def vina_dock(smi, protein_root, thread_name, output_path, exhaustiveness=8):
    mol = Chem.MolFromSmiles(smi)

    #print(mol, protein_root)
    vina_task = VinaDockingTask.from_2d_mol(mol,protein_root)
    
    
    docking_results = vina_task.run(mode='dock', exhaustiveness=exhaustiveness)
    pose = docking_results[0]["pose"]


    with open(f"{protein_root}/{thread_name}.pdbqt", "w") as f:
        f.write(pose)

    os.system(f"obabel {protein_root}/{thread_name}.pdbqt -O {output_path}")
    
    
    
    command = ["obabel", f"{protein_root}/{thread_name}.pdbqt", "-O", output_path]
    result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=10)
    print(result.stdout)  # 命令的标准输出
    

    return docking_results[0]["affinity"]

def vina_dock_crossdocked(smi, protein_root, protein_path, ligand_path, thread_name, output_path, exhaustiveness=16):
    mol = Chem.MolFromSmiles(smi)

    #print(mol, protein_root)
    vina_task = VinaDockingTask.from_2d_mol_crossdocked(mol,protein_path, ligand_path)
    
    
    docking_results = vina_task.run(mode='dock', exhaustiveness=8)
    pose = docking_results[0]["pose"]

    print(f"{protein_root}/{thread_name}.pdbqt")
    with open(f"{protein_root}/{thread_name}.pdbqt", "w") as f:
        f.write(pose)

    os.system(f"obabel {protein_root}/{thread_name}.pdbqt -O {output_path}")
    
    
    
    command = ["obabel", f"{protein_root}/{thread_name}.pdbqt", "-O", output_path]
    result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=10)
    print(result.stdout)  # 命令的标准输出
    

    return docking_results[0]["affinity"]




if __name__ == '__main__':
    
    print("Starting Vina Docking Evaluation")



            


