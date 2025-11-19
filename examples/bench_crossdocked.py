#!/usr/bin/env python3
"""
CIDD Benchmark Script for CrossDocked Dataset

This script runs CIDD molecular generation on a batch of molecules
from the CrossDocked2020 dataset, demonstrating the full pipeline
on realistic protein-ligand pairs.

Usage:
    python bench_crossdocked.py --data_path <path_to_pt_file> --output_dir <output_directory>

Requirements:
    - CrossDocked2020 dataset
    - PyTorch for loading data files
    - All CIDD dependencies
"""

import sys
import os
import argparse
import pickle
import json
import time
import random
import subprocess
import multiprocessing as mp
from typing import List, Dict, Any, Optional

# Add the src directory to Python path  
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
src_dir = os.path.join(parent_dir, 'src')
sys.path.append(src_dir)

try:
    import torch
    from rdkit import Chem
    from tqdm import tqdm
    import psutil
except ImportError as e:
    print(f"Error: Missing required dependencies: {e}")
    print("Please install: torch, rdkit, tqdm, psutil")
    sys.exit(1)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run CIDD benchmark on CrossDocked dataset')
    parser.add_argument('--data_path', type=str, required=True,
                       help='Path to CrossDocked data file (.pt format)')
    parser.add_argument('--output_dir', type=str, default='./benchmark_results',
                       help='Output directory for results')
    parser.add_argument('--max_samples', type=int, default=100,
                       help='Maximum number of samples to process')
    parser.add_argument('--samples_per_target', type=int, default=10,
                       help='Number of samples per protein target')
    parser.add_argument('--num_workers', type=int, default=10,
                       help='Number of parallel workers')
    parser.add_argument('--timeout', type=int, default=10800,
                       help='Timeout per task in seconds')
    parser.add_argument('--seed', type=int, default=2025,
                       help='Random seed for reproducibility')
    
    return parser.parse_args()


def main():
    """Main benchmark execution."""
    args = parse_args()
    
    print(f"🔬 CIDD Benchmark on CrossDocked Dataset")
    print(f"Data: {args.data_path}")
    print(f"Output: {args.output_dir}")
    print(f"Max samples: {args.max_samples}")
    
    # Load molecular data
    try:
        mol_data_list = torch.load(args.data_path)
        print(f"✓ Loaded {len(mol_data_list)} molecules")
    except Exception as e:
        print(f"✗ Failed to load data: {e}")
        sys.exit(1)

    # new_mol_data_list = []
    # for mol in mol_data_list:
    #     new_mol_data_list.extend(mol)
    # mol_data_list = new_mol_data_list
    # mol_data_list = mol_data_list["all_results"]
    
    
    
    
    import subprocess

    def process_one_pair_new(ligand_filename, ligand_smi, idx, cur_time):
        """Process one protein-ligand pair with subprocess timeout control."""
        script_path = os.path.join(current_dir, "task_script_smiles_crossdocked.py")
        command = ["python", script_path, ligand_filename, ligand_smi, str(idx), cur_time]

        print(" ".join(command))

        try:
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            stdout, stderr = process.communicate(timeout=10800)
            return stderr
        except subprocess.TimeoutExpired:
            import psutil
            parent = psutil.Process(process.pid)
            print(f"Killing parent process {parent.pid}")
            for child in parent.children(recursive=True):  # 递归杀死子进程
                print(f"Killing child process {child.pid}")
                child.kill()
            parent.kill()
            print(f"Task {ligand_filename} {idx} timed out.")
            return "Timeout"
        except Exception as e:
            print(f"Task {ligand_filename} {idx} failed: {e}")
            return "Failed"

        
        



    import multiprocessing as mp
    from tqdm import tqdm
    import random


   

    

    args = []

    
   
    
    results = []
    failed_tasks = []

    import time

    start_time = time.time()

    cur_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    
    
    data_dic = {}

    for mol in mol_data_list:
        ligand_filename = mol["ligand_filename"]
        if ligand_filename in data_dic:
            data_dic[ligand_filename].append(mol)
        else:
            data_dic[ligand_filename] = [mol]
    
    random.seed(2025)

    data_dic_keys = list(data_dic.keys())

    # shuffle

    random.shuffle(data_dic_keys)

    data_dic_keys = data_dic_keys[:10]

    finished_dic = {}

    with mp.Pool(100) as pool:
        for ligand_filename in tqdm(data_dic_keys):
            # random select 20 mols
            # if ligand_filename.split("/")[0] == "4x0p":
            #     continue
            mols = data_dic[ligand_filename]
            #print(len(mols))
            mols = random.sample(mols, min(100, len(mols)))
            
            for idx, mol in enumerate(mols):
                mol_object = mol["mol"]
                
                try:
                    ligand_smi = Chem.MolToSmiles(mol_object)
                except:
                    continue
                ligand_smi = mol["smiles"]
                if "." in ligand_smi:
                    continue
                
                tmp_dir = os.path.join(parent_dir, "tmp")
                os.makedirs(tmp_dir, exist_ok=True)
                
                results.append(pool.apply_async(process_one_pair_new, args=(ligand_filename, ligand_smi, idx, cur_time)))
                time.sleep(10)
                # target_name = ligand_filename.split("/")[0]
                # #mol["smiles"] = ligand_smi
                # if target_name in finished_dic:
                #     finished_dic[target_name].append(mol)
                # else:
                #     finished_dic[target_name] = [mol]

        pool.close()
        pool.join()

    # save finished_dic

    # with open("decomp_iter1_finished_dic.pkl", "wb") as f:
    #     pickle.dump(finished_dic, f)

    for i, res in enumerate(results):
        print(f"Task {i}: {res.get()}")
    
    end_time = time.time()

    duration = end_time - start_time


    print(f"Duration: {duration}")

    # 
        
    

    

