#!/usr/bin/env python3
"""
Task script for running CIDD molecular generation on CrossDocked data.

This script processes individual protein-ligand pairs and generates
optimized molecules using the CIDD framework.

Usage:
    python task_script_smiles_crossdocked.py <pdb_path> <ligand_file> <smiles> <index> <timestamp>
"""

import sys
import os
import pickle
import json
import numpy as np
from typing import Optional, Tuple

# Add the src directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
src_dir = os.path.join(parent_dir, 'src')
sys.path.append(src_dir)

try:
    from generation.cidd_generation import Generation
except ImportError as e:
    print(f"Error: Could not import CIDD Generation module: {e}")
    print("Please ensure all dependencies are installed and CIDD is properly set up.")
    sys.exit(1)


def main():
    """Main execution function."""
    if len(sys.argv) != 5:
        print("Usage: python task_script_smiles_crossdocked.py <ligand_file> <smiles> <index> <timestamp>")
        sys.exit(1)
    
    ligand_filename, ligand_smi, idx, save_name = sys.argv[1:]
    
    # Set up output directory
    output_dir = os.environ.get('CIDD_OUTPUT_DIR', './results')
    root_dir = os.path.dirname(ligand_filename)
    save_dir = os.path.join(output_dir, root_dir, save_name)
    os.makedirs(save_dir, exist_ok=True)
    
    # Initialize CIDD Generation
    data_root = os.environ.get('CIDD_DATA_ROOT', './data')
    generation = Generation("Generation Module", save_dir, data_root)
    
    print(f"Processing: {ligand_filename}")
    print(f"Starting SMILES: {ligand_smi}")
    print(f"Output directory: {save_dir}")
    
    # Run CIDD generation
    try:
        res = generation.new_execute_all(ligand_filename, ligand_smi, idx)
        if res is None:
            print("Generation failed - no results returned")
            sys.exit(1)
        
        # Save results
        output_file = os.path.join(save_dir, f"result_{idx}.pkl")
        
        try:
            # Try to unpack full result set (with molecule objects)
            chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, chosen_failed_count, ori_mol, chosen_gen_mol = res
            result_data = (chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, chosen_failed_count, ori_mol, chosen_gen_mol)
        except ValueError:
            # Fallback to basic result set (without molecule objects)
            chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, chosen_failed_count = res
            result_data = (chosen_smi, chosen_report, chosen_score, ori_smi, ori_score, chosen_failed_count)
        
        with open(output_file, "wb") as f:
            pickle.dump(result_data, f)
        
        print(f"Results saved to: {output_file}")
        print(f"Generated SMILES: {chosen_smi}")
        print(f"Score improvement: {ori_score} -> {chosen_score}")
        
    except Exception as e:
        print(f"Generation failed with error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

        
        
        



            
           
        
        

        
        



    

    #

