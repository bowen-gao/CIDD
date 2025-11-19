from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
import json
import os
from tqdm import tqdm
import multiprocessing as mp
from collections import defaultdict


def get_csk(mol):
    """Generate the Canonical Simplified Scaffold (CSK) for a molecule."""
    Chem.RemoveStereochemistry(mol)  # Important for canonization of CSK!
    scaff = MurckoScaffold.GetScaffoldForMol(mol)
    scaff = MurckoScaffold.MakeScaffoldGeneric(scaff)
    return scaff


def draw_top_scaffolds(scaffold_smiles_list, scaffold_counts, top_n=100, output_file="top_scaffolds_100.png"):
    """Draw the top N scaffolds based on their frequency."""
    scaff_mols = [Chem.MolFromSmiles(smi) for smi in scaffold_smiles_list[:top_n]]
    scaff_legends = [f"Scaffold {i+1} ({scaffold_counts[i]} times)" for i in range(top_n)]
    img = Draw.MolsToGridImage(scaff_mols, legends=scaff_legends, molsPerRow=100)
    img.save(output_file)
    print(f"Top {top_n} scaffolds saved to {output_file}")


def get_smi(ligand_file_path):
    """Extract the scaffold SMILES from a ligand file."""
    try:
        suppl = Chem.SDMolSupplier(ligand_file_path, removeHs=True)
        ligand = [mol for mol in suppl if mol is not None][0]
        scaff = get_csk(ligand)
        return Chem.MolToSmiles(scaff)
    except Exception as e:
        print(f"Error processing file {ligand_file_path}: {e}")
        return None


def process_jobs(jobfile, num_processes=128):
    """Process ligand files to extract scaffolds using multiprocessing."""
    with mp.Manager() as manager:
        scaff_list = manager.list()
        tbar = tqdm(total=len(jobfile))
        finished = open('finished.jobs', 'w')

        def call_back(result):
            if result:
                scaff_list.append(result)
                finished.write(result + '\n')
                finished.flush()
            tbar.update(1)

        with mp.Pool(num_processes) as pool:
            for job in jobfile:
                pool.apply_async(func=get_smi, args=(job.strip(),), callback=call_back)
            pool.close()
            pool.join()

        finished.close()
        tbar.close()

        return list(scaff_list)


if __name__ == '__main__':
    # split_file_path = '/mnt/nfs-ssd/data/huangyanwen/InterNet/extract_inter_for_crossdock/crossdock_split_by_name.json'
    # dir = '/nfs/data/crossdock/crossdocked_pocket10'
    #
    # jobfile = []
    #
    # # Load job file paths
    # try:
    #     with open(split_file_path, 'r') as split_file:
    #         data = json.load(split_file)
    #         for k, v in data.items():
    #             print(f"Processing key: {k}, {len(v)} items")
    #             for t in v:
    #                 _, ligand_dir = t
    #                 ligand_file_path = os.path.join(dir, ligand_dir)
    #                 jobfile.append(ligand_file_path)
    # except Exception as e:
    #     print(f"Error loading split file: {e}")
    #     exit(1)
    #
    # # Process jobs and extract scaffolds
    # scaff_list = process_jobs(jobfile, num_processes=128)
    #
    #
    # # Count occurrences of each scaffold
    # scaff_dict = defaultdict(int)
    # for s in scaff_list:
    #     scaff_dict[s] += 1
    #
    # # Sort scaffolds by frequency
    # sorted_scaffolds = sorted(scaff_dict.items(), key=lambda x: x[1], reverse=True)
    #
    # with open('crossdock_scaffold_counts.json', 'w') as jf:
    #     json.dump(sorted_scaffolds, jf, indent=6)

    with open('crossdock_scaffold_counts.json', 'r') as jf:
        sorted_scaffolds = json.load(jf)

    sorted_scaffold_smiles = [item[0] for item in sorted_scaffolds]
    sorted_scaffold_counts = [item[1] for item in sorted_scaffolds]

    # Draw the top 10 scaffolds
    draw_top_scaffolds(sorted_scaffold_smiles, sorted_scaffold_counts, top_n=100)


