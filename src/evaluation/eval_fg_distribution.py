from rdkit import Chem
from rdkit.Chem import Fragments
import multiprocessing as mp
from typing import List
import os
import json
from tqdm import tqdm


def identify_rdkit_fragments(ligand_file_path: str, error_counts: dict) -> List[str] or None:
    """Identify fragments in the ligand from the given file path."""
    try:
        # Load molecule from file
        suppl = Chem.SDMolSupplier(ligand_file_path, removeHs=True)
        ligand = [mol for mol in suppl if mol is not None][0]

        # If ligand is None, return None
        if ligand is None:
            raise ValueError(f"Failed to process ligand: {ligand_file_path}")

        fragment_names = [name for name in dir(Fragments) if name.startswith('fr_')]

        identified_groups = []

        for fragment_name in fragment_names:
            fragment_func = getattr(Fragments, fragment_name)
            fragment_count = fragment_func(ligand)

            if fragment_count > 0:
                identified_groups.append(fragment_name)

        return identified_groups

    except Exception as e:
        # Increment error count for this type of exception
        error_type = str(type(e)).split("'")[1]  # Extract the exception type
        error_counts[error_type] = error_counts.get(error_type, 0) + 1
        print(f"Error processing {ligand_file_path}: {e}")
        return None


def process_jobs(jobfile: List[str], num_processes=128) -> List[List[str]]:
    """Process ligand files to extract fragments using multiprocessing."""
    with mp.Manager() as manager:
        fg_list = manager.list()
        error_counts = manager.dict()  # Shared dictionary for error counts
        tbar = tqdm(total=len(jobfile))
        finished = open('finished.jobs', 'w')

        def call_back(result):
            if result:
                fg_list.append(result)
                finished.write(" ".join(result) + '\n')  # Write the result in the file
                finished.flush()
            tbar.update(1)

        with mp.Pool(num_processes) as pool:
            for job in jobfile:
                pool.apply_async(func=identify_rdkit_fragments, args=(job, error_counts), callback=call_back)
            pool.close()
            pool.join()

        finished.close()
        tbar.close()

        return list(fg_list), dict(error_counts)


def count_fragment_occurrences(fg_list: List[List[str]]) -> dict:
    """Count the occurrences of each functional group."""
    fragment_counts = {}

    for fg in fg_list:
        if fg:
            for fragment in fg:
                fragment_counts[fragment] = fragment_counts.get(fragment, 0) + 1

    return fragment_counts


def rank_fragments(fragment_counts: dict) -> dict:
    """Sort the fragment counts from most to least frequent."""
    # Sort the dictionary by count in descending order
    sorted_fragments = dict(sorted(fragment_counts.items(), key=lambda item: item[1], reverse=True))
    return sorted_fragments


if __name__ == '__main__':
    split_file_path = '/mnt/nfs-ssd/data/huangyanwen/InterNet/extract_inter_for_crossdock/crossdock_split_by_name.json'
    dir_path = '/nfs/data/crossdock/crossdocked_pocket10'

    jobfile = []

    # Load job file paths
    try:
        with open(split_file_path, 'r') as split_file:
            data = json.load(split_file)
            for k, v in data.items():
                print(f"Processing key: {k}, {len(v)} items")
                for t in v:
                    _, ligand_dir = t
                    ligand_file_path = os.path.join(dir_path, ligand_dir)
                    jobfile.append(ligand_file_path)
    except Exception as e:
        print(f"Error loading split file: {e}")
        exit(1)

    # Process jobs and extract fragments
    fg_list, error_counts = process_jobs(jobfile, num_processes=128)

    # Count occurrences of each fragment
    fragment_counts = count_fragment_occurrences(fg_list)

    # Rank the fragments from most to least frequent
    ranked_fragments = rank_fragments(fragment_counts)

    # Print the ranked fragments
    print("Ranked Fragment Occurrences:")
    for fragment, count in ranked_fragments.items():
        print(f"{fragment}: {count}")

    # Print the error counts
    print("\nError Counts:")
    for error_type, count in error_counts.items():
        print(f"{error_type}: {count}")

    # Optionally, save the ranked counts to a file
    with open("ranked_fragment_counts.json", "w") as count_file:
        json.dump(ranked_fragments, count_file, indent=4)

    # Optionally, save the error counts to a file
    with open("error_counts.json", "w") as error_file:
        json.dump(error_counts, error_file, indent=4)

