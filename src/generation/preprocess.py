import os
from tqdm import tqdm
from pymol import cmd
import copy


def load_sdf_and_add_polar_hs(input_ligand, input_protein, output_combined):
    # Load the ligand .sdf file
    cmd.load(input_ligand, "ligand")

    # Load the protein .pdb file
    cmd.load(input_protein, "protein")

    # Add hydrogens to both ligand and protein (PyMOL adds both polar and non-polar hydrogens)
    cmd.h_add("ligand")
    cmd.h_add("protein")

    # # Select and remove non-polar hydrogens from the ligand
    # cmd.select("non_polar_hydrogens_ligand", "ligand and hydro and (neighbor elem C)")
    # cmd.remove("non_polar_hydrogens_ligand")
    #
    # # Select and remove non-polar hydrogens from the protein
    # cmd.select("non_polar_hydrogens_protein", "protein and hydro and (neighbor elem C)")
    # cmd.remove("non_polar_hydrogens_protein")

    # Combine ligand and protein into one object
    cmd.create("combined", "ligand or protein")

    # Save the combined object as a .pdb file
    cmd.save(output_combined, "combined")

    # Clean up the selections and objects
    cmd.delete("non_polar_hydrogens_ligand")
    cmd.delete("ligand")
    cmd.delete("non_polar_hydrogens_protein")
    cmd.delete("protein")
    cmd.delete("combined")


def process_one_dir(indir):
    print(indir)

    ligand = None
    protein = None
    for f in os.listdir(indir):
        if '_ligand.sdf' in f:
            ligand = copy.deepcopy(os.path.join(indir, f))
            print(f)
        if '_protein.pdb' in f:
            protein = copy.deepcopy(os.path.join(indir, f))
            print(f)

    if ligand and protein:
        combined_out_filepath = os.path.join(indir, "combined_structures.pdb")
        load_sdf_and_add_polar_hs(ligand, protein, combined_out_filepath)

    return indir