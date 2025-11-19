import logging
from rdkit import Chem
from collections import defaultdict
import os
from rdkit.Chem import Descriptors


# Configure logging
logging.basicConfig(filename="unreasonable_molecules.log", level=logging.INFO, format="%(message)s")


def get_co_id(mol):
    """
    Detect C=O or C=N bonds in the molecule and return the atom indices of the carbon atoms
    involved in these bonds.

    Parameters:
        mol (rdkit.Chem.Mol): RDKit molecule object.

    Returns:
        List[int]: A list of atom indices for carbon atoms in C=O bonds.
    """
    co_ids = set()
    for bond in mol.GetBonds():
        # Check if the bond is a double bond
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()

            # Check for C=O configuration
            if atom1.GetSymbol() == "C" and atom2.GetSymbol() == "N" or "O":
                co_ids.add(atom1.GetIdx())
            elif atom2.GetSymbol() == "C" and atom1.GetSymbol() == "N" or "O":
                co_ids.add(atom2.GetIdx())
    return co_ids


def get_ring_groups(mol):
    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()  # List of tuples of atom indices

    # Initialize dictionary to store groups
    groups_dict = defaultdict(tuple)

    # Create a list to keep track of assigned groups
    assigned_groups = [-1] * len(ring_atoms)

    # Function to merge rings with common atoms into a single group
    def merge_groups(ring_idx, group_idx):
        # Start with the current ring
        stack = [ring_idx]
        groups_dict[group_idx] = (ring_atoms[ring_idx],)  # Initialize with current ring
        assigned_groups[ring_idx] = group_idx

        while stack:
            current_ring_idx = stack.pop()

            for idx, other_ring in enumerate(ring_atoms):
                # Skip already assigned rings
                if assigned_groups[idx] != -1:
                    continue

                # Convert both rings to sets to perform set operations
                current_ring_set = set(ring_atoms[current_ring_idx])
                other_ring_set = set(other_ring)

                # Check for common atoms between the current ring and other rings
                if not current_ring_set.isdisjoint(other_ring_set):
                    # Add the new ring as a tuple to the group's list of tuples
                    assigned_groups[idx] = group_idx
                    groups_dict[group_idx] += (other_ring,)
                    stack.append(idx)

    # Iterate over all rings
    group_counter = 0
    for ring_idx in range(len(ring_atoms)):
        # If the ring hasn't been assigned a group yet, assign it to a new group
        if assigned_groups[ring_idx] == -1:
            merge_groups(ring_idx, group_counter)
            group_counter += 1

    return groups_dict


def check_groups(groups_dict, mol, smiles, co_ids):
    reasonable = True
    unreasonable_atom_ratio = 0
    ring_atom_idex = set()
    mw = Descriptors.MolWt(mol)

    for ring_info in groups_dict.values():
        for t in ring_info:
            ring_atom_idex.update(set(t))

    ring_atom_num = len(ring_atom_idex)
    ring_atom_ratio = len(ring_atom_idex)/mol.GetNumAtoms()

    for group_idx, group in groups_dict.items():
        verified_sp2 = set()
        verified_non_sp2 = set()
        remaining_rings = []

        for ring in group:
            sp2_atoms = []
            non_sp2_atoms = []

            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == "C" and atom_idx not in co_ids:  # Check only carbon atoms
                    if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                        sp2_atoms.append(atom_idx)
                    else:
                        non_sp2_atoms.append(atom_idx)

            # Check if the entire ring is sp2 or not
            if len(sp2_atoms) == len(ring):
                verified_sp2.update(ring)
            elif len(non_sp2_atoms) == len(ring):
                verified_non_sp2.update(ring)
            else:
                remaining_rings.append(ring)

        # Iteratively handle remaining mixed rings
        while remaining_rings:
            new_sp2 = set()
            new_non_sp2 = set()

            for ring in remaining_rings[:]:
                remaining_atoms = set(ring) - verified_sp2 - verified_non_sp2

                # Recheck remaining atoms for consistency
                is_consistently_sp2 = all(
                    mol.GetAtomWithIdx(idx).GetHybridization() == Chem.rdchem.HybridizationType.SP2
                    for idx in remaining_atoms
                    if mol.GetAtomWithIdx(idx).GetSymbol() == "C" and idx not in co_ids  # Check only carbon atoms
                )
                is_consistently_non_sp2 = all(
                    mol.GetAtomWithIdx(idx).GetHybridization() != Chem.rdchem.HybridizationType.SP2
                    for idx in remaining_atoms
                    if mol.GetAtomWithIdx(idx).GetSymbol() == "C" and idx not in co_ids  # Check only carbon atoms
                )

                if is_consistently_sp2:
                    new_sp2.update(remaining_atoms)
                    remaining_rings.remove(ring)
                elif is_consistently_non_sp2:
                    new_non_sp2.update(remaining_atoms)
                    remaining_rings.remove(ring)

            # Update verified sets with newly determined atoms
            verified_sp2.update(new_sp2)
            verified_non_sp2.update(new_non_sp2)

            # Stop iteration if no new sp2 or non-sp2 atoms are identified
            if not new_sp2 and not new_non_sp2:
                break

        # If remaining_rings is not empty after the iteration, mark molecule as unreasonable
        if remaining_rings:
            unreasonable_atom_num = 0
            reasonable = False
            logging.info(smiles)
            for ring in remaining_rings:
                unreasonable_atom_num += len(ring)
            # unreasonable_atom_ratio = unreasonable_atom_num/mol.GetNumAtoms()
            unreasonable_atom_ratio = unreasonable_atom_num/ring_atom_num


            break

    return reasonable, unreasonable_atom_ratio, ring_atom_ratio, mw


def get_fda_smis():
    fda_list = []
    fda_path = '/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/get_fragment/FDA_Approved_structures.csv'
    count = 0
    with open(fda_path, 'r') as f:
        lines = f.readlines()[1:]

        for line in lines:
            mol = None
            smi = line.strip().split(',')[-1].strip()
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                count += 1
                fda_list.append(smi)

    print(count)
    return fda_list


def get_data(path):
    data_list = []
    count = 0
    for f in os.listdir(path):
        if f.endswith('.sdf'):
            ligand_file_path = os.path.join(path, f)
            # Load molecule from file
            suppl = Chem.SDMolSupplier(ligand_file_path, removeHs=True)
            ligands = [mol for mol in suppl if mol is not None]

            # If ligand is None, return None
            for ligand in ligands:
                if ligand is None:
                    count += 1
                else:
                    data_list.append(Chem.MolToSmiles(ligand))

    print(count)
    return data_list


def main(process='fda', model='bfn'):
    print('Processing', process)

    if process=='fda':
        smiles_list = get_fda_smis()
    elif process=='ori':
        # ori_path = '/mnt/nfs-ssd/data/gaobowen/pa_mols/2024-12-26-22-27-57_ori_mols'
        # ori_path = '/mnt/nfs-ssd/data/gaobowen/pa_mols/2024-12-27-18-59-09_ori_mols'
        ori_path = '/mnt/nfs-ssd/data/gaobowen/pa_mols/2024-12-28-16-36-19_ori_mols'
        smiles_list = get_data(ori_path)
    elif process=='gen':
        # gen_path = '/mnt/nfs-ssd/data/gaobowen/pa_mols/2024-12-26-22-27-57_gen_mols'
        # gen_path = '/mnt/nfs-ssd/data/gaobowen/pa_mols/2024-12-27-18-59-09_gen_mols'
        gen_path = '/mnt/nfs-ssd/data/gaobowen/pa_mols/2024-12-28-16-36-19_gen_mols'
        smiles_list = get_data(gen_path)
    else:
        print('Using example molecules with reasonable and unreasonable structures!')
        smiles_list = [
            "C1CCN(CC1)C2(CCN(CC2)CCCN3C4=CC=CC=C4CCC5=C3C=C(C=C5)Cl)C(=O)N",
            "CC1=C(C(C(=C(N1)C)C(=O)OC)C2=CC=CC=C2[N+](=O)[O-])C(=O)OC",
            'C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O'
        ]

    reasonable_count = 0
    unreasonable_count = 0
    group_num_list = []
    max_ring_num_list = []
    unreasonable_atom_ratio = []
    ring_atom_ratio_list = []
    mw_list = []

    for smiles in smiles_list:
        molecule = Chem.MolFromSmiles(smiles)
        groups_dict = get_ring_groups(molecule)

        # print(groups_dict)
        # print(len(groups_dict))

        group_num_list.append(len(groups_dict))

        ring_num_list = [len(ring_info) for ring_info in groups_dict.values()]
        try:
            max_ring_num_list.append(max(ring_num_list))
        except ValueError:
            max_ring_num_list.append(0)

        co_ids = get_co_id(molecule)

        reasonable, unreasonable_ratio, ring_atom_ratio, mw = check_groups(groups_dict, molecule, smiles, co_ids)
        ring_atom_ratio_list.append(ring_atom_ratio)
        unreasonable_atom_ratio.append(unreasonable_ratio)
        mw_list.append(mw)
        if reasonable:
            reasonable_count += 1
        else:
            unreasonable_count += 1

    total = (unreasonable_count+reasonable_count)
    print('For', process)
    print(f"    Reasonable molecules: {reasonable_count}")
    print(f"    Unreasonable molecules: {unreasonable_count}")
    print(f"    Total molecule count: {total}")
    print(f"    Reasonable rate: {reasonable_count/total}")
    # print(f"    Average group number: {sum(group_num_list)/total}")
    print(f"    Averaged max ring number in groups: {sum(max_ring_num_list)/total}")
    print(f"    Averaged unreasonable atom ratio: {sum(unreasonable_atom_ratio)/total}")
    # print(f"    Averaged ring atom ratio: {sum(ring_atom_ratio_list)/total}")
    print(f"    Averaged wm: {sum(mw_list)/total}")
    print("=" * 60)



if __name__ == '__main__':
    main('fda')
    main('ori')
    main('gen')

