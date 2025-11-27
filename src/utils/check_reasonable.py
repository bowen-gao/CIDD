from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict
import logging

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
    ring_atom_index = set()
    mw = Descriptors.MolWt(mol)

    # Collect ring atom indices
    for ring_info in groups_dict.values():
        for ring in ring_info:
            ring_atom_index.update(ring)

    ring_atom_num = len(ring_atom_index)
    ring_atom_ratio = ring_atom_num / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0

    # Identify shared atoms
    atom_to_rings = defaultdict(set)
    for group in groups_dict.values():
        for ring in group:
            for atom_idx in ring:
                atom_to_rings[atom_idx].add(tuple(ring))

    shared_atoms = [atom for atom, rings in atom_to_rings.items() if len(rings) > 1]
    shared_atom_count = len(shared_atoms)

    # Early SP2 check
    if any(len(group) > 1 for group in groups_dict.values()):  # 判断是否需要进行检查
        no_sp2 = True
        unreasonable_atom_count = 0

        for group in groups_dict.values():
            if len(group) > 1:  # 仅检查包含多个 ring 的 group
                for ring in group:
                    for atom_idx in ring:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == "C" and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                            no_sp2 = False  # 找到 sp2 碳，立即标记为 False
                            break
                    if not no_sp2:
                        break

                # 如果未提前退出，计算不合理的原子数量
                unreasonable_atom_count += sum(len(ring) for ring in group)

        # 如果所有 group 中都没有 sp2 碳，标记为 unreasonable
        if no_sp2 and unreasonable_atom_count > 0:
            reasonable = False
            unreasonable_atom_ratio = unreasonable_atom_count / ring_atom_num if ring_atom_num > 0 else 0
            return reasonable, unreasonable_atom_ratio, ring_atom_ratio, mw, shared_atom_count

    # SP2/non-SP2 classification
    for group in groups_dict.values():
        verified_sp2 = set()
        verified_non_sp2 = set()
        remaining_rings = []

        for ring in group:
            sp2_atoms = []
            non_sp2_atoms = []

            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == "C" and atom_idx not in co_ids:
                    if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                        sp2_atoms.append(atom_idx)
                    else:
                        non_sp2_atoms.append(atom_idx)

            if len(sp2_atoms) == len(ring):
                verified_sp2.update(ring)
            elif len(non_sp2_atoms) == len(ring):
                verified_non_sp2.update(ring)
            else:
                remaining_rings.append(ring)

        # Resolve mixed rings iteratively
        while remaining_rings:
            new_sp2 = set()
            new_non_sp2 = set()

            for ring in remaining_rings[:]:
                remaining_atoms = set(ring) - verified_sp2 - verified_non_sp2
                is_consistently_sp2 = all(
                    mol.GetAtomWithIdx(idx).GetHybridization() == Chem.rdchem.HybridizationType.SP2
                    for idx in remaining_atoms if mol.GetAtomWithIdx(idx).GetSymbol() == "C" and idx not in co_ids
                )
                is_consistently_non_sp2 = all(
                    mol.GetAtomWithIdx(idx).GetHybridization() != Chem.rdchem.HybridizationType.SP2
                    for idx in remaining_atoms if mol.GetAtomWithIdx(idx).GetSymbol() == "C" and idx not in co_ids
                )

                if is_consistently_sp2:
                    new_sp2.update(remaining_atoms)
                    remaining_rings.remove(ring)
                elif is_consistently_non_sp2:
                    new_non_sp2.update(remaining_atoms)
                    remaining_rings.remove(ring)

            verified_sp2.update(new_sp2)
            verified_non_sp2.update(new_non_sp2)

            if not new_sp2 and not new_non_sp2:
                break

        # Handle unresolved rings
        if remaining_rings:
            unreasonable_atom_count = sum(len(ring) for ring in remaining_rings)
            unreasonable_atom_ratio = unreasonable_atom_count / ring_atom_num if ring_atom_num > 0 else 0
            logging.info(f"Unresolved rings in molecule: {smiles}")
            reasonable = False
            break

    return reasonable, unreasonable_atom_ratio, ring_atom_ratio, mw, shared_atom_count




def check_aro_cp(molecule, smiles):
    if molecule is not None:
        groups_dict = get_ring_groups(molecule)
        co_ids = get_co_id(molecule)
        reasonable, unreasonable_ratio, ring_atom_ratio, mw, shared_atom_count = check_groups(groups_dict, molecule,
                                                                                              smiles, co_ids)
        return reasonable, unreasonable_ratio
    else:
        return None


if __name__ == "__main__":

    smi = "CN1CCC(c2cccc(C(=O)c3nc4c5c(ncc4s3)NC3=NCc4c(N)c(O)cc6c4C3=C5CC6)n2)CC1=O"

    print(check_aro_cp(Chem.MolFromSmiles(smi), smi))