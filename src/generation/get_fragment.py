import pickle
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Contrib.IFG.ifg import identify_functional_groups
#from extract_feat_labels import process_fg, rdkit_functional_group_label_features_generator
#from frag_graph_smiles import fragment_graph_cutbonds
from tqdm import tqdm
from collections import defaultdict



def brics_decomp(mol, addition_rule=False, return_all_bonds=True):
    """
    return break bonds, use additional rule or not
    """
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 1:
        return [[0]], []

    cliques = []
    breaks = []
    all_bonds = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        cliques.append([a1, a2])
        all_bonds.append([a1, a2])

    res = list(BRICS.FindBRICSBonds(mol))

    cut_bonds_set = [list(ele[0]) for ele in res]

    if not addition_rule:
        return cut_bonds_set

    if len(res) == 0:
        return [list(range(n_atoms))], []
    else:
        for bond in res:
            if [bond[0][0], bond[0][1]] in cliques:
                cliques.remove([bond[0][0], bond[0][1]])
            else:
                cliques.remove([bond[0][1], bond[0][0]])
            cliques.append([bond[0][0]])
            cliques.append([bond[0][1]])

    # break bonds between rings and non-ring atoms
    for c in cliques:
        if len(c) > 1:
            if mol.GetAtomWithIdx(c[0]).IsInRing() and not mol.GetAtomWithIdx(c[1]).IsInRing():
                cliques.remove(c)
                cliques.append([c[1]])
                breaks.append(c)
            if mol.GetAtomWithIdx(c[1]).IsInRing() and not mol.GetAtomWithIdx(c[0]).IsInRing():
                cliques.remove(c)
                cliques.append([c[0]])
                breaks.append(c)

    # select atoms at intersections as motif
    for atom in mol.GetAtoms():
        if len(atom.GetNeighbors()) > 2 and not atom.IsInRing():
            cliques.append([atom.GetIdx()])
            for nei in atom.GetNeighbors():
                if [nei.GetIdx(), atom.GetIdx()] in cliques:
                    cliques.remove([nei.GetIdx(), atom.GetIdx()])
                    breaks.append([nei.GetIdx(), atom.GetIdx()])
                elif [atom.GetIdx(), nei.GetIdx()] in cliques:
                    cliques.remove([atom.GetIdx(), nei.GetIdx()])
                    breaks.append([atom.GetIdx(), nei.GetIdx()])
                cliques.append([nei.GetIdx()])

    # merge breaks
    cut_bonds_set.extend(breaks)
    if return_all_bonds:
        return cut_bonds_set, all_bonds
    else:
        return cut_bonds_set

    # merge cliques
    for c in range(len(cliques) - 1):
        if c >= len(cliques):
            break
        for k in range(c + 1, len(cliques)):
            if k >= len(cliques):
                break
            if len(set(cliques[c]) & set(cliques[k])) > 0:
                cliques[c] = list(set(cliques[c]) | set(cliques[k]))
                cliques[k] = []
        cliques = [c for c in cliques if len(c) > 0]
    cliques = [c for c in cliques if len(c) > 0]

    # edges
    edges = []
    for bond in res:
        for c in range(len(cliques)):
            if bond[0][0] in cliques[c]:
                c1 = c
            if bond[0][1] in cliques[c]:
                c2 = c
        edges.append((c1, c2))
    for bond in breaks:
        for c in range(len(cliques)):
            if bond[0] in cliques[c]:
                c1 = c
            if bond[1] in cliques[c]:
                c2 = c
        edges.append((c1, c2))

    return cliques, edges

def fragment_graph_cutbonds(mol, cut_bonds_set):
    mol_num = len(list(mol.GetAtoms()))
    bond_set = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bond_set.append([a1, a2])

    left_bond_set = []

    for ele in bond_set:
        if [ele[0], ele[1]] not in cut_bonds_set and \
            [ele[1], ele[0]] not in cut_bonds_set:
                left_bond_set.append(ele)            

    left_bond_set = [list(ele) for ele in list(left_bond_set)]
    graph = defaultdict(list)

    for x, y in left_bond_set:
        graph[x].append(y)
        graph[y].append(x)

    visited = set()

    labels = [-1 for _ in range(mol_num)]

    def dfs(i, lb=-1):
        visited.add(i)
        labels[i] = lb
        for j in graph[i]:
            if j not in visited:
                dfs(j, lb)

    lb = 0
    for i in range(mol_num):
        if i not in visited:
            dfs(i, lb)
            lb += 1

    return labels

def frag_mol_brics(mol, addition_rule=False):
    #mol = Chem.MolFromSmiles(smi)
    cut_bonds_set = brics_decomp(mol, addition_rule=addition_rule)
    g_labels = fragment_graph_cutbonds(mol, cut_bonds_set)
    # get max label and get subsmiles
    group_num = max(g_labels) + 1


    frag_smiles_lst = []
    for g_id in range(group_num):
        g_mol_lst = []
        for i, l in enumerate(g_labels):
            if l == g_id:
                g_mol_lst.append(i)

        frag_smiles = Chem.MolFragmentToSmiles(mol, g_mol_lst, kekuleSmiles=True)
        frag_smiles_lst.append(frag_smiles)

    try:
    
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        coordiantes = mol.GetConformer().GetPositions()


        # the returned text should contain following info. each lines is
        # atom_index, atom_type, x, y, z, fragment_type

        text = ""

        for i in range(len(atoms)):
            text += f"{i} {atoms[i]} {coordiantes[i][0]} {coordiantes[i][1]} {coordiantes[i][2]} {frag_smiles_lst[g_labels[i]]}\n"
        return frag_smiles_lst, g_labels, text
    except:
        return frag_smiles_lst, g_labels, "Unknown"

    #print(text)


    #print(coordiantes)


    

def ifg_detect_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fgs = identify_functional_groups(mol)
    frag_smiles_lst = []
    for fg in fgs:
        atomdids = list(fg.atomIds)
        frag_smiles = Chem.MolFragmentToSmiles(mol, atomdids, kekuleSmiles=True)
        frag_smiles_lst.append(frag_smiles)
    return frag_smiles_lst


def frag_smi_files(smi_file_path, detect_frag_fg=True, save_file_prefix='', save_folder='', use_addtion_rules=False, mix_fg_brics=False):
    save_info_lst = []
    vocab = {}
    with open(smi_file_path, 'r') as sr:
        lines = sr.readlines()
        for line in tqdm(lines):
            smiles = line.strip()
            fg_label = process_fg(smiles)

            if detect_frag_fg:
                if mix_fg_brics:
                    try:
                        fragment_lst = ifg_detect_mol_brics(smiles)
                    except:
                        print(f'fail on mol {smiles}\n')
                        fragment_lst = frag_mol_brics(smiles, use_addtion_rules)
                else:
                    fragment_lst = frag_mol_brics(smiles, use_addtion_rules)
                frag_dict = {}
            else:
                fragment_lst = ifg_detect_mol(smiles)

            for frag_smi in fragment_lst:
                if frag_smi in vocab:
                    vocab[frag_smi] += 1
                else:
                    vocab[frag_smi] = 1
                if detect_frag_fg:
                    try:
                        frag_fg_label = process_fg(frag_smi)
                    except:
                        print(f'detect frag smile fail f{frag_smi}\n')
                        frag_fg_label = np.zeros(85, dtype=np.int64)
                    frag_dict[frag_smi] = frag_fg_label

            if detect_frag_fg:
                save_info_lst.append([smiles, fg_label, frag_dict])
            else:
                save_info_lst.append([smiles, fg_label])


    save_file = os.path.join(save_folder, f"{save_file_prefix}_frags.pickle")
    with open(save_file, 'wb') as handle:
        pickle.dump(save_info_lst, handle, protocol=pickle.HIGHEST_PROTOCOL)

    save_file = os.path.join(save_folder, f"{save_file_prefix}_vocab.pickle")
    with open(save_file, 'wb') as handle:
        pickle.dump(vocab, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # Load data (deserialize)
    # with open('filename.pickle', 'rb') as handle:
    #     unserialized_data = pickle.load(handle)


def ifg_detect_mol_brics(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fgs = identify_functional_groups(mol)
    fg_mol_lst = []
    frag_smiles_lst = []
    for fg in fgs:
        atomdids = list(fg.atomIds)
        fg_mol_lst.append(atomdids)
        frag_smiles = Chem.MolFragmentToSmiles(mol, atomdids, kekuleSmiles=True)
        frag_smiles_lst.append(frag_smiles)

    # print(f'fg list is {fg_mol_lst}\n')

    cut_bonds_set, all_bonds = brics_decomp(mol, addition_rule=True, return_all_bonds=True)


    g_labels_brics = fragment_graph_cutbonds(mol, cut_bonds_set)

    gl_dict = defaultdict(list)
    brics_group_num = 0
    for i, gl in enumerate(g_labels_brics):
        gl_dict[gl].append(i)

    brics_group_num = len(gl_dict)
    # print(f'orgin brics dict is {gl_dict}\n')
    # function group number



    # erase cut_bonds_set, if two node of an edge both occur in same function group
    remove_lst = []
    for cut_bond in cut_bonds_set:
        for fg in fg_mol_lst:
            if cut_bond[0] in fg and cut_bond[1] in fg:
                remove_lst.append(cut_bond)
    for cut_bond in remove_lst:
        cut_bonds_set.remove(cut_bond)
    # print(f'cut_bonds_set is {cut_bonds_set}\n')
    # add to cut_bonds_set new cut_bond
    for bond in all_bonds:
        for fg in fg_mol_lst:
            if (bond[0] in fg and bond[1] not in fg) or (bond[1] in fg and bond[0] not in fg): # keep such bonds
                cut_bonds_set.append(bond)
    # print(f'cut_bonds_set is {cut_bonds_set}\n')

    g_labels = fragment_graph_cutbonds(mol, cut_bonds_set)

    gl_dict = defaultdict(list)
    # brics_group_num = 0
    for i, gl in enumerate(g_labels):
        gl_dict[gl].append(i)

    # print(f'mix bircs and fg is {gl_dict}\n')
    group_num = max(g_labels) + 1

    frag_smiles_lst = []
    for g_id in range(group_num):
        g_mol_lst = []
        for i, l in enumerate(g_labels):
            if l == g_id:
                g_mol_lst.append(i)

        frag_smiles = Chem.MolFragmentToSmiles(mol, g_mol_lst, kekuleSmiles=True)
        frag_smiles_lst.append(frag_smiles)


    return frag_smiles_lst