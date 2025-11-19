"""
AutoDock Vina docking interface for molecular evaluation.

This module provides classes and functions for performing molecular docking
using AutoDock Vina to evaluate binding affinity of generated molecules.
"""

import os
import contextlib
import tempfile
import subprocess
from typing import Optional, Union

try:
    from openbabel import pybel
    from meeko import MoleculePreparation, obutils
    from vina import Vina
    import rdkit.Chem as Chem
    from rdkit.Chem import AllChem
    import AutoDockTools
    DEPENDENCIES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Some docking dependencies not available: {e}")
    DEPENDENCIES_AVAILABLE = False

from evaluation.docking_qvina import get_random_id, BaseDockingTask


def suppress_stdout(func):
    """Decorator to suppress stdout during function execution."""
    def wrapper(*args, **kwargs):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*args, **kwargs)
    return wrapper


class PrepLig:
    """
    Ligand preparation class for molecular docking.
    
    Handles conversion between molecular formats and preparation
    of ligand structures for docking calculations.
    """
    
    def __init__(self, input_mol: Union[str, bytes], mol_format: str):
        """
        Initialize ligand preparation.
        
        Args:
            input_mol: Molecular input (SMILES string or file path)
            mol_format: Format of input molecule ('smi' or 'sdf')
        """
        if not DEPENDENCIES_AVAILABLE:
            raise ImportError("Required dependencies not available for ligand preparation")
            
        if mol_format == 'smi':
            self.ob_mol = pybel.readstring('smi', input_mol)
        elif mol_format == 'sdf': 
            self.ob_mol = next(pybel.readfile(mol_format, input_mol))
        else:
            raise ValueError(f'Molecular format {mol_format} not supported')
        
    def addH(self, polaronly: bool = False, correctforph: bool = True, PH: float = 7.0):
        """Add hydrogen atoms to the molecule."""
        self.ob_mol.OBMol.AddHydrogens(polaronly, correctforph, PH)
        obutils.writeMolecule(self.ob_mol.OBMol, 'tmp_h.sdf')

    def gen_conf(self):
        """Generate 3D conformer using RDKit."""
        sdf_block = self.ob_mol.write('sdf')
        rdkit_mol = Chem.MolFromMolBlock(sdf_block, removeHs=False)
        
        if rdkit_mol is None:
            raise ValueError("Failed to parse molecule")
            
        AllChem.EmbedMolecule(rdkit_mol, Chem.rdDistGeom.ETKDGv3())
        self.ob_mol = pybel.readstring('sdf', Chem.MolToMolBlock(rdkit_mol))
        obutils.writeMolecule(self.ob_mol.OBMol, 'conf_h.sdf')

    @suppress_stdout
    def get_pdbqt(self, lig_pdbqt: Optional[str] = None):
        """Convert molecule to PDBQT format for docking."""
        preparator = MoleculePreparation()
        preparator.prepare(self.ob_mol.OBMol)
        if lig_pdbqt is not None: 
            preparator.write_pdbqt_file(lig_pdbqt)
            return 
        else: 
            return preparator.write_pdbqt_string()
        

class PrepProt(object): 
    def __init__(self, pdb_file): 
        self.prot = pdb_file
        self.prot_pqr = "None"
    
    def del_water(self, dry_pdb_file): # optional
        with open(self.prot) as f: 
            lines = [l for l in f.readlines() if l.startswith('ATOM') or l.startswith('HETATM')] 
            dry_lines = [l for l in lines if not 'HOH' in l]
        
        with open(dry_pdb_file, 'w') as f:
            f.write(''.join(dry_lines))
        self.prot = dry_pdb_file
        
    def addH(self, prot_pqr):  # call pdb2pqr
        self.prot_pqr = prot_pqr
        #subprocess.call(['pdb2pqr30','--ff=AMBER',self.prot, self.prot_pqr])
        subprocess.Popen(['pdb2pqr30','--ff=AMBER',self.prot, self.prot_pqr],
                         stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL).communicate()

    def get_pdbqt(self, prot_pdbqt):
        prepare_receptor = os.path.join(AutoDockTools.__path__[0], 'Utilities24/prepare_receptor4.py')
        print(prepare_receptor)
        # print(prepare_receptor)
        # print(prot_pdbqt)

        # # use subproces
        print(self.prot_pqr)

        # subprocess.call(['python3', prepare_receptor, '-r', self.prot_pqr, '-o', prot_pdbqt])

        subprocess.Popen(['python3', prepare_receptor, '-r', self.prot_pqr, '-o', prot_pdbqt],
                         stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL).communicate()
        # check whether success




class VinaDock(object): 
    def __init__(self, lig_pdbqt, prot_pdbqt): 
        self.lig_pdbqt = lig_pdbqt
        self.prot_pdbqt = prot_pdbqt
    
    def _max_min_pdb(self, pdb, buffer):
        with open(pdb, 'r') as f: 
            lines = [l for l in f.readlines() if l.startswith('ATOM') or l.startswith('HEATATM')]
            xs = [float(l[31:39]) for l in lines]
            ys = [float(l[39:47]) for l in lines]
            zs = [float(l[47:55]) for l in lines]
            print(max(xs), min(xs))
            print(max(ys), min(ys))
            print(max(zs), min(zs))
            pocket_center = [(max(xs) + min(xs))/2, (max(ys) + min(ys))/2, (max(zs) + min(zs))/2]
            box_size = [(max(xs) - min(xs)) + buffer, (max(ys) - min(ys)) + buffer, (max(zs) - min(zs)) + buffer]
            return pocket_center, box_size
    
    def get_box(self, ref=None, buffer=0):
        '''
        ref: reference pdb to define pocket. 
        buffer: buffer size to add 

        if ref is not None: 
            get the max and min on x, y, z axis in ref pdb and add buffer to each dimension 
        else: 
            use the entire protein to define pocket 
        '''
        if ref is None: 
            ref = self.prot_pdbqt
        self.pocket_center, self.box_size = self._max_min_pdb(ref, buffer)
        print(self.pocket_center, self.box_size)

    def dock(self, score_func='vina', seed=0, mode='dock', exhaustiveness=8, save_pose=False, **kwargs):  # seed=0 mean random seed
        v = Vina(sf_name=score_func, seed=seed, verbosity=0, **kwargs)
        v.set_receptor(self.prot_pdbqt)
        v.set_ligand_from_file(self.lig_pdbqt)
        v.compute_vina_maps(center=self.pocket_center, box_size=self.box_size)
        #print(v)
        if mode == 'score_only': 
            score = v.score()[0]
        elif mode == 'minimize':
            score = v.optimize()[0]
        elif mode == 'dock':
            v.dock(exhaustiveness=exhaustiveness, n_poses=20)
            #print(v)
            #print(v.energies(n_poses=9))
            score = v.energies(n_poses=1)[0][0]
        else:
            raise ValueError
        
        if not save_pose: 
            return score
        else: 
            if mode == 'score_only': 
                pose = None 
            elif mode == 'minimize': 
                tmp = tempfile.NamedTemporaryFile()
                with open(tmp.name, 'w') as f: 
                    v.write_pose(tmp.name, overwrite=True)             
                with open(tmp.name, 'r') as f: 
                    pose = f.read()
   
            elif mode == 'dock': 
                pose = v.poses(n_poses=1)
            else:
                raise ValueError
            return score, pose


class VinaDockingTask(BaseDockingTask):

    
    @classmethod
    def from_original_data(cls, data, ligand_root='./data/crossdocked_pocket10', protein_root='./data/crossdocked',
                           **kwargs):
        protein_fn = os.path.join(
            os.path.dirname(data.ligand_filename),
            os.path.basename(data.ligand_filename)[:10] + '.pdb'
        )
        protein_path = os.path.join(protein_root, protein_fn)

        ligand_path = os.path.join(ligand_root, data.ligand_filename)
        ligand_rdmol = next(iter(Chem.SDMolSupplier(ligand_path)))
        return cls(protein_path, ligand_rdmol, **kwargs)

    @classmethod
    def from_generated_mol(cls, ligand_rdmol, ligand_filename, protein_root='./data/crossdocked', **kwargs):
        # load original pdb
        protein_fn = os.path.join(
            os.path.dirname(ligand_filename),
            os.path.basename(ligand_filename)[:10] + '.pdb'  # PDBId_Chain_rec.pdb
        )
        protein_path = os.path.join(protein_root, protein_fn)
        return cls(protein_path, ligand_rdmol, **kwargs)
    
    @classmethod
    def from_2d_mol(cls, gen_rdmol, protein_root='./data/crossdocked', **kwargs):
        # load original pdb
        #protein_root = os.path.join(protein_root, target)

        protein_path = os.path.join(protein_root, "receptor.pdb")
        ligand_path = os.path.join(protein_root, "crystal_ligand.sdf")
        ligand_rdmol = next(iter(Chem.SDMolSupplier(ligand_path, sanitize=False)))
        pos = ligand_rdmol.GetConformer(0).GetPositions()

        center = (pos.max(0) + pos.min(0)) / 2
        
        # generate conformation for gen_rdmol
        # Generate conformers
        num_conformers = 3  # Number of conformers to generate
        # add H

        gen_rdmol = Chem.AddHs(gen_rdmol, addCoords=True)
        AllChem.EmbedMultipleConfs(gen_rdmol, numConfs=num_conformers, randomSeed=42)

        # Optimize conformers
        AllChem.MMFFOptimizeMoleculeConfs(gen_rdmol)

        # remove H

        gen_rdmol = Chem.RemoveHs(gen_rdmol)
        
        tmp_dir = os.path.join(protein_root, 'tmp')

        os.makedirs(tmp_dir, exist_ok=True)

        return cls(protein_path, gen_rdmol, tmp_dir=tmp_dir, center=center, size_factor = None, **kwargs)


    @classmethod
    def from_2d_mol_crossdocked(cls, gen_rdmol, protein_path, ligand_path, **kwargs):
        ligand_rdmol = next(iter(Chem.SDMolSupplier(ligand_path, sanitize=False)))
        ligand_rdmol = Chem.AddHs(ligand_rdmol, addCoords=True)
        pos = ligand_rdmol.GetConformer(0).GetPositions()

        center = (pos.max(0) + pos.min(0)) / 2
        
        num_conformers = 3
        gen_rdmol = Chem.AddHs(gen_rdmol, addCoords=True)
        AllChem.EmbedMultipleConfs(gen_rdmol, numConfs=num_conformers, randomSeed=42)
        AllChem.MMFFOptimizeMoleculeConfs(gen_rdmol)
        gen_rdmol = Chem.RemoveHs(gen_rdmol)
        
        protein_dir = os.path.dirname(protein_path)
        tmp_dir = os.path.join(protein_dir, 'tmp')
        os.makedirs(tmp_dir, exist_ok=True)

        return cls(protein_path, gen_rdmol, tmp_dir=tmp_dir, center=center, size_factor = None, **kwargs)

    def __init__(self, protein_path, ligand_rdmol, tmp_dir='./tmp', center=None,
                 size_factor=1., buffer=5.0):
        super().__init__(protein_path, ligand_rdmol)
        # self.conda_env = conda_env
        self.tmp_dir = os.path.realpath(tmp_dir)
        os.makedirs(tmp_dir, exist_ok=True)

        self.task_id = get_random_id()
        self.receptor_id = self.task_id + '_receptor'
        self.ligand_id = self.task_id + '_ligand'

        self.receptor_path = protein_path
        self.ligand_path = os.path.join(self.tmp_dir, self.ligand_id + '.sdf')

        self.recon_ligand_mol = ligand_rdmol
        ligand_rdmol = Chem.AddHs(ligand_rdmol, addCoords=True)

        sdf_writer = Chem.SDWriter(self.ligand_path)
        sdf_writer.write(ligand_rdmol)
        sdf_writer.close()
        self.ligand_rdmol = ligand_rdmol

        pos = ligand_rdmol.GetConformer(0).GetPositions()
        if center is None:
            self.center = (pos.max(0) + pos.min(0)) / 2
        else:
            self.center = center

        if size_factor is None:
            self.size_x, self.size_y, self.size_z = 20, 20, 20
        else:
            self.size_x, self.size_y, self.size_z = (pos.max(0) - pos.min(0)) * size_factor + buffer

        self.proc = None
        self.results = None
        self.output = None
        self.error_output = None
        self.docked_sdf_path = None

    def run(self, mode='dock', exhaustiveness=16, **kwargs):
        ligand_pdbqt = self.ligand_path[:-4] + '.pdbqt'
        protein_pqr = self.receptor_path[:-4] + '.pqr'
        protein_pdbqt = self.receptor_path[:-4] + '.pdbqt'

        lig = PrepLig(self.ligand_path, 'sdf')
        lig.get_pdbqt(ligand_pdbqt)

        prot = PrepProt(self.receptor_path)
        if not os.path.exists(protein_pqr):
            prot.addH(protein_pqr)
        if not os.path.exists(protein_pdbqt):
            prot.get_pdbqt(protein_pdbqt)

        dock = VinaDock(ligand_pdbqt, protein_pdbqt)
        dock.pocket_center, dock.box_size = self.center, [self.size_x, self.size_y, self.size_z]
        score, pose = dock.dock(score_func='vina', mode=mode, exhaustiveness=exhaustiveness,cpu=1, save_pose=True, **kwargs)

        return [{'affinity': score, 'pose': pose}]



    

