#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SMILES Structure Repair and Validation Utility

Official Research Code Repository - Chemical Structure Correction Module

This module provides comprehensive functionality for repairing and validating
problematic SMILES chemical structures in large-scale chemical databases.
It addresses common issues such as valence errors, aromatic system problems,
charge distribution inconsistencies, and structural anomalies that can occur
in chemical data processing pipelines.

Key Repair Functions:
- Valence error detection and automatic correction
- Aromatic system validation and repair
- Charge distribution optimization and balancing
- Radical and unpaired electron handling
- Stereochemistry preservation during repair
- Custom error pattern recognition and fixing

Technical Features:
- RDKit-based molecular structure manipulation
- Systematic error detection and classification
- Multi-step repair algorithms with fallback mechanisms
- Preservation of chemical information during repair
- Batch processing for large chemical datasets
- Comprehensive error logging and reporting

Research Applications:
- Large-scale chemical database cleaning and standardization
- Chemical structure validation for coordination complexes
- Data preprocessing for chemical informatics applications
- Quality control for chemical structure pipelines

Academic Citation:
If you use this code in your research, please cite the corresponding paper:
[Paper citation to be added upon publication]

Author: ChemDB Research Team
Institution: Shanghai Jiao Tong University
Created: 2025
Version: 1.0

License: Academic Open Source License

Copyright (c) 2025 ChemDB Research Project
For questions regarding this research or licensing, please contact:
jack111@sjtu.edu.cn
"""

try:
    import pandas as pd
except ImportError:
    # 创建一个假的pandas模块
    class FakePandas:
        def __getattr__(self, name):
            return lambda *args, **kwargs: None
    pd = FakePandas()
from rdkit import Chem
from rdkit.Chem import Draw
import warnings
from rdkit.Chem import rdmolops
from rdkit import RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
import math
from PIL import Image
import sys
import os
from contextlib import redirect_stdout
import io
from tqdm import tqdm
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Enable pandas tqdm progress bar
tqdm.pandas()
warnings.filterwarnings("ignore", category=UserWarning, module='rdkit')
RDLogger.DisableLog('rdApp.warning')

class MolecularFormulaRepairer:
    """
    Advanced molecular structure repair and validation system for chemical databases.
    
    This class provides comprehensive functionality for detecting and repairing
    common structural problems in SMILES chemical representations. It addresses
    valence errors, aromatic system inconsistencies, charge distribution problems,
    and other structural anomalies that can occur in large-scale chemical databases.
    
    Key Repair Capabilities:
    - Valence error detection and automatic correction
    - Aromatic system validation and standardization
    - Charge distribution optimization and balancing
    - Ring structure analysis and repair
    - Special structure pattern recognition (porphyrins, phthalocyanines)
    - Cumulated double bond system correction
    - Hydrogen atom addition/removal for valence correction
    
    Attributes:
        smiles (str): Input SMILES string to be repaired
        debug (bool): Enable debug output for troubleshooting
        draw (bool): Enable molecular visualization during repair
        typical_valences (dict): Element-specific valence rules
        error_found (set): Collection of error types encountered during repair
        need_kill (bool): Flag indicating if molecule should be discarded
        readable (bool): Flag indicating if molecule is parseable
    """
    
    def __init__(self, smiles,debug = False,draw = False):
        self.smiles = smiles
        self.exempt_atom_idx = set()
        self.debug = debug
        self.draw = draw
        self.typical_valences = {
        1: [1],            # H (Hydrogen)
        2: [0],            # He (Helium)
        3: [1],            # Li (Lithium)
        4: [2],            # Be (Beryllium)
        5: [3],            # B (Boron)
        6: [4],            # C (Carbon)
        7: [3,5],         # N (Nitrogen)
        8: [2],            # O (Oxygen)
        9: [1],            # F (Fluorine)
        10: [0],           # Ne (Neon)
        11: [1],           # Na (Sodium)
        12: [2],           # Mg (Magnesium)
        13: [3],           # Al (Aluminum)
        14: [4],           # Si (Silicon)
        15: [3, 5],        # P (Phosphorus)
        16: [2, 4, 6],     # S (Sulfur)
        17: [1],           # Cl (Chlorine)
        18: [0],           # Ar (Argon)
        19: [1],           # K (Potassium)
        20: [2],           # Ca (Calcium)
        21: [3],           # Sc (Scandium)
        22: [2, 3, 4],     # Ti (Titanium)
        23: [2, 3, 5],     # V (Vanadium)
        24: [2, 3, 6],     # Cr (Chromium)
        25: [2, 3, 4, 7],  # Mn (Manganese)
        26: [2, 3, 6],     # Fe (Iron)
        27: [2, 3],        # Co (Cobalt)
        28: [2, 3, 4],     # Ni (Nickel)
        29: [1, 2],        # Cu (Copper)
        30: [2],           # Zn (Zinc)
        31: [3],           # Ga (Gallium)
        32: [2, 4],        # Ge (Germanium)
        33: [3, 5],        # As (Arsenic)
        34: [2, 4, 6],     # Se (Selenium)
        35: [1, 3, 5, 7],  # Br (Bromine)
        36: [0],           # Kr (Krypton)
        37: [1],           # Rb (Rubidium)
        38: [2],           # Sr (Strontium)
        39: [3],           # Y (Yttrium)
        40: [4],           # Zr (Zirconium)
        41: [2, 3, 5],     # Nb (Niobium)
        42: [2, 3, 4, 6],  # Mo (Molybdenum)
        43: [2, 3, 4, 6, 7], # Tc (Technetium)
        44: [2, 3, 4, 6, 8], # Ru (Ruthenium)
        45: [2, 3, 4],     # Rh (Rhodium)
        46: [2, 4],        # Pd (Palladium)
        47: [1],           # Ag (Silver)
        48: [2],           # Cd (Cadmium)
        49: [3],           # In (Indium)
        50: [2, 4],        # Sn (Tin)
        51: [3, 5],        # Sb (Antimony)
        52: [2, 4, 6],     # Te (Tellurium)
        53: [1, 5, 7],     # I (Iodine)
        54: [0],           # Xe (Xenon)
        55: [1],           # Cs (Cesium)
        56: [2],           # Ba (Barium)
        57: [3],           # La (Lanthanum)
        58: [3, 4],        # Ce (Cerium)
        59: [3],           # Pr (Praseodymium)
        60: [3],           # Nd (Neodymium)
        61: [3],           # Pm (Promethium)
        62: [3],           # Sm (Samarium)
        63: [3],           # Eu (Europium)
        64: [3],           # Gd (Gadolinium)
        65: [3],           # Tb (Terbium)
        66: [3],           # Dy (Dysprosium)
        67: [3],           # Ho (Holmium)
        68: [3],           # Er (Erbium)
        69: [3],           # Tm (Thulium)
        70: [3],           # Yb (Ytterbium)
        71: [3],           # Lu (Lutetium)
        72: [4],           # Hf (Hafnium)
        73: [5],           # Ta (Tantalum)
        74: [6],           # W (Tungsten)
        75: [4, 7],        # Re (Rhenium)
        76: [2, 3, 4, 6, 8], # Os (Osmium)
        77: [2, 3, 4, 6],  # Ir (Iridium)
        78: [2, 4],        # Pt (Platinum)
        79: [1, 3],        # Au (Gold)
        80: [1, 2],        # Hg (Mercury)
        81: [1, 3],        # Tl (Thallium)
        82: [2, 4],        # Pb (Lead)
        83: [3, 5],        # Bi (Bismuth)
        84: [2, 4, 6],     # Po (Polonium)
        85: [1, 5, 7],     # At (Astatine)
        86: [0],           # Rn (Radon)
        87: [1],           # Fr (Francium)
        88: [2],           # Ra (Radium)
        89: [3],           # Ac (Actinium)
        90: [4],           # Th (Thorium)
        91: [4, 5],        # Pa (Protactinium)
        92: [4, 6],        # U (Uranium)
        93: [5, 6],        # Np (Neptunium)
        94: [3, 4, 6],     # Pu (Plutonium)
        95: [3, 4, 5, 6],  # Am (Americium)
        96: [3, 4],        # Cm (Curium)
        97: [3, 4],        # Bk (Berkelium)
        98: [3],           # Cf (Californium)
        99: [3],           # Es (Einsteinium)
        100: [3],          # Fm (Fermium)
        101: [3],          # Md (Mendelevium)
        102: [3],          # No (Nobelium)
        103: [3],          # Lr (Lawrencium)
        104: [4],          # Rf (Rutherfordium)
        105: [5],          # Db (Dubnium)
        106: [6],          # Sg (Seaborgium)
        107: [7],          # Bh (Bohrium)
        108: [8],          # Hs (Hassium)
        109: [8],          # Mt (Meitnerium)
        110: [9],          # Ds (Darmstadtium)
        111: [9],          # Rg (Roentgenium)
        112: [10],         # Cn (Copernicium)
        113: [1, 3],       # Nh (Nihonium)
        114: [4],          # Fl (Flerovium)
        115: [1, 3],       # Mc (Moscovium)
        116: [2, 4],       # Lv (Livermorium)
        117: [1, 3, 5, 7], # Ts (Tennessine)
        118: [0],          # Og (Oganesson)
    }
        self.anion = ['O']
        self.atoms_to_modify = []
        self.error_found = set()
        self.stable_isotopes = {
            1: 1,   # H: Hydrogen-1
            6: 12,  # C: Carbon-12
            7: 14,  # N: Nitrogen-14
            8: 16,  # O: Oxygen-16
            16: 32, # S: Sulfur-32
            9: 19,  # F: Fluorine-19
            17: 35, # Cl: Chlorine-35
            35: 79, # Br: Bromine-79
            53: 127 # I: Iodine-127
        }
        self.has_isotope = False
        self.culmul_map = {'C[N+]1=C=[N+](C=C1)C':'CN1CN(C)C=C1','C1=S=CCC=1':'C1=CC=CS1','C1CN=S=N1':'C1CNSN1',\
              'C1=CCCC=1':'C1=CC=C[C-]1','C1=C=CC=C=1':'C1=CC=C[C-]1','C1=CCSC=1':'C1=CC=CS1',\
                'C1=CCOC=1':'C1=COC=C1','C1=CC=C[S+]=1':'C1=CC=CS1','C1=CSCN=1':'C1=CSC=N1',\
                'C1=CCC=S=1':'C1=CC=CS1','C1=CCOC=1':'C1=COC=C1','C1=CC=CCN=1':'C1=CC=NC=C1',\
                'C1=CNOC=1':'C1=CON=C1','C1=CCNC=1':'C1=CNC=C1','C1=CC=CP=1':'C1=CC=CP1',\
                'C1=CNNC=1':'C1=CNN=C1','C1=CC=C[NH+]=1':'C1=CNC=C1','C1=Ncc[NH+]=1':'C1NC=CN1',\
                'C1=COCN=1':'C1=COC=N1','C1=CN=[N+]=C1':'C1=CNN=C1','C1=NC=CCN=1':'C1=NC=CC=N1',\
                'C1=CNCNC=1':'C1=CN=CN=C1','C1CN=S=N1':'C1CNSN1','C1=CNNN=1':'C1=NNN=C1','C1=[N+]=NNN1':'C1=NN=NN1',\
                'C1=CCN=CN=1':'C1=NC=NC=C1','C1=CC=[S+]C=1':'C1=CC=CS1'
                }
        self.culmul_map2 = {'C1=CC=CCN=1':'C1=CC=NC=C1','C1=CC=CC=1':'C1=CC=C[C-]1','C1=CC=NCC=1':'C1=CC=NC=C1',\
               'C1=CC=CNC=1':'C1=CC=NC=C1','C1=CC=COC=1':'C1=CC=C[O+]=C1','C1=CN=CNC=1':'C1=CN=CN=C1',\
                'C1=CC=NNC=1':'C1=CC=NN=C1','C[N+]1=C=[N+](C=C1)C':'CN1CN(C)C=C1'}
        self.culmul_bond_flag = False

    def is_halogen(self,atom):
        """
        Determine if an atom is a halogen element.
        
        This method checks if the given atom belongs to the halogen group
        (F, Cl, Br, I, At) based on its atomic number.
        
        Args:
            atom (rdkit.Chem.Atom): RDKit atom object to be checked
            
        Returns:
            bool: True if atom is a halogen, False otherwise
        """
        halogens = {9, 17, 35, 53, 85}
        # Get atomic number of the atom
        atomic_num = atom.GetAtomicNum()
        # Check if it's a halogen atom
        return atomic_num in halogens
    def is_CH(self,atom):
        """
        Identify carbon atoms with exactly one hydrogen (CH groups).
        
        This method determines if a carbon atom has exactly one hydrogen
        neighbor and three total connections, which is characteristic of
        unsaturated carbon centers that may require repair.
        
        Args:
            atom (rdkit.Chem.Atom): RDKit atom object to be analyzed
            
        Returns:
            bool: True if atom is a CH carbon with degree 3, False otherwise
        """
        if atom.GetSymbol() == 'C':  # Atom is carbon
            # Get neighbor atoms connected to this atom
            neighbors = atom.GetNeighbors()
            # Count hydrogen atoms connected to carbon
            hydrogen_count = sum(1 for nbr in neighbors if nbr.GetSymbol() == 'H')
            
            return hydrogen_count == 1 and atom.GetDegree() == 3 # If carbon atom is connected to only 1 hydrogen and has exactly 3 bonds
        return False
    def is_neutral_carbon(self,atom):
        """
        Check if an atom is a neutral carbon atom.
        
        Args:
            atom (rdkit.Chem.Atom): RDKit atom object to be checked
            
        Returns:
            bool: True if atom is carbon with formal charge 0, False otherwise
        """
        return atom.GetSymbol() == 'C' and atom.GetFormalCharge() == 0
    def reset_mol(self,mol):
        """
        Reset and normalize molecule object through SMILES conversion.
        
        This method performs a complete molecule reset by converting to SMILES
        and back to molecule object, which can help resolve certain structural
        inconsistencies and aromaticity issues.
        
        Args:
            mol (rdkit.Chem.RWMol): Editable molecule object to be reset
            
        Returns:
            rdkit.Chem.RWMol: Reset editable molecule object
        """
        mol = mol.GetMol()
        new_smiles = Chem.MolToSmiles(mol,allHsExplicit=False)
        
        mol_ori = Chem.MolFromSmiles(new_smiles)

        try:
            mol = Chem.AddHs(mol_ori)
            mol = Chem.RWMol(mol)  # Convert to editable molecule
        except:
            mol = Chem.RWMol(mol)  # Convert to editable molecule
        return mol
    # Assuming mol is an RDKit Mol object
    def is_halogen_and_inorganic(self,mol):
        """
        Determine if molecule contains halogens but no carbon (inorganic halide).
        
        This method identifies inorganic halogen compounds that should be
        excluded from organic structure repair algorithms.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object to be analyzed
            
        Returns:
            bool: True if molecule contains halogens but no carbon, False otherwise
        """
        # Define halogen element symbols
        halogens = ['F', 'Cl', 'Br', 'I', 'At']
        
        # Check if molecule contains halogen elements
        contains_halogen = any(atom.GetSymbol() in halogens for atom in mol.GetAtoms())
        
        # Check if molecule contains carbon to determine if inorganic
        contains_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
        
        # Return whether contains halogens and no carbon (inorganic)
        return contains_halogen and not contains_carbon
    def is_legitimate(self,mol):
        """
        Validate if a molecule object represents a legitimate chemical structure.
        
        This method performs comprehensive validation including SMILES
        round-trip conversion and Kekulization checks to ensure the
        molecule represents a valid chemical structure.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object to be validated
            
        Returns:
            bool: True if molecule is valid, False if problematic
        """
        try:
            with open(os.devnull, 'w') as fnull:
                with redirect_stdout(fnull):
                    test_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        except  Chem.KekulizeError:
            return False
        if test_mol is None:
            return False
        else:
            return True
    # Function to determine if a ring contains double bonds
    def ring_has_double_bond(self,ring, mol):
        """
        Check if a ring structure contains any double bonds or aromatic bonds.
        
        This method examines all bonds within a ring to determine if any
        are double bonds or aromatic bonds, which affects repair strategies.
        
        Args:
            ring (tuple): Ring atom indices
            mol (rdkit.Chem.Mol): Molecule containing the ring
            
        Returns:
            bool: True if ring contains double/aromatic bonds, False otherwise
        """
        for i in range(len(ring)):
            atom1_idx = ring[i]
            atom2_idx = ring[(i + 1) % len(ring)]  # Ensure last atom connects to first atom
            
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE or bond.GetIsAromatic():  # Check if double bond
                return True
        return False
    def has_two_double_bonds(self,atom):
        """
        Determine if an atom is connected to two double bonds.
        
        This method identifies atoms involved in cumulated double bond
        systems (e.g., allenes, ketenes) which require special handling.
        
        Args:
            atom (rdkit.Chem.Atom): RDKit Atom object to be analyzed
            
        Returns:
            bool: True if atom has two double bonds, False otherwise
        """
        # Get all bonds around the atom
        bonds = atom.GetBonds()
        
        # Count the number of double bonds
        double_bond_count = sum(1 for bond in bonds if bond.GetBondType() == Chem.BondType.DOUBLE)
        
        # Check if there are two double bonds
        return double_bond_count == 2
    def contains_metal(self,mol):
        """
        Check if molecule contains any metallic elements.
        
        This method scans the molecule for metal atoms which may require
        special handling or exclusion from certain repair algorithms.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object to be checked
            
        Returns:
            bool: True if molecule contains metals, False otherwise
        """
        # Define common metal element symbols
        metals = {'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra', 
                'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 
                'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 
                 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 
                'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 
                'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 
                'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 
                'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'}
        # Check if metal atoms exist in the molecule
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in metals:
                return True
        return False
    def convert_to_benzene(self,ring, mol_edit):
        """
        Convert a ring structure to a benzene-like aromatic structure.
        
        This method modifies the bond types in a ring to convert it to a
        benzene-like structure with alternating double bonds. It also marks
        the atoms in the ring as aromatic.
        
        Args:
            ring (tuple): Ring atom indices
            mol_edit (rdkit.Chem.RWMol): Editable molecule object containing the ring
            
        Note:
            This method modifies the molecule in-place by updating bond types
            and atomic properties to reflect the aromatic nature of benzene.
        """
        for i in range(len(ring)):
            atom1 = ring[i]
            atom2 = ring[(i + 1) % len(ring)]  # Next atom in ring (cyclic connection)
            bond = mol_edit.GetBondBetweenAtoms(atom1, atom2)
            if bond:
                bond.SetBondType(Chem.BondType.AROMATIC)  # Set as aromatic bond
                bond.SetIsAromatic(True)  # Mark bond as aromatic
        # Mark atoms in ring as aromatic atoms
        for atom_idx in ring:
            atom = mol_edit.GetAtomWithIdx(atom_idx)
            atom.SetIsAromatic(True)
    def merge_chains(self,matches):
        """
        Merge fragmented chain structures into unified chains.
        
        This method takes matched fragments of chain structures and merges
        them into unified chains based on shared atom indices. It's used to
        reconstruct continuous chains from fragmented substructure matches.
        
        Args:
            matches (list): List of matched fragment sets
            
        Returns:
            list: List of merged chains, each represented by a sorted list of atom indices
        """
        chains = []
        for match in matches:
            added_to_chain = False
            for chain in chains:
                # Check if current fragment has common atoms with existing chain
                if any(atom in chain for atom in match):
                    # Merge into existing chain
                    chain.update(match)
                    added_to_chain = True
                    break
            if not added_to_chain:
                # If current fragment has no common atoms with any chain, create new chain
                chains.append(set(match))
        # Convert results to list and sort
        return [sorted(list(set(chain))) for chain in chains]
    def get_atom_valence(self,atom):
        """
        Calculate the valence of an atom based on its bonds.
        
        This method calculates the total bond count for an atom, considering
        double and triple bonds as multiple counts. It also adds the number
        of explicit hydrogens for carbon atoms to determine the complete
        valence state of the atom.
        
        Args:
            atom (rdkit.Chem.Atom): RDKit atom object for which to calculate valence
            
        Returns:
            float: Calculated valence of the atom
        """
        # Calculate total bond count for atom (double bonds count as 2, triple bonds count as 3)
        valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if atom.GetSymbol() == 'C':
            valence += atom.GetTotalNumHs()
        return valence
    def find_matching_fragments(self,mol):
        """
        Identify and validate [CH][CH] fragment matches in a molecule.
        
        This method searches for substructures matching the [CH][CH] pattern
        and validates them based on the actual bond counts of the carbon atoms.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object to be searched
            
        Returns:
            list: List of valid matching fragment atom index tuples
        """
        # Define SMARTS pattern for [CH][CH] structure
        pattern = Chem.MolFromSmarts('[CH;R0][CH;R0]')
        matches = mol.GetSubstructMatches(pattern)
        
        valid_matches = []
        for match in matches:
            atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            # Check actual bond count of the two C atoms
            if all(self.get_atom_valence(atom) < 4 for atom in atoms):
                valid_matches.append(match)
        valid_matches = self.merge_chains(valid_matches)
        return valid_matches
    def find_closest_element(self,lst, target):
        """ Find the list element closest to the target value """
        if not lst:
            raise ValueError("The list is empty.")
        closest_element = min(lst, key=lambda x: abs(x - target))
        return closest_element
    def exam_readable(self):
        mol_ori = Chem.MolFromSmiles(self.smiles)
        try:
            self.mol = Chem.AddHs(mol_ori)
        except:
            self.need_kill = True
            self.readable = False
            return self.need_kill,self.readable
        self.mol_edit = Chem.RWMol(self.mol)
        self.need_kill = False
        self.readable = True

        if mol_ori is None or self.contains_metal(self.mol):
            if self.debug:
                print('There is wrong in smiles: ',self.smiles)
            self.need_kill = True
    def draw_exempt_atom(self):
        # Clone molecule object to avoid modifying original structure
        mol = Chem.RWMol(self.mol_edit)
        
        # Mark valence (bond count + absolute value of charge)
        for atom in mol.GetAtoms():
            valence = math.floor(sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]))
            valence += abs(atom.GetFormalCharge())
            #atom.SetProp("atomNote", str(valence))
            atom.SetProp("atomNote", str(atom.GetIdx()))
        
        # Get atom indices that need highlighting
        highlight_atoms = [
            atom.GetIdx() 
            for atom in mol.GetAtoms() 
            if atom.GetIntProp('fixed_already') == 1
        ]

        # Use vector drawing engine to ensure annotations are displayed
        drawer = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
        drawer.FinishDrawing()

        # Convert to PIL image object
        img = Image.open(io.BytesIO(drawer.GetDrawingText()))
        return img
    def update_exemption(self,atom_idx):
        self.exempt_atom_idx.update(atom_idx)
    def get_culmul_smart(self,smiles):
        # Convert SMILES to SMART, keeping only double bonds connected to the same atom
        # Input molecule
        mol = Chem.MolFromSmiles(smiles)

        # Count the number of double bonds each atom participates in
        atom_double_bonds = [0] * mol.GetNumAtoms()
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                begin = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                atom_double_bonds[begin] += 1
                atom_double_bonds[end] += 1

        # Modify bond query conditions
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            # Handle double bonds
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Check if connected to an atom with at least two double bonds
                if atom_double_bonds[begin_idx] >= 2 or atom_double_bonds[end_idx] >= 2:
                    continue
                else:
                    bond.SetBondType(Chem.BondType.UNSPECIFIED)
            # Handle other types of bonds
            else:
                bond.SetBondType(Chem.BondType.UNSPECIFIED)

        # Generate SMARTS string
        smarts = Chem.MolToSmarts(mol)
        return smarts
    
    def smiles_to_smarts_with_wildcard_bonds(self,special_sub):
        """
        Convert SMILES to SMARTS ignoring bond types and aromaticity (all bonds replaced with wildcard bonds ~)
        """
        mol = Chem.MolFromSmiles(special_sub)
        if not mol:
            raise ValueError("Invalid SMILES string")

        # Set all atoms to non-aromatic
        for atom in mol.GetAtoms():
            atom.SetIsAromatic(False)

        # Convert molecule to SMARTS
        smarts = Chem.MolToSmarts(mol)

        # Replace all bond symbols with wildcard bonds ~
        smarts_with_wildcard_bonds = smarts.replace('-', '~').replace('=', '~').replace(':', '~')
        return smarts_with_wildcard_bonds
    def repair_specical_structure(self):
        special_sub = {
            'porphyrin': r'C12=C/C(C=C/3)=NC3=C/C4=CC=C(/C=C5C=CC(/C=C(C=C/2)\N1)=N/5)N4',
            'Phthalocyanine': r'C12=N\C(C=C/3)=NC3=N/C4=CC=C(N4)/N=C(C=C5)\N=C5/N=C(N/2)/C=C1'
        }

        # Iterate through special substructures
        for name in special_sub:

            # Get SMILES of special substructure
            special_smiles = special_sub[name]
            
            # Generate SMARTS for special substructure (ignoring bond types)
            special_smarts = self.smiles_to_smarts_with_wildcard_bonds(special_sub[name])
            
            special_template = Chem.MolFromSmarts(special_smarts)
            if special_template is None:
                raise ValueError("repair_porphyrin invalid special substructure SMARTS!")

            # Match substructures in molecule
            matches = self.mol.GetSubstructMatches(special_template)
            if not matches:
                continue  # If no matching substructure found, skip this template
            self.error_found.add(name+'_error')
            # Because directly reading mol from SMILES will cause kekulization errors, first convert to SMARTS then to mol
            special_mol = Chem.MolFromSmiles(special_smiles)
            special_smarts = Chem.MolToSmarts(special_mol)
            special_mol = Chem.MolFromSmarts(special_smarts)
            # Fix all matching substructures
            for match in matches:
                self.update_exemption(match)
                # Iterate through bonds in template and update bonds in matching part of target molecule
                for bond in special_mol.GetBonds():
                    start_atom_idx = bond.GetBeginAtomIdx()
                    end_atom_idx = bond.GetEndAtomIdx()
                    bond_type = bond.GetBondType()

                    # Get corresponding atom indices in target molecule
                    mol_start_idx = match[start_atom_idx]
                    mol_end_idx = match[end_atom_idx]

                    # Get bond in molecule and modify its bond type
                    target_bond = self.mol_edit.GetBondBetweenAtoms(mol_start_idx, mol_end_idx)
                    if target_bond:  # Ensure bond exists
                        target_bond.SetBondType(bond_type)
                # Iterate through atoms in matching part and neutralize charged atoms
                for atom_idx in match:
                    atom = self.mol_edit.GetAtomWithIdx(atom_idx)
                    atom.SetIntProp('fixed_already', 1)
                    if atom.GetFormalCharge() != 0:  # If atom is charged
                        atom.SetFormalCharge(0)      # Set to neutral
                        atom.SetNumExplicitHs(atom.GetTotalNumHs())  # Add explicit hydrogen atoms
        if self.mol_edit is None:
            return self.smiles
        # Update molecule object
        self.mol = self.mol_edit.GetMol()
        # The following two lines help avoid can't Kekulize errors
        # try:
        #     smiles_tmp = Chem.MolToSmiles(self.mol)
        # except:
        #     smiles_tmp = Chem.MolToSmiles(self.mol,isomericSmiles=False)
        # self.mol = Chem.MolFromSmiles(smiles_tmp,sanitize=False)
        # Output repaired molecule SMILES
        #res_smiles = Chem.MolToSmiles(self.mol, isomericSmiles=True)
        #return res_smiles
    def single_haevy_atom(self):
        halogen_flag = False
        for atom in self.mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            expected_valence = self.typical_valences.get(atomic_num, "Unknown")
            total_valence = atom.GetTotalValence()
            formal_charge = atom.GetFormalCharge()
            
            if atom.GetSymbol() == 'O' and formal_charge == 1:
            # Find hydrogen atoms connected to oxygen
            # For [OH+] cases
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H':
                        atom.SetFormalCharge(0)
                        break
            # About charge correction valence
            if halogen_flag:
                expected_valence = [-1,1,3,5,7]
                
            elif expected_valence[0] <= formal_charge or atom.GetSymbol() in self.anion:
                expected_valence = [val + formal_charge for val in expected_valence]

            else:
                expected_valence = [val - formal_charge for val in expected_valence]
            if self.debug:
                print('atom:', atom.GetSymbol(), 'original_charge:', total_valence, 'expected_valence:', expected_valence,'atom.GetFormalCharge():',atom.GetFormalCharge())
            if total_valence in expected_valence and formal_charge != 0:
                    
                # Valence is correct, but charge is not 0
                label_set_elec_zero = False
                # Check if there are hydrogen atoms around this atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H':
                        label_set_elec_zero = True
                        break
                if label_set_elec_zero:
                    
                    atom.SetFormalCharge(0)
                    expected_valence = self.typical_valences.get(atomic_num, "Unknown")
                    total_valence = atom.GetTotalValence()
                    needed_hydrogens = self.find_closest_element(expected_valence, total_valence) - total_valence
                    self.atoms_to_modify.append((atom.GetIdx(), needed_hydrogens))
    def repair_CH_error(self):
        # Get all rings
        rings = rdmolops.GetSymmSSSR(self.mol_edit)

        # Iterate through each ring, check if CH count >= 2, contains only C,N atoms, and has no double bonds
        flag_aromaticity = [] # Mark whether this ring is a [CH] ring
        flag_turn_benzene = []
        # Check if it's a benzene ring, criteria: has 6 C atoms and 3 double bonds
        # Handle cyclic [CH] structures
        for i, ring in enumerate(rings):
            if self.need_kill:
                break
            ring_mols = [self.mol_edit.GetAtomWithIdx(atom_idx) for atom_idx in ring]
            ring_atoms = [mol.GetSymbol() for mol in ring_mols]
            len_ring = len(ring_atoms)
            if len_ring == 5:
                for atom in ring_mols:
                    atom.SetIntProp("in_5_ring", 1)  # Mark as five-membered ring

            all_C_N = all(atom.GetSymbol() in ['C', 'N'] for atom in ring_mols)
            C_count = sum(1 for atom in ring_mols if atom.GetSymbol() == 'C')
            ch_count = sum(1 for atom in ring_mols if self.is_CH(atom))
            
            is_full_occupy = False
            has_two_double_bond = False
            for atom in ring_mols:
                total_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                if self.has_two_double_bonds(atom):
                    has_two_double_bond = True
                    self.culmul_bond_flag = True
                    if len_ring in {3,4} or len_ring >= 7:

                        self.need_kill = True
                if total_valence >= 4:
                    is_full_occupy = True
                if has_two_double_bond and is_full_occupy:
                    break

            flag_Benzene = C_count == 6 and all_C_N and has_two_double_bond
            flag = (ch_count>=2 and (not self.ring_has_double_bond(ring, self.mol_edit)) and all_C_N and (not is_full_occupy))
            # Record rings that need modification
            flag_aromaticity = flag_aromaticity + [flag]
            flag_turn_benzene = flag_turn_benzene + [flag_Benzene]

        if any(flag_turn_benzene):
            # Record error type
            self.error_found.add('Benzene')
            for i, ring in enumerate(rings):
                if flag_turn_benzene[i]:
                    self.convert_to_benzene(ring, self.mol_edit)
        if any(flag_aromaticity):
            self.error_found.add('ring_CH_error')
            has_flag_CH_error = True   
            if self.debug:
                print('ring:',flag_aromaticity)
                for i, ring in enumerate(rings):
                    print('ring_atoms:',[self.mol_edit.GetAtomWithIdx(atom_idx).GetSymbol() for atom_idx in ring])
            for i, ring in enumerate(rings):
                
                if flag_aromaticity[i]:
                    HAC_ring = 0
                    self.update_exemption(ring)
                    atom_idx_list = list(ring)
                    for atom_idx in atom_idx_list:
                        
                        atom = self.mol_edit.GetAtomWithIdx(atom_idx)
                        # Set all atom charges in [CH] ring to 0
                        if atom.GetSymbol() == 'C':
                            atom.SetFormalCharge(0)
                        HAC_ring += (atom.GetSymbol() != 'H')
                        atom.SetIntProp('fixed_already', 1)
                    #print(Chem.MolToSmiles(mol))
                    if HAC_ring % 4 == 0:
                        for k in range(len(list(ring))):
                            atom_1 = self.mol_edit.GetAtomWithIdx(atom_idx_list[k])
                            if atom_1.GetSymbol() == 'C':
                                break
                        atom_1.SetFormalCharge(1)
                        atom_2 = self.mol_edit.GetAtomWithIdx(atom_idx_list[-1])
                        atom_2.SetFormalCharge(-1)
                    if HAC_ring % 4 == 1:
                        atom_1 = self.mol_edit.GetAtomWithIdx(atom_idx_list[-1])
                        atom_1.SetFormalCharge(-1)
                    if HAC_ring % 4 == 3:
                        for k in range(len(list(ring))):  
                            atom_1 = self.mol_edit.GetAtomWithIdx(atom_idx_list[k])
                        atom_1.SetFormalCharge(1)
                    if self.is_neutral_carbon(self.mol_edit.GetAtomWithIdx(atom_idx_list[0])):

                        atom_start_idx = 0
                        atom_end_idx = 1
                    else:

                        atom_start_idx = 1
                        atom_end_idx = 2
                    while atom_end_idx < len(atom_idx_list):
                        bond = self.mol_edit.GetBondBetweenAtoms(atom_idx_list[atom_start_idx], atom_idx_list[atom_end_idx])
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
                        atom_start_idx += 2
                        atom_end_idx += 2

        if self.debug:
            print('after ring:',Chem.MolToSmiles(self.mol_edit))
        # Process chain-type [CH] structures
        matches = self.find_matching_fragments(self.mol_edit)
        # Find connected unsaturated carbon atoms
        if len(matches) != 0: 
            self.error_found.add('chain_CH_error')
            for match in matches:
                unsaturated_carbons = []
                HAC_ring = 0
                for atom_idx in match:
                    atom = self.mol_edit.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        # Confirm if it's an unsaturated carbon
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetNumRadicalElectrons() == 1:
                            unsaturated_carbons.append(neighbor.GetIdx())


                atom_idx_list = sorted(list(set(unsaturated_carbons)))
                if len(atom_idx_list) < 3:
                    continue
                for atom_idx in atom_idx_list:
                    
                    atom = self.mol_edit.GetAtomWithIdx(atom_idx)
                    # Set all atom charges in [CH] chain to 0
                    if atom.GetSymbol() == 'C':
                        atom.SetFormalCharge(0)
                    HAC_ring += (atom.GetSymbol() != 'H')
                    atom.SetIntProp('fixed_already', 1)
                if HAC_ring % 4 == 0:
                    for k in range(len(list(atom_idx_list))):
                        atom_1 = self.mol_edit.GetAtomWithIdx(atom_idx_list[k])
                        if atom_1.GetSymbol() == 'C':
                            break
                    atom_1.SetFormalCharge(1)
                    atom_2 = self.mol_edit.GetAtomWithIdx(atom_idx_list[-1])
                    atom_2.SetFormalCharge(-1)
                if HAC_ring % 4 == 1:
                    atom_1 = self.mol_edit.GetAtomWithIdx(atom_idx_list[-1])
                    atom_1.SetFormalCharge(-1)
                if HAC_ring % 4 == 3:
                    for k in range(len(list(atom_idx_list))):  
                        atom_1 = self.mol_edit.GetAtomWithIdx(atom_idx_list[k])
                    atom_1.SetFormalCharge(1)
                if self.is_neutral_carbon(self.mol_edit.GetAtomWithIdx(atom_idx_list[0])):

                    atom_start_idx = 0
                    atom_end_idx = 1
                else:

                    atom_start_idx = 1
                    atom_end_idx = 2
                while atom_end_idx < len(atom_idx_list):
 
                    bond = self.mol_edit.GetBondBetweenAtoms(atom_idx_list[atom_start_idx], atom_idx_list[atom_end_idx])
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
                    atom_start_idx += 2
                    atom_end_idx += 2
        # Update molecule object
        self.mol = self.mol_edit.GetMol()

    def repair_charge_separation(self):
        # Record information of atoms that need modification
        
        for atom in self.mol_edit.GetAtoms():
            # Traverse heavy atoms in the molecule
            if atom.GetAtomicNum() != 1:
                atomic_num = atom.GetAtomicNum()
                if_trans_elec = False
                formal_charge = atom.GetFormalCharge()
                    
                if atom.GetSymbol() == 'O' and formal_charge == 1:
                # Find hydrogen atoms connected to oxygen
                # For [OH+] cases
                    for neighbor in atom.GetNeighbors():
                        
                        if neighbor.GetSymbol() == 'H':
                            atom.SetFormalCharge(0)
                            self.error_found.add('OH+_error')
                            break
                total_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                if (atom.GetSymbol() == 'N' or atom.GetSymbol() == 'P') and formal_charge == 1:
                    # Find hydrogen atoms connected to oxygen
                    # For [NH+] and [PH+] cases

                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'H':
                            self.error_found.add('NH+_and_PH+_error')
                            total_valence -= 2
                            break

                # Traverse neighboring atoms
                if formal_charge > 0:

                    for neighbor in atom.GetNeighbors():
                    
                        # If current atom is positively charged and neighboring atom is negatively charged, try to transfer charge

                        if neighbor.GetFormalCharge() < 0 and formal_charge + neighbor.GetFormalCharge() == 0 and neighbor.GetIntProp("charge_se") < 2:
                            if_trans_elec = True
                            # If C- and C+ are together, directly form a double bond and remove charges
                            if atom.GetSymbol() == 'C' and neighbor.GetSymbol() == 'C':
                                self.error_found.add('C-_C+_error')
                                atom.SetFormalCharge(0)
                                neighbor.SetFormalCharge(0)
                                bond = self.mol_edit.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                bond.SetBondType(Chem.BondType.DOUBLE)
                                break
                            self.error_found.add('charge_se_error')
                            atom.SetIntProp("charge_se",atom.GetIntProp("charge_se") + 1)
                            atom.SetIntProp('fixed_already', 1)
                            break
                    
                    
                    # If current atom is negatively charged and neighboring atom is positively charged, try to transfer charge
                if formal_charge < 0:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetFormalCharge() > 0 and formal_charge + neighbor.GetFormalCharge() == 0 and neighbor.GetIntProp("charge_se") < 2:
                            if_trans_elec = True
                            if atom.GetSymbol() == 'C' and neighbor.GetSymbol() == 'C':
                                self.error_found.add('C-_C+_error')
                                atom.SetFormalCharge(0)
                                neighbor.SetFormalCharge(0)
                                bond = self.mol_edit.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                bond.SetBondType(Chem.BondType.DOUBLE)
                                break                         
                            atom.SetIntProp("charge_se",atom.GetIntProp("charge_se") + 1)
                            self.error_found.add('charge_se_error')
                            atom.SetIntProp('fixed_already', 1)
                            break
        self.mol = self.mol_edit.GetMol()

    def repair_valance_error(self):
        for atom in self.mol_edit.GetAtoms():

            atomic_num = atom.GetAtomicNum()

            # If this element is in the dictionary, set its isotope to the most stable isotope
            if atomic_num in self.stable_isotopes:
                self.has_isotope = True
                stable_isotope = self.stable_isotopes[atomic_num]
                atom.SetIsotope(stable_isotope)

            if atom.GetAtomicNum() == 1:
                # Skip hydrogen atoms
                continue
            expected_valence = self.typical_valences.get(atomic_num, "Unknown")
            total_valence = math.floor(sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]))# Round down to prevent aromatic bonds from being read as 1.5 instead of 1
            formal_charge = atom.GetFormalCharge()
            
            # About charge-corrected valence
            if self.is_halogen(atom):
                expected_valence = [-1,1,3,5,7]
            elif atom.GetSymbol() == 'C' and atom.GetIntProp("charge_se") > 0:
                expected_valence = [4]
            elif expected_valence[0] <= formal_charge or atom.GetSymbol() in self.anion:
                # This handles cases with more charges, like Cl2-, where Cl valence should be -1
                total_valence += formal_charge
            else:
                total_valence -= formal_charge
            
            if self.debug:
                print('atom:', atom.GetSymbol(), 'original_charge:', total_valence, 'expected_valence:', expected_valence,'atom.GetFormalCharge():',atom.GetFormalCharge())
            if total_valence in expected_valence and formal_charge == 0:
                
                # Best case scenario
                continue
            elif atom.GetIntProp("charge_se") > 0:
                # Cases with charge transfer and already fixed cases
                continue
            elif total_valence in expected_valence and formal_charge != 0:
                if atom.GetIntProp("fixed_already") == 1:
                    continue
                # Valence is correct, but charge is not 0
                label_set_elec_zero = False
                # Check if there are hydrogen atoms around this atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H' and formal_charge > 0:
                        self.error_found.add('elec_error')
                        label_set_elec_zero = True
                        break
                    if formal_charge < 0:
                        label_set_elec_zero = True
                        self.error_found.add('elec_error')
                        break
                
                if label_set_elec_zero:
                    
                    atom.SetFormalCharge(0)
                    expected_valence = self.typical_valences.get(atomic_num, "Unknown")
                    total_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                    needed_hydrogens = self.find_closest_element(expected_valence, total_valence) - total_valence
                    self.atoms_to_modify.append((atom.GetIdx(), needed_hydrogens))

            else:

                self.error_found.add('elec_error')
                # If none of the above conditions are met, reset expected_valence
                expected_valence = self.typical_valences.get(atomic_num, "Unknown")
                total_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                if atom.GetIntProp('in_5_ring') == 1 and atom.GetSymbol() == 'C' and formal_charge == -1:
                    continue
                    
                if (atom.GetSymbol() == 'N' or atom.GetSymbol() == 'P') and all(bond.GetBondType() == Chem.BondType.AROMATIC for bond in atom.GetBonds()) and formal_charge == -1:
                    # Special handling for N atoms, directly add hydrogen at this step
                    total_valence -= 1
                    needed_hydrogens = self.find_closest_element(expected_valence, total_valence) - total_valence
                    self.atoms_to_modify.append((atom.GetIdx(), needed_hydrogens))
                    continue
                if total_valence != 0:
                    atom.SetFormalCharge(0)
                    needed_hydrogens = self.find_closest_element(expected_valence, total_valence) - total_valence
                    
                    self.atoms_to_modify.append((atom.GetIdx(), needed_hydrogens))
        for atom_idx, needed_hydrogens in self.atoms_to_modify:
            
            try:
                atom = self.mol_edit.GetAtomWithIdx(atom_idx)
            except:
                print('atom',atom.GetSymbol(),'fixed_smiles',Chem.MolToSmiles(self.mol_edit))
            atom.SetFormalCharge(0)  # Remove charge
            
            # Get indices of all hydrogen atoms connected to this atom
            bonds = atom.GetBonds()
            try:
                hydrogen_indices = [bond.GetOtherAtomIdx(atom_idx) for bond in bonds if bond.GetOtherAtom(atom).GetAtomicNum() == 1]
            except:
                need_kill = True
                break
                    
            needed_hydrogens = int(needed_hydrogens)
            # Add hydrogen atoms (if needed_hydrogens is positive)
            if needed_hydrogens > 0:
                for _ in range(needed_hydrogens):
                    temp_mol =  Chem.RWMol(self.mol_edit) 
                    hydrogen_idx = temp_mol.AddAtom(Chem.Atom(1))  # Add hydrogen atom
                    temp_mol.AddBond(atom_idx, hydrogen_idx, Chem.BondType.SINGLE)  # Add single bond
                    temp_mol.GetAtomWithIdx(hydrogen_idx).SetIntProp('fixed_already', 0)
                    if self.is_legitimate(temp_mol):
                        self.mol_edit = temp_mol
            # Remove hydrogen atoms (if needed_hydrogens is negative)
            elif needed_hydrogens < 0:
                # Calculate number of hydrogen atoms to remove
                hydrogens_to_remove = min(len(hydrogen_indices), -needed_hydrogens)
                # Remove excess hydrogen atoms
                for hydrogen_idx in hydrogen_indices[:hydrogens_to_remove]:
                    temp_mol =  Chem.RWMol(self.mol_edit) 
                    temp_mol.RemoveAtom(hydrogen_idx)
                    if self.is_legitimate(temp_mol):
                        self.mol_edit = temp_mol
        self.mol = self.mol_edit.GetMol()
    def repair_H_error(self):

        for atom_idx, needed_hydrogens in self.atoms_to_modify:
            try:
                atom = self.mol_edit.GetAtomWithIdx(atom_idx)
            except:
                print('atom',atom.GetSymbol(),'fixed_smiles',Chem.MolToSmiles(self.mol_edit))
            atom.SetFormalCharge(0)  # Remove charge
            
            # Get indices of all hydrogen atoms connected to this atom
            bonds = atom.GetBonds()
            try:
                hydrogen_indices = [bond.GetOtherAtomIdx(atom_idx) for bond in bonds if bond.GetOtherAtom(atom).GetAtomicNum() == 1]
            except:
                need_kill = True
                break
                    
            needed_hydrogens = int(needed_hydrogens)
            # Add hydrogen atoms (if needed_hydrogens is positive)
            if needed_hydrogens > 0:
                for _ in range(needed_hydrogens):
                    hydrogen_idx = self.mol_edit.AddAtom(Chem.Atom(1))  # Add hydrogen atom
                    self.mol_edit.AddBond(atom_idx, hydrogen_idx, Chem.BondType.SINGLE)  # Add single bond
            
            # Remove hydrogen atoms (if needed_hydrogens is negative)
            elif needed_hydrogens < 0:
                # Calculate number of hydrogen atoms to remove
                hydrogens_to_remove = min(len(hydrogen_indices), -needed_hydrogens)
                # Remove excess hydrogen atoms
                for hydrogen_idx in hydrogen_indices[:hydrogens_to_remove]:
                    self.mol_edit.RemoveAtom(hydrogen_idx)
        self.mol = self.mol_edit.GetMol()
    def match_substitute(self,smiles,special_template,culmul_map):
        """
        Perform template-based substructure matching and replacement.
        
        This method finds substructures matching a template pattern and
        replaces them with corrected structures from a mapping dictionary.
        It's used for fixing known problematic structural patterns.
        
        Args:
            smiles (str): SMILES pattern to search for
            special_template (rdkit.Chem.Mol): Template molecule for matching
            culmul_map (dict): Mapping of problematic to corrected structures
            
        Note:
            This method modifies the molecule in-place by updating bond types
            and atomic properties based on the template mapping.
        """
        matches = self.mol_edit.GetSubstructMatches(special_template)
        target_mol = Chem.MolFromSmiles(culmul_map[smiles])
        for match in matches:

            # Iterate through bonds in template and update bonds in matching part of target molecule
            for bond in target_mol.GetBonds():
                start_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()
                bond_type = bond.GetBondType()

                # Get corresponding atom indices in target molecule
                mol_start_idx = match[start_atom_idx]
                mol_end_idx = match[end_atom_idx]

                # Get bond in molecule and modify its bond type
                target_bond = self.mol_edit.GetBondBetweenAtoms(mol_start_idx, mol_end_idx)
                if target_bond:  # Ensure bond exists
                    target_bond.SetBondType(bond_type)
            # Iterate through atoms in matching part and neutralize charged atoms
            i = 0
            for atom_idx in match:
                atom = self.mol_edit.GetAtomWithIdx(atom_idx)
                atom.SetIntProp('fixed_already', 1)

                atom.SetFormalCharge(target_mol.GetAtoms()[i].GetFormalCharge())     
                atom.SetNumExplicitHs(atom.GetTotalNumHs())  # Add explicit hydrogen atoms
                i = i+1

    def repair_culmul_ring(self):
        for smiles in self.culmul_map:
            special_template = Chem.MolFromSmarts(self.get_culmul_smart(smiles))
            self.match_substitute(smiles,special_template,self.culmul_map)
        for smiles in self.culmul_map2:
            special_template = Chem.MolFromSmiles(smiles)
            self.match_substitute(smiles,special_template,self.culmul_map2)
        for atom in self.mol_edit.GetAtoms():
            if self.has_two_double_bonds(atom):
                self.need_kill = True
                break
    def draw_comparison(self,repaired_mol):
        """Generate comparison image"""
        original_mol = Chem.MolFromSmiles(self.smiles)
        try:
            img = Draw.MolsToGridImage(
                [original_mol, repaired_mol],
                molsPerRow=2,
                subImgSize=(400, 300),
                legends=['Original', 'Repaired']
            )
            return img
        except Exception as e:
            print(f"Drawing error: {str(e)}")
            return None
    def repair_smiles(self):
        """
        Execute comprehensive SMILES structure repair and validation pipeline.
        
        This is the main repair method that orchestrates all repair algorithms
        in the proper sequence. It handles molecule initialization, applies
        various repair strategies, and returns the corrected structure or
        appropriate error flags.
        
        Returns:
            tuple: (repaired_smiles, need_kill, readable) where:
                - repaired_smiles (str): Corrected SMILES string
                - need_kill (bool): True if molecule should be discarded
                - readable (bool): True if molecule was parseable
                
        Algorithm:
            1. Initial molecule parsing and validation
            2. Special structure repair (porphyrins, phthalocyanines)
            3. Ring structure analysis and CH error correction
            4. Cumulated bond system repair
            5. Charge separation correction
            6. Valence error repair and hydrogen adjustment
            7. Final canonicalization and validation
            
        Note:
            The method includes extensive error handling and can gracefully
            handle molecules that cannot be repaired by marking them for
            exclusion rather than causing pipeline failures.
        """

        self.exam_readable()

        if self.need_kill or not self.readable:
            return self.smiles,self.need_kill,self.readable
        # Initialize molecule object
        self.mol_ori = Chem.MolFromSmiles(self.smiles)
        if self.is_halogen_and_inorganic(self.mol_ori):
            return self.smiles,self.need_kill,self.readable
        self.mol = Chem.AddHs(self.mol_ori)
        self.mol_edit = Chem.RWMol(self.mol)  # Convert to editable molecule object
        self.HAC = self.mol.GetNumHeavyAtoms()
        if self.HAC == 1:
            self.need_kill = True
            return self.smiles,self.need_kill,self.readable
        for atom in self.mol_edit.GetAtoms():
            atom.SetIntProp("charge_se",0) # Create charge_separation property to mark charge separation cases, set "charge_se" property for each atom
            atom.SetIntProp("fixed_already", 0)  # Set boolean flag
            atom.SetIntProp("in_5_ring", 0)  # Whether in five-membered ring
        self.repair_specical_structure()

        self.repair_CH_error()

        if self.need_kill:
            return self.smiles,self.need_kill,self.readable
        if self.culmul_bond_flag:
            self.repair_culmul_ring()
            if self.need_kill:
                return self.smiles,self.need_kill,self.readable
        self.repair_charge_separation()
        self.repair_valance_error()
        if self.has_isotope:
            res_mol = self.mol
            res_smiles = Chem.MolToSmiles(res_mol,isomericSmiles=False)
            self.mol = Chem.MolFromSmiles(res_smiles,sanitize=False)
        res_mol = Chem.RemoveHs(self.mol, sanitize=False)
        # Canonicalization
        res_mol = Chem.MolFromSmiles(Chem.MolToSmiles(res_mol))
        if self.draw:
            img = self.draw_comparison(res_mol)
            img.show()
        
        if res_mol == None:
            self.need_kill = True
            return self.smiles,self.need_kill,self.readable
        return Chem.MolToSmiles(res_mol,isomericSmiles=False),self.need_kill,self.readable
if __name__ == '__main__':
    # Test case
    test_smiles = 'c1ccc(-c2cccc(N=c3ncc[n-]3)n2)nc1'
    
    test_repairer = MolecularFormulaRepairer(test_smiles,debug=False,draw=True)
    res_smiles,need_kill,readable = test_repairer.repair_smiles()
    print(res_smiles,need_kill,readable)