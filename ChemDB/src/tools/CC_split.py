#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chemical Complex Splitter and Analyzer Module

Official Research Code Repository for Chemical Coordination Complex Analysis

This module is part of the official codebase accompanying the research paper on 
automated analysis and decomposition of chemical coordination complexes. The code 
provides comprehensive functionality for identifying coordination compounds, 
extracting metal centers, ligands, and other molecular components from chemical 
structures represented as SMILES strings.

Research Focus:
- Automated chemical coordination complex decomposition and analysis
- Metal center identification and coordination environment characterization
- Ligand extraction and structural fingerprinting
- Database storage optimization for large-scale chemical datasets
- High-throughput processing of PubChem coordination compounds

Key Features:
- Chemical complex decomposition and analysis
- Metal center identification and coordination environment analysis
- Ligand extraction and characterization
- Molecular fingerprinting for unique identification
- Error handling for invalid molecular structures

Academic Citation:
If you use this code in your research, please cite the corresponding paper:
[Paper citation to be added upon publication]

Author: ChemDB Research Team
Institution: Shanghai Jiao Tong University
Created: 2025
Version: 1.0

License: Academic Open Source License

Copyright (c) 2025 ChemDB Research Project

This code is released under an open source license for academic and research purposes.
Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction for academic and research purposes, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software for non-commercial purposes, subject to the following conditions:

1. The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.
2. Any publications or research that uses this code must cite the corresponding paper.
3. Commercial use requires separate licensing agreement.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For questions regarding this research or licensing, please contact:
jack111@sjtu.edu.cn
"""

import re
import os
import sys
import csv
import hashlib
from collections import deque
from rdkit import Chem
from rdkit.Chem import AllChem      
from rdkit.Chem import Draw
from rdkit.Chem import GetPeriodicTable
from utils import load_elements_list
from tools.DID_calculate import calculate_canonical_did


class ChemicalComplex:
    """
    A comprehensive class for analyzing and decomposing chemical coordination complexes.
    
    This class provides functionality to parse, analyze, and extract structural information
    from chemical coordination complexes represented as SMILES strings. It identifies
    metal centers, ligands, coordination environments, and generates unique molecular
    fingerprints for database storage.
    
    Key Features:
    - Molecular structure parsing and validation
    - Metal center identification and analysis
    - Ligand extraction and characterization
    - Coordination environment determination
    - Molecular fingerprinting for unique identification
    - Error handling for invalid structures
    
    Attributes:
        cid (str): PubChem Compound ID
        smiles (str): SMILES representation of the chemical complex
        corcomp (bool): Flag indicating if this is a coordination compound
        did (str): Unique molecular fingerprint identifier
        metal_element_list (list): List of metal elements for analysis
        inactive (int): Flag for inactive compounds
    """
    
    def __init__(self, pubchem_info_dict, metal_element_list, error_file="error.out"):
        """
        Initialize a ChemicalComplex instance.
        
        Args:
            pubchem_info_dict (dict): Dictionary containing PubChem compound information
                                     with keys: 'cid', 'smiles', 'corcomp'
            metal_element_list (list): List of metal element symbols for analysis
            error_file (str): Path to error log file (default: "error.out")
        """
        self.cid = pubchem_info_dict['cid']
        # print(self.cid)
        self.smiles = pubchem_info_dict['smiles']
        #self.gibberish = pubchem_info_dict['gibberish']
        self.corcomp = pubchem_info_dict['corcomp']
        # mol = Chem.MolFromSmiles(smiles)
        # self.smiles = Chem.MolToSmiles(mol,canonical = True)
        self.inactive = 0
        #self.illegitimate = 0
        #self.metalcompound = 0
        #self.questionable = 0
        self.did = self.calculate_DID_from_fingerprint(self.smiles)

        # print(self.did)
        self.metal_element_list = metal_element_list
        self.error_file = error_file
        self.complex_info = {
            "complex_info": {
                "DID": self.did,  
                "complex_smiles": self.smiles,
                #"gibberish": self.gibberish,
                "corcomp": self.corcomp,
                "inactive": 0
            },
            "central_metal_info": [],
            "neighbor_info": []
        }
        if self.did is None:
            print('SB SMILES')
            return
        
        #print(self.questionable)
        #self.did = self.generate_did_from_cid(self.cid)  # Generate DID


        self.metal_info_list = self.extract_metal_content()
        self.metal_did_list = [self.generate_metal_did(symbol) for symbol in self.metal_info_list]

        self.process_metal_structure()
        self.process_complex()
        # Questionable boundary conditions can be changed by modifying parameters
        self.questionable_detection()
        self.is_metal_only(self.smiles)
        self.complex_info["complex_info"]["inactive"] = self.inactive
        # print(aaa)
        #self.merge_neighbor_info()
        
    # @staticmethod
    # def generate_did_from_cid(cid):
    #     """
    #     Generate corresponding DID based on given CID.

    #     Args:
    #         cid (str): Compound CID.

    #     Returns:
    #         str: Generated DID.
    #     """
    #     return 'D' + cid

        
    def generate_metal_did(self, symbol):
        """
        Generate a unique metal identifier (DID) from element symbol.
        
        This method creates a standardized identifier for metal elements
        using their atomic number. The format is 'D{atomic_number}E'.
        
        Args:
            symbol (str): Chemical element symbol (e.g., 'Fe', 'Cu', 'Zn')
            
        Returns:
            str: Metal DID in format 'D{atomic_number}E' or None if invalid
            
        Example:
            generate_metal_did('Fe') returns 'D26E'
        """
        if not symbol:
            print("Element symbol is not provided.")
            return None
        
        # Create PeriodicTable instance
        pt = Chem.GetPeriodicTable()

        # Use PeriodicTable instance to get atomic number corresponding to element symbol
        atomic_number = pt.GetAtomicNumber(symbol)

        # Check if atomic number is valid
        if atomic_number is None:
            print("Element symbol is invalid.")
            return None

        # Add 'D' before atomic number to generate metal_did
        metal_did = f"D{atomic_number}E"
        return metal_did

    
    def mark_smiles(self):
        """
        Mark atoms in the input molecule with their indices for tracking.
        
        This method creates an RDKit molecule object from the SMILES string and
        adds index annotations to each atom. This is essential for tracking
        atom positions during coordination complex analysis and ligand extraction.
        
        Returns:
            rdkit.Mol: RDKit molecule object with atom indices marked as properties,
                      or an empty molecule if SMILES parsing fails
                      
        Error Handling:
            - Logs parsing errors to the specified error file
            - Returns empty molecule object to prevent pipeline crashes
            - Continues execution even with invalid SMILES input
        """
        input_smiles = self.smiles
        input_mol = Chem.MolFromSmiles(input_smiles)
        if input_mol is None:
            with open(self.error_file, "a") as f:
                f.write(f"Error creating molecule from SMILES '{input_smiles}'\n")
            # Return a default empty molecule object so code can continue execution
            return Chem.MolFromSmiles('')
        for atom in input_mol.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIdx()))
        return input_mol


    def extract_metal_content(self):
        """
        Extract metal elements from the SMILES string using pattern matching.
        
        This method identifies metal centers in coordination complexes by parsing
        bracketed atomic symbols in SMILES notation. It handles both simple and
        complex metal specifications, including charged states and multiple metals.
        
        Algorithm:
        1. Use regex to find all bracketed atomic symbols [...]
        2. Extract alphabetic characters representing element symbols
        3. Split compound symbols into individual elements using uppercase delimiters
        4. Match against the provided metal element list
        5. Return list of identified metal symbols
        
        Returns:
            list: List of metal element symbols found in the structure (e.g., ['Fe', 'Cu'])
            
        Note:
            - Handles bimetallic and multimetallic complexes
            - Breaks on first metal match per bracket to avoid duplicates
            - Critical function for coordination complex identification
            
        Example:
            For SMILES '[Fe+2].CCO.[Zn]' returns ['Fe', 'Zn']
        """
        # Regular expression to match bracketed content in SMILES
        pattern = re.compile(r'\[([^\]]+)\]')
        matches = pattern.findall(self.smiles)
        metal_matches = []
        
        for match in matches:
            # Remove non-alphabetic characters, keep only letter parts
            clean_match = ''.join([char for char in match if char.isalpha()])
            
            # Split into individual element symbols using uppercase letter boundaries
            metal_symbols = re.findall('[A-Z][a-z]*', clean_match)
            
            # Check if each extracted symbol is in the metal element list
            for metal_symbol in metal_symbols:
                if metal_symbol in self.metal_element_list:
                    metal_matches.append(metal_symbol)
                    # Break loop once a metal element is matched to avoid duplicates
                    break
        
        return metal_matches
    
    def questionable_detection(self, ligand_num_boundary=11):
        """
        Detect potentially questionable coordination complexes based on ligand count.
        
        This method identifies coordination complexes that may be questionable for
        analysis based on the number of molecular fragments (ligands) in the SMILES
        string. Complexes with excessive numbers of fragments may indicate data
        quality issues or may not represent genuine coordination compounds.
        
        Args:
            ligand_num_boundary (int): Maximum acceptable number of ligands before
                                     marking as questionable (default: 11)
                                     
        Side Effects:
            Sets self.inactive to 3 if the ligand count exceeds the boundary
            
        Detection Logic:
            - Counts '.' separators in SMILES string as ligand boundaries
            - Compares count against configurable threshold
            - Marks complex as inactive (status 3) if threshold exceeded
            
        Note:
            The boundary parameter can be adjusted based on research requirements
            and the specific characteristics of the dataset being analyzed.
        """
        cc_smiles = self.smiles
        ligand_num = cc_smiles.count('.')
        if ligand_num > ligand_num_boundary:
            # Mark as questionable due to excessive ligand count
            self.inactive = 3

    def calculate_DID_from_fingerprint(self, smiles):
        """
        Calculate a canonical Database Identifier (DID) from SMILES string.
        
        This method generates a unique molecular identifier using advanced
        cryptographic fingerprinting algorithms. It delegates to the external
        calculate_canonical_did function for the actual computation.
        
        Args:
            smiles (str): SMILES string representation of the molecule
            
        Returns:
            str: Canonical database identifier (DID) for the molecule
            
        Note:
            This method serves as a wrapper for the external DID calculation
            function, providing consistent interface within the ChemicalComplex class.
        """
        result = calculate_canonical_did(smiles)
        return result

    # def calculate_DID_from_fingerprint(self, smiles):
    #     try:
    #         # Use RDKit to validate SMILES string
    #         mol = Chem.MolFromSmiles(smiles)
    #         if mol is None:
    #             print(f"Invalid SMILES: '{smiles}'")
    #             self.illegitimate = 1  # Set self.illegitimate to 1
    #             return None
            
    #         # Get canonical SMILES
    #         canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    #         # Generate molecular fingerprint, parameters need modification
    #         canonical_mol = Chem.MolFromSmiles(canonical_smiles)
    #         fingerprint = AllChem.GetMorganFingerprintAsBitVect(canonical_mol, radius=2, nBits=2048)
    #         fingerprint_bytes = fingerprint.ToBitString().encode()
    #         # Use hashlib to generate hash value

    #         hash_value = hashlib.sha256(fingerprint_bytes).hexdigest()
    #         hash_decimal = int(hash_value, 16)
    #         truncated_hash_decimal = hash_decimal % (10 ** 10)
    #         truncated_hash_decimal_str = 'D' + str(truncated_hash_decimal)
            
    #         return truncated_hash_decimal_str
    #     except Exception as e:
    #         print(f"Error calculating DID for SMILES '{smiles}': {e}")
    #         self.illegitimate = 1  # Set self.illegitimate to 1
    #         return None

    def generate_ligand_smiles(self, coordinate_atom_index, central_metal_index, original_mol, central_metals_set):
        """
        Generate SMILES for ligand structure connected to metal center
        
        Args:
            coordinate_atom_index: Index of coordinating atom
            central_metal_index: Index of central metal atom
            original_mol: Original molecule object
            central_metals_set: Set of central metal atom indices
            
        Returns:
            tuple: (ligand_smiles, coordinating_atom_new_index, ligand_did)
        """
        bonds = original_mol.GetBonds()
        connected_atoms = set()
        visited = set()
        queue = deque([coordinate_atom_index])
        while queue:
            current_atom_idx = queue.popleft()
            visited.add(current_atom_idx)
            if current_atom_idx == central_metal_index:
                continue
            connected_atoms.add(current_atom_idx)
            for bond in bonds:
                if current_atom_idx == bond.GetBeginAtomIdx():
                    neighbor_idx = bond.GetEndAtomIdx()
                elif current_atom_idx == bond.GetEndAtomIdx():
                    neighbor_idx = bond.GetBeginAtomIdx()
                else:
                    continue
                if neighbor_idx in central_metals_set:
                    continue
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)
        ligand_index_list = list(connected_atoms)
        new_mol = self.get_mol_from_indices(original_mol, ligand_index_list)
        coordinating_atom_new_index = ligand_index_list.index(coordinate_atom_index)
        ligand_smiles = Chem.MolToSmiles(new_mol, canonical=True)
        did = self.calculate_DID_from_fingerprint(ligand_smiles)
        return ligand_smiles, coordinating_atom_new_index, did

    def find_metals(self, smiles):
        """
        Find metal elements in SMILES structure
        
        Args:
            smiles: SMILES string to analyze
            
        Returns:
            list: List of metal symbols found in the structure
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Invalid SMILES input!")
            return
        metals_found = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in self.metal_element_list:
                metals_found.append(symbol)
        return metals_found

    def get_mol_from_indices(self, original_mol, indices):
        """
        Create new molecule from specified atom indices
        
        Args:
            original_mol: Original RDKit molecule object
            indices: List of atom indices to include in new molecule
            
        Returns:
            RDKit molecule object containing only specified atoms and their bonds
        """
        new_mol = Chem.RWMol()
        atom_mapping = {}
        for idx in indices:
            atom = original_mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)
            atom_mapping[idx] = new_idx
        for bond in original_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in indices and end_idx in indices:
                new_mol.AddBond(atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType())
        new_mol.UpdatePropertyCache()
        return new_mol.GetMol()

    def get_neighbors_info(self, central_metals):
        """
        Get information about atoms neighboring central metals
        
        Args:
            central_metals: List of central metal substructures
            
        This method analyzes the coordination environment around metal centers,
        identifying ligands and their binding modes.
        """
        marked_rdkit_mol = self.mark_smiles()
        bonds = marked_rdkit_mol.GetBonds()
        metal_atom_indices = set()
        for i in range(len(central_metals)):
            central_metal = central_metals[i]
            atom_index_tuples = marked_rdkit_mol.GetSubstructMatches(central_metal)
            metal_atom_indices.update([atom_tuple[0] for atom_tuple in atom_index_tuples])
        #print(metal_atom_indices)
        for atom_index in metal_atom_indices:
            central_metal_elements = marked_rdkit_mol.GetAtomWithIdx(atom_index).GetSymbol()
            central_metal_elements_did = self.generate_metal_did(central_metal_elements)
            #print(central_metal_elements)
            for bond in bonds:
                if bond.GetBeginAtomIdx() == atom_index:
                    coordinating_atom_index = bond.GetEndAtomIdx()
                elif bond.GetEndAtomIdx() == atom_index:
                    coordinating_atom_index = bond.GetBeginAtomIdx()
                else:
                    continue
                ligand_smiles, coordinating_atom_new_index, ligand_did = self.generate_ligand_smiles(
                    coordinating_atom_index, atom_index, marked_rdkit_mol, metal_atom_indices)
                bond_type = bond.GetBondType()
                self.complex_info["neighbor_info"].append({
                    "coordinating_atom_index": coordinating_atom_new_index,
                    "central_atom": central_metal_elements,
                    "central_atom_did": central_metal_elements_did,
                    "bond_type": bond_type,
                    "ligand_smiles": ligand_smiles,
                    "ligand_did": ligand_did
                })

    def is_metal_only(self, smiles):
        """
        Detect if structure contains only metal atoms (intermetallic compound)
        
        Args:
            smiles: SMILES string to analyze
            
        Sets self.inactive to 2 when all atoms are metals, 1 for invalid SMILES
        """
        try:
            # Parse SMILES string to generate molecule object
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Invalid SMILES: '{smiles}'")
                self.inactive = 1  # Set self.inactive to 1
                return 
            
            # Get all atoms in the molecule
            atoms = mol.GetAtoms()

            # Check if all atoms are metal elements
            for atom in atoms:
                symbol = atom.GetSymbol()
                if symbol not in self.metal_element_list:
                    return  # If there are non-metal elements, exit directly without any setting

            # If all atoms are metal elements
            self.inactive = 2

        except Exception as e:
            print(f"Error processing SMILES '{smiles}': {e}")

    
    def calculate_valence(self, bonds, metal_atom_index):
        """
        Calculate valence state of metal atom (needs modification, currently not very useful)
        
        Args:
            bonds: List of bonds in the molecule
            metal_atom_index: Index of the metal atom
            
        Returns:
            float: Total bond order around the metal atom
        """
        total_bond_state = 0
        for bond in bonds:
            atom1_index = bond.GetBeginAtom().GetIdx()
            atom2_index = bond.GetEndAtom().GetIdx()
            if atom1_index == metal_atom_index or atom2_index == metal_atom_index:
                bond_order = bond.GetBondTypeAsDouble()
                total_bond_state += bond_order
        return total_bond_state


    def get_metal_info(self, central_metals):
        """
        Get detailed information about metal centers in the complex
        
        Args:
            central_metals: List of central metal substructures
            
        This function needs comments and optimization.
        Extracts metal atom information including position, symbol, and valence state.
        """
        marked_rdkit_mol = self.mark_smiles()
        bonds = marked_rdkit_mol.GetBonds()
        metal_atom_indices = set()
        #print(metal_atom_indices)
        for central_metal in central_metals:
            atom_index_tuples = marked_rdkit_mol.GetSubstructMatches(central_metal)
            metal_atom_indices.update([atom_tuple[0] for atom_tuple in atom_index_tuples])
        #print(metal_atom_indices)
        for atom_index in metal_atom_indices:
            #print(atom_index)
            #central_metal_elements = self.find_metals(Chem.MolToSmiles(central_metals[0]))
            central_metal_elements = marked_rdkit_mol.GetAtomWithIdx(atom_index).GetSymbol()
            #atom_number = Chem.PeriodicTable.GetAtomicNumber(central_metal_elements)
            # print(type(central_metal_elements))
            # print(central_metal_elements)
            #print(atom_number)
            #print(type(central_metal_elements))
            metal_did = self.generate_metal_did(central_metal_elements)
            central_metal_info = {
                "atom_index": atom_index,
                "central_metal": central_metal_elements,
                "valence": 0,
                "metal_did": metal_did
            }
            total_bond_state = self.calculate_valence(bonds, atom_index)
            central_metal_info["valence"] = total_bond_state
            self.complex_info["central_metal_info"].append(central_metal_info)

    
    def process_metal_structure(self):
        """
        Process and analyze metal center structures in the coordination complex.
        
        This method initiates the analysis of metal centers by creating SMARTS
        patterns for each identified metal and extracting detailed information
        about their coordination environments.
        
        Process:
        1. Retrieves metal symbols from the extracted metal information
        2. Creates SMARTS patterns for metal center identification
        3. Delegates to get_metal_info for detailed metal analysis
        
        Note:
            This method only processes complexes that contain identified metals.
            Empty metal lists will result in no processing.
        """
        metal_symbols = self.metal_info_list
        if metal_symbols:
            metal_substructures = [Chem.MolFromSmarts('[' + symbol + ']') for symbol in metal_symbols]
            self.get_metal_info(metal_substructures)

    def process_covalent_structure(self):
        """
        Process covalent coordination structures and ligand environments.
        
        This method analyzes covalent bonds between metal centers and their
        coordinating ligands. It creates SMARTS patterns for metal identification
        and extracts detailed information about the coordination environment.
        
        Process:
        1. Creates SMARTS patterns for each identified metal
        2. Analyzes coordination bonds and ligand structures
        3. Extracts neighbor information for database storage
        
        Note:
            This method focuses on direct covalent coordination bonds
            between metals and ligands.
        """
        metal_symbols = self.metal_info_list
        if metal_symbols:
            metal_substructures = [Chem.MolFromSmarts('[' + symbol + ']') for symbol in metal_symbols]
            self.get_neighbors_info(metal_substructures)

    def process_ion_structure(self):
        """
        Process ionic coordination structures and separated ligand components.
        
        This method handles coordination complexes where metals and ligands
        are present as separate ionic components rather than covalently bonded.
        It identifies metal-containing and ligand-containing substructures
        and processes them accordingly.
        
        Process:
        1. Splits SMILES into individual molecular components
        2. Identifies which components contain metals
        3. Processes metal-containing parts as covalent structures
        4. Treats non-metal components as ionic ligands
        5. Creates ionic ligand information with appropriate metadata
        
        Note:
            This method is essential for handling coordination complexes
            with ionic rather than covalent metal-ligand interactions.
        """
        metal_symbols = self.metal_info_list
        metal_symbols_did = self.metal_did_list
        
        if metal_symbols:
            # Split SMILES into individual molecular components
            substructure_pattern = re.compile(r'[^\.]+')
            substructures = substructure_pattern.findall(self.smiles)
            
            for substructure in substructures:
                if any(metal in substructure for metal in metal_symbols):
                    # Process metal-containing components as covalent structures
                    self.process_covalent_structure()
                else:
                    # Process non-metal components as ionic ligands
                    did = self.calculate_DID_from_fingerprint(substructure)
                    ligand_info = {
                        'coordinating_atom_index': None, 
                        'ligand_smiles': substructure, 
                        'bond_type': 'IONIC',
                        'central_atom': metal_symbols, 
                        'central_atom_did': metal_symbols_did, 
                        'ligand_did': did
                    }
                    self.complex_info['neighbor_info'].append(ligand_info)

    def merge_neighbor_info(self):
        merged_info = self.complex_info.copy()  # Copy input dictionary to avoid modifying original data
        merged_neighbor_info_list = []
        seen = set()  # Track already seen ligand information
        try:
            for ligand_info in self.complex_info["neighbor_info"]:
                # Handle cases where central_atom_did might be a list
                central_atom_did = ligand_info.get("central_atom_did", "")
                if isinstance(central_atom_did, list):
                    central_atom_did = ",".join(map(str, central_atom_did))  # Convert list to comma-separated string
                tuple_key = (ligand_info.get("coordinating_atom_index", ""), central_atom_did, ligand_info.get("ligand_did", ""))
                
                # If ligand information is the same, don't add duplicates
                if tuple_key not in seen:
                    merged_neighbor_info_list.append(ligand_info)
                    seen.add(tuple_key)
                else:
                    # If ligand information already exists, merge dictionaries
                    for existing_ligand_info in merged_neighbor_info_list:
                        if (existing_ligand_info.get("coordinating_atom_index") == ligand_info.get("coordinating_atom_index") and
                            existing_ligand_info.get("central_atom_did") == central_atom_did and
                            existing_ligand_info.get("ligand_did") == ligand_info.get("ligand_did")):
                            # Merge dictionaries
                            existing_ligand_info.update(ligand_info)
                            break
            merged_info["neighbor_info"] = merged_neighbor_info_list
        except Exception as e:
            print(f"Error processing SMILES: {e}")
        return merged_info


    def is_multi_molecule(self):
        return '.' in self.smiles


    def process_complex(self):
        if self.is_multi_molecule():
            return self.process_ion_structure()
        else:
            return self.process_covalent_structure()


# Main function call also needs to be updated
# Note: There's still an unresolved issue where ligands might be inserted into the dictionary repeatedly
if __name__ == "__main__":
    metal_path = '/Users/weilinzou/Code/DatabaseProject/ChemDB/Resources/metal_list.txt'
    metal_list = load_elements_list(metal_path)
    # print(metal_list)
    # print(len(metal_list))
    test_smiles = 'CC[Sn+2]CC.C(=O)[O-].C(=O)[O-].[Sn]'
    test_smiles2 = 'CCCC[Sn](CCCC)CCC(C)COC(=C)C(=O)OC'
    test_smiles3 = '[2H][C-]([2H])[2H].[2H]C([2H])([2H])C(C1=CC=NC=C1)(C2=NN(C=C2)CC(F)F)O.[Li+].CC(C1=CC=NC=C1)C2=NN(C=C2)CC(F)F.C1=CN=CC=C1C(C2=NN(C=C2)CC(F)F)O.C1=CN=CC=C1C(=O)C2=NN(C=C2)CC(F)F.C1=CN=CC=C1I.C1=CN(N=C1C=O)CC(F)F.C1=C(NN=C1)C=O.C(C(F)F)OS(=O)(=O)C(F)(F)F.O=[Mn]=O.[Zn]'
    test_smiles4 = 'C1=CC=C2C(=C1)C=CC=C2O.C1=CC=C2C(=C1)C=CC=C2O.[Mn]'
    cid_file_name = "/home/weilin/ScientificProject/scientific/PubChemData/processed_data/Mn.csv"  # Please replace with your CID file name
    # Create ChemicalComplex instance
    print('ok')
    test_info_dict = {'cid': '111111', 'smiles': 'CC1=C2CC(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(CC1OC(=O)[C@@H]([C@@H](C)C=C(C)C)O)O)OC(=O)C)(CO4)O)O)C.[Ac]','corcomp':True}
    complex_instance = ChemicalComplex(test_info_dict, metal_list)
    #print(complex_instance.complex_info)
    merged_dict = complex_instance.merge_neighbor_info()
    print(merged_dict)
    print(complex_instance.inactive)
    # Use add_did_info function to process DID information
    #complex_instance_with_did = complex_instance.add_did_info(cid_file_name)
    # print(complex_instance_with_did)

