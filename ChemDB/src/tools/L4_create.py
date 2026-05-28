#!/usr/bin/env python3
"""
ChemDB: Chemical Database Analysis Pipeline
==========================================

L4 Molecular Fragment Creator Module

This module provides comprehensive functionality for creating Level-4 (L4) molecular fragments
from marked ligand data. It enables advanced molecular structure analysis by decomposing
complex molecular structures into smaller, analyzable fragments based on GA (Guest Agent)
and IRL (In-situ Reactive Ligand) marking patterns.

This is part of the official code release for the ChemDB research project.

Academic Citation:
    If you use this code in your research, please cite:
    [Citation information to be added upon publication]

Author: ChemDB Research Team
Institution: Shanghai Jiao Tong University
Contact: jack111@sjtu.edu.cn

License:
    This code is released under an academic open source license for research and non-commercial use.
    Commercial use requires a separate licensing agreement.

Copyright (c) 2025 ChemDB Research Team. All rights reserved.

Module Features:
    - L4 fragment generation from marked ligand structures
    - IRL and GA pattern-based molecular decomposition
    - Fragment visualization and structural analysis
    - Progress tracking for batch processing
    - Molecular fingerprinting for unique identification

Research Context:
    Level-4 fragments represent chemically meaningful substructures that preserve
    important coordination chemistry information while enabling efficient database
    storage and retrieval for large-scale chemical datasets.

Module Classes:
    MoleculeSplitter: Core class for molecular fragmentation and analysis

Dependencies:
    - RDKit: Chemical informatics toolkit
    - matplotlib: Visualization capabilities
    - PIL: Image processing
    - tqdm: Progress tracking
"""

import io
import csv
import sys
import random
import os
import json
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from rdkit import Chem
from PIL import Image
from rdkit import Chem
from collections import deque
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from collections import deque
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from tools.DID_calculate import calculate_canonical_did

class MoleculeSplitter:
    """
    Advanced molecular splitter for creating Level-4 (L4) fragments from marked ligand data.
    
    This class provides sophisticated functionality for decomposing complex molecular structures
    into chemically meaningful fragments based on IRL (In-situ Reactive Ligand) and GA (Guest Agent)
    marking patterns. The resulting L4 fragments preserve important coordination chemistry
    information while enabling efficient database storage and analysis.
    
    Attributes:
        smiles (str): SMILES representation of the input molecule
        ga_data (list): List of Guest Agent reference data
        marked_ligand_data (dict): Atom-level marking information with GA and IRL annotations
        irl_list (list): List of IRL reference structures
        mol (rdkit.Chem.Mol): RDKit molecule object
        marked_molecule_data (dict): Initialized atom status and marking data
        IRL_in_ligand (list): List of all IRL DIDs found in the ligand
        source_did (str): Source molecule DID for database traceability
    """
    
    def __init__(self, smiles, marked_ligand_data, irl_list, ga_list=None, source_did=None):
        """
        Initialize MoleculeSplitter for L4 fragment generation.
        
        Args:
            smiles (str): SMILES string of the target molecule
            marked_ligand_data (dict): Atom-level marking information containing GA and IRL annotations
            irl_list (list): List of IRL reference structures for pattern matching
            ga_list (list, optional): List of Guest Agent reference structures
            source_did (str, optional): Source molecule DID for traceability
        """
        self.smiles = smiles
        self.ga_data = ga_list or []
        self.marked_ligand_data = marked_ligand_data
        self.irl_list = irl_list
        self.source_did = source_did
        
        # Initialize RDKit molecule
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Initialize atom status and marking data
        self.marked_molecule_data = self.initialize_atom_status()
        
        # Get all IRL DIDs found in the ligand
        self.IRL_in_ligand = self.get_all_IRL_DIDs()

    def initialize_atom_status(self):
        """
        Initialize marking data for each atom in the molecule with default active status.
        
        Creates a comprehensive dictionary containing atom-level information including
        IRL IDs, GA IDs, atomic symbols, activity status, and atom indices. All atoms
        are initially marked as active for fragment generation.

        Returns:
            dict: Dictionary mapping atom identifiers to their marking status and properties
                 Keys: 'atom_{index}' format
                 Values: Dict with 'IRL_ids', 'GA_ids', 'atom_symbol', 'active', 'atom_idx'
        """
        atom_status = {}
        # Iterate through each atom in the molecule, initialize atom marking status
        for atom_idx in range(self.mol.GetNumAtoms()):
            atom_status[f'atom_{atom_idx}'] = {
                'IRL_ids': [],  # Initialize IRL_ids as empty list
                'GA_ids': [],   # Initialize GA_ids as empty list
                'atom_symbol': self.mol.GetAtomWithIdx(atom_idx).GetSymbol(),  # Get atomic symbol
                'active': True,   # All atoms are active by default
                'atom_idx': atom_idx  # Store atom index for reference
            }

        # Populate each atom's IRL_ids and GA_ids based on marked_ligand_data
        for idx, atom in enumerate(self.mol.GetAtoms()):
            atom_idx = f'atom_{idx}'  # Generate atom index identifier
            if atom_idx in self.marked_ligand_data:
                atom_data = self.marked_ligand_data[atom_idx]
                atom_status[atom_idx]['IRL_ids'] = atom_data.get('IRL_ids', [])  # Populate IRL_ids
                atom_status[atom_idx]['GA_ids'] = atom_data.get('ga_ids', [])  # Populate GA_ids
                atom_status[atom_idx]['atom_symbol'] = atom_data.get('atom_symbol', '')  # Update atomic symbol

        return atom_status

    def get_all_IRL_DIDs(self):
        """
        Extract all IRL DIDs present in the ligand molecule.
        
        Scans through all marked atom data to identify unique IRL identifiers
        that are present in the current ligand structure.
        
        Returns:
            list: Unique list containing all IRL DIDs found in the ligand
        """
        irl_dids = set()  # Use set to avoid duplicate DIDs
        marked_ligand_data = self.marked_ligand_data  # Get marked atom data

        # Iterate through marked atom data to extract DIDs from IRL_ids
        for atom_data in marked_ligand_data.values():
            irl_ids = atom_data.get('IRL_ids', [])  # Get current atom's IRL_ids
            for irl_id in irl_ids:  # Add each ID from IRL_ids to irl_dids set
                irl_dids.add(irl_id)

        return list(irl_dids)  # Return list of IRL DIDs

    def bfs_find_substructure(self, start_atom_indices):
        """
        Use BFS algorithm to find all first-order neighbors of a group of atoms,
        including atoms in their GA groups or regular atoms.
        
        This method performs breadth-first search starting from IRL atoms to identify
        the complete molecular fragment that should be grouped together, considering
        GA boundaries and molecular connectivity.

        Args:
            start_atom_indices (tuple): Set of starting atom indices representing an IRL
            
        Returns:
            set: Set containing all relevant atom indices representing the complete molecular fragment
        """
        from collections import deque

        fragment_atoms = set()  # Store all atom indices in the fragment
        visited = set()         # Record visited atoms
        queue = deque(start_atom_indices)  # Initialize queue with starting atoms

        while queue:
            current_idx = queue.popleft()  # Get atom index from queue
            if current_idx in visited:     # Skip if already visited
                continue
            visited.add(current_idx)       # Mark as visited
            fragment_atoms.add(current_idx)  # Add current atom to fragment set

            # Get neighbors of current atom
            if isinstance(current_idx, str):  # Check if index is string type
                atom_idx = int(current_idx.split('_')[1])  # Extract atom number
            else:
                atom_idx = current_idx  # If index is already integer, use it directly

            current_atom = self.mol.GetAtomWithIdx(atom_idx)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = f'atom_{neighbor.GetIdx()}'

                # Get neighbor atom information
                neighbor_data = self.marked_molecule_data.get(neighbor_idx, {})
                neighbor_ga_ids = neighbor_data.get('GA_ids', [])
                
                if neighbor_ga_ids:  # If neighbor belongs to GA
                    # Find all atoms of the entire GA group
                    ga_id = neighbor_ga_ids[0]
                    connected_ga_atoms = self._find_connected_ga_atoms(neighbor_idx, ga_id)
                    for ga_atom in connected_ga_atoms:
                        if self.marked_molecule_data[ga_atom]['active']:
                            fragment_atoms.add(ga_atom)
                else:  # If neighbor does not belong to GA
                    if self.marked_molecule_data[neighbor_idx]['active']:
                        fragment_atoms.add(neighbor_idx)

        return fragment_atoms



    def _find_connected_ga_atoms(self, start_atom_idx, ga_id):
        """
        Find all atoms connected to the specified atom that belong to the same GA instance.

        :param start_atom_idx: Starting atom index
        :param ga_id: Current GA ID
        :return: A list containing indices of all atoms related to this GA instance
        """
        connected_ga_atoms = set()  # Store atoms belonging to the same GA instance
        queue = deque([start_atom_idx])  # Initialize queue
        visited = set()  # Record visited atoms

        while queue:
            current_idx = queue.popleft()
            if current_idx in visited:
                continue
            visited.add(current_idx)

            # Currently, by default, the same GA only appears once in a molecule
            # Check if current atom belongs to the specified GA instance
            current_data = self.marked_molecule_data.get(current_idx, {})
            if ga_id in current_data.get('GA_ids', []):
                connected_ga_atoms.add(current_idx)

                # Add adjacent atoms of this atom to the queue
                current_atom = self.mol.GetAtomWithIdx(int(current_idx.split('_')[1]))
                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = f'atom_{neighbor.GetIdx()}'
                    if neighbor_idx not in visited:
                        queue.append(neighbor_idx)

        return connected_ga_atoms

    
    def match_substructure_and_get_indices(self, substructure_smiles):
        """
        Find all matching substructures in the target molecule and return the list of matching atom indices.

        :param substructure_smiles: SMILES string of the substructure to match
        :return: List of matching atom indices (each element is a tuple of matching atom indices)
        """
        substructure_mol = Chem.MolFromSmiles(substructure_smiles)  # Generate molecule object for substructure
        if substructure_mol is None:
            print(f"Unable to parse substructure: {substructure_smiles}")
            return []

        # Perform substructure matching to find all matching substructures
        matches = self.mol.GetSubstructMatches(substructure_mol)
        if matches:
            # print(f"Matched atom indices: {matches}")
            return matches  # Return list of all matching atom indices
        else:
            print(f"No substructure found: {substructure_smiles}")
            return []
        
    def split_molecule_based_on_IRL(self):
        """
        Split the molecule based on IRL markings and return a list of split molecules, each containing atom indices and marking information.
        """
        split_molecules = []  # Store split molecular fragments
        IRL_unprocessed_list = self.IRL_in_ligand

        # Iterate through IRL list, process each IRL in order
        for IRL_info in self.irl_list:
            IRL_did = IRL_info['DID']  # Get IRL DID
            IRL_smiles = IRL_info['complex_smiles']  # Get IRL SMILES

            if IRL_did not in IRL_unprocessed_list:
                continue
            # print(f"Processing IRL: {IRL_did}")
            # print(f"IRL SMILES: {IRL_smiles}")

            # Find all atom indices containing this IRL
            atom_index_related_to_IRL = [
                atom_idx
                for atom_idx, atom_data in self.marked_molecule_data.items()
                if IRL_did in atom_data['IRL_ids']
            ]

            # If any atom in the IRL is in inactive state, skip this IRL
            if any(not self.marked_molecule_data[atom_idx]['active'] for atom_idx in atom_index_related_to_IRL):
                continue

            # Get specific matching atom indices of IRL in the molecule
            IRL_index_list = self.match_substructure_and_get_indices(IRL_smiles)
            for IRL_index in IRL_index_list:
                # Ensure all atoms in IRL are in active state
                if any(not self.marked_molecule_data[f'atom_{idx}']['active'] for idx in IRL_index):
                    continue
                formatted_IRL_index = [f'atom_{idx}' for idx in IRL_index]
                # Use BFS to find related fragments of this IRL
                fragment_atoms = self.bfs_find_substructure(formatted_IRL_index)
                # Mark all atoms in the fragment as inactive
                for atom in fragment_atoms:
                    self.mark_inactive(atom)
                
                # Add split result to molecule list
                # Save IRL metadata together with fragment atom indices
                split_molecules.append({
                    'IRL_did': IRL_did,  # Add IRL DID
                    'IRL_smiles': IRL_smiles,  # Add IRL SMILES
                    'fragment_atoms': fragment_atoms  # Fragment atom indices
                })

        return split_molecules
    
    def get_substructure_smiles_from_custom_indices(self, atom_indices):
        """
        Extract the complete topological structure of molecular substructure based on custom format atom index list (like 'atom_13'),
        and return its SMILES representation.

        :param atom_indices: List containing custom format atom indices (like 'atom_13')
        :return: SMILES string of the substructure
        """
        if not atom_indices:
            return ""  # If index list is empty, return empty string

        # Extract numeric parts of atom indices
        numeric_indices = [int(idx.split('_')[1]) for idx in atom_indices]

        # Create a new molecule object
        new_mol = Chem.RWMol()
        atom_mapping = {}

        # Add atoms to new molecule
        for idx in numeric_indices:
            atom = self.mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)  # Create new atom
            atom_mapping[idx] = new_idx  # Establish atom index mapping

        # Add bonds to new molecule
        for bond in self.mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            # If both atoms of the bond are in the submolecule, add this bond
            if begin_idx in atom_mapping and end_idx in atom_mapping:
                new_mol.AddBond(
                    atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType()
                )

        # Update molecule property cache and return SMILES
        new_mol.UpdatePropertyCache()
        substructure = new_mol.GetMol()
        return Chem.MolToSmiles(substructure)

    def mark_inactive(self, atom_idx):
        """
        Mark the given atom as inactive, meaning it has been processed and will no longer serve as an IRL starting point.
        """
        atom_data = self.marked_molecule_data.get(atom_idx)  # Get atom marking data
        if atom_data:
            atom_data['active'] = False  # Mark this atom as inactive

    def get_active_atoms(self):
        """
        Get all atom indices that are in active state.

        :return: List of active atom indices
        """
        return [atom_idx for atom_idx, atom_data in self.marked_molecule_data.items() if atom_data['active']]  # Return list of active atom indices

    def get_inactive_atoms(self):
        """
        Get all atom indices that are in inactive state.

        :return: List of inactive atom indices
        """
        return [atom_idx for atom_idx, atom_data in self.marked_molecule_data.items() if not atom_data['active']]  # Return list of inactive atom indices

    @property
    def fragments(self):
        """
        Property that returns the split fragments with their metadata.
        
        Returns:
            List of dictionaries, each containing:
            - fragment_smiles: SMILES representation of the fragment
            - fragment_DID: Canonical DID of the fragment
            - fragment_IRL_did: IRL DID associated with the fragment
            - fragment_IRL_smiles: IRL SMILES associated with the fragment
            - fragment_atoms: List of atom indices in the fragment
        """
        return self.split_and_get_smiles()

    def split_and_get_smiles(self):
        """
        Split the molecule and get SMILES representation and DID for each fragment.
        :return: Each fragment contains fragment_smiles, fragment_DID, fragment_IRL_did, fragment_IRL_smiles
        """
        split_fragments = self.split_molecule_based_on_IRL()
        fragment_smiles_list = []

        for fragment in split_fragments:
            fragment_atoms = fragment['fragment_atoms']
            fragment_smiles = self.get_substructure_smiles_from_custom_indices(fragment_atoms)
            fragment_DID = calculate_canonical_did(fragment_smiles)

            fragment_smiles_list.append({
                'fragment_IRL_did': fragment['IRL_did'],
                'fragment_IRL_smiles': fragment['IRL_smiles'],
                'fragment_smiles': fragment_smiles,
                'fragment_DID': fragment_DID,
                'fragment_atoms': fragment_atoms
            })

        return fragment_smiles_list