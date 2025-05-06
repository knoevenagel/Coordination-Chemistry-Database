import json
import yaml
from neo4j import GraphDatabase

class Neo4jWriter:
    def __init__(self, yaml_config):
        neo4j_config_path = yaml_config['database']['neo4j_config_path']
        with open(neo4j_config_path,'r') as f:
            config = yaml.safe_load(f)
        self.url = config['url']
        self.user = config['user']
        self.password = config['password']
        self.driver = GraphDatabase.driver(self.url, auth=(self.user, self.password))

    def close(self):
        self.driver.close()

    def insert_complex_and_ligands(self, complex_info_dict):
        with self.driver.session() as session:
            session.write_transaction(self._create_nodes_and_relationships, complex_info_dict)

    @staticmethod
    def _create_nodes_and_relationships(tx, complex_info_dict):
        complex_data = complex_info_dict['complex_info']
        complex_did = complex_data['DID']
        complex_smiles = complex_data['complex_smiles']
        complex_inactive = complex_data['inactive']

        # 1. 创建 Complex 节点（唯一）
        tx.run("""
            MERGE (c:Complex {DID: $DID})
            SET c.smiles = $smiles, c.inactive = $inactive
        """, DID=complex_did, smiles=complex_smiles, inactive=complex_inactive)

        # 2. 遍历 Ligand（neighbor_info）创建 Ligand 节点和关系
        for neighbor in complex_info_dict.get('neighbor_info', []):
            ligand_did = neighbor['ligand_did']
            ligand_smiles = neighbor['ligand_smiles']
            bond_type = neighbor.get('bond_type', 'UNKNOWN')
            coordinating_atom_index = neighbor.get('coordinating_atom_index')
            central_atom = neighbor.get('central_atom', [])

            tx.run("""
                MERGE (l:Ligand {DID: $ligand_did})
                SET l.smiles = $ligand_smiles

                MERGE (c:Complex {DID: $complex_did})
                MERGE (c)-[r:HAS_LIGAND]->(l)
                SET r.bond_type = $bond_type,
                    r.coordinating_atom_index = $coordinating_atom_index,
                    r.central_atom = $central_atom
            """, ligand_did=ligand_did,
                 ligand_smiles=ligand_smiles,
                 complex_did=complex_did,
                 bond_type=bond_type,
                 coordinating_atom_index=coordinating_atom_index,
                 central_atom=central_atom)