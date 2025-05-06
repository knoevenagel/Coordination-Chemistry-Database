import os
import csv
import json
import mysql.connector

def load_mysql_config(mysql_connect_config_path):

    with open(mysql_connect_config_path,'r') as f:
        config = json.load(f)
        
    mysql_config =  config['database']

    return mysql_config

def load_neo4j_config(neo4j_connect_config_path):

    pass

