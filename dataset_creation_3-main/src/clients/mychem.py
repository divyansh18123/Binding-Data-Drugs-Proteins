import requests
from ..utils import fetch_json

MYCHEM_URL = "https://mychem.info/v1/query"

class MyChemClient:
    def __init__(self):
        pass

    def fetch_data(self, chembl_id):
        """
        Fetches aggregated drug data from MyChem.info using ChEMBL ID.
        
        Returns a dictionary containing:
        - _id: MyChem ID
        - drugbank: {id, ...}
        - pubchem: {cid, ...}
        - pharmgkb: {id, ...}
        - sider: [{side_effect, ...}]
        - admet: { ... } (if available)
        """
        # We assume chembl_id is like 'CHEMBL123'
        # MyChem query syntax: q=chembl.molecule_chembl_id:CHEMBL123
        params = {
            "q": f"chembl.molecule_chembl_id:{chembl_id}",
            "fields": "chembl,drugbank,pubchem,pharmgkb,sider,admet,aeolus"
        }
        
        try:
            data = fetch_json(MYCHEM_URL, params=params)
            if data and "hits" in data and len(data["hits"]) > 0:
                return data["hits"][0]
        except Exception as e:
            print(f"Error fetching MyChem data for {chembl_id}: {e}")
            
        return None

    def fetch_data_by_name(self, name):
        """
        Fetches data by drug name (fallback).
        """
        params = {
            "q": f'"{name}"', # Exact phrase match
            "fields": "chembl,drugbank,pubchem,pharmgkb,sider,admet,aeolus"
        }
        try:
            data = fetch_json(MYCHEM_URL, params=params)
             # We might get multiple hits, take the best one? 
             # Usually the first one is best match by score.
            if data and "hits" in data and len(data["hits"]) > 0:
                return data["hits"][0]
        except Exception as e:
            print(f"Error fetching MyChem data for {name}: {e}")
        return None
