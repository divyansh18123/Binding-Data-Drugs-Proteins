from ..config import PDB_ENTRY_URL
from ..utils import fetch_json

def fetch_pdb_entry(pdb_id):
    return fetch_json(PDB_ENTRY_URL.format(pdb_id))
