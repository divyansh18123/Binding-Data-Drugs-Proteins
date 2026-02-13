from config import UNIPROT_BASE_URL, UNIPROT_JSON_URL
from utils import fetch_json

def fetch_uniprot_data(uniprot_id):
    url = UNIPROT_BASE_URL.format(uniprot_id)
    data = fetch_json(url)
    if data:
        return data
    # fallback: try explicit .json URL
    return fetch_json(UNIPROT_JSON_URL.format(uniprot_id))