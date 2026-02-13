from ..config import QUICKGO_TERM_URL, QUICKGO_ANCESTOR_URL, QUICKGO_MULTI_TERM_URL
from ..utils import fetch_json

def fetch_term(go_id):
    data = fetch_json(QUICKGO_TERM_URL.format(go_id))
    if data and data.get("results"):
        return data["results"][0]
    return None

def fetch_ancestors(go_id):
    data = fetch_json(QUICKGO_ANCESTOR_URL.format(go_id))
    if data and data.get("results"):
        return data["results"][0].get("ancestors", [])
    return []

def fetch_multiple_terms(ids_list):
    if not ids_list:
        return []
    ids_str = ",".join(ids_list)
    data = fetch_json(QUICKGO_MULTI_TERM_URL.format(ids_str))
    if data:
        return data.get("results", [])
    return []
