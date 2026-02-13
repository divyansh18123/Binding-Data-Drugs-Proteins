from chembl_webresource_client.new_client import new_client

def fetch_molecule(chembl_id):
    try:
        return new_client.molecule.get(chembl_id)
    except Exception:
        return None

def fetch_mechanisms(chembl_id):
    try:
        return new_client.mechanism.filter(molecule_chembl_id=chembl_id)
    except Exception:
        return []

def fetch_activities(chembl_id):
    try:
        return new_client.activity.filter(molecule_chembl_id=chembl_id)
    except Exception:
        return []

def fetch_indications(chembl_id):
    try:
        return new_client.drug_indication.filter(molecule_chembl_id=chembl_id)
    except Exception:
        return []


def fetch_target(target_chembl_id):
    try:
        return new_client.target.get(target_chembl_id)
    except Exception:
        return None

def fetch_metabolism(chembl_id):
    try:
        return new_client.metabolism.filter(drug_chembl_id=chembl_id)
    except Exception:
        return []

def fetch_atc_details(atc_code):
    try:
        return new_client.atc_class.filter(level5=atc_code)
    except Exception:
        return []

def fetch_adme_data(chembl_id):
    """
    Fetches ADME specific activities (Assay Type 'A').
    """
    try:
        # Filter for ADME assays
        return new_client.activity.filter(molecule_chembl_id=chembl_id, assay_type='A')
    except Exception:
        return []
