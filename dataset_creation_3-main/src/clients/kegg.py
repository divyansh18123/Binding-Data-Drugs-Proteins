from ..config import KEGG_LINK_PATHWAY, KEGG_LINK_DISEASE, KEGG_LIST_PATHWAY, KEGG_LIST_DISEASE, KEGG_GET_PATHWAY
from ..utils import fetch_text

def fetch_kegg_lookup(url, prefix=""):
    lookup = {}
    text = fetch_text(url)
    if text:
        for line in text.strip().split("\n"):
            if "\t" in line:
                k, v = line.split("\t")
                lookup[k] = v
                if k.startswith(prefix):
                    lookup[k[len(prefix):]] = v
    return lookup

def fetch_pathways(kegg_gene_id):
    rows = []
    text = fetch_text(KEGG_LINK_PATHWAY.format(kegg_gene_id))
    if text:
        for line in text.strip().split("\n"):
            if "\t" in line:
                _, pathway = line.split("\t")
                rows.append(pathway)
    return rows

def fetch_diseases(kegg_gene_id):
    rows = []
    text = fetch_text(KEGG_LINK_DISEASE.format(kegg_gene_id))
    if text:
        for line in text.strip().split("\n"):
            if "\t" in line:
                _, disease = line.split("\t")
                rows.append(disease)
    return rows

def fetch_kegg_entry(kegg_id):
    """Fetches key-value pairs from a KEGG entry, focusing on DESCRIPTION."""
    entry_data = {}
    text = fetch_text(KEGG_GET_PATHWAY.format(kegg_id)) # KEGG_GET uses /get/ which works for both
    if not text:
        return entry_data
    
    current_key = None
    buffer = []

    for line in text.split("\n"):
        if not line.startswith(" "):
            # New Key
            if current_key:
                entry_data[current_key] = " ".join(buffer).strip()
            
            parts = line.split(maxsplit=1)
            if not parts: continue
            current_key = parts[0]
            if len(parts) > 1:
                buffer = [parts[1]]
            else:
                buffer = []
        else:
            # Continuation
            buffer.append(line.strip())
    
    if current_key:
        entry_data[current_key] = " ".join(buffer).strip()

    return entry_data
