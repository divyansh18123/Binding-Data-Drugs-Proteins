import pandas as pd

name = ["protein","protein_features","protein_claims","protein_claims","protein_structures","protein_kegg_pathways","protein_kegg_diseases","protein_interactions","protein_go","ontology_edges","ontology_terms","ontology_pathways", "ontology_diseases","protein_tissue","protein_variants","tissue_nodes","tissue_nodes","tissue_hierarchy"]


for i in name:

    try:
        df = pd.read_csv(f"{i}.csv", engine = 'python')
    except Exception as e:
        print(f"Couldn't read {i}", e)
        continue


    df = df.drop_duplicates()
    df = df.fillna('None')
    df.to_csv(f"{i}.csv", index = False)

