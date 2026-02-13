import pandas as pd

name = ["drug","drug_protein_edges","drug_indication","drug_atc","drug_metabolism","drug_similarity","drug_side_effects","drug_pharmacogenomics","drug_clinical_labels","enzyme","atc_nodes","atc_edges","drug_atc_hierarchy","adme_nodes","drug_adme_edges"]

for i in name:


    if i == "drug_side_effects":
        df = pd.read_csv(f"{i}.csv", engine = 'python')
        df = df.fillna('None')
        df.to_csv(f"{i}.csv", index = False)
        continue


    try:
        df = pd.read_csv(f"{i}.csv", engine = 'python')
    except Exception as e:
        print(f"Couldn't read {i}", e)
        continue


    df = df.drop_duplicates()
    df = df.fillna('None')
    df.to_csv(f"{i}.csv", index = False)

