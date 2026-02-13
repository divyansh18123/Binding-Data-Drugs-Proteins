import os
import pandas as pd

ROOT_DIR = "pdbs"

# fpocket HEADER index → column name
HEADER_MAP = {
    "0": "pocket_score",
    "1": "drug_score",
    "2": "num_alpha_spheres",
    "3": "mean_alpha_radius",
    "4": "mean_alpha_solvent_acc",
    "5": "mean_bfactor",
    "6": "hydrophobicity_score",
    "7": "polarity_score",
    "8": "aa_volume_score",
    "9": "pocket_volume_mc",
    "10": "pocket_volume_ch",
    "11": "charge_score",
    "12": "local_hydrophobic_density",
    "13": "num_apolar_alpha",
    "14": "prop_apolar_alpha",
}

pocket_rows = []
atom_rows = []

for pdb_id in os.listdir(ROOT_DIR):
    pocket_dir = os.path.join(ROOT_DIR, pdb_id, "pockets")
    if not os.path.isdir(pocket_dir):
        continue

    for fname in os.listdir(pocket_dir):
        if not fname.endswith("_atm.pdb"):
            continue

        pocket_path = os.path.join(pocket_dir, fname)

        pocket_data = {
            "pdb_id": pdb_id[:len(pdb_id) - 4],
            "pocket_file": fname
        }

        with open(pocket_path) as f:
            for line in f:

                # ---- POCKET-LEVEL DATA ----
                import re

                if line.startswith("HEADER") and ":" in line:
                    left, value = line.split(":", 1)

                    # extract the LAST number in the HEADER line
                    match = re.findall(r"\d+", left)
                    if not match:
                        continue

                    idx = match[-1]  # correct fpocket index

                    val = value.strip()
                    if idx in HEADER_MAP and val != "":
                        pocket_data[HEADER_MAP[idx]] = float(val)



                # ---- ATOM-LEVEL DATA ----
                elif line.startswith("ATOM"):
                    atom_rows.append({
                        "pdb_id": pdb_id[:len(pdb_id) - 4],
                        "pocket_file": fname,
                        "atom_serial": int(line[6:11].strip()),
                        "atom_name": line[12:16].strip(),
                        "res_name": line[17:20].strip(),
                        "chain": line[21].strip(),
                        "res_num": int(line[22:26].strip()),
                        "x": float(line[30:38].strip()),
                        "y": float(line[38:46].strip()),
                        "z": float(line[46:54].strip()),
                        "occupancy": float(line[54:60].strip()),
                        "b_factor": float(line[60:66].strip()),
                        "element": line[76:78].strip()
                    })

        pocket_rows.append(pocket_data)

# ---- WRITE CSVs ----
pd.DataFrame(pocket_rows).to_csv("pockets.csv", index=False)
pd.DataFrame(atom_rows).to_csv("pocket_atoms.csv", index=False)

print("Saved: pockets.csv, pocket_atoms.csv")
