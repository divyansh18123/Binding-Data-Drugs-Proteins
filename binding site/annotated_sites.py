import os
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# ================= CONFIG =================

CIF_DIR = "p2rank_2.5.1/pdb_ids_cifs"
OUTPUT_CSV = "annotated_sites.csv"



rows = []

for fname in os.listdir(CIF_DIR):
    if not fname.endswith(".cif"):
        continue

    pdb_id = os.path.splitext(fname)[0]
    cif_path = os.path.join(CIF_DIR, fname)


    cif = MMCIF2Dict(cif_path)

    # ---------------- struct_site ----------------

    struct_site = {}
    if "_struct_site.id" in cif:
        for i, sid in enumerate(cif["_struct_site.id"]):
            struct_site[sid] = {
                "_struct_site.id": sid,
                "_struct_site.pdbx_evidence_code": cif.get("_struct_site.pdbx_evidence_code", [None])[i],
                "_struct_site.pdbx_auth_asym_id": cif.get("_struct_site.pdbx_auth_asym_id", [None])[i],
                "_struct_site.pdbx_auth_comp_id": cif.get("_struct_site.pdbx_auth_comp_id", [None])[i],
                "_struct_site.pdbx_auth_seq_id": cif.get("_struct_site.pdbx_auth_seq_id", [None])[i],
                "_struct_site.pdbx_auth_ins_code": cif.get("_struct_site.pdbx_auth_ins_code", [None])[i],
                "_struct_site.pdbx_num_residues": cif.get("_struct_site.pdbx_num_residues", [None])[i],
                "_struct_site.details": cif.get("_struct_site.details", [None])[i],
            }

    # ---------------- safety check ----------------

    if "_struct_site_gen.site_id" not in cif:
        print(f"  → No struct_site_gen in {pdb_id}, skipping")
        continue

    # ---------------- atom_site columns ----------------

    atom_site_cols = [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ]

    atom_count = len(cif["_atom_site.id"])

    # ---------------- residue → atom join ----------------

    for i in range(len(cif["_struct_site_gen.site_id"])):

        site_id = cif["_struct_site_gen.site_id"][i]
        auth_chain = cif["_struct_site_gen.auth_asym_id"][i]
        auth_resnum = cif["_struct_site_gen.auth_seq_id"][i]
        auth_comp = cif["_struct_site_gen.auth_comp_id"][i]

        for j in range(atom_count):
            if (
                cif["_atom_site.auth_asym_id"][j] == auth_chain
                and cif["_atom_site.auth_seq_id"][j] == auth_resnum
                and cif["_atom_site.auth_comp_id"][j] == auth_comp
                and cif["_atom_site.pdbx_PDB_ins_code"][j] == cif["_struct_site_gen.pdbx_auth_ins_code"][i]
            ):
                row = {
                    "pdb_id": pdb_id,

                    # -------- struct_site --------
                    **struct_site.get(site_id, {}),

                    # -------- struct_site_gen --------
                    "_struct_site_gen.id": cif["_struct_site_gen.id"][i],
                    "_struct_site_gen.site_id": site_id,
                    "_struct_site_gen.pdbx_num_res": cif["_struct_site_gen.pdbx_num_res"][i],
                    "_struct_site_gen.label_comp_id": cif["_struct_site_gen.label_comp_id"][i],
                    "_struct_site_gen.label_asym_id": cif["_struct_site_gen.label_asym_id"][i],
                    "_struct_site_gen.label_seq_id": cif["_struct_site_gen.label_seq_id"][i],
                    "_struct_site_gen.pdbx_auth_ins_code": cif["_struct_site_gen.pdbx_auth_ins_code"][i],
                    "_struct_site_gen.auth_comp_id": auth_comp,
                    "_struct_site_gen.auth_asym_id": auth_chain,
                    "_struct_site_gen.auth_seq_id": auth_resnum,
                    "_struct_site_gen.label_atom_id": cif["_struct_site_gen.label_atom_id"][i],
                    "_struct_site_gen.label_alt_id": cif["_struct_site_gen.label_alt_id"][i],
                    "_struct_site_gen.symmetry": cif["_struct_site_gen.symmetry"][i],
                    "_struct_site_gen.details": cif["_struct_site_gen.details"][i],
                }

                # -------- atom_site --------
                for col in atom_site_cols:
                    row[col] = cif[col][j] if col in cif else None

                rows.append(row)

# ---------------- write CSV ----------------

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)

print(f"\nSaved {OUTPUT_CSV}")

