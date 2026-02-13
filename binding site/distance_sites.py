import os
import pandas as pd
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# ================= CONFIG =================

CIF_DIR = "p2rank_2.5.1/pdb_ids_cifs"
DISTANCE_CUTOFF = 5.0  # Å
OUTPUT_CSV = "distance_sites.csv"

EXCLUDE = {"HOH", "WAT", "DOD", "SO4", "PO4", "PEG", "MPD", "GOL"}



all_rows = []

for fname in os.listdir(CIF_DIR):
    if not fname.endswith(".cif"):
        continue


    pdb_id = os.path.splitext(fname)[0]
    cif_path = os.path.join(CIF_DIR, fname)

    cif = MMCIF2Dict(cif_path)

    # ---------------- atom_site ----------------

    atom_cols = [
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

    atom_df = pd.DataFrame({c: cif[c] for c in atom_cols})

    atom_df[["_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z"]] = atom_df[
        ["_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z"]
    ].astype(float)

    # ---------------- identify ligands ----------------

    ligand_atoms = atom_df[
        (atom_df["_atom_site.group_PDB"] == "HETATM") &
        (~atom_df["_atom_site.auth_comp_id"].isin(EXCLUDE))
    ]

    # AlphaFold safety check
    if ligand_atoms.empty:
        print(f"No ligands found in {pdb_id}, skipping")
        continue

    protein_atoms = atom_df[atom_df["_atom_site.group_PDB"] == "ATOM"]

    lig_coords = ligand_atoms[[
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z"
    ]].to_numpy()

    prot_coords = protein_atoms[[
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z"
    ]].to_numpy()

    # ---------------- distance calculation ----------------

    close_residues = set()

    for i, pcoord in enumerate(prot_coords):
        dists = np.linalg.norm(lig_coords - pcoord, axis=1)
        if np.any(dists <= DISTANCE_CUTOFF):
            row = protein_atoms.iloc[i]
            close_residues.add((
                row["_atom_site.auth_asym_id"],
                row["_atom_site.auth_seq_id"],
                row["_atom_site.pdbx_PDB_ins_code"],
            ))

    # ---------------- expand to full residues ----------------

    mask = atom_df.apply(
        lambda r: (
            r["_atom_site.auth_asym_id"],
            r["_atom_site.auth_seq_id"],
            r["_atom_site.pdbx_PDB_ins_code"],
        ) in close_residues,
        axis=1
    )

    binding_site_atoms = atom_df.loc[mask].copy()
    binding_site_atoms["pdb_id"] = pdb_id

    # ---------------- metadata ----------------

    binding_site_atoms["_reflns.d_resolution_high"] = (
        cif.get("_reflns.d_resolution_high", [None])[0]
    )

    binding_site_atoms["_exptl.method"] = (
        cif.get("_exptl.method", [None])[0]
    )

    binding_site_atoms["_entity_poly.pdbx_seq_one_letter_code"] = (
        cif.get("_entity_poly.pdbx_seq_one_letter_code", [None])[0]
    )

    if "_entity_poly.pdbx_seq_one_letter_code" in binding_site_atoms:
        binding_site_atoms["_entity_poly.pdbx_seq_one_letter_code"] = (
            binding_site_atoms["_entity_poly.pdbx_seq_one_letter_code"]
            .astype(str)
            .str.replace("\n", "", regex=False)
        )

    # ---------------- struct_conn ----------------

    struct_conn_cols = [
        "_struct_conn.id",
        "_struct_conn.conn_type_id",
        "_struct_conn.pdbx_leaving_atom_flag",
        "_struct_conn.pdbx_PDB_id",
        "_struct_conn.ptnr1_label_asym_id",
        "_struct_conn.ptnr1_label_comp_id",
        "_struct_conn.ptnr1_label_seq_id",
        "_struct_conn.ptnr1_label_atom_id",
        "_struct_conn.pdbx_ptnr1_label_alt_id",
        "_struct_conn.pdbx_ptnr1_PDB_ins_code",
        "_struct_conn.pdbx_ptnr1_standard_comp_id",
        "_struct_conn.ptnr1_symmetry",
        "_struct_conn.ptnr2_label_asym_id",
        "_struct_conn.ptnr2_label_comp_id",
        "_struct_conn.ptnr2_label_seq_id",
        "_struct_conn.ptnr2_label_atom_id",
        "_struct_conn.pdbx_ptnr2_label_alt_id",
        "_struct_conn.pdbx_ptnr2_PDB_ins_code",
        "_struct_conn.ptnr1_auth_asym_id",
        "_struct_conn.ptnr1_auth_comp_id",
        "_struct_conn.ptnr1_auth_seq_id",
        "_struct_conn.ptnr2_auth_asym_id",
        "_struct_conn.ptnr2_auth_comp_id",
        "_struct_conn.ptnr2_auth_seq_id",
        "_struct_conn.ptnr2_symmetry",
        "_struct_conn.pdbx_ptnr3_label_atom_id",
        "_struct_conn.pdbx_ptnr3_label_seq_id",
        "_struct_conn.pdbx_ptnr3_label_comp_id",
        "_struct_conn.pdbx_ptnr3_label_asym_id",
        "_struct_conn.pdbx_ptnr3_label_alt_id",
        "_struct_conn.pdbx_ptnr3_PDB_ins_code",
        "_struct_conn.details",
        "_struct_conn.pdbx_dist_value",
        "_struct_conn.pdbx_value_order",
        "_struct_conn.pdbx_role",
    ]

    if "_struct_conn.id" in cif:
        for col in struct_conn_cols:
            binding_site_atoms[col] = cif.get(col, [None])[0]
    else:
        for col in struct_conn_cols:
            binding_site_atoms[col] = None

    binding_site_atoms["_struct_conn_type.id"] = (
        cif.get("_struct_conn_type.id", [None])[0]
    )
    binding_site_atoms["_struct_conn_type.criteria"] = (
        cif.get("_struct_conn_type.criteria", [None])[0]
    )
    binding_site_atoms["_struct_conn_type.reference"] = (
        cif.get("_struct_conn_type.reference", [None])[0]
    )

    all_rows.append(binding_site_atoms)

# ---------------- save ----------------

if not all_rows:
    print("No binding sites detected.")

else:
    final_df = pd.concat(all_rows, ignore_index=True)

    cols = ["pdb_id"] + [c for c in final_df.columns if c != "pdb_id"]
    final_df = final_df.loc[:, cols]

    final_df.to_csv(OUTPUT_CSV, index=False)

    print(f"\nSaved {OUTPUT_CSV}")
