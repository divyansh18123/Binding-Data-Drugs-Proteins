import os
import re
import csv

ROOT_DIR = "pdbs"          # change if needed
OUTPUT_CSV = "fpocket_info.csv"

# Regex patterns
pocket_header_re = re.compile(r"Pocket\s+(\d+)\s*:")
key_value_re = re.compile(r"\s*([^:]+?)\s*:\s*([-\d\.]+)")

rows = []
all_keys = set()

for root, dirs, files in os.walk(ROOT_DIR):
    for file in files:
        if file.endswith("_info.txt"):
            pdb_id = file.replace("_info.txt", "")
            info_path = os.path.join(root, file)

            with open(info_path) as f:
                current_pocket = None
                pocket_data = {}

                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    # Pocket start
                    pocket_match = pocket_header_re.match(line)
                    if pocket_match:
                        # Save previous pocket
                        if current_pocket is not None:
                            rows.append(pocket_data)
                            all_keys.update(pocket_data.keys())

                        pocket_number = int(pocket_match.group(1))
                        current_pocket = pocket_number
                        pocket_data = {
                            "pdb_id": pdb_id,
                            "pocket_id": pocket_number
                        }
                        continue

                    # Key-value lines
                    kv_match = key_value_re.match(line)
                    if kv_match and current_pocket is not None:
                        key = kv_match.group(1).strip()
                        value = kv_match.group(2).strip()

                        # Normalize column name
                        key = key.lower()
                        key = key.replace(" ", "_")
                        key = key.replace(".", "")
                        key = key.replace("-", "")
                        key = key.replace("__", "_")

                        pocket_data[key] = float(value)

                # Save last pocket
                if current_pocket is not None:
                    rows.append(pocket_data)
                    all_keys.update(pocket_data.keys())

# Write CSV
priority_cols = ["pdb_id", "pocket_id"]
other_cols = sorted(k for k in all_keys if k not in priority_cols)
fieldnames = priority_cols + other_cols

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"Saved {len(rows)} pockets to {OUTPUT_CSV}")
