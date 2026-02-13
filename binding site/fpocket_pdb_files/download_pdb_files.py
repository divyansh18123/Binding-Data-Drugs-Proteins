# Instructions
# path should be fpocket_pdb_files
# change the pdb_ids in this file to include all the pdb_ids of the dataset
# make sure that there are no files after the header in fpocket_rescore.ds so that you dont get repeated data
# run this file using python download_pdb_files.py
# run fpocket using (after installing it through docker)
# docker run --platform linux/amd64 --rm \
#   -v "$(pwd)/pdbs:/data" \
#   fpocket/fpocket \
#   sh -c 'for f in /data/*.pdb; do fpocket -f "$f"; done'
# merge all
# run info_extract.py

import os
import requests


# pdb_ids = ['5TE0', '5OAZ', '4EIH', '7MRZ', '3SOC', '4PED', '1UNQ', '8ZXW', '6C9H', '6E2N', '6GR8', '5VXZ', '4W9W', '3MDY', '2HLQ', '8X2A', '8C7X', '6CZ4', '6T29', '6Y3O', '5EFQ', '6Q4G', '7XQK', '6P8E', '7SJ3', '3MI9', '4AAA', '3VA4', '6FT8', '6FYR', '4HNI', '2CMW', '2IZR', '3US4', '3BHY', '1MG4', '2WUH', '7QUE', '6ZJF', '8C2Z', '8A27', '3HIL', '7OFV', '2R2P', '3KUL', '5L6O', '7BHF', '4QTB', '8AOJ', '7AQB', '5BYZ', '6W39', '6YOJ', '5EW8', '2PVF', '8UDU', '1RJB', '1FLT', '4BSJ', '4O38', '7E2K', '9H8F', '1P4O', '3BRT', '3BRV', '6ZIW', '4HCU', '6N7A', '8BXH', '5LWM', '2XRW', '7N8T', '7ORF', '2VIF', '7LWH', '4ZRI', '7QHG', '2WTK', '6EIM', '5MY9', '7NX1', '6LDW', '2C60', '2OKR', '8TXY', '5ES1', '7B9L', '4H3Q', '3ALN', '2NPT', '2Y8O', '4R1V', '5Y90', '5JQA', '3DTC', '2RF0', '5K26', '5OTE', '4L0N', '3A7I', '5YF4', '4DRI', '8S9P', '3KF9', '2X4F', '6BXI', '5XQZ', '3ZKE', '2FST', '3GP0', '8X23', '1CM8', '3FXZ', '3PCS', '6FD3', '8C12', '9ASF', '1E8Y', '8Q6F', '3A99', '3H1Z', '3X01', '7QPN', '4WB8', '5VCZ', '3ORK', '4I5P', '4B6L', '3COK', '1WMH', '5F9E', '5JIX', '8BI7', '6CNH', '3CC6', '8R4V', '6FDM', '8X92', '7Z5X', '1VZO', '7P74', '4NW6', '6G77', '5WBH', '6HVD', '5YKS', '2BUJ', '7NTH', '6BDN', '8XTT', '3BEA', '7LO0', '4YFI', '5JFW', '4ASZ', '6KZD', '4H7Y', '2DM0', '3LXP', '6HYO', '2XIR', '2V62', '5VDK', '4FR4', '2HDA', '7Z4V', '5X5O', '7SIY'] #Davis ids

# pdb_ids = ['4X9R', '8A27', '1MFG', '7SXF', '1O6L', '7ORF', '3GP0', '3BRV', '3BRT', '4DNL', '4HNI', '7SIY', '7Z87', '3O0G', '4WB8', '4AU8', '6VPM', '2UZP', '1SHA', '4QTB', '8BJU', '4BSJ', '2XIR', '1WMH', '4EIH', '6GZD', '6FT8', '6BFN', '3BEA', '6S14', '2FST', '6GU2', '2WH0', '6TGU', '6Q38', '3WAR', '1BLX', '3BHH', '2V7O', '3GP2', '6X5G', '3S95', '2D9Z', '5OAZ', '2I0E', '2FK9', '5F9E', '1YRK', '4YJR', '3BU3', '2OKR', '8VCF', '4DRI', '6P8E', '6Q4G', '2IIM', '7A2P', '6NMW', '9ASF', '2Q8G', '5EW8', '2PVF', '8KH9', '8UDU', '6YOJ', '6PYR', '8Q6F', '2R5T', '1UNQ', '8AOJ', '2XRW', '4R1V', '1FLT', '6T29', '7DUA', '6FNK', '1E8Y', '5H0B', '8ZXW', '5GRN', '7B7R', '3A99', '1RHF', '6SRH', '4U3Y', '6FYL', '2W5A', '6C9H', '3BHY', '8JOT', '6KZD', '1P4O', '5LWM', '5P9J', '4HZR', '7MZY', '4ASZ', '2J0I', '5XS2', '8BXH', '4FWW', '8C2Z', '4H3Q', '8P79', '3IQU', '5WBH', '2Y8O', '3MI9', '8C7X', '1RJB', '4NW6', '5JFW', '7NX1', '7N8T', '1VZO', '4OTH', '4AF3', '5VXZ', '6RCG', '2C47', '2CMW', '4TND', '5LVO', '2YEX', '4B6L', '5MY9', '5EP6', '8P0S', '6FYV', '3VA4', '8TXY', '5V62', '4Y85', '8X23', '1CM8', '2VIF', '4HCU', '2HDA', '6Q7D', '6CZ4', '2UV4', '8BIK', '6EGE', '4X7Q', '5YF4', '3EAZ', '3COK', '8X92', '7T1K', '7Z5X', '7JT9', '2Y7J', '6KC4', '8X2A', '3DTC', '3PCS', '4L0N', '7MBJ', '3CC6', '2AHX', '7NRB', '3DLS', '6E2N', '7QUE', '5J0A', '2W4O', '4FG8', '3FXZ', '2H6D', '6Y8A', '6N7A', '5X5O', '3LXP', '8EBL', '5K26', '8UOJ', '1MG4', '2F15', '6HVD', '2AC3', '3BGM', '7AB0', '4RRV', '6P5S', '2RF0', '6YA6', '8GM5', '2Y9Q', '2IZR', '5K00', '5ES1', '5NNG', '3US4', '5Y86', '5JIX'] #Kiba

# pdb_ids = ['4H3Q', '6C9H', '4WB8', '6VPM', '4RRV', '2UZP', '3FXZ', '4BSJ', '5J0A', '3BHY', '8C2Z', '3MI9', '4TND', '5LWM', '4ASZ', '2IZR', '5EW8', '2D9Z', '7NX1', '7B7R', '3A99', '1XJD', '6Q4G', '3A7I', '6HVD', '5JIX', '1WMH', '4AF3', '7Z5X', '2OKR', '7QUE', '2CMW', '4HCU', '4HZR', '1UNQ', '8X92', '8ZXW', '6T29', '5X5O', '8TXY', '5VXZ', '3COK', '7A2P', '6CZ4', '1MG4', '3VA4', '4B6L', '2VIF', '1P4O', '7O86', '3BRT', '5ES1', '8GM4', '6Y3O'] #metz

# Output directory
OUTPUT_DIR = "pdbs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASE_URL = "https://files.rcsb.org/download"

def download_pdb(pdb_id):
    pdb_id = pdb_id.lower()
    url = f"{BASE_URL}/{pdb_id}.pdb"
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}.pdb")

    try:
        response = requests.get(url, timeout=20)
        if response.status_code == 200:
            with open(out_path, "wb") as f:
                f.write(response.content)
            print(f"Downloaded {pdb_id}.pdb")

            with open("fpocket_rescore.ds", "a") as f2:
                f2.write(f"pdbs/{pdb_id}_out/{pdb_id}_out.pdb  pdbs/{pdb_id}.pdb\n")

    
        else:
            print(f"Failed ({response.status_code}) for {pdb_id}")
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")

if __name__ == "__main__":
    for pdb_id in pdb_ids:
        download_pdb(pdb_id)
