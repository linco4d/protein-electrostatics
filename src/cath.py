import os
import pandas as pd
from collections import defaultdict

#Download all cath classifications, domains
os.system('wget -c ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt')
os.system('wget -c ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-superfamily-list.txt')

def parse_cath_domain_list(filepath):
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split()
            if len(tokens) < 7:
                continue
            domain_id = tokens[0]
            c, a, t, h = tokens[1:5]
            #resolution = float(tokens[6]) if tokens[6] != '-' else float('inf')
            data.append((domain_id, f"{c}.{a}.{t}", f'{h}'))#,  resolution))
    return data

#Do not use this function -- we care less about resolution.
def get_best_resolution_per_superfamily(domain_data):
    best_structures = {}
    for domain_id, cat_key, homologues, resolution in domain_data:
        if cat_key not in best_structures or resolution < best_structures[(cat_key, homologues)][1]:
            best_structures[(cat_key, homologues)] = (domain_id, resolution)
    return best_structures

domain_file = "cath-domain-list.txt"

domain_data = parse_cath_domain_list(domain_file)
#This returns a list of tuples (i.e. domain_id, CAT, H)

#df_best = pd.DataFrame([(k[0], k[1], v[0], v[1]) for k, v in best_per_cat.items()],
#                      columns=['C.A.T.', 'H',  'Domain_ID', 'Resolution'])
#rossman_set = df_best[df_best['C.A.T.'] == '3.40.50']

#@TODO: After you have a list of all the domains as well as the CAT, H classifications
# extract PDB ids, and move on install the fastas, and cif/pdb files

