import json

with open("/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/get_brics/drugbank_fgs.json", "r") as f:
    drugbank_fgs = json.load(f)


frags = []

for key in drugbank_fgs:
    print(key)
    for key2 in drugbank_fgs[key]:
        print(key2)
        for frag in drugbank_fgs[key][key2]:
            frags.append(frag)
    


print(len(frags))

# save to smi file 

with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.smi", "w") as f:
    for frag in frags:
        if frag in ["[H]N=C(N)N", "[111InH3]"]:
            continue
        f.write(frag+"\n")



import lmdb


# open(/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/data.mdb)


import lmdb
import pickle

# Open the LMDB environment
env = lmdb.open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/", readonly=True)

with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.smi", "r") as f:
    lines = f.readlines()
    

smiles = [line.strip() for line in lines]


dic = {}

# Begin a read-only transaction
with env.begin() as txn:
    # Create a cursor to iterate through the database
    cursor = txn.cursor()
    
    keys = []
    for key, value in cursor:
        keys.append(key)
    num_keys = len(keys)
    for i in range(num_keys):
        key = str(i).encode()
        value = pickle.loads(txn.get(key))
        smi = smiles[i]
        dic[smi] = value

# save to json

import pickle

with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.pkl", "wb") as f:
    pickle.dump(dic, f)
        
        

    

# Close the environment when done
env.close()


    
