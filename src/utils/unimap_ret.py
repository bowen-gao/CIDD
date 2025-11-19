# import json

# with open("/mnt/nfs-ssd/data/huangyanwen/SubstructureDetection/get_brics/drugbank_fgs.json", "r") as f:
#     drugbank_fgs = json.load(f)


# frags = []

# for key in drugbank_fgs:
#     print(key)
#     for key2 in drugbank_fgs[key]:
#         print(key2)
#         for frag in drugbank_fgs[key][key2]:
#             frags.append(frag)
    


# print(len(frags))

# # save to smi file 

# with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.smi", "w") as f:
#     for frag in frags:
#         if frag in ["[H]N=C(N)N", "[111InH3]"]:
#             continue
#         f.write(frag+"\n")



# import lmdb


# # open(/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/data.mdb)


# import lmdb
# import pickle

# # Open the LMDB environment
# env = lmdb.open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/", readonly=True)

# with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.smi", "r") as f:
#     lines = f.readlines()
    

# smiles = [line.strip() for line in lines]


# dic = {}

# # Begin a read-only transaction
# with env.begin() as txn:
#     # Create a cursor to iterate through the database
#     cursor = txn.cursor()
    
#     keys = []
#     for key, value in cursor:
#         keys.append(key)
#     num_keys = len(keys)
#     for i in range(num_keys):
#         key = str(i).encode()
#         value = pickle.loads(txn.get(key))
#         smi = smiles[i]
#         dic[smi] = value

# # save to json

# import pickle


# home_dir = "/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/"

# with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.pkl", "wb") as f:
#     pickle.dump(dic, f)
        
        

    

# # Close the environment when done
# env.close()


home_dir = "/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/"


import os

import sys

import subprocess

import numpy as np
import lmdb

import pickle

def ret_fragments(query_frag_dic):

    query_frags = list(query_frag_dic.keys())

    
    # save to .smi

    # remove .lmdb file in home_dir

    if os.path.exists(f"{home_dir}/data.mdb"):
        os.remove(f"{home_dir}/data.mdb")


    with open(f"{home_dir}/temp_frags.smi", "w") as f:
        for frag in query_frags:

            f.write(frag+"\n")
    
    # get cur_dir of the file

    cur_dir = os.path.dirname(os.path.realpath(__file__))

    # run bash cur_dir/extract_smi.sh using subprocess

    result = subprocess.run([f"{cur_dir}/extract_smi.sh"], stdout=subprocess.PIPE)

    #os.system(f"{cur_dir}/extract_smi.sh")

    # check job done

    if result.returncode != 0:
        print("Error")
        sys.exit(1)
    
    # open lmdb

    env = lmdb.open(f"{home_dir}/", readonly=True)

    with open(f"{home_dir}/temp_frags.smi_new", "r") as f:
        lines = f.readlines()
    
    frags = [line.strip() for line in lines]

    dic = {}
    with env.begin() as txn:
        cursor = txn.cursor()
        
        keys = []
        for key, value in cursor:
            keys.append(key)
        num_keys = len(keys)
        for i in range(num_keys):
            key = str(i).encode()
            value = pickle.loads(txn.get(key))
            smi = frags[i]
            dic[smi] = value
    env.close()

    with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.pkl", "rb") as f:
        all_frags = pickle.load(f)
    
    # get max ten cos sim frags for each frag in dic
    res_dic = {}
    for frag in dic:
        emb = dic[frag]
        sims = []
        visited = set()
        for frag2 in query_frag_dic[frag]:
            if frag2 not in all_frags:
                continue
            if frag2 in visited:
                continue
            visited.add(frag2)
            emb2 = all_frags[frag2]
            cos_sim = np.dot(emb, emb2)/(np.linalg.norm(emb)*np.linalg.norm(emb2))
            sims.append((frag2, cos_sim))
        sims = sorted(sims, key=lambda x:x[1], reverse=True)
        res_dic[frag] = [x[0] for x in sims[:10]]
    print(res_dic)
    return res_dic


def ret_fragments_debug(query_frags):

   

    
    # save to .smi

    # remove .lmdb file in home_dir

    if os.path.exists(f"{home_dir}/data.mdb"):
        os.remove(f"{home_dir}/data.mdb")


    with open(f"{home_dir}/temp_frags.smi", "w") as f:
        for frag in query_frags:

            f.write(frag+"\n")
    
    # get cur_dir of the file

    cur_dir = os.path.dirname(os.path.realpath(__file__))

    # run bash cur_dir/extract_smi.sh using subprocess

    result = subprocess.run([f"{cur_dir}/extract_smi.sh"], stdout=subprocess.PIPE)

    #os.system(f"{cur_dir}/extract_smi.sh")

    # check job done

    if result.returncode != 0:
        print("Error")
        sys.exit(1)
    
    # open lmdb

    env = lmdb.open(f"{home_dir}/", readonly=True)

    with open(f"{home_dir}/temp_frags.smi_new", "r") as f:
        lines = f.readlines()
    
    frags = [line.strip() for line in lines]

    dic = {}
    with env.begin() as txn:
        cursor = txn.cursor()
        
        keys = []
        for key, value in cursor:
            keys.append(key)
        num_keys = len(keys)
        for i in range(num_keys):
            key = str(i).encode()
            value = pickle.loads(txn.get(key))
            smi = frags[i]
            dic[smi] = value
    env.close()

    with open("/mnt/nfs-ssd/data/fengshikun/UniMAP/iupac-pretrain/all_feats/all_frags.pkl", "rb") as f:
        all_frags = pickle.load(f)
    
    # get max ten cos sim frags for each frag in dic

    print(dic)
    res_dic = {}
    for frag in dic:
        emb = dic[frag]
        sims = []
        visited = set()
        for frag2 in all_frags:
            if frag2 not in all_frags:
                continue
            if frag2 in visited:
                continue
            visited.add(frag2)
            emb2 = all_frags[frag2]
            cos_sim = np.dot(emb, emb2)/(np.linalg.norm(emb)*np.linalg.norm(emb2))
            sims.append((frag2, cos_sim))
        sims = sorted(sims, key=lambda x:x[1], reverse=True)
        res_dic[frag] = [x[0] for x in sims[:10]]
    print(res_dic)
    return res_dic


if __name__ == "__main__":
    frags = ["C1=CC=CC=C1"]
    ret_fragments_debug(frags)
            
        
        









        
        




    
