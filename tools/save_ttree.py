import os, sys, re
from coffea.util import load
import numpy as np
import uproot3
from tqdm import tqdm
from pprint import pprint
import gc

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def load_acc(file_list):
    acc = load(file_list[0])
    for i in tqdm(file_list[1:], desc="Loading", unit=" files", total=(len(file_list)-1)):
        acc += load(i)
    return acc

def create_ttree(acc, out, n=None):
    if n is None:
        n = ""
    else:
        n = "_" + str(n)

    Upsilon = acc['Dimu']
    Dstar = acc['Dstar']
    UpsilonDstar = acc['DimuDstar']

    filename = f"{out}/Upsilon{n}.root"
    print(f"Saving Upsilon data to {filename}")
    with uproot3.recreate(filename) as f:
        tree_dict = {}
        data_dict = {}
        for key in Upsilon.keys():
            tree_dict[key] = "float64"
            data_dict[key] = Upsilon[key].value
        f['Upsilon'] = uproot3.newtree(tree_dict)
        f['Upsilon'].extend(data_dict)

    filename = f"{out}/Dstar{n}.root"
    print(f"Saving Dstar data to {filename}")
    with uproot3.recreate(filename) as f:
        tree_dict = {}
        data_dict = {}
        for key in Dstar.keys():
            tree_dict[key] = "float64"
            data_dict[key] = Dstar[key].value
        f['Dstar'] = uproot3.newtree(tree_dict)
        f['Dstar'].extend(data_dict)

    filename = f"{out}/UpsilonDstar{n}.root"
    print(f"Saving UpsilonDstar data to {filename}")
    with uproot3.recreate(f"{filename}") as f:
        tree_dict = {}
        data_dict = {}
        for key in UpsilonDstar.keys():
            tree_dict[key] = "float64"
            data_dict[key] = UpsilonDstar[key].value
        f['UpsilonDstar'] = uproot3.newtree(tree_dict)
        f['UpsilonDstar'].extend(data_dict)
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Create TTree from coffea files")
    parser.add_argument("-p", "--path", help="Path of the coffea files", type=str, required=True)
    parser.add_argument("-o", "--out", help="Output file to be written", type=str, required=True)

    args = parser.parse_args()

    with os.scandir(args.path) as it:
        file_list = []
        for f in it:
            if f.name.endswith('.coffea') and f.name.find('hists') == -1:
                file_list.append(f.path)
    
    if len(file_list) == 0:
        print("No coffea files found in the given directory. Check the path and try again!")
        sys.exit()
    file_list.sort(key=natural_keys)

    chunks = [file_list[x:x+50] for x in range(0, len(file_list), 50)]

    if len(chunks) > 1:
        for idx, c in enumerate(chunks):
            acc = load_acc(c)
            create_ttree(acc, args.out, idx)
            del acc
            gc.collect()
    else:
        acc = load_acc(chunks[0])
        create_ttree(acc, args.out)