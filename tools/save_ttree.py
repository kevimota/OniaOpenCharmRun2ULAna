import os, sys, re
from coffea.util import load
import numpy as np
import uproot
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
    if len(file_list) > 1:
        for i in tqdm(file_list[1:], desc="Loading", unit=" files", total=(len(file_list)-1)):
            acc += load(i)
    return acc

def create_ttree(acc, out, n=None, association_only=False):
    if n is None:
        n = ""
    else:
        n = "_" + str(n)

    if not association_only:
        Upsilon = acc['Dimu']
        Dstar = acc['Dstar']
        filename = f"{out}/Upsilon{n}.root"
        print(f"Saving Upsilon data to {filename}")
        with uproot.recreate(filename) as f:
            data_dict = {}
            for key in Upsilon.keys():
                data_dict[key] = Upsilon[key].value
            f['Upsilon'] = data_dict

        filename = f"{out}/Dstar{n}.root"
        print(f"Saving Dstar data to {filename}")
        with uproot.recreate(filename) as f:
            data_dict = {}
            for key in Dstar.keys():
                data_dict[key] = Dstar[key].value
            f['Dstar'] = data_dict
    
    UpsilonDstar = acc['DimuDstar']
    
    filename = f"{out}/UpsilonDstar{n}.root"
    print(f"Saving UpsilonDstar data to {filename}")
    with uproot.recreate(filename) as f:
            data_dict = {}
            for key in UpsilonDstar.keys():
                data_dict[key] = UpsilonDstar[key].value
            f['UpsilonDstar'] = data_dict
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Create TTree from coffea files")
    parser.add_argument("-p", "--path", help="Path of the coffea files", type=str, required=True)
    parser.add_argument("-o", "--out", help="Output file to be written", type=str, required=True)
    parser.add_argument("-a", "--association", help="Save association ttree only", action="store_true")

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
        if args.association: create_ttree(acc, args.out, association_only=True)
        else: create_ttree(acc, args.out)