import time
import os, sys, subprocess
from tqdm import tqdm

from coffea.util import save, load
import coffea.processor as processor
import numpy as np

def merger(name):
    tstart = time.time()

    # find the directory named as the analyzer name
    if (subprocess.run("find output/ -type d -name '" + name + "'", shell=True, stdout=subprocess.PIPE).stdout.decode("utf-8") == ''):
        raise Exception("Directory not found!")
    print("Merging files in output/" + name)

    # find the files in the directory
    files = subprocess.run("find output/" + name + "/ -type f -name '" + name + "*' -not -path *merged*", shell=True, stdout=subprocess.PIPE)
    file_list = files.stdout.decode("utf-8").splitlines()
    if len(file_list) == 0:
        raise Exception("no files in directory!")

    #loading files into the acumulator and merging then
    for idx, f in tqdm(enumerate(file_list), desc="Merging", unit=" files", total=len(file_list)):
        if (idx == 0): 
            acc = load(file_list[0])
        else:
            acc += load(f)
        os.system("rm -rf " + f)

    #finally saving the merged file into the folder merged/
    print("Saving as output/" + name + "/merged/" + name + "_merged.coffea")
    os.system("mkdir -p output/" + name + "/merged")
    save(acc, "output/" + name + "/merged/" + name + "_merged.coffea")

    elapsed = round(time.time() - tstart, 2)
    print(f"Merge finished in: {elapsed} s")