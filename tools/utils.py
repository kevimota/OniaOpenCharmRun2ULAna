import os, re
import awkward as ak
import numpy as np

D0_PDG_MASS = 1.864

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def build_p4(acc):
    p4 = ak.zip({'x': acc['x'].value,
                 'y': acc['y'].value,
                 'z': acc['z'].value,
                 't': acc['t'].value}, with_name="LorentzVector")

    return p4

def build_acc(obj, name=None):
    from coffea import processor
    acc = processor.dict_accumulator({})
    for var in obj.fields:
        acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(obj[var])))
    if name is not None:
        acc[f"n{name}"] = processor.column_accumulator(ak.to_numpy(ak.num(obj)))
    return acc


def association(cand1, cand2):
    ''' Function for association of the particles. The cuts that operates on all of them and 
    computation of quantities can go here. individual cuts can go on the main processing'''
    cand2 = cand2[cand2.associationIdx >= 0]
    asso = ak.zip({'0': cand1[cand2.associationIdx], '1': cand2})
    
    cand1 = ak.zip({
            'pt': asso.slot0.pt,
            'eta': asso.slot0.eta,
            'phi': asso.slot0.phi,
            'mass': asso.slot0.mass,
            'charge': asso.slot0.charge}, with_name="PtEtaPhiMCandidate")

    cand2 = ak.zip({
            'pt': asso.slot1.pt,
            'eta': asso.slot1.eta,
            'phi': asso.slot1.phi,
            'mass': asso.slot1.mass,
            'charge': asso.slot1.charge}, with_name="PtEtaPhiMCandidate")

    asso['deltarap'] = asso.slot0.rap - asso.slot1.rap
    asso['deltapt'] = asso.slot0.pt - asso.slot1.pt
    asso['deltaeta'] = asso.slot0.eta - asso.slot1.eta
    asso['deltaphi'] = np.remainder(np.abs(asso.slot0.phi - asso.slot1.phi), np.pi)
    asso['cand'] = cand1 + cand2
    
    return asso

def get_files(paths, pattern='.root', exclude=None):
    files = []
    for path in paths:
        for it in os.scandir(path):
            exclude_file = True
            if it.name.find(pattern) > -1 and (it.stat().st_size != 0):
                exclude_file = False
                if not exclude == None:
                    for e in exclude:
                        if it.name.find(e) > -1: exclude_file = True
            if not exclude_file: files.append(it.path)
    files.sort(key=natural_keys)
    return files

def save_kin_hists(hists, cand, gen=False, get_deltam=False):
    hists['pt'].fill(pt=ak.flatten(cand.pt))
    hists['eta'].fill(eta=ak.flatten(cand.eta))
    hists['phi'].fill(phi=ak.flatten(cand.phi))
    if not gen: 
        if get_deltam: hists['deltam'].fill(deltam=ak.flatten(cand.deltamr))
        else: hists['mass'].fill(mass=ak.flatten(cand.mass))

def get_n(array):
    array2 = array[ak.fill_none(array, -1) > -1]
    return ak.num(array2)

def remove_none(array):
    return ak.fill_none(array, -99) > -1