import re, os, math
import numpy as np


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]
            
def get_root_files(paths):
    files = []
    for path in paths:
        with os.scandir(path) as it:
            for file in it:
                if file.name.endswith('.root') and (file.stat().st_size != 0):
                    files.append(file.path)
    files.sort(key=natural_keys)
    return files

def get_vars_dict(events, col_list):
    dict = {}
    col = ''
    for c in col_list:
        col = c[c.find("_")+1:]
        if len(events[c]) == 0:
            dict[col] = np.array([])
        else:
            dict[col] = events[c]

        if not c[:c.find("_")] == 'PVtx':
            if col == 'x' or col == 'y' or col == 'z':
                col = 'v' + col

        if len(events[c]) == 0:
            dict[col] = np.array([])
        else:
            dict[col] = events[c]
    return dict

def build_4mom_string(candidate):
    return f'({candidate.pt:.2f}, {candidate.eta:.2f}, {candidate.phi:.2f}, {candidate.mass:.2f})'

def to_cartesian(cand):
    x = cand.pt*math.cos(cand.phi)
    y = cand.pt*math.sin(cand.phi)
    z = cand.pt*math.sinh(cand.eta)
    t = math.sqrt(x*x + y*y + z*z + cand.mass*cand.mass)

    return t, x, y, z

def sum_cand(cand1, cand2):
    t1, x1, y1, z1 = to_cartesian(cand1)
    t2, x2, y2, z2 = to_cartesian(cand2)

    tr = t1+t2
    xr = x1+x2
    yr = y1+y2
    zr = z1+z2

    r = math.sqrt(xr*xr + yr*yr + zr*zr)
    pt = math.sqrt(xr*xr + yr*yr)
    eta = math.asinh(zr / r)
    phi = math.atan2(yr, xr)
    mass = math.sqrt(tr*tr - xr*xr - yr*yr - zr*zr)

    return f'({pt:.2f},{eta:.2f},{phi:.2f},{mass:.2f})'


def print_gen_candidate(candidate):
    quadrimom = build_4mom_string(candidate)
    print(f"PDGId: {candidate.pdgId}, Id: {candidate.Id}, motherPDGId: {candidate.parpdgId}, remoteMotherPDGId: {candidate.sparpdgId}, motherId: {candidate.genPartIdxMother}, n_daughters: {candidate.numberOfDaughters}, status: {candidate.status}, 4-momentum (pt, eta, phi, mass): {quadrimom}, vtx: ({candidate.vx:.2f}, {candidate.vy:.2f}, {candidate.vz:.2f})")

def get_hlt(events, cols):
    dict = {}
    for col in cols:
        if not col in events.fields:
            if len(events['run']) == 0:
                dict[col] = np.array([])
            else:
                dict[col] = np.zeros(len(events), dtype=bool)
        else:
            if len(events[col]) == 0:
                dict[col] = np.array([])
            else:
                dict[col] = events[col]
    return dict

def get_n_from_hist(hist):
    n_cand = np.sum(hist.axes[0].centers * hist.counts())
    n_evt = hist.sum() - hist.values()[0]
    return n_cand, n_evt

cuts_category = {
    0: "No cuts",
    1: "Dimu charge = 0",
    2: "8.5 < Dimu mass < 11.5",
    3: "Muon pt > 3 GeV",
    4: "Muon |eta| < 2.5",
    5: "Muon Soft Id",
    6: "D* tracks are not muons",
    7: "K and pi charge = 0 ",
    8: "D* |eta| < 2.5",
    9: "D0 from D* |eta| < 2.5",
    10: "D* pt > 2 GeV",
    11: "D* D0 tracks reduced chi2 < 2.5",
    12: "D* D0 tracks hits >4 in tracker >1 pix",
    13: "D* D0 tracks dxy < 0.1",
    14: "D* D0 tracks dz < 1",
    15: "D* pi slow track reduced chi2 < 3",
    16: "D* pi slow # hits in the tracker > 2",
    17: "D0 of D* |mass - D0_PDG_mass| < 0.028",
    18: "D0 of D* cossine of pointing angle > 0.99",
    19: "D0 of D* decay length significance > 3",
    20: "Apply the Triggers",
}

""" cuts_category = {
    0: "Muon > 1 GeV",
    1: "Muon > 2 GeV",
    2: "Muon > 3 GeV",
    3: "Muon > 4 GeV",
    4: "Muon > 5 GeV",
    5: "Muon > 6 GeV",
    6: "Muon > 7 GeV",
    7: "Muon > 8 GeV",
    8: "Muon > 9 GeV",
} """