import numpy as np

muon_cols = ['Muon_charge', 'Muon_dxy', 'Muon_dxyErr', 'Muon_dz', 'Muon_dzErr', 'Muon_eta', 'Muon_isGlobal', 'Muon_mass',
             'Muon_phi', 'Muon_pt', 'Muon_ptErr', 'Muon_softId', 'Muon_vtxIdx', 'Muon_vtxFlag',]

dimu_cols = ['Dimu_pt', 'Dimu_eta', 'Dimu_phi', 'Dimu_rap', 'Dimu_mass', 'Dimu_charge', 'Dimu_vtxIdx', 'Dimu_chi2', 'Dimu_dl',
             'Dimu_dlErr', 'Dimu_dlSig', 'Dimu_cosphi', 'Dimu_x', 'Dimu_y', 'Dimu_z', 'Dimut1_muIdx', 'Dimut2_muIdx',]

d0_cols = ['D0_pt', 'D0_eta', 'D0_phi', 'D0_rap', 'D0_mass12', 'D0_mass21', 'D0_vtxIdx', 'D0_chi2', 'D0_dl', 'D0_dlErr', 'D0_dlSig',
           'D0_cosphi', 'D0_x', 'D0_y', 'D0_z', 'D0_hasMuon',
           'D0t1_pt', 'D0t1_eta', 'D0t1_phi', 'D0t1_chindof', 'D0t1_nValid', 'D0t1_nPix', 'D0t1_dxy', 'D0t1_dz', 
           'D0t2_pt', 'D0t2_eta', 'D0t2_phi', 'D0t2_chindof', 'D0t2_nValid', 'D0t2_nPix', 'D0t2_dxy', 'D0t2_dz',]

dstar_cols = ['Dstar_pt', 'Dstar_eta', 'Dstar_phi', 'Dstar_rap', 'Dstar_deltam', 'Dstar_deltamr', 'Dstar_vtxIdx', 'Dstar_hasMuon',
              'DstarD0_pt', 'DstarD0_eta', 'DstarD0_phi', 'DstarD0_mass', 'DstarD0_chi2', 'DstarD0_dl', 'DstarD0_dlErr',
              'DstarD0_dlSig', 'DstarD0_cosphi', 'DstarD0_x', 'DstarD0_y', 'DstarD0_z',
              'DstarK_pt', 'DstarK_eta', 'DstarK_phi', 'DstarK_vtxIdx', 'DstarK_chindof', 'DstarK_nValid', 'DstarK_nPix', 'DstarK_dxy',
              'DstarK_dz', 'DstarK_chg',
              'Dstarpi_pt', 'Dstarpi_eta', 'Dstarpi_phi', 'Dstarpi_vtxIdx', 'Dstarpi_chindof', 'Dstarpi_nValid', 'Dstarpi_nPix',
              'Dstarpi_dxy', 'Dstarpi_dz', 'Dstarpi_chg',
              'Dstarpis_pt', 'Dstarpis_eta', 'Dstarpis_phi', 'Dstarpis_vtxIdx', 'Dstarpis_chindof', 'Dstarpis_nValid', 'Dstarpis_nPix',
              'Dstarpis_dxy', 'Dstarpis_dz',]


def get_vars_dict(events, col_list):
    dict = {}
    col = ''
    for c in col_list:
        if c.startswith('Muon'):
            col = c[5:]
        elif c.startswith('Dimu'):
            col = c[4:]
            if col.startswith('_'): col = col[1:]
        elif c.startswith('D0'):
            col = c[2:]
            if col.startswith('_'): col = col[1:]
        elif c.startswith('Dstar'):
            col = c[5:]
            if col.startswith('_'): col = col[1:]
        else:
            Exception('Not good!')

        if len(events[c]) == 0:
            dict[col] = np.array([])
        else:
            dict[col] = events[c]
    return dict
