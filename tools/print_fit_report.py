import os
import yaml
from uncertainties import ufloat

years = ['2016APV', '2016', '2017', '2018', 'all']

for year in years:
    with open(f"output/fit/{year}/UpsilonDstar_fit_params.yaml", 'r') as f:
        fit_params = yaml.load(f, Loader=yaml.FullLoader)
    
    n_evt = fit_params['n_evt']
    chi2_upsilon = fit_params['chi2_upsilon']
    chi2_dstar = fit_params['chi2_dstar']
    mscale = ufloat(fit_params['mscale']['value'], fit_params['mscale']['error_hi'])
    dstar_mean = ufloat(fit_params['dstar_mean']['value'], fit_params['dstar_mean']['error_hi'])
    signal_frac = ufloat(fit_params['signal_frac']['value'], fit_params['signal_frac']['error_hi'])
    upsilon1S_frac = ufloat(fit_params['upsilon1S_frac']['value'], fit_params['upsilon1S_frac']['error_hi'])
    upsilon2S_frac = ufloat(fit_params['upsilon2S_frac']['value'], fit_params['upsilon2S_frac']['error_hi'])

    N_signal = n_evt * signal_frac
    N_upsilon1S = N_signal*upsilon1S_frac
    N_upsilon2S = N_signal*(1-upsilon1S_frac)*(upsilon2S_frac)
    N_upsilon3S = N_signal*(1-upsilon1S_frac)*(1-upsilon2S_frac)


    print(f"------------- Report for year {year:^10} -------------")
    print(f"n_evt             = {n_evt:.0f}")
    print(f"N_upsilon1S       = {N_upsilon1S:.0f}")
    print(f"N_upsilon2S       = {N_upsilon2S:.0f}")
    print(f"N_upsilon3S       = {N_upsilon3S:.0f}")
    print(f"mscale            = {mscale:.4f}")
    print(f"dstar_mean        = {dstar_mean*1e3:.2f}")
    print(f"chi2_upsilon      = {chi2_upsilon:.2f}")
    print(f"chi2_dstar        = {chi2_dstar:.2f}")

    if not os.path.exists(f"output/fit/{year}/Upsilon_SPS_fit_params.yaml"): 
        print('\n')
        continue
    if not os.path.exists(f"output/fit/{year}/Dstar_SPS_fit_params.yaml"): 
        print('\n')   
        continue
    with open(f"output/fit/{year}/Upsilon_SPS_fit_params.yaml", 'r') as f:
        fit_params_y_sps = yaml.load(f, Loader=yaml.FullLoader)
    with open(f"output/fit/{year}/Dstar_SPS_fit_params.yaml", 'r') as f:
        fit_params_dstar_sps = yaml.load(f, Loader=yaml.FullLoader)

    n_evt_upsilon_sps = fit_params_y_sps['n_evt']
    #chi2_upsilon_sps = fit_params_y_sps['chi2_upsilon']
    mscale_sps = ufloat(fit_params_y_sps['mscale']['value'], fit_params_y_sps['mscale']['error_hi'])
    upsilon1S_frac_sps = ufloat(fit_params_y_sps['upsilon1S_frac']['value'], fit_params_y_sps['upsilon1S_frac']['error_hi'])
    upsilon2S_frac_sps = ufloat(fit_params_y_sps['upsilon2S_frac']['value'], fit_params_y_sps['upsilon2S_frac']['error_hi'])
    upsilon3S_frac_sps = ufloat(fit_params_y_sps['upsilon3S_frac']['value'], fit_params_y_sps['upsilon3S_frac']['error_hi'])

    N_upsilon1S_sps = n_evt_upsilon_sps*upsilon1S_frac_sps
    N_upsilon2S_sps = n_evt_upsilon_sps*(1-upsilon1S_frac_sps)*(upsilon2S_frac_sps)
    N_upsilon3S_sps = n_evt_upsilon_sps*(1-upsilon1S_frac_sps)*(1-upsilon2S_frac_sps)*upsilon3S_frac_sps

    n_evt_dstar_sps = fit_params_dstar_sps['n_evt']
    #chi2_dstar_sps = fit_params_dstar_sps['chi2_dstar']
    dstar_mean_sps = ufloat(fit_params_dstar_sps['dstar_mean']['value'], fit_params_dstar_sps['dstar_mean']['error_hi'])
    signal_frac_sps = ufloat(fit_params_dstar_sps['signal_frac']['value'], fit_params_dstar_sps['signal_frac']['error_hi'])

    N_dstar_sps = n_evt_dstar_sps*signal_frac_sps

    print(f"N_evt_upsilon_SPS = {n_evt_upsilon_sps:.0f}")
    print(f"N_upsilon1S_SPS   = {N_upsilon1S_sps:.0f}")
    print(f"N_upsilon2S_SPS   = {N_upsilon2S_sps:.0f}")
    print(f"N_upsilon3S_SPS   = {N_upsilon3S_sps:.0f}")
    print(f"mscale_SPS        = {mscale_sps:.4f}")
    print(f"N_evt_dstar_sps   = {n_evt_dstar_sps:.0f}")
    print(f"N_dstar_SPS       = {N_dstar_sps:.0f}")
    print(f"dstar_mean_SPS    = {dstar_mean_sps*1e3:.2f}")
    #print(f"chi2_upsilon_SPS = {chi2_upsilon_sps:.2f}")
    #print(f"chi2_dstar_SPS   = {chi2_dstar_sps:.2f}")
    print("\n")
