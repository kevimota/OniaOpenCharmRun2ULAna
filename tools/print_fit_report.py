import yaml
from uncertainties import ufloat

years = ['2016APV', '2016', '2017', '2018']

for year in years:
    with open(f"output/RunII_trigger_processed_vtxfit/{year}/UpsilonDstar_fit_params.yaml", 'r') as f:
        fit_params = yaml.load(f, Loader=yaml.FullLoader)
    n_evt = fit_params['n_evt']
    chi2_upsilon = fit_params['chi2_upsilon']
    chi2_dstar = fit_params['chi2_dstar']
    mscale = ufloat(fit_params['mscale']['value'], fit_params['mscale']['error_hi'])
    dstar_mean = ufloat(fit_params['dstar_mean']['value'], fit_params['dstar_mean']['error_hi'])
    signal_frac = ufloat(fit_params['signal_frac']['value'], fit_params['signal_frac']['error_hi'])
    upsilon1S_frac = ufloat(fit_params['upsilon1S_frac']['value'], fit_params['upsilon1S_frac']['error_hi'])
    upsilon2S_frac = ufloat(fit_params['upsilon1S_frac']['value'], fit_params['upsilon2S_frac']['error_hi'])

    N_signal = n_evt * signal_frac
    N_upsilon1S = N_signal*upsilon1S_frac
    N_upsilon2S = N_signal*(1-upsilon1S_frac)*(upsilon2S_frac)
    N_upsilon3S = N_signal*(1-upsilon1S_frac)*(1-upsilon2S_frac)

    print(f"------------- Report for year {year:^10} -------------")
    print(f"n_evt        = {n_evt:.0f}")
    print(f"N_upsilon1S  = {N_upsilon1S:.0f}")
    print(f"N_upsilon2S  = {N_upsilon2S:.0f}")
    print(f"N_upsilon3S  = {N_upsilon3S:.0f}")
    print(f"mscale       = {mscale:.4f}")
    print(f"dstar_mean   = {dstar_mean*1e3:.2f}")
    print(f"chi2_upsilon = {chi2_upsilon:.2f}")
    print(f"chi2_dstar   = {chi2_dstar:.2f}")
    print("\n")
