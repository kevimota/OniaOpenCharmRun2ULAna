import yaml, sys, os
from fit.fit_fom import *

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run fit routines to data")
    parser.add_argument("-y", "--year", help="Year to fit", nargs='+', type=str, required=True)
    parser.add_argument("-p", "--plot", help="Create FOM plots", action="store_true")
    args = parser.parse_args()

    with open('config/fom.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    with open('config/fit.yaml', 'r') as f:
        fit = yaml.load(f, Loader=yaml.FullLoader)

    for year in args.year:
        if year not in config['path']:
            print(f"Year {year} not found in the config file!")
            sys.exit()

    y = max(args.year)
    with open(f'output/RunII_trigger_processed_vtxfit/{y}/Upsilon_fit_params.yaml') as f:
        fit_params = yaml.load(f, Loader=yaml.FullLoader)
        alpha_CB = fit_params['alpha_CB']['value']
        n_CB = fit_params['n_CB']['value']

    processed_lumi = 0
    with open("config/skim_trigger.yaml", 'r') as f:
        trigger = yaml.load(f, Loader=yaml.FullLoader)['trigger']
    with open("config/lumi.yaml", 'r') as f:
        lumis = yaml.load(f, Loader=yaml.FullLoader)
        for year in args.year:
            for era in lumis[year]:
                processed_lumi += lumis[year][era][trigger[year]]

    f_name = ""
    if len(args.year) > 1:
        for year in args.year:
            f_name += year+'-'
        f_name = f_name[:-1]
    else: 
        f_name = args.year[0]

    if args.plot:
        os.system(f'mkdir -p plots/fom_vtxfit/{f_name}')  

    failed = {}
    for param in config:
        if param == 'path': continue
        failed[param] = []
        if args.plot: 
            path = f'output/fom_vtxfit/{f_name}'
            plot_fom(param, config[param], path, f_name, processed_lumi)
        for value in config[param]:
            path = [p for y in args.year for p in config['path'][y]]
            if args.plot: 
                path = f'output/fom_vtxfit/{f_name}'
                r = plot_results(param, str(value).replace('.', 'p'), path, f_name, processed_lumi)
                if r != 0:
                    failed[param].append(value)
            else: fit_fom(param, path, str(value).replace('.', 'p'), fit, alpha_CB, n_CB, f_name)
    if args.plot:
        if len(failed.keys()) > 0:
            print("Fit failed for following values:")
            print(failed)
        else: print("Good values for chi2 for all the plots!!!")
            