import yaml, sys, os
from fit.fit_fom import *

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run fit routines to data")
    parser.add_argument("-y", "--year", help="Year to fit", type=str, required=True)
    parser.add_argument("-p", "--plot", help="Create FOM plots", action="store_true")
    args = parser.parse_args()

    with open('config/fom.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    with open('config/fit.yaml', 'r') as f:
        fit = yaml.load(f, Loader=yaml.FullLoader)

    if args.year not in config['path']:
        print("Year not found in the config file!")
        sys.exit()

    with open(f'output/RunII_trigger_processed/{args.year}/Upsilon_fit_params.yaml') as f:
        fit_params = yaml.load(f, Loader=yaml.FullLoader)
        alpha_CB = fit_params['alpha_CB']['value']
        n_CB = fit_params['n_CB']['value']

    with open("config/lumi.yaml", 'r') as f:
        lumi = yaml.load(f, Loader=yaml.FullLoader)[args.year]
    with open("config/skim_trigger.yaml", 'r') as f:
        trigger = yaml.load(f, Loader=yaml.FullLoader)['trigger'][args.year]

    if args.plot:
        os.system(f'mkdir -p plots/fom/{args.year}')

    processed_lumi = 0
    for era in lumi:
        processed_lumi += lumi[era][trigger]

    failed = {}
    for param in config:
        if param == 'path': continue
        if args.plot: 
            path = config['path'][args.year]
            path = path[0][:path[0].rfind('/')]
            plot_fom(param, config[param], path, args.year, processed_lumi)
        for value in config[param]:
            path = config['path'][args.year]
            if args.plot: 
                path = path[0][:path[0].rfind('/')]
                r = plot_results(param, str(value).replace('.', 'p'), path, args.year, processed_lumi)
                if r != 0:
                    if param in failed.keys(): failed[param] = [value]
                    else: failed[param].append(value)
            else: fit_fom(param, config['path'][args.year], str(value).replace('.', 'p'), fit, alpha_CB, n_CB, args.year)
    if args.plot:
        if len(failed.keys()) > 0:
            print("Fit failed for following values:")
            print(failed)
        else: print("Good values for chi2 for all the plots!!!")
            