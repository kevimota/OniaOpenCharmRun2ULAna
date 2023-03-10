# Caller for fitting Upsilon, Dstar and UpsilonDstar
# Upsilon must be fitted before fitting UpsilonDstar, n_CB and alpha_CB parameters
# are fixed in UpsilonDstar and are fetched from previous Upsilon fitting.

import yaml, sys
from fit.fit_functions import *

years = ['2016APV', '2016', '2017', '2018', 'all']

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run fit routines to data")
    parser.add_argument("-y", "--year", help="Year to fit", type=str, required=True, choices=years)
    parser.add_argument("-c", "--channel", help="Particle to fit", type=str, required=True, choices={"Upsilon", "Dstar", "UpsilonDstar"})
    parser.add_argument("-s", "--sps", action="store_true", help="Do full Run 2 fit on data")
    parser.add_argument("-p", "--plot", action="store_true", help="Plot the results of previous fit")
    parser.add_argument("-cm", "--check-matrix", action="store_true", help="Get Covariant matrix status")
    args = parser.parse_args()

    with open('config/fit.yaml') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    save_path = f'output/fit/{args.year}'
    if not os.path.exists(save_path): os.makedirs(save_path)

    if not args.plot:
        if not args.check_matrix:
            if args.sps:
                print(f"Fitting {args.channel} SPS for year {args.year}")
                if args.channel == "Upsilon":
                        fit_upsilon_sps(config, args.year)
                if args.channel == "Dstar":
                    fit_dstar_sps(config, args.year)
                if args.channel == "UpsilonDstar":
                    print('No fit for UpsilonDstar SPS')

            else:
                print(f"Fitting {args.channel} for year {args.year}")
                if args.channel == "Upsilon":
                    fit_upsilon(config, args.year)

                if args.channel == "Dstar":
                    fit_dstar(config, args.year)
                    
                if args.channel == "UpsilonDstar":
                    fit_upsilondstar(config, args.year)
        else:
            path = f'output/fit/{args.year}'
            check_cov_matrix(path, args.year)
    else:
        print(f"Ploting {args.channel} for year {args.year}")
        plot_fit(args.year, args.channel)
        

