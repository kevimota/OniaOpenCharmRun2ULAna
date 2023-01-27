import os
import sys

particles = ['Upsilon', 'Dstar', 'UpsilonDstar']
years = ['2016APV', '2016', '2017', '2018']

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Helper to create fits")
    parser.add_argument("-y", "--years", help="Year to fit", type=str, nargs='+', required=True)
    parser.add_argument("-c", "--channels", help="Particle to fit", type=str, nargs='+', required=True)
    parser.add_argument("-cf", "--create_files", action="store_true", help="create files for the fit", default=False)
    args = parser.parse_args()

    years_to_run = [year for year in args.years if year in years]
    particles_to_run = [particle for particle in args.channels if particle in particles]

    if (len(years_to_run) == 0) or (len(particles_to_run) == 0):
        sys.exit(1)

    print('Running for the following years')
    print(years_to_run)
    print('and the following channels')
    print(particles_to_run)

    if args.create_files:
        base_folder = '/Users/kevimota/cernbox'
        for year in years_to_run:
            for it in os.scandir(base_folder):
                if it.name.find('_output') < 0: continue
                if year == '2016APV':
                    if it.name.find('HIPM') < 0: continue
                else:
                    if it.name.find('HIPM') > -1: continue
                    if it.name.find(year) < 0: continue
                print(f"Running for {it.name}")
                if not it.is_dir(): continue
                os.system(f'python nanoAODplus_trigger.py -p {it.path} -y {year} -m')

        base_folder = 'output/RunII_trigger_processed_vtxfit'
        for year in years_to_run:
            for it in os.scandir(f'{base_folder}/{year}'):
                if not it.is_dir(): continue
                os.system(f'python tools/save_ttree.py -p {it.path} -o {it.path}')

    base_folder = 'output/fom_vtxfit'

    for year in years_to_run:
        for particle in particles_to_run:
            os.system(f'python nanoAODplus_fit.py -y {year} -c {particle}')
            os.system(f'python nanoAODplus_fit.py -y {year} -c {particle} -p')