import os

base_folder = '/Users/kevimota/cernbox'
years = ['2016APV', '2016', '2017', '2018']
particles = ['Upsilon', 'Dstar', 'UpsilonDstar']

for year in years:
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
        os.system(f'python nanoAODplus_trigger.py -p {it.path} -y {year} -f -m')

base_folder = 'output/RunII_trigger_processed_vtxfit'

for year in years:
    for it in os.scandir(f'{base_folder}/{year}'):
        if not it.is_dir(): continue
        os.system(f'python tools/save_ttree.py -p {it.path} -o {it.path}')

base_folder = 'output/fom_vtxfit'

for year in years:
    for it in os.scandir(f'{base_folder}/{year}'):
        for it2 in os.scandir(it.path):
            for it3 in os.scandir(it2.path):
                if not it3.is_dir(): continue
                os.system(f'python tools/save_ttree.py -p {it3.path} -o {it3.path} -a')

for year in years:
    for particle in particles:
        os.system(f'python nanoAODplus_fit.py -y {year} -c {particle}')
        os.system(f'python nanoAODplus_fit.py -y {year} -c {particle} -p')