import os

years = ['2016APV', '2016', '2017', '2018']

for year in years:
    print(f'Creating efficiency for year {year}')
    os.system(f'python nanoAODplus_efficiency.py -y {year} -p')