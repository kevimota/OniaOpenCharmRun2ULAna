import os

years = ['2016APV', '2016', '2017', '2018']
for year in years:
    os.system(f'python nanoAODplus_plotter.py -y {year} -s raw')
    os.system(f'python nanoAODplus_plotter.py -y {year} -s sel')
    os.system(f'python nanoAODplus_plotter.py -y {year} -m')

os.system('python tools/data_mc_comparison.py')