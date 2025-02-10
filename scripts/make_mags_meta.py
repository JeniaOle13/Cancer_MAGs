import argparse
import sys
import os.path
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'prepare bins')
    parser.add_argument('--work_dir', help = "full path to /path/to/data/analysis folder of your project")
    args = parser.parse_args(sys.argv[1:])
    work_dir = args.work_dir

# merge classification tables
    
if os.path.isfile(work_dir +'/gtdbtk/classify/gtdbtk.bac120.summary.tsv'):
    classify_bac = pd.read_csv(work_dir +'/gtdbtk/classify/gtdbtk.bac120.summary.tsv', sep='\t')
    classify_bac = classify_bac[['user_genome', 'classification']]
    classify_bac['species'] = [ann.split('s__')[-1] if ann.split('s__')[-1] != '' else ann.split('g__')[-1].split(';')[0] + '_spp'\
                            for ann in classify_bac.classification.to_list()]
    classify_bac.rename(columns = {'user_genome': 'Feature ID', 'classification': 'Taxon'}, inplace=True)
else:
    classify_bac = pd.DataFrame()
    print('No Bacteria have been identified.')

if os.path.isfile(work_dir +'/gtdbtk/classify/gtdbtk.ar53.summary.tsv'):
    classify_arch = pd.read_csv(work_dir +'/gtdbtk/classify/gtdbtk.ar53.summary.tsv', sep='\t')
    classify_arch = classify_arch[['user_genome', 'classification']]
    classify_arch['species'] = [ann.split('s__')[-1] if ann.split('s__')[-1] != '' else ann.split('g__')[-1].split(';')[0] + '_spp'\
                            for ann in classify_arch.classification.to_list()]
    classify_arch.rename(columns = {'user_genome': 'Feature ID', 'classification': 'Taxon'}, inplace=True)
else:
    classify_arch = pd.DataFrame()
    print('No Archaea have been identified.')

classify = pd.concat([classify_bac, classify_arch])
    
# pull genome quality

checkm_i = pd.read_csv(work_dir +'/drep/data_tables/genomeInfo.csv', 
                       sep = ',')
checkm_i['genome'] = checkm_i.loc[:, 'genome'].apply(
    lambda x: x.split('.fa')[0])
checkm_i.rename(columns = {'genome': 'Feature ID'}, inplace = True)
checkm_i

classify.merge(checkm_i, on = 'Feature ID').to_csv(work_dir +'/drep/mags_metadata.tsv', sep = '\t', index = False)