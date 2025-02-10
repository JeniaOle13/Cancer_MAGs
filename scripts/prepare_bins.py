import argparse
import sys
import os.path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'prepare bins')
    parser.add_argument('--work_dir', help = "full path to /path/to/data/analysis folder of your project")
    args = parser.parse_args(sys.argv[1:])
    work_dir = args.work_dir 

os.mkdir(work_dir + '/dastool_bins')

samples = [s for s in os.listdir(work_dir + '/dastool/') 
           if os.path.isdir(work_dir + '/dastool/' + s)]

for smp in samples:
    temp_dir = work_dir +'/dastool/'+ smp +'/'+ smp +'_DASTool_bins/'
    
    if os.path.isdir(temp_dir):
        files = os.listdir(temp_dir)
        for f in files:
            if f.startswith('maxbin_') or f.startswith('metabat_'):
                os.symlink(temp_dir + f, work_dir +'/dastool_bins/'+ f)
            elif f.find('.0') != -1:
                os.symlink(temp_dir + f, work_dir +'/dastool_bins/'+'maxbin_'+ f)
            elif f.find('SemiBin_') != -1:
                os.symlink(temp_dir + f, work_dir +'/dastool_bins/'+'semibin_'+ smp +'.'+ f.split('SemiBin_')[-1])
            else:
                os.symlink(temp_dir + f, work_dir +'/dastool_bins/'+'metabat_'+ f)
                
    else:
        continue