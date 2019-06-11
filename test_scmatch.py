# test results for scmatch...

import os
import subprocess
import shlex

scmatch_path = '/home/yjzhang/Grad_School/bioNLP/scMatch/scMatch.py'
toterms_path = '/home/yjzhang/Grad_School/bioNLP/scMatch/toTerms.py'
ref_path = '/home/yjzhang/Grad_School/bioNLP/scMatch/refDB/FANTOM5'
data_file = 'tm_droplet_data_subset.csv'

run_command = 'python {0} --refType mouse --testType mouse --refDS {1} --dFormat csv --testDS {2}'.format(scmatch_path, ref_path, data_file)
result = subprocess.call(shlex.split(run_command))

data_file = 'tm_droplet_data_means.csv'
run_command = 'python {0} --refType mouse --testType mouse --refDS {1} --dFormat csv --testDS {2}'.format(scmatch_path, ref_path, data_file)
result = subprocess.call(shlex.split(run_command))