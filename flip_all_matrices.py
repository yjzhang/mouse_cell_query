# TODO: flip all matrices from cells x genes to genes x cells

import os
import subprocess

import scipy.io

# folders:

path1 = 'tabula_muris/cell_type_matrices_tm_droplet'
path2 = 'tabula_muris/cell_type_matrices_tm_facs'
path3 = 'mca_microwell_seq/cell_type_matrices_mca'

for path in [path1, path2, path3]:
    for filename in os.listdir(path):
        if filename.endswith('.mtx.gz'):
            full_filename = os.path.join(path, filename)
            data = scipy.io.mmread(full_filename)
            data = data.T
            scipy.io.mmwrite(full_filename[:-3], data)
            subprocess.call(['gzip', '-f', full_filename[:-3]])


for path in [path1, path2, path3]:
    for filename in os.listdir(path):
        if filename.endswith('.mtx'):
            full_filename = os.path.join(path, filename)
            subprocess.call(['gzip', '-f', full_filename])
