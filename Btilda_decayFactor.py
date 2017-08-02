import sys

suffix = sys.argv[1]
if not suffix == '.dat':
    suffix = '_' + suffix

efit_file_name = sys.argv[2]

print suffix, efit_file_name
