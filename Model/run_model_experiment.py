import os
import sys

# This script takes in a field test location (SVALBARD or CMBRDGBY)
# and a purpose (SPINUP, CNTRL, FLOOD). An optional site name can be 
# provided for the CNTRL and FLOOD purposes. The script will create a
# new Icepack case, build it, copy in the appropriate namelist file, 
# and run the case.


# read in command line arguments
test_loc = sys.argv[1]
purpose = sys.argv[2]
if len(sys.argv) > 3:
    site = sys.argv[3]
    case_name = test_loc + '_' + site + '_' + purpose
else:
    case_name = test_loc + '_' + purpose

# set up directory paths
base = '/Users/mollyw/Desktop/ARISE-AIM/'
Model_dir = base + 'Model/namelists/'
Icepack_dir = base + 'Icepack/'
exp_dir = base + 'icepack-dirs/cases/'

# create the new Icepack case
os.chdir(Icepack_dir)
command = './icepack.setup -m conda -e macos -s ionetcdf,debug -c ~/icepack-dirs/cases/' + case_name
os.system(command)

# build the Icepack case
os.chdir(exp_dir + case_name)
command = './icepack.build'
os.system(command)

# copy in the ARISE-AIM tempate namelist file
if purpose == 'SPINUP':
    command = 'cp '+ Model_dir + 'icepack_in.' + test_loc + '_spinup ' + exp_dir + case_name + '/icepack_in'
elif purpose == 'CNTRL':
    command = 'cp '+ Model_dir + 'icepack_in.' + test_loc + '_' + site + '_cntrl ' + exp_dir + case_name + '/icepack_in'
elif purpose == 'FLOOD':
    command = 'cp '+ Model_dir + 'icepack_in.' + test_loc + '_' + site + '_flood ' + exp_dir + case_name + '/icepack_in'
else:
    print('Invalid purpose. Please use SPINUP, CNTRL, or FLOOD.')
    sys.exit(1)
os.system(command)

# run the Icepack case
command = './icepack.run'
os.system(command)

print('Icepack case run script complete. Please check the Icepack output files.')