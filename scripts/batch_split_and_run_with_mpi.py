#!/usr/bin/env python

# T. Carman, H. Genet
# Feb, March 2018

import errno            # nice errors
import os               # mkdir, etc
import shutil           # for copying files, removing dirs
import json             # for parsing the config file
import re               # for stripping comments from json files.
import textwrap         # for nicely formatting the inline slurm script
import subprocess       # for calling squeue
import numpy as np      # for generating numeric data
import netCDF4 as nc    # for handling netcdf files

# This script is used to split a dvmdostem run into "sub domains" that can be
# run individually (submitted to the queue manager) and then merged together
# at the end. In this case, the "full domain" is NOT the entire IEM domain, but
# is simply the full area that you are trying to run, i.e. a 10x10 region, or
# a 50x50 region.

# 1) Log on to atlas, cd into your dvmdostem directory.
#
# 2) Checkout desired git version, setup environment.
#
# 3) Compile (make).
#
# 4) Setup your run as though you were going to run serially (single 
#    processor). Configure as necessary the following things:
#      - paths in the config file to the input data and full-domain run mask
#      - adjust the output_spec.csv to your needs
#      - turn outputs on/off for various run-stages (eq, sp, etc)
#      - path to output data is not important - it will be overwritten
#
# 5) Figure out how many cells you want per batch, and set the constant below.
#
# 6) Run this script.
#
# This script will split your run into however many batches are necesseray to
# run all the cells and keep the max cells per batch in line with the constant
# you set below. The script will setup two directory hierarchies: one for the 
# outputs of the individual batch runs and one for the "staging" area for each
# batch run. The staging area allows each run to have a different config file
# and different run mask (which is what actually controls which cells are 
# in which batch. Then the script will submit a job to slurm for each batch.
# 
# To process the outputs, use the "batch_merge.sh" script.


# USER SHOULD SET THIS VALUE
IDEAL_CELLS_PER_BATCH = 15


# Look in the config file to figure out where the full-domain runmask is.
with open("config/config.js", 'r') as f:
    input_str = f.read()
j = json.loads(re.sub('//.*\n','\n', input_str))
BASE_RUNMASK = j['IO']['runmask_file']

# Figure out how many batches are necessary to complete the full run.
with nc.Dataset(BASE_RUNMASK, 'r') as runmask:
  TOTAL_CELLS_TO_RUN = np.count_nonzero(runmask.variables['run'])  

NBATCHES = TOTAL_CELLS_TO_RUN / IDEAL_CELLS_PER_BATCH
if (TOTAL_CELLS_TO_RUN % IDEAL_CELLS_PER_BATCH != 0) or (TOTAL_CELLS_TO_RUN < IDEAL_CELLS_PER_BATCH):
  print "Adding another batch to pick up stragglers!"
  NBATCHES += 1
   
print "NUMBER OF BATCHES: ", NBATCHES



# Utility function
def mkdir_p(path):
  '''Provides similar functionality to bash mkdir -p'''
  try:
     os.makedirs(path)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise  

#
# SETUP DIRECTORIES
#
print "Removing any existing staging or batch run directories"
if os.path.isdir('staging-batch-run'):
  shutil.rmtree('staging-batch-run')
if os.path.isdir('batch-run'):
  shutil.rmtree('batch-run')

for batch_id in range(0, NBATCHES):

  print "Making directories for batch {}".format(batch_id)
  mkdir_p('staging-batch-run/batch-{}'.format(batch_id))
  mkdir_p('batch-run/batch-{}'.format(batch_id))
  
  print "Copy run mask, config file, etc for batch {}".format(batch_id)
  shutil.copy(BASE_RUNMASK,  'staging-batch-run/batch-{}/'.format(batch_id))
  shutil.copy('config/config.js', 'staging-batch-run/batch-{}/'.format(batch_id))
  
  print "Reset the run mask for batch {}".format(batch_id)
  with nc.Dataset('staging-batch-run/batch-{}/run-mask.nc'.format(batch_id), 'a') as runmask:
    runmask.variables['run'][:] = np.zeros(runmask.variables['run'].shape)
  
#
# BUILD BATCH SPECIFIC RUN-MASKS
#
with nc.Dataset(BASE_RUNMASK, 'r') as runmask:
  nz_ycoords = runmask.variables['run'][:].nonzero()[0]
  nz_xcoords = runmask.variables['run'][:].nonzero()[1]

# For every cell that is turned on in the main run-mask, we assign this cell
# to a batch to be run, and turn on the corresponding cell in the batch's
# run mask.
print "Turning on pixels in each batch's run mask..."
batch = 0
cells_in_sublist = 0
coord_list = zip(nz_ycoords, nz_xcoords)
for i, cell in enumerate(coord_list):

  with nc.Dataset("staging-batch-run/batch-{}/run-mask.nc".format(batch), 'a') as grp_runmask:
    grp_runmask.variables['run'][cell] = True
    cells_in_sublist += 1

  if (cells_in_sublist == IDEAL_CELLS_PER_BATCH) or (i == len(coord_list)-1):
    print "Group {} will run {} cells...".format(batch, cells_in_sublist)
    batch += 1
    cells_in_sublist = 0 


#
# SUMMARIZE
#
number_batches = batch
assert (NBATCHES == number_batches), "PROBLEM: Somthing is wrong with the batch numbers."
print "Split cells into {} batches...".format(number_batches)

#
# MODIFY THE CONFIG FILE FOR EACH BATCH
#
print "Modifying each batch's config file; changing path to run mask and to output directory..."
for batch_num in range(0, number_batches):

  with open('staging-batch-run/batch-{}/config.js'.format(batch_num), 'r') as f:
    input_string = f.read()
  
  j = json.loads(re.sub('//.*\n','\n', input_string)) # Strip comments from json file
  j['IO']['runmask_file'] = 'staging-batch-run/batch-{}/run-mask.nc'.format(batch_num)
  j['IO']['output_dir'] = 'batch-run/batch-{}/output/'.format(batch_num)
  
  output_str = json.dumps(j, indent=2, sort_keys=True)

  with open('staging-batch-run/batch-{}/config.js'.format(batch_num), 'w') as f:
    f.write(output_str)



#
# SUBMIT SBATCH SCRIPT FOR EACH BATCH
#
for batch in range(0, number_batches):

  with nc.Dataset("staging-batch-run/batch-{}/run-mask.nc".format(batch), 'r') as runmask:
    cells_in_batch = np.count_nonzero(runmask.variables['run'])

  assert (cells_in_batch > 0), "PROBLEM! Groups with no cells activated to run!"
  
  slurm_runner_scriptlet = textwrap.dedent('''\
  #!/bin/bash -l

  # Job name, for clarity
  #SBATCH --job-name="ddt-batch-{0}"

  # Reservation
  #SBATCH --reservation=snap_8 

  # Partition specification
  #SBATCH -p main

  # Number of MPI tasks
  #SBATCH -n {1}

  echo $SBATCH_RESERVATION
  echo $SLURM_JOB_NODELIST

  # Load up my custom paths stuff
  module purge
  module load jsoncpp/1.8.1-foss-2016a netCDF/4.4.0-foss-2016a Boost/1.55.0-foss-2016a-Python-2.7.11
  module load netcdf4-python/1.2.2-foss-2016a-Python-2.7.11

  mpirun --mca btl self,tcp \
  --mca btl_tcp_if_include eth2 \
  --mca oob_tcp_if_include eth2 \
  -n {1} \
  ./dvmdostem -f staging-batch-run/batch-{0}/config.js -l disabled --max-output-volume 25GB -p 100 -e 1000 -s 250 -t 109 -n 91 

  '''.format(batch, cells_in_batch + 1))
  print "Writing sbatch script for batch {}".format(batch)
  with open("staging-batch-run/batch-{}/slurm_runner.sh".format(batch), 'w') as f:
    f.write(slurm_runner_scriptlet)
  
  print "STUB::: Submitting sbatch script to squeue for batch {}".format(batch)
  sbatch_output = subprocess.check_output(["sbatch", "staging-batch-run/batch-{}/slurm_runner.sh".format(batch)])
  print sbatch_output



