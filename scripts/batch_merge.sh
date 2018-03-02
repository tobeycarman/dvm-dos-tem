#!/bin/bash

# Merge dvm-dos-tem outputs along spatial dimensions

# H. Genet, T. Carman 
# Feb 2018

OUTPUT_SPEC_PATH="./config/output_spec.csv"
STAGES="eq sp tr sc"
TIMESTEPS="daily monthly yearly"
BATCH_DIR="batch-run"
FINAL_DIR="all-merged"

mkdir -p "all-merged"

# First handle all the normal outputs.
for variable in $(cat $OUTPUT_SPEC_PATH | cut -d, -f1)
do
  echo "Processing variable: $variable"
  if [ $variable != 'Name' ] ; then   # ignore the header
    # Loop through the various run modes
    for stage in $STAGES
    do
      echo "  --> stage: $stage"
      # Loop through the various time steps
      for timestep in $TIMESTEPS
      do
        echo "  --> timestep: $timestep"
        # Determine the file name of the outputs variable for the specific run mode and time step
        filename="${variable}_${timestep}_${stage}.nc" # Not sure why we need {}??
        # List all the output files for the variable in question in every output sub-directory (one directory = one sub-regional run)
        echo "  --> find $filename"
        filelist=$(find $BATCH_DIR -maxdepth 4 -type f -name $filename)
        
        echo "THE FILE LIST IS: $filelist"
        # Concatenate all these files together
        echo "merge files"
        if [ ! -z "$filelist" ] ; then
          # Something is messed up with my quoting, as this only works with 
          # the filelist variable **unquoted** which I think is bad practice.
          ncea -O -h -y avg $filelist "$FINAL_DIR/$filename"
        else
          echo "nothing to do"
        fi
      done
    done
  fi
done

# Next handle the run_status file
filelist=$(find $BATCH_DIR -maxdepth 4 -type f -name "run_status.nc")
if [ ! -z "$filelist" ] ; then
  # Something is messed up with my quoting, as this only works with 
  # the filelist variable **unquoted** which I think is bad practice.
  ncea -O -h -y avg $filelist "$FINAL_DIR/run_status.nc"
else
  echo "nothing to do - no run_status.nc files found?"
fi

# Finally, handle the fail log
filelist=$(find $BATCH_DIR -maxdepth 4 -type f -name "fail_log.txt")
if [ ! -z "$filelist" ] ; then
  for f in $filelist
  do
    cat $f >> "$FINAL_DIR/fail_log.txt"
  done
else
  echo "nothing to do - no fail_log.txt files found?"
fi




