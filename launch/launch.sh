#!/usr/bin/bash
# Launch script to launch DNS code on Richardson Cluster

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Run Settings -----------------------------------------------------------------------------------
runName=kz3_8_rk1_2  # must not exist yet in runFolder
runTime=1:00:00  # format: hh:mm:ss (Millikan/Richardson)
jobDependency="none"  # jobid of run (Millikan/Richardson) for dependency condition
# jobDependency is optional, put "none" if not needed
mpiProcessors=96  # must be equal to nprocs in ctes3D

inputFile="/home/yh/channel_rough_data/runs/kz3_7/kz3_7.001"

# Run Paramters
# Reynolds number
Re=8850
# write a restart file every nimag
nimag=15000
# total time steps, must be multiple of nimag + 1
nstep=3
# time step to update CFL and write to .cf
# no larger thatn 10
nhist=5

# Time Stepping prameters
#  Adaptive time stepping:
#    Set CFL to non zero ( no larger than 1.5 ), and FixTimeStep to zero
#  Constant time stepping:
#    Set CFL to zero, and FixTimeStep to desired value
#
CFL=1.0 
FixTimeStep=0

# Run Settings -----------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


# Directory Settings
runFolder="/home/yh/channel_rough/runs/$runName"
scratchFolder="/home/yh/channel_rough_data/runs/$runName"

# Strings to replace in the hre.dat file that sets the input/output file path
stringToReplace_output="output_filepath_set_by_launchscript"  # in hre.dat, the file output path
stringToReplace_input="input_filepath_set_by_launchscript"
stringToReplace_Re="Re_set_by_launchscript"
stringToReplace_nstep="nstep_set_by_launchscript"
stringToReplace_nimag="nimag_set_by_launchscript"
stringToReplace_nhist="nhist_set_by_launchscript"
stringToReplace_CFL="CFL_set_by_launchscript"
stringToReplace_FixTimeStep="FixTimeStep_set_by_launchscript"

if [[ ! -e $runFolder ]]  # only if run name does not exist yet (avoid data loss)
then
  # Change to the directory this script is located 
  # Since all the directories below are relative paths
  SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
  cd $SCRIPT_DIR
  
  # recompile script
  cd ../build/
  make
  cd ../launch/
  # create run directory where all run files are saved
  mkdir $runFolder
  # create scratch directory where all data are saved
  mkdir $scratchFolder
  # save files so that run can be repeated if needed
  cp ./launch.sh $runFolder
  cp ../build/channel $runFolder
  cp -r ../src $runFolder/src
  cp hre.dat $runFolder
  cd $runFolder

  # set missing file paths in hre.dat
  # which points to the file output directory
  sed -i "s|$stringToReplace_output|$scratchFolder/$runName|g" hre.dat
  sed -i "s|$stringToReplace_input|$inputFile|g" hre.dat
  sed -i "s|$stringToReplace_Re|$Re|g" hre.dat
  sed -i "s|$stringToReplace_nstep|$nstep|g" hre.dat
  sed -i "s|$stringToReplace_nimag|$nimag|g" hre.dat
  sed -i "s|$stringToReplace_nhist|$nhist|g" hre.dat
  sed -i "s|$stringToReplace_CFL|$CFL|g" hre.dat
  sed -i "s|$stringToReplace_FixTimeStep|$FixTimeStep|g" hre.dat

  echo "Submitting job to Richardson Cluster"
  # write a job.sh script for sbatch to submit job
  echo "#!/usr/bin/bash" > job.sh
  echo "#SBATCH --partition=normal" >> job.sh
  echo "#SBATCH --job-name=$runName" >> job.sh
  echo "#SBATCH --output=$runName.txt" >> job.sh
  echo "#SBATCH --mail-type=all" >> job.sh
  echo "#SBATCH --mail-user=yhuang1@caltech.edu" >> job.sh

  echo "#SBATCH -n $mpiProcessors --tasks-per-node=24" >> job.sh
  echo "#SBATCH -t $runTime" >> job.sh
  echo "#SBATCH --mem-per-cpu=5000M" >> job.sh

  if [[ $jobDependency != "none" ]]; then
    echo "#SBATCH --dependency=aftercorr:$jobDependency" >> job.sh
  fi

  echo "echo \"Running on hosts: \$SLURM_NODELIST\"" >> job.sh
  echo "echo \"Running on \$SLURM_NNODES nodes.\"" >> job.sh
  echo "echo \"Running on \$SLURM_NPROCS processors.\"" >> job.sh

  echo "mpirun ./channel" >> job.sh

  # submit job with sbatch
  sbatch ./job.sh

else  # don't do anything if run name already exists
  echo "Runname already exists. Exiting"
fi

