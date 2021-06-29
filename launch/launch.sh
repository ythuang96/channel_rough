#!/usr/bin/bash
# Launch script to launch DNS code on Richardson Cluster

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Run Settings -----------------------------------------------------------------------------------
runName=smooth_1  # must not exist yet in runFolder
runTime=22:00:00  # format: hh:mm (Euler), hh:mm:ss (Millikan) hh:mm:ss (Richardson)
#jobDependency=${3:-none}  # name of run (Euler), jobid of run (Millikan) for dependency condition, optional
mpiProcessors=96  # must be equal to nprocs in ctes3D

inputFile="/scratch/yh/channel_rough_data/runs/smooth_init9/smooth_init9.006"

# Run Paramters
# Reynolds number
Re=113000
# total time steps, must be multiple of 500 + 1
nstep=10001
# time step to update CFL and write to .cf
# no larger thatn 10
nhist=5
# CFL condition, no larger than 1.5
CFL=1.0 

# Run Settings -----------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


# Directory Settings
runFolder="/home/yh/channel_rough/runs/$runName"
scratchFolder="/scratch/yh/channel_rough_data/runs/$runName"

# Strings to replace in the hre.dat file that sets the input/output file path
stringToReplace_output="output_filepath_set_by_launchscript"  # in hre.dat, the file output path
stringToReplace_input="input_filepath_set_by_launchscript"
stringToReplace_Re="Re_set_by_launchscript"
stringToReplace_nstep="nstep_set_by_launchscript"
stringToReplace_nhist="nhist_set_by_launchscript"
stringToReplace_CFL="CFL_set_by_launchscript"

if [[ ! -e $runFolder ]]  # only if run name does not exist yet (avoid data loss)
then
  # recompile script
  cd ../build/
  make
  cd ../launch/
  # create run directory where all run files are saved
  mkdir $runFolder
  # create scratch directory where all data are saved
  mkdir $scratchFolder
  # save files so that run can be repeated if needed
  cp ../build/channel $runFolder
  cp ../src/ctes3D $runFolder
  cp hre.dat $runFolder
  cd $runFolder

  # set missing file paths in hre.dat
  # which points to the file output directory
  sed -i "s|$stringToReplace_output|$scratchFolder/$runName|g" hre.dat
  sed -i "s|$stringToReplace_input|$inputFile|g" hre.dat
  sed -i "s|$stringToReplace_Re|$Re|g" hre.dat
  sed -i "s|$stringToReplace_nstep|$nstep|g" hre.dat
  sed -i "s|$stringToReplace_nhist|$nhist|g" hre.dat
  sed -i "s|$stringToReplace_CFL|$CFL|g" hre.dat

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
  echo "#SBATCH --mem-per-cpu=5000" >> job.sh

  echo "echo \"Running on hosts: \$SLURM_NODELIST\"" >> job.sh
  echo "echo \"Running on \$SLURM_NNODES nodes.\"" >> job.sh
  echo "echo \"Running on \$SLURM_NPROCS processors.\"" >> job.sh

  echo "mpirun ./channel" >> job.sh

  # submit job with sbatch
  sbatch ./job.sh

else  # don't do anything if run name already exists
  echo "Runname already exists. Exiting"
fi

