#!/usr/bin/env bash
runName=$1  # name of the run to be removed
hostName=$(hostname | cut -b 1-8 -n)  # host system is identified by first 8 characters of hostname
if [[ "$hostName" == "eu-login" ]]; then
  runFolder="$HOME/channel_rough/runs/$runName"
  scratchFolder="$SCRATCH/$runName"
  [[ ! -z "$runName" && -d "$runFolder" ]] && rm -rf $runFolder $scratchFolder || echo "Wrong runname. Exiting."
# add other architectures as needed: elif ...
# For Richardson:
elif [[ "$hostName" == "richards" ]]; then
  runFolder="/home/yh/channel_rough/runs/$runName"
  scratchFolder="/scratch/yh/channel_rough_data/runs/$runName"
  [[ ! -z "$runName" && -d "$runFolder" ]] && rm -rf $runFolder $scratchFolder || echo "Wrong runname. Exiting."
fi

