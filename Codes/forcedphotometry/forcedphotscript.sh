#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=7GB
#SBATCH --job-name=Sumedha

#SBATCH --account=b19-meerlicht-ag
#SBATCH --partition=Main
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=p.vreeswijk@astro.ru.nl
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "date              = $(date)"
echo "hostname          = $(hostname -s)"
echo "working directory = $(pwd)"

echo "#task allocated   = $SLURM_NTASKS"
echo "#cores/task       = $SLURM_CPUS_PER_TASK"

#singularity exec --env PYTHONPATH="/users/pmv/BBtest:\$PYTHONPATH" /idia/projects/meerlicht/Containers/ML_latest.sif python ~/BBtest/force_phot.py ~/4Sumedha/all_gaia_data_mjd_pmtest.dat ~/4Sumedha/all_gaia_data_mjd_MLforcephot_30min_3sigma_20230313_pmtest.fits --radecs_file_format csv --ra_col RA --dec_col Dec --date_col ObsTime --date_format mjd --dtime_max 0.5 --thumbnails True --logfile log_Sumedha_20230313 --fullsource True --input_cols2copy Name,SourceID,ObsTime,RA,Dec,GMag,GMagErr,PMRA,PMDEC --include_fluxes True
singularity exec /idia/projects/meerlicht/Containers/ML_latest.sif python ~/BBtest/force_phot.py ~/4Sumedha/all_gaia_data_mjd.dat ~/4Sumedha/all_gaia_data_mjd_MLforcephot_notimelimit_3sigma_20230525.fits --radecs_file_format csv --ra_col RA --dec_col Dec --thumbnails True --logfile log_Sumedha_20230525 --fullsource True --input_cols2copy Name,SourceID,ObsTime,RA,Dec,GMag,GMagErr --include_fluxes True
