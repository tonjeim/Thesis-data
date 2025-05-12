#!/bin/bash 

#SBATCH --nodes=1 
#SBATCH --ntasks=10 
#SBATCH --mem=90G 
#SBATCH --time=01:00:00 
#SBATCH --partition=RHEL9,smallmem,orion
#SBATCH --job-name=R_analysis  
#SBATCH --output=R_analysis_bdis_bsyl%j.log 

# Navigate to the scratch directory 
cd /mnt/users/tonjeei/RScripts/bdis_bsyl_r
 
# Load R module 
module load R 

# Set the RSTUDIO_PANDOC environment variable 
export RSTUDIO_PANDOC="$HOME/bin/pandoc"

# Render the R Markdown file 
Rscript -e "rmarkdown::render('bdis_bsyl_r.Rmd')"