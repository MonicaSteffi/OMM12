#!/bin/bash
#SBATCH --output=pul_marbouty.out       # Standard output and error log
#SBATCH --error=error_pul_marbouty.err            # Standard output and error log
#SBATCH --job-name=pul_marbouty                 # Job name
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --nodes=1                               # Number of nodes for this job   
#SBATCH --cpus-per-task=10                     # Number of CPU cores per task
#SBATCH --time=48:00:00                         # Time limit hrs:min:sec
#SBATCH --mail-type=END                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Durai@mvp.lmu.de     # Where to send mail   
#SBATCH --export=NONE

source ~/anaconda3/etc/profile.d/conda.sh
conda activate kofamscan
# conda install pfam_scan


# Downloaded files from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/
# wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.dat.gz
# wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/active_site.dat.gz
# wget wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
# wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/database_files/pfamA_HMM.txt.gz

export PATH=/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/Tools/hmmer-3.3.2/bin/:${PATH}

pfamDir=/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/Tools/pfam35.0/

#### hmmpress
## cd ${pfamDir}
## hmmpress Pfam-A.hmm

cd /dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/Marbouty_annotations/PULs_new/
prtfaa=$1
 spNm=($(echo $prtfaa | sed 's|.*Marbouty_annotations_proteins_integrated/||g' | sed 's|_proteins_integrated.faa||g'))
 rm -f ${prtfaa}Edited.faa
 cat ${prtfaa} | sed 's|-||g' | sed 's|*||g' > ${prtfaa}Edited.faa
 rm -f ${spNm}.pfam
 pfam_scan.pl -fasta ${prtfaa}Edited.faa -outfile ${spNm}.pfam  -as -cpu 50 -dir ${pfamDir}

