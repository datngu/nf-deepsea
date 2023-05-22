#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=DS-cattle   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm

git pull

#cd /mnt/SCRATCH/ngda

# git clone https://github.com/datngu/nf-deepsea.git


# removing chr

# for fi in /mnt/SCRATCH/ngda/cattle/data/*
# do
#     cat $fi | sed 's/chr//g' > ${fi}.removed_chr.bed
# done


genome='/mnt/users/ngda/genomes/cattle/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa'

export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity

# nextflow_res_dir=/mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/results/${tis}
# nextflow_trace_dir=/mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/results/trace_dir_${tis}
# nextflow_work_dir=/mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/work_dir/${tis}


nextflow run main.nf -resume -w work_dir \
    --genome ${genome} \
    --chrom 29 \
    --window 200 \
    --seqlen 1000 \
    --peaks '/mnt/SCRATCH/ngda/cattle/data/*.removed_chr.bed'