#!/bin/bash
## Cistanthe RNA seq project
## Anri Chomentowska

###### Download illumina data with sbatch slurm_getdata_rnaseq.sh
```bash
#!/bin/bash
#SBATCH --partition transfer
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --job-name download_rnaseq
#SBATCH --output download_rnaseq_output-%j.txt
#SBATCH --error download_rnaseq_error-%j.txt
#SBATCH --mail-type ALL

cd /SAY/standard2/eje5-CC0522-FASEEB-2/shared/montiaceae/cistanthe_rnaseq

## get the data
wget -r -nH --cut-dirs=2 --no-parent --reject "index.html*" http://fcb.ycga.yale.edu:3010/3hDWQCO0dTIVXW7WPQPhGngs6kqupnn/sample_dir_000014407/

```


###### Merge files from sample that were reran:
srun --pty -t 20:00 -p devel bash

cat *_R1_001.fastq.gz > L3I10PM_R1.fastq.gz
cat *_R2_001.fastq.gz > L3I10PM_R2.fastq.gz


###### Run fastQC, then Trimmomatic to trim, then fastQC again

### fastqc
sbatch slurm_fastqc.sh
```
#!/bin/bash
#SBATCH --job-name fastqc
#SBATCH --output fastqc_output-%j.log
#SBATCH --error fastqc_error-%j.txt
#SBATCH --time 24:00:00
#SBATCH --partition pi_edwards,day
#SBATCH --nodes 1                    # number of cores and nodes
#SBATCH --cpus-per-task 12           # number of cores
#SBATCH --mem-per-cpu 4G             # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/palmer_scratch/rnaseq

# Set up modules
module purge # Unload any existing modules that might conflict
module load FastQC

fastqc -t 12 --outdir ./fastqc_initial ./merged_reads/*/*.fastq.gz

```

### MultiQC
sbatch slurm_multiqc.sh
```
#!/bin/bash
#SBATCH --job-name multiqc
#SBATCH --output multiqc_output-%j.log
#SBATCH --error multiqc_error-%j.txt
#SBATCH --time 24:00:00
#SBATCH --partition pi_edwards,day
#SBATCH --nodes 1                    # number of cores and nodes
#SBATCH --cpus-per-task 1            # number of cores
#SBATCH --mem-per-cpu 8G             # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/palmer_scratch/rnaseq/fastqc_initial

# Set up modules
module purge # Unload any existing modules that might conflict
module load MultiQC

multiqc .

```

### Now, trim samples with Trimmomatic using Snakemake
#No need to use Snakemake, but just wanted to learn how to use it
#Snakefile (Snakefile_rnaseq_trim) and config file uploaded in trim_workflow
srun --pty -t 30:00 -p devel bash

module load miniconda
conda create -n snakemake_env -c conda-forge -c bioconda snakemake

conda activate snakemake_env

cd project/rnaseq/trim_workflow

## In trim_workflow/ directory, put the config.yaml and Snakefile_rnaseq_trim files:
# Dry run, make sure everything will run correctly
snakemake -n -p --cores 1

# Submit following batch script
sbatch slurm_snakemake_rnaseq.sh
```
#!/bin/bash
#SBATCH --job-name snakemake
#SBATCH --output snakemake_output-%j.log
#SBATCH --error snakemake_error-%j.txt
#SBATCH --requeue
#SBATCH --time 72:00:00
#SBATCH --partition pi_edwards,week
#SBATCH --nodes 1                    # number of cores and nodes
#SBATCH --cpus-per-task 6            # number of cores
#SBATCH --mem-per-cpu 4G             # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/project/rnaseq/trim_workflow

# Set up modules
module purge # Unload any existing modules that might conflict
module load Trimmomatic
module load miniconda

conda activate snakemake_env

snakemake --scheduler greedy --verbose --rerun-incomplete --cores $SLURM_CPUS_PER_TASK --latency-wait 90
```

# Run FastQC again for trimmed reads in ./trimmed_reads/*.fastq.gz, and multiqc again


###### Run STAR
# Make sure we have our genome (cislon_corrected.fasta) and annotation (cislon_genes.gff3) files
# Uploaded in github.com/anriiam/Cistanthe-cachinalensis-genome

### Convert our gff3 file to gtf
ssh ac2767@mccleary.ycrc.yale.edu
srun --pty -t 1:00:00 -p devel bash

module load miniconda
conda create -n gffread_env

conda activate gffread_env
conda install gffread -c conda-forge -c bioconda

gffread cislon_genes.gff3 -T -o cislon_genes.gtf

### Index genome
sbatch slurm_star_index.sh
```
#!/bin/bash
#SBATCH --job-name STAR_index
#SBATCH --output STAR_index_output-%j.log
#SBATCH --error STAR_index_error-%j.txt
#SBATCH --time 24:00:00
#SBATCH --partition pi_edwards,day
#SBATCH --nodes 1                    # number of cores and nodes
#SBATCH --cpus-per-task 4            # number of cores
#SBATCH --mem-per-cpu 16G            # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/palmer_scratch/rnaseq/

module purge
module load STAR/2.7.11a-GCC-12.2.0

STAR --runMode genomeGenerate \
--runThreadN 4 \
--genomeDir ./star_index_new \
--genomeFastaFiles ./star_ref_new/cislon_corrected.fasta \
--sjdbGTFfile ./star_ref_new/cislon_genes.gtf \
--sjdbOverhang 149

```

### Save the following script as run_star.sh and make it executable, then run it on our trimmed Reads

## run_star.sh script
```
#!/bin/bash

# Path to the sample list file
SAMPLE_LIST="./rna_samples.txt"

# Path to the genome index directory
GENOME_DIR="./star_index_new"

# Path to the input and output directories
INPUT_DIR="./trimmed_reads"
OUTPUT_DIR="./star_results_new"

# Number of threads
THREADS=8

# Maximum number of mismatches
MISMATCHES=2

# Loop through each sample in the sample list
while IFS= read -r sample; do
    echo "Processing sample: $sample"

    STAR --genomeDir $GENOME_DIR \
        --readFilesCommand zcat \
        --readFilesIn $INPUT_DIR/${sample}_R1_trimmed.fastq.gz $INPUT_DIR/${sample}_R2_trimmed.fastq.gz \
        --outFilterMismatchNmax $MISMATCHES \
        --runThreadN $THREADS \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix $OUTPUT_DIR/${sample}

    echo "Finished processing sample: $sample"
done < "$SAMPLE_LIST"
```
chmod +x run_star.sh

## Run script
sbatch slurm_star_map.sh
```
#!/bin/bash
#SBATCH --job-name STAR_map
#SBATCH --output STAR_map_output-%j.log
#SBATCH --error STAR_map_error-%j.txt
#SBATCH --time 96:00:00
#SBATCH --partition pi_edwards,week
#SBATCH --ntasks 1                   # Number of tasks (1 since we're running a single script)
#SBATCH --cpus-per-task 8            # number of cores
#SBATCH --mem-per-cpu 16G            # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/palmer_scratch/rnaseq/

module purge
module load STAR/2.7.11a-GCC-12.2.0

# Run the run_star.sh script
./run_star.sh
```


###### Count features with featureCounts

### Save the following batch script and submit
sbatch slurm_featurecounts.sh
```
#!/bin/bash
#SBATCH --job-name featurecounts
#SBATCH --output featurecounts_output-%j.log
#SBATCH --error featurecounts_error-%j.txt
#SBATCH --time 96:00:00
#SBATCH --partition pi_edwards,week
#SBATCH --ntasks 1                   # Number of tasks (1 since we're running a single script)
#SBATCH --cpus-per-task 8            # number of cores
#SBATCH --mem-per-cpu 16G            # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/palmer_scratch/rnaseq/

module load Subread/2.0.3-GCC-10.2.0

# Store list of files as a variable
dirlist=$(ls -t ./star_results_new/*.bam | tr '\n' ' ')
echo $dirlist

# paired-end reads while counting multimapping reads fractionally with featureCount
featureCounts -p -g gene_id -s 2 -T 8 -M --fraction -F GTF -a ./star_ref_new/cislon_genes.gtf -o ./feature_counts_new/final_counts_stranded_geneid_NEW.txt $dirlist

# if you want to work with transcripts rather than genes
#featureCounts -p -g transcript_id -s 2 -T 8 -M --fraction -F GTF -a ./star_ref_new/cislon_genes.gtf -o ./feature_counts_new/final_counts_stranded_transcript_NEW.txt $dirlist
```

###### MultiQC at the end

### Now on all the previous folders
```
#!/bin/bash
#SBATCH --job-name multiqc
#SBATCH --output multiqc_output-%j.log
#SBATCH --error multiqc_error-%j.txt
#SBATCH --time 2:00:00
#SBATCH --partition pi_edwards,day
#SBATCH --nodes 1                    # number of cores and nodes
#SBATCH --cpus-per-task 1            # number of cores
#SBATCH --mem-per-cpu 8G             # shared memory, scaling with CPU request
#SBATCH --mail-type ALL
#SBATCH --chdir /home/ac2767/palmer_scratch/rnaseq

# Set up modules
module purge # Unload any existing modules that might conflict
module load MultiQC

multiqc ./fastqc_trimmed ./feature_counts_new ./star_results_new -n multiqc_new
```

###fin
