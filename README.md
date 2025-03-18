# Cistanthe-RNASeq

Data, scripts, and outputs from Chomentowska et al (in press): "High-quality genome of the Atacama Desert plant Cistanthe cachinalensis and its photosynthetic behavior related to drought and life history."

Raw RNA sequencing data will be available at NCBI under BioProject PRJNA1181828.


### In this repository, you will find:

- Pipeline to process RNASeq reads: Cistanthe_RNAseq_analysis.sh

- RNASeq Metadata: RNA_metadata.xlsx

- Snakemake files for trimming and QC to work within the above pipeline: trim_workflow/

- Raw and normalized gene counts, list of photosynthetic orthologs, and R scripts for differential gene expression analyses + TPM calculations: DGE_analysis/

- Data and R script for visualizing gas exchange: gas_exchange_measurements/

- Data and R script for visualizing titratable acidity: titratable_acidity/

- Output from XSTREME analyses (https://meme-suite.org/meme/doc/xstreme.html?man_type=web): XSTREME_outputs/
