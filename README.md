# cfDNA_CHIP_MLM
This repository contains the code and data used for the cfDNA CHIP detection paper by Canzoniero et al.

### Pipeline Overview:
1. Sample Alignment and pre-processing
2. Somatic variant calling
3. Identify reads/fragments with somatic variants
4. Calculate summary statistics for each position
5. Run model to predict the origin of the variant

### Data
The JHU-based cohorts used for the initial training and testing can be found here:

Data used for the validation cohort was acquired from [Ravazi et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7061455/)

Variant level annotation data was accessed via [Open-Cravat](https://github.com/KarchinLab/open-cravat)
