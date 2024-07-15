# cfDNA Variant Origin Classifier

This repository contains the code and analysis from the paper "Machine Learning Approaches to Determine Variant Origin in Liquid Biopsies" by Canzoniero et al.

### Pipeline Overview:
1. **Sample Alignment and Pre-processing**
   - **Tools:** BWA (or Bowtie2), Samtools
   - **Description:** Raw sequencing data is aligned to the reference genome using BWA. The resulting alignments are then sorted and indexed with Samtools. Depending on samples, quality control and trimming may be necessary. 

2. **Somatic Variant Calling**
   - **Tools:** Mutect2, Strelka2
   - **Description:** Any variant caller optimized for detecting low-frequency variants in cfDNA samples is acceptable. 
3. **Identify Reads/Fragments with Somatic Variants**
   - **Tools:** [Python Script](cfDNA_CHIP_MLM/Data_Preprocessing/Validation_Extract_Fragments.py) and R code (need to upload script here!)
   - **Description:** This script will identify reads corresponding to the fragments containing the somatic variants. Before running must configure path to variant list and path to Bam files. 

4. **Calculate Summary Statistics for Each Position**
   - **Tools:** Script [here](cfDNA_CHIP_MLM/Data_Preprocessing/CHIP_EGA_fragmentLevel_sumStats.R) 
   - **Description:** Script to calculate various statistical tests on the variant-containing fragments and the other fragments at that positions.  

5. **Run Model to Predict Origin of Variant**
   - **Tools:** [Prediction Model](cfDNA_CHIP_MLM/Model/)
   - **Description:** R Code and weights needed to generate the origin prediction based on the calculations in the previous step. 

## Data

### Training and Testing 
The cohorts used for the initial training and testing are associated with the following studies:

1. Parikh AR, Leshchiner I, Elagina L, et al. "Liquid versus tissue biopsy for detecting acquired resistance and tumor heterogeneity in gastrointestinal cancers." *Nat Med*. 2019;25(9):1415-1421. doi: [10.1038/s41591-019-0561-9](https://doi.org/10.1038/s41591-019-0561-9)
   
2. Sivapalan L, Murray JC, Canzoniero JV, et al. "Liquid biopsy approaches to capture tumor evolution and clinical outcomes during cancer immunotherapy." *J Immunother Cancer*. 2023;11(1). doi: [10.1136/jitc-2022-005924](https://doi.org/10.1136/jitc-2022-005924). PMID: 36657818; PMCID: PMC9853269.
   
3. Liu J, Chen X, Wang J, et al. "Biological background of the genomic variations of cf-DNA in healthy individuals." *Ann Oncol*. 2019;30(3):464-470. doi: [S0923-7534(19)31073-7](https://doi.org/10.1093/annonc/mdz041)
   
4. Barbany G, Arthur C, Liedén A, et al. "Cell-free tumour DNA testing for early detection of cancer – a potential future tool." *J Intern Med*. 2019;286:118–136. doi: [10.1111/joim.12958](https://doi.org/10.1111/joim.12958)
   
5. Hu Z, Chen H, Long Y, et al. "The main sources of circulating cell-free DNA: Apoptosis, necrosis and active secretion." *Crit Rev Oncol Hematol*. 2021;157:103166. doi: [10.1016/j.critrevonc.2020.103166](https://doi.org/10.1016/j.critrevonc.2020.103166)
   
6. Brannon AR, Jayakumaran G, Diosdado M, et al. "Enhanced specificity of clinical high-sensitivity tumor mutation profiling in cell-free DNA via paired normal sequencing using MSK-ACCESS." *Nat Commun*. 2021;12:3770. doi: [10.1038/s41467-021-24109-5](https://doi.org/10.1038/s41467-021-24109-5)

### Independent Validation 
Data used for the validation cohort was acquired from MSK-IMPACT trial and can be found here: 

Razavi P, Li BT, Brown DN, et al. High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. Nat Med. 2019;25(12):1928-1937. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7061455/. Accessed Jun 24, 2023. doi: 10.1038/s41591-019-0652-7


## R Dependencies
The following R packages are required for the analysis:

- `tidyr`
- `dplyr`
- `magrittr`
- `caret`
- `naivebayes`
- `xgboost`


## Annotation 

Variant level annotation data was accessed via [Open-Cravat](https://github.com/KarchinLab/open-cravat)

Cancer hotspots were identified with [COSMIC](https://cancer.sanger.ac.uk/cosmic)


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact
For any questions or further information, please contact Jenna Canzoniero at [jcanzon1@jhmi.edu](mailto:jcanzon1@jhmi.edu).
