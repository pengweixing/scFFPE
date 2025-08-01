This repository is for single cell FFPE project


### software requrements
    - Python 3.9.16
    - computeMatrix 3.5.1
    - bamCoverage 3.5.1
    - plotHeatmap 3.5.1
    - Rscript (R) version 4.3.1 
    - plotProfile 3.5.1
    - sambamba 0.7.0
    - BWA Version: 0.7.17-r1188
    - samtools Version: 1.17
    - macs2 2.1.2
    - snakename 7.18.2
    
### package requrements
    python:
        pysam = '0.22.0'
        fuzzywuzzy = '0.18.0'
        Levenshtein = '0.20.8'
        pandas = '1.4.0'
        matplotlib = '3.8.2'
    R version 4.1.3
        ggplot2 = '3.3.6'
        tidyverse = '1.3.2'
        
### installation
        
      - download the pipeline to your owndirectory  
      - modify config.hg38.yaml and config.mm10.yaml
      - install all required softwares and packages
      
### usage
    
    Please check the script in the example folder
    
#### If you want to extract cells based on FRiP or FRiT mannully
    please run the `extract_cell.sh` script


#### If you want to get the qc file for all cells
    please run `get_single_cell_qc.sh` script

    | barcode                   | total | Dups  | LowQ | Unpaired | Final | FRiP  | FRiT  |
    |---------------------------|-------|-------|------|----------|-------|-------|-------|
    | TGTCGTCTTGGTATG_Tn5_17    | 3     | 0.00  | 0    | 1        | 2     | 0.00  | 0.00  |
    | CCCTTCCTTTGGATA_Tn5_35    | 98    | 8.67  | 13   | 2        | 80    | 13.75 | 20.00 |
    | CCCGAACGTCGGTCA_Tn5_22    | 264   | 12.12 | 36   | 3        | 210   | 14.52 | 17.86 |
    | CTATTCTCTCGTCAC_Tn5_20    | 414   | 13.04 | 77   | 8        | 313   | 13.74 | 13.74 |
    | TAGTCTGGGATAAGA_Tn5_31    | 334   | 13.77 | 61   | 3        | 254   | 12.80 | 16.54 |
    | TCTGACCCTTTCGGT_Tn5_27    | 551   | 12.89 | 136  | 8        | 404   | 10.52 | 13.24 |
    | ATTATTCACCCAAGT_Tn5_8     | 2     | 0.00  | 0    | 0        | 2     | 0.00  | 0.00  |
    | ATAAGGTAAATCACC_Tn5_21    | 2725  | 10.42 | 695  | 32       | 2061  | 6.91  | 11.52 |
    | AATCTAACACGAGGT_Tn5_1     | 1266  | 12.56 | 207  | 16       | 987   | 10.89 | 15.50 |
