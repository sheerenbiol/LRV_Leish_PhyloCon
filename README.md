# Code and data files related to co-phylogenetic analysis of Leishmania parasitesthe Leishmania-LRV symbiosis

This repo contains scripts and data files related to the study of co-evolutionary history of *Leishmania* parasites and their LRV viruses.
Analyses are based on bulk RNA seq data from cultivated isolates which were positive for LRV.
RNA seq allowed for the simultaneous extraction of the parasite's transcriptome and the viral dsRNA genome.

Seperate analyses were done for both LRV and *Leishmania* and subsequently a joint evolutionary analysis was conducted using both global-fit and event-based co-phylogenetic methods.


### Overview of the repo:
`$ tree -d .`
```.
├── 1_LRV_analysis
│   ├── 1_Reference_Mapping
│   │   ├── 1_HPC_Scripts
│   │   ├── 2_FASTQs
│   │   └── 3_REFERENCE
│   ├── 2_Viralmetagenome
│   │   └── 1_HPC_Scripts
│   ├── 3_Assembly_QC
│   │   ├── 1_HPC_Scripts
│   │   └── 2_ASSEMBLIES
│   ├── 4_LRV_MSA
│   ├── 5_ML_Phylogenies
│   │   ├── 1_IQtree
│   │   │   ├── LRV_CP_w_REFS
│   │   │   ├── LRV_RDRP_w_REFS
│   │   │   ├── LRV_WGS42
│   │   │   └── LRV_WGS_w_REFS
│   │   ├── 2_Genetic_Distance
│   │   │   └── 1_R_Scripts
│   │   └── 3_Gene_Cophylo
│   │       └── 1_R_Scripts
│   ├── 6_Ancestral_State_Reconstruction
│   └── 7_Maps
├── 2_Leish_analysis
│   ├── 1_Mapping
│   │   └── 1_HPC_Scripts
│   ├── 2_Depths
│   │   ├── 1_HPC_Scripts
│   │   └── 2_Python_Scripts
│   ├── 3_VariantCalling
│   │   └── 1_HPC_Scripts
│   ├── 4_SNPbased_Fasta_and_Frequencies
│   │   ├── 1_VCF2Fasta
│   │   └── 2_VCF2Freq
│   │       ├── FreqHistograms
│   │       └── FreqScatterplots
│   │           ├── IOCL3539
│   │           └── IOCL3574
│   └── 5_SNP_Counts
│       ├── 1_Conversion_012_Format
│       └── 2_R_Scripts
└── 3_Cophylo_analysis
    ├── 1_Preparing_MSA
    │   ├── 1_Subgenus
    │   │   ├── 1_LRV
    │   │   └── 2_LEISH
    │   │       ├── 1_HPC_Scripts
    │   │       ├── 2_Reference_orthologs
    │   │       ├── 3_Ortholog_specific_alignments
    │   │       └── 4_Concatenated_alignment_MLtree
    │   ├── 2_Species
    │   │   ├── 1_Viannia_SNP_alignment
    │   │   └── 2_Viannia_phylo
    │   └── 3_Population
    ├── 3_Genetic_Dist_Regression
    │   └── R_Scripts
    ├── 4_Global_fit_cophylo
    └── 5_Event_based_cophylo
        ├── 1_Subgenus
        ├── 2_Species
        ├── 3_Population
        │   ├── 1_Lbraziliensis
        │   └── 2_Lguyanensis
        └── 4_Cost_Distribution_plots

65 directories, 249 files
```

### Citation:
> **Note** Pre-print on bioRxiv:
> Heeren, S., Motta Cantanhêde, L., Chourabi, K., de Oliveira Santana, M. C., Klaps, J., Kostygov, A. Y., Yurchenko, V., Lemey, P., Dujardin, J.-C., Van den Broeck, F., & Cupolillo, E. (2025). (Genomic evidence for co-evolution and sporadic host shifts in leishmaniaviruses)[https://doi.org/10.1101/2025.11.20.685321]. *bioRxiv*.
