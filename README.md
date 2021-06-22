# Heterotrophic Bacteria Dominate Catalase Expression During Microcystis Blooms

This repository contains all of the R analysis code and input data used to generate figures for the manuscript "Heterotrophic Bacteria Dominate Catalase Expression During Microcystis Blooms."  

All perl scripts mentioned in README documents or methods sections are available at this GitHub repositroy: /Geo-omics/scripts  

Description of README files:  
Describes metagenomic assembly and binning workflow on /geomicro computation servers: Erie2014_coAssembly_readme.txt   
Describes ROS read mapping on /geomicro computation servers: WLE_2014_ROS_gene_blasts_readme.txt  
Describes Microcystis mcy and 16S rRNA read mapping and non-specific mapping troubleshooting on /geomicro computation servers: Mcy_BLAST_README.txt  

Description of Mardownfiles:  
Catalase_Peroxidase_Mapping_Analysis.Rmd contains code used to generate plots of gene abundances in metagenomes and metatranscriptomes from Lake Erie.  
Pyruvate_Light_Experiment.Rmd contains code and analysis for data collected from cultivation experiments with Microcystis isolates.  
Mcy_mapped_reads_vs_URDB.Rmd contains code and analysis used to evaluate Microcystis mcyD and 16S rRNA mapping cutoffs and databases.  

The metaG and metaT folders contain the results of ROS read mapping used in the Catalase_Peroxidase_Mapping_Analysis markdown file.  
Raw chemiluminescent data from Felume H2O2 measurements are located in Raw_Felume_Files.  
