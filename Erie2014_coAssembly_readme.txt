#DJS 5/10/18:

#Run QC of Illumina reads for all samples to be considered using bbmap pipeline:
nohup sh bbmap_qc.sh &> qc.log & ##Bash script is annotated to provide info on what each step does

#DJS 5/11/18:

#Convert interleaved QC reads to fasta format for stations 53600 and 49615 (working directories are the respective sample folders):
nohup reformat in=Sample_49615_adtrim_clean_qtrim_derep_int.fastq out=Sample_49615_adtrim_clean_qtrim_derep_int.fasta &> fastq_to_fasta.log &
nohup reformat in=Sample_53600_adtrim_clean_qtrim_derep_int.fastq out=Sample_53600_adtrim_clean_qtrim_derep_int.fasta &> fastq_to_fasta.log &

#Check that adapters found by bbmerge were removed in Sample 49615 and 53600 with blast (wd= repsective sample directories):
nohup blastn -query Sample_49615_adtrim_clean_qtrim_derep_int.fasta -db /omics/HABs/Microcystis/Downsampling/Finished_Run_1345_1351/Sample_49615/bbmerge_adapters.fa -outfmt "7 std qlen slen qcovs" -out bbmerge_adapters_v_Sample_49615_qc_reads.blastn &> bbmerge_adapters_v_Sample_49615_qc_reads.log &
nohup blastn -query Sample_53600_adtrim_clean_qtrim_derep_int.fasta -db /omics/HABs/Microcystis/Lake_Erie_Meta-omes/Run_1423/ruddle.brcf.med.umich.edu/Run_1423/denef/Finished_Assemblies/Full_Assemblies/Sample_53600/bbmerge_adapters.fa -outfmt "7 std qlen slen qcovs" -out bbmerge_adapters_v_Sample_53600_qc_reads.blastn &> bbmerge_adapters_v_Sample_53600_qc_reads.log &

#DJS 5/13/18:

#Remove nonsignificant blast hits (bit score must be at least 50 and e-value at most 1e-5; wd= respective sample directories):
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b bbmerge_adapters_v_Sample_53600_qc_reads.blastn -o bbmerge_adapters_v_Sample_53600_qc_reads.blastn.filtered
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b bbmerge_adapters_v_Sample_49615_qc_reads.blastn -o bbmerge_adapters_v_Sample_49615_qc_reads.blastn.filtered
#Only one hit found in Sample_49615 at 93% identity.

#Fastqc report of quality checked interleaved reads:
nohup sh batch_fastqc.sh &> batch_fastqc.log & ##Bash script is annotated to provide info on what each step does

#Fastqc report was flagging a sequence duplication warning in some samples.
#10-50% of the reads were often present in 10-50 copies
#Unlikely to be due to uneven trimming of mate pairs, which only accounted for 6% of the sample.
#Likely due to running analysis on interleaved reads (fwd reads could have matched with reverse reads from another pair)
#Another possibility: over sequencing over library or biological redundancy (not high enough to be serious technical duplication)

#Reran fastqc on the forward and reverse read files separately:
#Rewrote old batch_fastqc.sh script
nohup sh batch_fastqc.sh &> batch_fastqc.log &

#DJS 5/16/18:

#Try Dereplicating reads using omics.derep script:
nohup sh bbmap_qc.sh &> qc.omics_derep.log &

#DJS 5/17/18:
#Run fastqc on the output of the pipeline using omics derep:
nohup sh batch_fastqc_omics_derep.sh &> fastqc_omics_derep.log &

#Dereplication script written by Robert seemed to out perform the dedupe script.
#Sequencing length distributions suggest very few reads are keeping full read length, trimming too agressively?
#Rewrote script to trim to Q10 and min kmer on adapter trim to 8 run the fastqc automatically:
nohup sh bbmap_qc.sh &> qc.log &

#DJS 5/20/18:

#Trimming to Q10 had some reads approaching poor quality in the last few bases, bump to Q15:
nohup sh bbmap_qc.sh &> qc.log &

#Some samples have a fair amount of the library present between 10-50 copies according to fastqc report.
#Fastqc report only considers the first 50 bases of the reads when checking for duplicates, and there were no sequences flagged as overrepresesented, so continue with these QC'd reads for the assembly of one sample and check the read mapping to the assembly for areas that could be impacted by duplicates later.

#DJS 5/24/18:

#Check that all adapters were removed with new trimming settings:
#Convert interleaved QC reads to fasta format for stations 53600 and 49615 (working directories are the respective sample folders):
nohup reformat in=Sample_49615_adtrim_clean_qtrim_derep_int.fastq out=Sample_49615_adtrim_clean_qtrim_derep_int.fasta &> fastq_to_fasta.log &
nohup reformat in=Sample_53600_adtrim_clean_qtrim_derep_int.fastq out=Sample_53600_adtrim_clean_qtrim_derep_int.fasta &> fastq_to_fasta.log &

nohup reformat in=Sample_49615_adtrim_clean_qtrim_fwd.derep.fastq out=Sample_49615_adtrim_clean_qtrim_fwd.derep.fasta &> fwd_fastq_to_fasta.log &
nohup reformat in=Sample_49615_adtrim_clean_qtrim_rev.derep.fastq out=Sample_49615_adtrim_clean_qtrim_rev.derep.fasta &> rev_fastq_to_fasta.log &

#Check that adapters found by bbmerge were removed in Sample 49615 and 53600 with blast (wd= repsective sample directories):
nohup blastn -query Sample_49615_adtrim_clean_qtrim_fwd.derep.fasta -db /omics/HABs/Microcystis/Downsampling/Finished_Run_1345_1351/Sample_49615/bbmerge_adapters.fa -outfmt "7 std qlen slen qcovs" -out bbmerge_adapters_v_Sample_49615_fwd_reads.blastn &> bbmerge_adapters_v_Sample_49615_fwd_reads.log &
nohup blastn -query Sample_49615_adtrim_clean_qtrim_rev.derep.fasta -db /omics/HABs/Microcystis/Downsampling/Finished_Run_1345_1351/Sample_49615/bbmerge_adapters.fa -outfmt "7 std qlen slen qcovs" -out bbmerge_adapters_v_Sample_49615_rev_reads.blastn &> bbmerge_adapters_v_Sample_49615_rev_reads.log &
nohup blastn -query Sample_53600_adtrim_clean_qtrim_fwd.derep.fasta -db /omics/HABs/Microcystis/Lake_Erie_Meta-omes/Run_1423/ruddle.brcf.med.umich.edu/Run_1423/denef/Finished_Assemblies/Full_Assemblies/Sample_53600/bbmerge_adapters.fa -outfmt "7 std qlen slen qcovs" -out bbmerge_adapters_v_Sample_53600_fwd_reads.blastn &> bbmerge_adapters_v_Sample_53600_fwd_reads.log &
nohup blastn -query Sample_53600_adtrim_clean_qtrim_rev.derep.fasta -db /omics/HABs/Microcystis/Lake_Erie_Meta-omes/Run_1423/ruddle.brcf.med.umich.edu/Run_1423/denef/Finished_Assemblies/Full_Assemblies/Sample_53600/bbmerge_adapters.fa -outfmt "7 std qlen slen qcovs" -out bbmerge_adapters_v_Sample_53600_rev_reads.blastn &> bbmerge_adapters_v_Sample_53600_rev_reads.log &

#Reassembled Sample 42896 using megahit (previously IDBA-UD) in an attempted to improve Bryobacter bin.
#IDBA-UD assembly and ESOM binning yielded genome with 70% completeness and 5% contamination.
#Megahit command, run in comics container:
#working directory = Sample_42896/
comics -- nohup megahit -1 Sample_42896_adtrim_clean_qtrim_fwd.derep.fastq -2 Sample_42896_adtrim_clean_qtrim_rev.derep.fastq --k-min 21 --k-max 141 --k-step 12 -t 30 -o Sample_42896_megahit_assembly --verbose &> Sample_42896_megahit_assembly.log &

#DJS 5/27/18:
#Filter blast output to only include blast hits with bit score >50 and e-value 1e-5:
#wd: respective Sample directories
nohup perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b bbmerge_adapters_v_Sample_49615_fwd_reads.blastn -o bbmerge_adapters_v_Sample_49615_fwd_reads.blastn.filtered &
nohup perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b bbmerge_adapters_v_Sample_49615_rev_reads.blastn -o bbmerge_adapters_v_Sample_49615_rev_reads.blastn.filtered &
nohup perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b bbmerge_adapters_v_Sample_53600_fwd_reads.blastn -o bbmerge_adapters_v_Sample_53600_fwd_reads.blastn.filtered &
nohup perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b bbmerge_adapters_v_Sample_53600_rev_reads.blastn -o bbmerge_adapters_v_Sample_53600_rev_reads.blastn.filtered &

#Only one significant hit in the middle of a read for one sample 49615, likely not an adapter.

#Coassembly of samples that had bryobacter bins:
comics -- nohup megahit -1 Sample_42896/Sample_42896_adtrim_clean_qtrim_fwd.derep.fastq,Sample_53600/Sample_53600_adtrim_clean_qtrim_fwd.derep.fastq -2 Sample_42896/Sample_42896_adtrim_clean_qtrim_rev.derep.fastq,Sample_53600/Sample_53600_adtrim_clean_qtrim_rev.derep.fastq --k-min 21 --k-max 141 --k-step 12 -t 30 -o Sample_42896_53600_megahit_coassembly --verbose &> Sample_42896_53600_megahit_coassembly.log &

#Indexed and chopped Sample 42896 assembly:
#wd: Sample_42896/
comics -- nohup omics mapping -a Sample_42896_megahit_assembly/final.contigs.fa --chop --index-only &> Sample_42896_megahit_assembly_index.log &

#move chopped contigs to the assembly directory:
mv final.contigs.chop.fa Sample_42896_megahit_assembly/

#Mapped short reads to Sample 42896 assembly:
#wd: Sample_42896/
comics -- nohup omics mapping --index-dir Sample_42896_assembly_bowtie2-index/ -a Sample_42896_megahit_assembly/final.contigs.chop.fa -f Sample_42896_adtrim_clean_qtrim_fwd.derep.fastq  -r Sample_42896_adtrim_clean_qtrim_rev.derep.fastq --cpus 30 -o Mapping &> Sample_42896_mapping.log &

#Indexed Sample 42896+53600 coassembly:
#wd:Sample_42896_53600_megahit_coassembly/
comics -- nohup omics mapping -a final.contigs.fa --chop --index-only &> Sample_42896_53600_megahit_coassembly_index.log &

#Mapped short reads to Sample 42896-53600 coassembly:
comics -- nohup sh Mapping_sample_42896_53600_mapping.sh &> Mapping_sample_42896_53600_coassembly.log &

#Coassembled all samples that contained Bryobacter OTUs in the corresponding 16S tag sequence samples:
comics -- nohup megahit -1 Sample_42896/Sample_42896_adtrim_clean_qtrim_fwd.derep.fastq,Sample_49615/Sample_49615_adtrim_clean_qtrim_fwd.derep.fastq,Sample_49622/Sample_49622_adtrim_clean_qtrim_fwd.derep.fastq,Sample_49625/Sample_49625_adtrim_clean_qtrim_fwd.derep.fastq,Sample_49629/Sample_49629_adtrim_clean_qtrim_fwd.derep.fastq,Sample_53599/Sample_53599_adtrim_clean_qtrim_fwd.derep.fastq,Sample_53600/Sample_53600_adtrim_clean_qtrim_fwd.derep.fastq,Sample_53601/Sample_53601_adtrim_clean_qtrim_fwd.derep.fastq -2 Sample_42896/Sample_42896_adtrim_clean_qtrim_rev.derep.fastq,Sample_49615/Sample_49615_adtrim_clean_qtrim_rev.derep.fastq,Sample_49622/Sample_49622_adtrim_clean_qtrim_rev.derep.fastq,Sample_49625/Sample_49625_adtrim_clean_qtrim_rev.derep.fastq,Sample_49629/Sample_49629_adtrim_clean_qtrim_rev.derep.fastq,Sample_53599/Sample_53599_adtrim_clean_qtrim_rev.derep.fastq,Sample_53600/Sample_53600_adtrim_clean_qtrim_rev.derep.fastq,Sample_53601/Sample_53601_adtrim_clean_qtrim_rev.derep.fastq --k-min 21 --k-max 141 --k-step 12 -t 30 -o Sample_Bryobacter_all_megahit_coassembly --verbose &> Sample_Bryobacter_all_megahit_coassembly.log &

#DJS 5/28/18

#Merged coverage files from each sample mapped to the 42896_53600 coassembly into one summary file for binning:
#wd: Sample_42896_53600_megahit_coassembly
comics -- merge-coverage -a final.contigs.chop.fa Mapping_42896/final.contigs.chop.genomeCovBed.tsv Mapping_53600/final.contigs.chop.genomeCovBed.tsv -o final.contigs.chop.cov.merged

#Bin the assembly of Sample 42896 using CONCOCT:
#wd: Sample_42896/Sample_42896_megahit_assembly/
nohup omics binning -a final.contigs.chop.fa --coverage-file ../Mapping/final.contigs.chop.cov --cpus 30 -o CONCOCT/ &> Sample_42896_concoct.log &

#Get CheckM tree placement of Sample 42896 CONCOCT bins:
#wd: Sample_42896/Sample_42896_megahit_assembly/CONCOCT/ 
nohup checkm tree_qa -f checkm/checkM_tree_qa.txt -o 2 --tab_table checkm/ &> checkm_tree_qa.log &

#Bin the coassembly of Sample 42896 and 53600 using CONCOCT:
#wd: Sample_42896_53600_megahit_coassembly/
comics -- nohup omics binning -a final.contigs.chop.fa --coverage-file final.contigs.chop.cov.merged --cpus 30 -o CONCOCT/ &> Sample_42896_53600_concoct.log &

#DJS 5/29/18

#CheckM of CONCOCT binning of 42896 and 53600 coassembly failed because bins 79 and 123 were empty.
#removed this bin and reran checkm:
#wd: Sample_42896_53600_megahit_coassembly/CONCOCT/
comics -- nohup checkm lineage_wf final.contigs.chop.bins/ checkM/ -f checkM/Sample_42896_53600_coassembly.chop.bin_stats.tsv --tab_table -x fa -t 30 &> checkM.log &

#Get CheckM tree placement of Sample 42896_53600 CONCOCT bins:
#wd: Sample_42896_53600_megahit_coassembly/CONCOCT/ 
nohup checkm tree_qa -f checkM/checkM_tree_qa.txt -o 2 --tab_table checkM/ &> checkM_tree_qa.log &

#Indexed Bryobacter_all coassembly:
#wd: Sample_Bryobacter_all_megahit_coassembly
comics -- nohup omics mapping -a final.contigs.fa --chop --index-only &> Sample_Bryobacter_all_coassembly_index.log &

#Mapped short reads to Sample Bryobacter all coassembly:
nohup sh Mapping_sample_Bryobacter_all.sh &> Mapping_sample_Bryobacter_all_coassembly.log &

#DJS 6/1/18

#Merged coverage files from each sample mapped to the Bryobacter_all coassembly into one summary file for binning:
#wd: Sample_Bryobacter_all_megahit_coassembly
comics -- merge-coverage -a final.contigs.chop.fa Mapping_42896/final.contigs.chop.genomeCovBed.tsv Mapping_49615/final.contigs.chop.genomeCovBed.tsv Mapping_49622/final.contigs.chop.genomeCovBed.tsv Mapping_49625/final.contigs.chop.genomeCovBed.tsv Mapping_49629/final.contigs.chop.genomeCovBed.tsv Mapping_53599/final.contigs.chop.genomeCovBed.tsv Mapping_53600/final.contigs.chop.genomeCovBed.tsv Mapping_53601/final.contigs.chop.genomeCovBed.tsv -o final.contigs.chop.cov.merged

#Binning of Sample_Bryobacter_all_coassembly using CONCOCT:
#wd: Sample_Bryobacter_all_megahit_coassembly
comics -- nohup omics binning -a final.contigs.chop.fa --coverage-file final.contigs.chop.cov.merged --cpus 30 -o CONCOCT/ &> Sample_Bryobacter_all_concoct.log &

#DJS 6/2/18

#Get assembly graph of Sample 42896 assembly:
#wd: Sample_42896/Sample_42896_megahit_assembly/intermediate_contigs/
comics -- megahit_toolkit contig2fastg 141 k141.contigs.fa > k141.fastg

#Tested bandage on Bin 33 of Sample_42896 assembly:
#See notes and protocol in the following directory: Sample_42896/Sample_42896_megahit_assembly/Bandage
#Bandage is up and running well, but will only be useful when looking at the graph of each bin individually while including the coloring of all the bin contigs.
#Should be run on targeted bins of interest to look for assembly and binning issues rather than whole datasets.

#DJS 6/3/18

#Some bins appear to be either spuriously split or not separated enough by CONCOCT. Probably due to differential coverage in strain heterogeneity for the first issue and organisms co-occurring at similar abundances for the second issue.
#Trying to bin based solely on tetramer signal with ESOM, then dereplicate and concactenate the bins from either method using DASTool.

#Generated ESOM files for Sample_42896 assembly:
#wd: Sample_42896/Sample_42896_megahit_assembly
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 2000 -max 5000 -DIR ESOM_2-5kb &> esomwrapper.log &

#Opened lrn file (Tetra_esom_2000.lrn) in the GUI and normalized the data using Robust ZT transform. Saved normalized data as Tetra_esom_2000.rzt.lrn

#Run the esom training:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_2-5kb/
#Run this perl scrip to get the train command, using the recommended cols and rows from the esomwrapper results:
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_2000.rzt.lrn -cls Tetra_esom_2000.cls -names Tetra_esom_2000.names -rows 297 -cols 594 -epochs 50 &> esomTrain.log &
#Excute the command (changed lrn file to the normalized rzt file):
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_2-5kb/
#First you need to give JAVA 5GB of heap memory with the following command:
export JAVA_TOOL_OPTIONS=-Xmx200g
#Then excute the training command generated in output of esomTrain.pl (modified so that it doesn't run the empty no_norm.lrn file):
nohup /opt/packages/ESOM/1.1/bin/esomtrn2 --permute --out Tetra_esom_2000.rzt.297x594e50.wts -b Tetra_esom_2000.rzt.297x594e50.bm --cls Tetra_esom_2000.cls --lrn Tetra_esom_2000.rzt.lrn --algorithm kbatch --rows 297 --columns 594 -bmc 37 --start-radius 149 --epochs 50 -k 0.15 --bmsearch standard --dist euc &>> esomTrain.log &

#DJS 6/8/18

#Get CheckM tree placement of Sample Bryobacter all coassembly CONCOCT bins:
#wd: Sample_Bryobacter_all_megahit_coassembly/CONCOCT/
nohup checkm tree_qa -f checkM/checkM_tree_qa.txt -o 2 --tab_table checkm/ &> checkM_tree_qa.log &

#Run checkm of bins using remerged contigs for Sample_42896:
#wd: Sample_42896/Sample_42896_megahit_assembly/CONCOCT
comics -- nohup checkm lineage_wf final.contigs.bins/ checkm_merged_contigs/ -f checkm_merged_contigs/Sample_42896.merged.bin_stats.tsv --tab_table -x fa -t 20 &> checkm_merged_contigs.log &

#DJS 6/12/18

#Got placement of Sample_42896 bins for Sample_42896 (using merged contigs) into checkm reference tree:
nohup checkm tree_qa -f checkm_merged_contigs/checkm_merged_contigs_tree_qa.txt -o 2 --tab_table checkm_merged_contigs/ &> checkM_merged_contigs_tree_qa.log &

#No difference between CheckM results using the bins containing the chopped or merged contigs
#While it appeared that some ESOM bins are clearly clustered on the ESOM map, other regions look unresolved.

#DJS 6/13/18

#Make new ESOM map with contig size window of 4-10kb:
#wd: Sample_42896/Sample_42896_megahit_assembly
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 4000 -max 10000 -DIR ESOM_4-10kb &> esomwrapper_4-10kb.log &

#Opened lrn file (ESOM_4-10kb/Tetra_esom_4000.lrn) in the GUI and normalized the data using Robust ZT transform. Saved normalized data as Tetra_esom_4000.rzt.lrn

#Run the esom training:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb/
#Run this perl scrip to get the train command, using the recommended cols and rows from the esomwrapper results:
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_4000.rzt.lrn -cls Tetra_esom_4000.cls -names Tetra_esom_4000.names -rows 174 -cols 348 -epochs 50 &> esomTrain.log &
#Run the training command below generated from the esomTrain.pl script:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb/
nohup /opt/packages/ESOM/1.1/bin/esomtrn --permute --out Tetra_esom_4000.rzt.174x348e50.wts -b Tetra_esom_4000.rzt.174x348e50.bm --cls Tetra_esom_4000.cls --lrn Tetra_esom_4000.rzt.lrn --algorithm kbatch --rows 174 --columns 348 -bmc 22 --start-radius 87 --epochs 50 -k 0.15 --bmsearch standard --dist euc &>> esomTrain.log &

#DJS 6/14/18:

#Clustering was much better with contig size window of 4-10 kb
#Load map:

#Clipped bins then saved class file as Tetra_esom_4000_binned.cls
#Get the bin fasta files:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb
mkdir Bin_fastas/
cd Bin_fastas/
for i in `seq 4 30`; do perl getClassFasta.pl -cls Tetra_esom_4000_binned.cls -names Tetra_esom_4000.names -fasta esom.fa -num $i -loyal 60 -id k141_*; done

#Run CheckM on ESOM bins:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb
comics -- nohup checkm lineage_wf Bin_fastas/ checkm/ -f checkm/Sample_42896.ESOM.bin_stats.tsv --tab_table -x fasta -t 20 &> checkm_ESOM.log &
nohup checkm tree_qa -f checkm/Sample_42896_ESOM_checkm_tree_qa.txt -o 2 --tab_table checkm/ &> checkm_ESOM_tree_qa.log &

#DJS 6/15/18:

#Removed ESOM_2-5kb directory, because the 4-10kb map had more clear boundaries between bins:

###Combine and dereplicate bins generated through CONCOCT and ESOM using DASTOOL:
#Convert the CONCOCT comma-separated list of scaffoldIDs and BinIDs to a tab-delimited one readable by DASTool:
#wd: Sample_42896/Sample_42896_megahit_assembly/CONCOCT
perl -pe 's/,/\tbin./g' clustering_gt1000.csv > concoct_scaffolds2bin.tsv

#Generate the scaffold2bin file for ESOM bins by entering these one-liners in succession:
#wd:
cat *.conf > esom_scaffolds2bin.tsv
sed '/^#/ d' esom_scaffolds2bin.tsv > esom_scaffolds2bin.cleaned.tsv
awk '{print $2,$1}' esom_scaffolds2bin.cleaned.tsv > esom_scaffolds2bin.tsv
rm esom_scaffolds2bin.cleaned.tsv

#Realized that split and merged contigs using concoct in Robert's pipeline do not match up with the original bin headers, which is not good for using DAS-Tool.
#Furthermore, chopping scaffolds before mapping and index is going to cause many downstream problems...

#Remove old mapping files, then re-indexed each assembly using bow-tie:
#First rename the contigs so they don't have complex spaces and symbols (which will create problems for anvio later):
#done for each assembly in their respective directories:
sed -e 's/\s/_/g' final.contigs.fa > final.contigs.fixed.fa
sed -e 's/=/_/g' final.contigs.fixed.fa > final.contigs.fixed2.fa
mv final.contigs.fixed2.fa final.contigs.fa
rm final.contigs.fixed.fa
#Index the assemblies:
comics -- nohup omics mapping -a final.contigs.fa --index-only &> Sample_42896_megahit_assembly_index.log & #wd: Sample_42896/Sample_42896_megahit_assembly/
comics -- nohup omics mapping -a final.contigs.fa --index-only &> Sample_42896_53600_megahit_assembly_index.log & #wd: Sample_42896_53600_megahit_coassembly/
comics -- nohup omics mapping -a final.contigs.fa --index-only &> Sample_Bryobacter_all_megahit_coassembly_index.log & #wd: Sample_Bryobacter_all_megahit_coassembly/

#Redo mapping -- map short reads to unchopped contigs (with simple headers):
#wd: Sample_42896/
comics -- nohup omics mapping --index-dir Sample_42896_megahit_assembly/bowtie2-index/ -a Sample_42896_megahit_assembly/final.contigs.fa -f Sample_42896_adtrim_clean_qtrim_fwd.derep.fastq  -r Sample_42896_adtrim_clean_qtrim_rev.derep.fastq --cpus 15 -o Mapping &> Sample_42896_mapping.log &
#Now for the 42896 and 53600 coassemlby (wd: ./):
nohup sh Mapping_sample_42896_53600_mapping.sh &> Mapping_sample_42896_53600_coassembly.log &
#Merge coverage files:
#wd: Sample_42896_53600_megahit_coassembly/
comics -- merge-coverage -a final.contigs.fa Mapping_42896/final.contigs.genomeCovBed.tsv Mapping_53600/final.contigs.genomeCovBed.tsv -o final.contigs.merged.cov
#Lastly, the Bryobacter_all coassembly (wd: ./):
nohup sh Mapping_sample_Bryobacter_all.sh &> Mapping_sample_Bryobacter_all_coassembly.log &
#Merge coverage files:
#wd: Sample_Bryobacter_all_megahit_coassembly
comics -- nohup merge-coverage -a final.contigs.fa Mapping_42896/final.contigs.genomeCovBed.tsv Mapping_49615/final.contigs.genomeCovBed.tsv Mapping_49622/final.contigs.genomeCovBed.tsv Mapping_49625/final.contigs.genomeCovBed.tsv Mapping_49629/final.contigs.genomeCovBed.tsv Mapping_53599/final.contigs.genomeCovBed.tsv Mapping_53600/final.contigs.genomeCovBed.tsv Mapping_53601/final.contigs.genomeCovBed.tsv -o final.contigs.merged.cov &> merge-coverage.log &

#DJS 10/7/18:
#Mapped the short reads from Sample 53600 to the single sample assembly of Sample 42896 in order to use differential coverage binning:
#wd: Sample_42896/
comics -- nohup omics mapping --index-dir Sample_42896_megahit_assembly/bowtie2-index/ -a Sample_42896_megahit_assembly/final.contigs.fa -f ../Sample_53600/Sample_53600_adtrim_clean_qtrim_fwd.derep.fastq -r ../Sample_53600/Sample_53600_adtrim_clean_qtrim_rev.derep.fastq --cpus 15 -o Sample_53600_reads_mapped_to_Sample_42896_contigs &> Sample_53600_reads_mapped_to_Sample_42896_contigs.log &

######################################################################
########## Start anvio binning of Sample_42896 assembly ##############
######################################################################
#First generate contigs database
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-gen-contigs-database -f ../final.contigs.fa -o anvio_contigs.db -n Sample_42896_assembly --split-length 10000 --kmer-size 4 &> anvi-gen-contigs-database.log &

#Look for gene start and stop sites and rRNA and single-copy house-keeping genes in the contigs database:
nohup anvi-run-hmms -c anvio_contigs.db --num-threads 15 &> anvi-run-hmms.log &

#DJS 7/1/18:

#Annotate Anvio gene calls with NCBI COGS:
nohup anvi-run-ncbi-cogs -c anvio_contigs.db --cog-data-dir /geomicro/data2/COMMON/anviCOGs/ -T 16 --search-with diamond --sensitive &> anvi-run-ncbi-cogs.log &

##Next annotate Anvio gene calls with functions via GhostKoala:
##The parser scripts and files were obtained from https://github.com/edgraham/GhostKoalaParser.git and made by Elaina Graham at the University of Southern California (parser endorsed by Murat).
#First, get the amino acid sequences of the genes found in the contigs database:
nohup anvi-get-aa-sequences-for-gene-calls -c anvio_contigs.db -o Sample_42896_gene_calls.faa &> anvi-get-aa-sequences-for-gene-calls.log

#Make gene sequences headers start with the prefix 'genecall' so the GhostKoala parser scripts won't receive an error (doesn't like seqIDs that start with numbers):
sed 's/>/>genecall_/g' Sample_42896_gene_calls.faa > Sample_42896_gene_calls_edited.faa
mv Sample_42896_gene_calls_edited.faa Sample_42896_gene_calls.faa

#Next upload the gene file to GhostKoala for annotation in the KEGG database.

#DJS 7/3/18:
#Download the GhostKOALA annotations when it is finished: filename user_ko.txt

#Run this command to add the necessary header line:
echo -e "contig\taccession_id" > .temp && cat user_ko.txt >> .temp && mv .temp user_ko.txt

#Run the script to parse the GhostKoala output into a format readable by ANVIO:
python KEGG-to-anvio.edited --KeggDB KO_Orthology_ko00001.txt -i user_ko.txt -o KeggAnnotations-AnviImportable.txt

#The script treated the column headers of the Ghostkoala output as a gene, had to remove this row from KeggAnnotations-AnviImportable.txt before importing functions. Done manually in nano.

#Next import the funcitons:
comics -i bio.img -- nohup anvi-import-functions -c anvio_contigs.db -i KeggAnnotations-AnviImportable.txt &> anvi-import-functions.log &

#DJS 7/7/18:

#Imported KEGG taxnomic annotations into anvio:
#wd: /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_42896/Sample_42896_megahit_assembly/ANVIO

#First parse the output from GhostKoala into a format readable by ANVIO (also written by Elaina Graham):
python GhostKOALA-taxonomy-to-anvio user.out.top KeggTaxonomy.txt

#Import to anvio contigs database:
comics -i bio.img -- nohup anvi-import-taxonomy -c anvio_contigs.db -p default_matrix -i KeggTaxonomy.txt &> anvi-import-taxonomy.log &

#DJS 7/11/18:

#Profile the anvio contigs database for Sample_42896 assembly and cluster contigs by euclidean distance:
comics -i bio.img -- nohup anvi-profile -i /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_42896/Mapping/final.contigs.sorted.bam -c anvio_contigs.db -o Sample_42896_anvio_profile --sample-name Sample_42896 --min-contig-length 2500 --cluster-contigs --num-threads 16 &> anvi-profile-42896.log &

#DJS 7/12/18:

#Tried exporting the splits to run CONCOCT, as well as profiling the sample individually. However, most of the desired functions in anvio require a merged profile database. Single assemblies will have to forgo the combined binning approach with differential coverage for now...

#DJS 10/7/2018
#Normalizing kmers of the coassembly did not yield a complete Microcystis bin.
#Now trying to map the short reads from Sample 53600 to the contigs from the single assembly of Sample 42896 to use differential coverage binning and anvio-profiling
#See above for the mapping commands

#Profile the Sample 42896 data with the Sample 53600 mapping data:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-profile -i ../../Sample_53600_reads_mapped_to_Sample_42896_contigs/final.contigs.sorted.bam -c anvio_contigs.db -o Sample_42896_anvio_profile_with_Sample_53600_reads/ --sample-name Sample_53600 --min-contig-length 2500 --num-threads 16 &> anvi-profile-53600.log &

#Merge the profiles created from reads in Samples 42896 and 53600, then bin using CONCOCT:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-merge Sample_42896_anvio_profile/PROFILE.db Sample_42896_anvio_profile_with_Sample_53600_reads/PROFILE.db -c anvio_contigs.db -o Sample_42896_merged_profile --sample-name Sample_42896_megahit_assembly &> anvi-merge.log &

#DJS 10/8/2018

#Export CONCOCT binning results:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-export-collection -p Sample_42896_merged_profile/PROFILE.db -C CONCOCT -O CONCOCT_binning_results &> anvi-export-collection-concoct.log &

#Export the splits and the coverage information for third-party binning software:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-export-splits-and-coverages -p Sample_42896_merged_profile/PROFILE.db -o Sample_42896_merged_profile/ -c anvio_contigs.db -O Sample_42896_megahit_assembly --splits-mode &> anvi-export-splits-and-coverages.log &

#Bin using metabat:
#wd: Sample_42896/Sample_42896_megahit_assembly/METABAT
comics -- nohup metabat -i Sample_42896_megahit_assembly-SPLITS.fa -o Bin -a Sample_42896_megahit_assembly-COVs.txt --minContig 2500 --cvExt -t 16 --saveCls --onlyLabel -v &> metabat.log &

#Parse the metabat output to get the binning results table to import the collection in anvio:
#wd: Sample_42896/Sample_42896_megahit_assembly/METABAT
sh Metabat_to_anvio_parser.sh

#Import the metabat bins as a collection in Anvio:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/

#Rerun ESOM using splits from ANVIO contigs database (kept contig size window of 4-10 kb):
#First generate the files:
#wd: Sample_42896/Sample_42896_megahit_assembly/
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 4000 -max 10000 -DIR ESOM_4-10kb &> esomwrapper_4-10kb.log &

#Opened esom files in GUI, normalized the lrn file using Robust ZT. Saved new data as Tetra_esom_4000.rzt.lrn

#Run training script:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb/
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_4000.rzt.lrn -cls Tetra_esom_4000.cls -names Tetra_esom_4000.names -rows 174 -cols 348 -epochs 50 &> esomTrain.log &

#Run the training command generated from the wrapper script:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb/
nohup esomtrn --permute --out Tetra_esom_4000.rzt.174x348e50.wts -b Tetra_esom_4000.rzt.174x348e50.bm --cls Tetra_esom_4000.cls --lrn Tetra_esom_4000.rzt.lrn --algorithm kbatch --rows 174 --columns 348 -bmc 22 --start-radius 87 --epochs 50 -k 0.15 --bmsearch standard --dist euc &>> esomTrain.log &

#ESOM map finished, binned contigs and saved class file as Sample_42896/Sample_42896_megahit_assembly/ESOM/ESOM_4-10kb/Tetra_esom_4000_binned.cls

#Get ESOM bin fastas:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb/
nohup sh batch_getClassFasta.sh &> batch_getClassFasta.log &

#Generate ESOM_binning_results_file for Anvio and DAStool:
#wd: Sample_42896/Sample_42896_megahit_assembly/ESOM_4-10kb/
sh ESOM_binning_results_parser.sh

#Import the ESOM bins into Anvio as a collection:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-import-collection -c anvio_contigs.db -p Sample_42896_merged_profile/PROFILE.db -C ESOM ../ESOM_4-10kb/ESOM_binning_results.txt &> anvi-import-collection-esom.log &

#Refined large Bin12 containing Microcystis and Anabaena contigs.
#Exported the update collection info into a binning_results_file for DASTool

#DJS 10/10/18:

#Binned the Sample 42896 assembly (anvio splits) using VizBin on personal laptop using tetramer frequency 
#Imported the VizBin bins as a collection into the anvio database:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-import-collection -c anvio_contigs.db -p Sample_42896_merged_profile/PROFILE.db -C VizBin VizBin_binning_results.txt &> anvi-import-collection-vizbin.log &

#Summarized the binning results and dereplicated bins using DASTool:
#wd: Sample_42896/Sample_42896_megahit_assembly/DASTOOL/
nohup das-tool -i CONCOCT_binning_results.txt,ESOM_binning_results_updated.txt,Metabat_binning_results.txt,VizBin_binning_results.txt -c Sample_42896_megahit_assembly-SPLITS.fa -o DASTool_bins --search_engine blast --threads 16 &> das-tool.log &

#Edit the DASTool bin ids so they are compatible with Anvio:
#wd: Sample_42896/Sample_42896_megahit_assembly/DASTOOL/
sed 's/results\./results_/g' DASTool_bins_DASTool_scaffolds2bin.txt > DASTool_binning_results.txt
sed 's/updated\./updated_/g' DASTool_binning_results.txt > DASTool_binning_results2.txt
mv DASTool_binning_results2.txt DASTool_binning_results.txt

#Import the DASTool collection into anvio:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-import-collection -c anvio_contigs.db -p Sample_42896_merged_profile/PROFILE.db -C DASTool ../DASTOOL/DASTool_binning_results.txt &> anvi-import-collection-dastool.log &

#Started refining DASTool bins in anvio using anvi-refine. Saved updated bins in DASTool collection

#DJS 10/11/2018:

#Summarize the refined DASTool bins:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup anvi-summarize -p Sample_42896_merged_profile/PROFILE.db -c anvio_contigs.db --init-gene-coverages -C DASTool -o DASTool_refined_anvi_summarize &> anvi-summarize.log &

#DJS 10/13/2018:

#Microcystis bin is only ~60% complete. Try to recruit smaller contigs to the bin using bandage.
#First generate the bin-colors csv file to label nodes in the graph based on which bin they belong to:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
python anvi-summary-to-bandage-csv.py
#Manually removed "partial" flag from contigs that were partially binned.

#Move the results to the bandage working directy:
mv bin-colors.csv ../Bandage/

#Generate the assembly graph:
#wd: Sample_42896/Sample_42896_megahit_assembly/intermediate_contigs/
comics -- megahit_toolkit contig2fastg 141 k141.contigs.fa > k141.fastg

##Next, we need a list of the nodes (contigs) in the Microcystis bin on one line separated by commas. 
##Generate this by manually selecting and copying the contigs in the Microcystis bin (listed in the bin-colors.csv file) to a separate text file.
##Using text wrangler, remove all text in the contig headers but the node/contig number, and the node names onto one line (separated by commas).

#Use bandage to generate an assembly graph around the binned contigs (range of 50 nodes around each listed node in the bin):
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
bandage reduce k141.fastg CONCOCT_binning_results_012.gfa --scope aroundnodes --distance 50 --nodes `cat CONCOCT_binning_results_012_nodes.txt`

#DJS 8/15/18:

#Visualize the graph for the bin:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
bandage load CONCOCT_binning_results_012.gfa 

#Assembly graph shows many small contigs connected to the binned Microcystis contigs
#Spot checked ~100-200 nodes in NCBI Blast nr database do confirm Microcystis
#Top 5 hits were all to Microcystis at 92-99% identity.
#Sequences saved in selected_sequences_test.fasta files

#Saved all node seqeunces in graph as CONCOCT_binning_results_012_plus_bandage_nodes.fa
#Manually removed a few nodes that were linked to the Pseudanabaena bin

#Next for the Phenylobacterium bin (ESOM_binning_results_023):
#Manually get a list of contigs/nodes in the bin

#Reduce the graph to be 50 nodes around the binned contigs:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
bandage reduce k141.fastg ESOM_binning_results_023.gfa --scope aroundnodes --distance 50 --nodes `cat ESOM_binning_results_updated_023_nodes.txt`

#Saved all node sequences in graph as ESOM_binning_results_023_plus_bandage_nodes.fa
#Manually removed some nodes that were linked to a Methylobacteracea bin (confirmed by checking with blast against NCBI nr database)

#DJS 10/16/18:

#Repeat again for the incomplete Acetobacteraceae bin:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/

#Generated nodes file in textwrangler using the contigs listed in bin-color.csv file for the Roseomonas bin (Vizbin_binning_results_012)
#Reduce the assembly graph around the Acetobacteraceae nodes:
bandage reduce k141.fastg VizBin_binning_results_012.gfa --scope aroundnodes --distance 50 --nodes `cat VizBin_binning_results_012_nodes.txt`

#Saved all nodes in graph as VizBin_binning_results_012_plus_bandage_nodes.fa

#DJS 10/19/2018:

#Double check that only contigs/nodes from the desired genomes are captured with bandage:
#Do this by blasting the fasta files generated through bandage against NCBI's nt database:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
##First run for the Microcystis bin and the Phenylobacterium bin
nohup blastn -query CONCOCT_binning_results_012_plus_bandage_nodes.fa -db /geomicro/data9/flux/reference-data/blast/nt -outfmt "7 std qlen slen qcovs" -out CONCOCT_binning_results_012_bandage_v_nt.blastn &> CONCOCT_binning_results_012_bandage_vs_nt.log &
nohup blastn -query ESOM_binning_results_023_plus_bandage_contigs.fa -db /geomicro/data9/flux/reference-data/blast/nt -outfmt "7 std qlen slen qcovs" -out ESOM_binning_results_023_bandage_v_nt.blastn &> ESOM_binning_results_023_bandage_vs_nt.log &
#Next for the Rosemonas bin:
nohup blastn -query VizBin_binning_results_012_plus_bandage_nodes.fa -db /geomicro/data9/flux/reference-data/blast/nt -outfmt "7 std qlen slen qcovs" -out VizBin_binning_results_012_bandage_vs_nt.blastn &> VizBin_binning_results_012_bandage_vs_nt.log &

##Manual curation of the Microcystis bin obtained from bandage:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
#First get all top 10 hits from the blast search against nt:
perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -b CONCOCT_binning_results_012_bandage_v_nt.blastn -t 10 -out CONCOCT_binning_results_012_bandage_v_nt.top10.blastn

#Then removed nodes on edges of graph that did not have significant top hits to Microcystis:
#Edge is defined as a path of edgeds and nodes ending.

#DJS 10/23/18:

#Ran checkm on the initial final bin fastas (before incorporating Bandage results):
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup checkm lineage_wf Final_bin_fastas/ Final_bin_fasta_initial_checkm/ -f Final_bin_fasta_initial_checkm/Sample_42896.initial.checkm.bin_stats.tsv --tab_table -x fa -t 16 &> checkm_initial.log &

#Get placement of initial bins in Checkm reference tree:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup checkm tree_qa -f Final_bin_fasta_initial_checkm/Sample_42896_final_bin_initial_tree_qa.txt -o 2 --tab_table Final_bin_fasta_initial_checkm/ &> checkm_initial_tree_qa.log &

#Manual curation of Phenylobacterium bin:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
#Get all top 5 hits from the blast search against nt:
perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -b ESOM_binning_results_023_bandage_v_nt.blastn -o ESOM_binning_results_023_bandage_v_nt.top5.blastn -t 5

#Manual curation of Acetobacteraceae bin:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
#Get all top 5 hits from the blast search against nt:
perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -b VizBin_binning_results_012_bandage_vs_nt.blastn -o VizBin_binning_results_012_bandage_vs_nt.top5.blastn -t 5
#No highly similar hits to nt; all within 70-82% ID range. All nodes had best hits to proteobacteria, so no nodes were removed from fasta.

#DJS 10/24/18:

#CheckM and Bandage results suggest that ESOM bin 23 and Metabat bin 27 could be merged.
#Get a combined list of nodes for ESOM bin 23 and Metabat bin 27 using the bin-colors spreadsheet.
#reduce the assembly graph around those nodes:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
bandage reduce k141.fastg Merged_Phenylobacterium_bin.gfa --scope aroundnodes --distance 50 --nodes `cat Merged_Phenylobacterium_bin_nodes.txt` 

#Saved all nodes in graph as a fasta file "Merged_Phenylobacterium_bin.fa"

#Blast the merged fasta file against nt to validate that the small contigs belong to the bin:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
nohup blastn -query Merged_Phenylobacterium_bin.fa -db /geomicro/data9/flux/reference-data/blast/nt -outfmt "7 std qlen slen qcovs" -out Merged_Phenylobacterium_bin_vs_nt.blastn &

#Get top 5 hits for each query sequence:
#wd: Sample_42896/Sample_42896_megahit_assembly/Bandage/
perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -b Merged_Phenylobacterium_bin_vs_nt.blastn -o Merged_Phenylobacterium_bin_vs_nt.top5.blastn -t 5

#Manually removed nodes that were on edges of graph connections that did not have significant hits to Caulobacteraceae.

#DJS 10/26/2018:

#Copied the refined bins including the bandage contigs over to the directory for checkM (Sample_42896/Sample_42896_megahit_assembly/ANVIO/Final_bin_fastas)

#Ran checkm of all final bins, including bins where bandage was used to recruit small contigs:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO
nohup checkm lineage_wf Final_bin_fastas/ Final_bin_fasta_post_bandage_checkm/ -f Final_bin_fasta_post_bandage_checkm/Sample_42896.postBandage.checkm.bin_stats.tsv --tab_table -x fa -t 16 &> checkm_postBandage.log &

#Get placement of bins in checkm reference tree:
#wd: Sample_42896/Sample_42896_megahit_assembly/ANVIO/
nohup checkm tree_qa -f Final_bin_fasta_post_bandage_checkm/Sample_42896_postBandage_tree_qa.txt -o 2 --tab_table Final_bin_fasta_post_bandage_checkm/ &> checkm_post_bandage_tree_qa.log &

########################################################################################################
                     Start binning process with Sample 42896_53600 coassembly
########################################################################################################

#DJS 7/12/18:

#First generate the contigs database:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- nohup anvi-gen-contigs-database -f ../final.contigs.fa -o Sample_42896_53600_anvio_contigs.db -n Sample_42896_53600_coassembly --split-length 10000 --kmer-size 4 &> anvi-gen-contigs-database.log &

#Look for rRNA and single-copy housekeeping genes in the contigs database:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- nohup anvi-run-hmms -c Sample_42896_53600_anvio_contigs.db --num-threads 15 &> anvi-run-hmms.log &

#DJS 7/13/18:

#Next start KEGG GhostKoala annotation of anvio gene calls:
#Get fasta file of gene calls found by ANVIO:
comics -i bio.img -- nohup anvi-get-aa-sequences-for-gene-calls -c Sample_42896_53600_anvio_contigs.db -o Sample_42896_53600_coassembly_gene_calls.faa &> anvi-get-aa-sequences-for-gene-calls.log &

#Edit the gene sequence file so that the headers start with the prefix "genecall_". This is necessary for GhostKoala:
sed 's/>/>genecall_/g' Sample_42896_53600_coassembly_gene_calls.faa > Sample_42896_53600_coassembly_gene_calls.edited.faa
mv Sample_42896_53600_coassembly_gene_calls.edited.faa Sample_42896_53600_coassembly_gene_calls.faa

#Download and submit the file to Kegg Ghostkoala for functional and taxonomic annotation.

#Download the GhostKOALA annotations when it is finished: filename user_ko.txt

#Run the script to parse the GhostKoala output into a format readable by ANVIO:
python /omics/HABs/Microcystis/Erie2014_coAssembly/GhostKoalaParser/KEGG-to-anvio --KeggDB /omics/HABs/Microcystis/Erie2014_coAssembly/GhostKoalaParser/samples/KO_Orthology_ko00001.txt -i user_ko.txt -o KeggAnnotations-AnviImportable.txt

#Import the functional annotations into ANVIO contigs database:
comics -i bio.img -- nohup anvi-import-functions -c Sample_42896_53600_anvio_contigs.db -i KeggAnnotations-AnviImportable.txt &> anvi-import-functions.log &

#Next parse the GhostKoala taxonomic annotations into a format readable by ANVIO:
#First parse the output from GhostKoala into a format readable by ANVIO (also written by Elaina Graham):
python /omics/HABs/Microcystis/Erie2014_coAssembly/GhostKoalaParser/GhostKOALA-taxonomy-to-anvio user.out.top KeggTaxonomy.txt

#Then import the taxonomy into the anvio contigs database:
nohup anvi-import-taxonomy -c Sample_42896_53600_anvio_contigs.db -p default_matrix -i KeggTaxonomy.txt &> anvi-import-taxonomy.log &

#DJS 7/14/18:

#Profile the contigs database with the mapping files:
#wd: Sample_42896_53600_megahit_coassembly/
nohup sh batch_anvio_profile.sh &> batch_anvio_profile.log &

#Merge the profiles from the different samples used in the coassembly, then bin with CONCOCT:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- nohup anvi-merge Sample_42896_profile/PROFILE.db Sample_53600_profile/PROFILE.db -c Sample_42896_53600_anvio_contigs.db -o Sample_42896_53600_merged_profile --sample-name Coassembly_42896_53600 &> anvi-merge.log &

#Get list of contigs and their assigned CONCOCT bins:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- anvi-export-collection -p Sample_42896_53600_merged_profile/PROFILE.db -C CONCOCT -O CONCOCT_binning_results 

#Export the splits and coverages from the contigs database for the 42896 and 53600 coassembly:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- nohup anvi-export-splits-and-coverages -p Sample_42896_53600_merged_profile/PROFILE.db -o Sample_42896_53600_merged_profile/ -c Sample_42896_53600_anvio_contigs.db -O Sample_42896_53600_coassembly --splits-mode &> anvi-export-splits-and-coverages.log &

#Run binning of splits from contigs database using metabat:
#wd: /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_42896_53600_megahit_coassembly/METABAT/
comics -- nohup metabat -i Sample_42896_53600_coassembly-SPLITS.fa -o Bin -a Sample_42896_53600_coassembly-COVs.txt --minContig 2500 --cvExt -t 16 --saveCls --onlyLabel -v &> metabat.log &

#Parse metabat output into format that can be import as a collection into the anvio profile:
#wd: /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_42896_53600_megahit_coassembly/METABAT/
#This command will add the Bin ID in a column after the split name
for i in Bin.*; do sed -i "s/$/\t$i/" $i; done
#This command concatenates all the files into one binning results file for anvio:
cat Bin.* > Metabat_binning_results.txt
#Replace the dots with underscores to make anvio happy...
sed -i 's/Bin./Bin_/g' Metabat_binning_results.txt

#Import the metabat binning results into the anvio profile:
comics -i bio.img -- nohup anvi-import-collection -c Sample_42896_53600_anvio_contigs.db -p Sample_42896_53600_merged_profile/PROFILE.db -C METABAT ../METABAT/Metabat_binning_results.txt &> anvi-import-collection_metabat.log &

#DJS 7/15/18:

##Start binning the anvio splits from the contigs database using tetraESOM:

#First generate the ESOM files, using the same contig size window as the splits in anvio:
#wd: Sample_42896_53600_megahit_coassembly/ESOM/
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 2500 -max 10000 -DIR ESOM_2.5-10kb &> esomwrapper_2.5-10kb.log &

#Opened lrn file (Tetra_esom_2500.lrn) in the GUI and normalized the data using Robust ZT transform. Saved normalized data as Tetra_esom_2500.rzt.lrn

#Run the ESOM training:
#wd: Sample_42896_53600_megahit_coassembly/ESOM/ESOM_2.5-10kb/
#Run this perl scrip to get the train command, using the recommended cols and rows from the esomwrapper results:
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_2500.rzt.lrn -cls Tetra_esom_2500.cls -names Tetra_esom_2500.names -rows 286 -cols 572 -epochs 50 &> esomTrain.log &

#Next execute the command generated by the script. Need to manually remove the .no_norm flag on the lrn file and add in the path to the command:
#wd: Sample_42896_53600_megahit_coassembly/ESOM/ESOM_2.5-10kb/ 
nohup /opt/packages/ESOM/1.1/bin/esomtrn --permute --out Tetra_esom_2500.rzt.286x572e50.wts -b Tetra_esom_2500.rzt.286x572e50.bm --cls Tetra_esom_2500.cls --lrn Tetra_esom_2500.rzt.lrn --algorithm kbatch --rows 286 --columns 572 -bmc 36 --start-radius 143 --epochs 50 -k 0.15 --bmsearch standard --dist euc &> esomTrain.log &

#DJS 7/22/18:

Opened finished map, drew bins:
esomana -c Tetra_esom_2500.cls -l Tetra_esom_2500.rzt.lrn -n Tetra_esom_2500.names -w Tetra_esom_2500.rzt.286x572e50.wts -b Tetra_esom_2500.rzt.286x572e50.bm
Saved binned class file as: Tetra_esom_2500_binned.cls

#Generate binning_results file for ESOM bins using these one-liners:
#wd: Sample_42896_53600_megahit_coassembly/ESOM/ESOM_2.5-10kb/Bin_fastas/
cat *.conf > esom_scaffolds2bin.tsv
sed '/^#/ d' esom_scaffolds2bin.tsv > esom_scaffolds2bin.cleaned.tsv
awk 'BEGIN{OFS="\t"}{print $2,$1}' esom_scaffolds2bin.cleaned.tsv > esom_scaffolds2bin.tsv
awk 'BEGIN{OFS="\t"}{$2="Bin_"$2; print}' esom_scaffolds2bin.tsv > ESOM_binning_results.txt
perl -pe 's/(?<=\d)_(?=\d)/./g' ESOM_binning_results.txt > ESOM_binning_results.txt.fixed
sed 's/k141\./k141_/g' ESOM_binning_results.txt.fixed > ESOM_binning_results.txt
rm esom_scaffolds2bin.cleaned.tsv
rm esom_scaffolds2bin.tsv
rm ESOM_binning_results.txt.fixed

#Import ESOM bins into ANVIO:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- nohup anvi-import-collection -c Sample_42896_53600_anvio_contigs.db -p Sample_42896_53600_merged_profile/PROFILE.db -C ESOM ../ESOM/ESOM_2.5-10kb/Bin_fastas/ESOM_binning_results.txt &> anvi-import-collection_esom.log &

#Ran DAStool to summarize and dereplicate the multiple binning results:
#wd: Sample_42896_53600_megahit_coassembly/DASTOOL
nohup das-tool -i CONCOCT_binning_results.txt,ESOM_binning_results.txt,Metabat_binning_results.txt -c Sample_42896_53600_coassembly-SPLITS.fa -o DASTool_bins --search_engine blast --threads 16 &> das-tool.log &

#DJS 7/26/18

#Edit DASTool bin names to be compatible with ANVIO:
#wd: Sample_42896_53600_megahit_coassembly/DASTOOL/
sed 's/results\./results_/g' DASTool_bins_DASTool_scaffolds2bin.txt > DASTool_binning_results.txt

#Import the DASTool bins into ANVIO contigs database:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -i bio.img -- nohup anvi-import-collection -c Sample_42896_53600_anvio_contigs.db -p Sample_42896_53600_merged_profile/PROFILE.db -C DASTool ../DASTOOL/DASTool_binning_results.txt &> anvi-import-collection_DASTool.log &

#DJS 7/27/18

#Refined DASTool bins in anvio using differential coverage and sequence composition in anvi-interactive interface.
##This was done for CONCOCT seeded bins 42, 14, 25, 44, 48, 29, 16 (DASTool bins with high completeness and reduncancy above 10%)
##Commandline (wd=Sample_42896_53600_megahit_coassembly/ANVIO/):
anvi-refine -c Sample_42896_53600_anvio_contigs.db -p Sample_42896_53600_merged_profile/PROFILE.db -C DASTool -b CONCOCT_binning_results_016 --server-only -P 8000

##Note: -b flag above will be replaced with the appropriate bin id for each DASTool bin.

#DJS 7/29/18

#Summarize refined DASTool collection and get bin fasta files:
comics -i bio.img -- nohup anvi-summarize -p Sample_42896_53600_merged_profile/PROFILE.db -c Sample_42896_53600_anvio_contigs.db --init-gene-coverages -C DASTool -o DASTool_anvi_refine_summarize &> anvi-summarize.log &

#Made softlinks to all of the high quality DASTool bins (completeness >70% and redundancy <5%) in the directory Sample_42896_53600_megahit_coassembly/ANVIO/HiQal_DASTool_refined_bin_fastas/ for CheckM and Phylosift analysis.

#Run CheckM to compare completeness and contamination metrics with those provided by anvio and to infer taxonomy of bins:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -- nohup checkm lineage_wf HiQal_DASTool_refined_bin_fastas/ HiQal_DASTool_refined_checkm/ -f HiQal_DASTool_refined_checkm/Sample_42896_53600_coassembly.checkm.bin_stats.tsv --tab_table -x fa -t 30 &> checkM.log &

#DJS 7/30/2018:

#Get placement of bins in checkm genome tree:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -- nohup checkm tree_qa -f HiQal_DASTool_refined_checkm/Sample_42896_53600_coassembly_checkm_tree_qa.txt -o 2 --tab_table HiQal_DASTool_refined_checkm/ &> checkm_tree_qa.log &

#Run phylosift on bins using bacterial/archaeal/eukaryotic markers:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/HiQal_DASTool_refined_bin_fastas/
nohup sh ../phylosift_wrapper.sh &> phylosift.log &

#DJS 8/16/18:

#Update anvio contigs database to version 12 in anvio v5:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup anvi-migrate-db Sample_42896_53600_anvio_contigs.db &> anvi-migrate-db.log &

#DJS 8/22/18:

#Annotation anvio gene calls with NCBI COGs:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup anvi-run-ncbi-cogs -c Sample_42896_53600_anvio_contigs.db --cog-data-dir /geomicro/data2/COMMON/anviCOGs/ -T 16 --search-with diamond --sensitive &> anvi-run-ncbi-cogs.log &

#Refined all bins (even low contamination bins) with anvi-refine, removed outlier contigs using coverage and taxonomy annotation where applicable.

#Resummarize the bins:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup anvi-summarize -p Sample_42896_53600_merged_profile/PROFILE.db -c Sample_42896_53600_anvio_contigs.db --init-gene-coverages -C DASTool -o DASTool_anvi_refine_summarize &> anvi-summarize.log &

#Update the softlinks to all of the high quality DASTool bins (completeness >70% and redundancy <5%) in the directory Sample_42896_53600_megahit_coassembly/ANVIO/HiQal_DASTool_refined_bin_fastas/ for CheckM and Phylosift analysis.

#Rerun checkm on the updated bins:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -- nohup checkm lineage_wf HiQal_DASTool_refined_bin_fastas/ HiQal_DASTool_refined_checkm/ -f HiQal_DASTool_refined_checkm/Sample_42896_53600_coassembly.checkm.bin_stats.tsv --tab_table -x fa -t 30 &> checkM.log &

#Rerun checkm tree qa:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
comics -- nohup checkm tree_qa -f HiQal_DASTool_refined_checkm/Sample_42896_53600_coassembly_checkm_tree_qa.txt -o 2 --tab_table HiQal_DASTool_refined_checkm/ &> checkm_tree_qa.log &

#DJS 9/4/18:

#Tried binning with VizBin on personal computer.
##Used 4mer with default settings and no coverage.

#Import vizbin bins as a collection into anvio:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup anvi-import-collection -c Sample_42896_53600_anvio_contigs.db -p Sample_42896_53600_merged_profile/PROFILE.db -C VizBin VizBin_binning_results.txt &> anvi-import-collection_VizBin.log &

#Run DASTool including the Vizbin results:
#wd: Sample_42896_53600_megahit_coassembly/DASTOOL_with_VizBin/
nohup das-tool -i CONCOCT_binning_results.txt,ESOM_binning_results.txt,Metabat_binning_results.txt,VizBin_binning_results.txt -c Sample_42896_53600_coassembly-SPLITS.fa -o DASTool_bins --search_engine blast --threads 16 &> das-tool.log &

#Edit DASTool bin names into format suitable for anvio:
#wd: Sample_42896_53600_megahit_coassembly/DASTOOL_with_VizBin/
sed 's/results\./results_/g' DASTool_bins_DASTool_scaffolds2bin.txt > DASTool_binning_results.txt

#Import the DASTool bins including Vizbin into the Anvio profile as a collection:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup anvi-import-collection -c Sample_42896_53600_anvio_contigs.db -p Sample_42896_53600_merged_profile/PROFILE.db -C DASTool_with_VizBin ../DASTOOL_with_VizBin/DASTool_binning_results.txt &> anvi-import-collection-DASTool-with-VizBin.log &

#DJS 9/6/2018:

#Refined the bins with anvi-refine. Updated anvio collection.

#Summarize the refined DASTool_with_VizBin collection, which had better bin stats than the DASTool collection without VizBin:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup anvi-summarize -p Sample_42896_53600_merged_profile/PROFILE.db -c Sample_42896_53600_anvio_contigs.db --init-gene-coverages -C DASTool_with_VizBin -o DASTool_with_VizBin_anvi_summarize &> anvi-summarize.log &

#Run CheckM on the new DASTool bins that included VizBin results:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup checkm lineage_wf HiQal_DASTool_refined_bin_fastas/ HiQal_DASTool_refined_checkm/ -f HiQal_DASTool_refined_checkm/Sample_42896_53600_coassembly.checkm.bin_stats.tsv --tab_table -x fa -t 30 &> checkM.log &

#Get placement of bins in CheckM reference genome tree:
#wd: Sample_42896_53600_megahit_coassembly/ANVIO/
nohup checkm tree_qa -f HiQal_DASTool_refined_checkm/Sample_42896_53600_coassembly_checkm_tree_qa.txt -o 2 --tab_table HiQal_DASTool_refined_checkm/ &> checkm_tree_qa.log &

########################################################################################################
                     Start binning process with Sample Bryobacter all coassembly
########################################################################################################

#DJS 8/9/2018:

#Generate contigs database in anvio:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-gen-contigs-database -f ../final.contigs.fa -o Bryobacter_all_anvio_contigs.db -n Bryobacter_all_coassembly --split-length 10000 --kmer-size 4 &> anvi-gen-contigs-database.log &

#Look for rRNA genes and single-copy marker genes in the contigs database:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-run-hmms -c Bryobacter_all_anvio_contigs.db --num-threads 15 &> anvi-run-hmms.log &

#DJS 8/11/2018:

#Get amino acid sequences of gene calls for annotation in KEGG GhostKoala:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-get-aa-sequences-for-gene-calls -c Bryobacter_all_anvio_contigs.db -o Bryobacter_all_coassembly_gene_calls.faa &> anvi-get-aa-sequences-for-gene-calls.log &

#Add the prefix "genecall" to each fasta header of the gene call file:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
sed 's/>/>genecall_/g' Bryobacter_all_coassembly_gene_calls.faa > Bryobacter_all_coassembly_gene_calls.edited.faa
mv Bryobacter_all_coassembly_gene_calls.edited.faa Bryobacter_all_coassembly_gene_calls.faa

#The file is too large for GhostKoala online submission, so split it into two equal size parts:
split -b 244m Bryobacter_all_coassembly_gene_calls.faa Bryobacter_all_coassembly_gene_calls.faa.

#One gene was split between the two file parts, manually pasted back together. Then submit to GhostKoala for annotation.

#DJS 8/15/2018:

#Run the script to parse the GhostKoala output into a format readable by ANVIO:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
python /omics/HABs/Microcystis/Erie2014_coAssembly/GhostKoalaParser/KEGG-to-anvio --KeggDB /omics/HABs/Microcystis/Erie2014_coAssembly/GhostKoalaParser/samples/KO_Orthology_ko00001.txt -i user_ko.txt -o KeggAnnotations-AnviImportable.txt

#Import the functions into anvio:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-import-functions -c Bryobacter_all_anvio_contigs.db -i KeggAnnotations-AnviImportable.txt --drop-previous-annotations &> anvi-import-functions.log &

#Next parse the GhostKoala taxonomic annotations into a format readable by ANVIO:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
python /omics/HABs/Microcystis/Erie2014_coAssembly/GhostKoalaParser/GhostKOALA-taxonomy-to-anvio user.out.top KeggTaxonomy.txt

#Import GhostKoala taxonomy into anvio contigs db:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-import-taxonomy -c Bryobacter_all_anvio_contigs.db -p default_matrix -i KeggTaxonomy.txt &> anvi-import-taxonomy.log &

#Annotate gene calls with COG functions:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
#Does not work with v4 due to a bug! Bug fixed in v5 'Margaret'
nohup anvi-run-ncbi-cogs -c Bryobacter_all_anvio_contigs.db --cog-data-dir /geomicro/data2/COMMON/anviCOGs/ -T 16 --search-with diamond --sensitive &> anvi-run-ncbi-cogs.log &

#DJS 8/16/18:

#Profile the contigs database with mapping data from each sample:
#wd: Sample_Bryobacter_all_megahit_coassembly/
nohup sh batch_anvio_profile.sh &> batch_anvio_profile.log &

#DJS 8/22/18:

#Merge the anvio profiles, then bin using CONCOCT:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-merge Sample_42896_profile/PROFILE.db Sample_49615_profile/PROFILE.db Sample_49622_profile/PROFILE.db Sample_49625_profile/PROFILE.db Sample_49629_profile/PROFILE.db Sample_53599_profile/PROFILE.db Sample_53600_profile/PROFILE.db Sample_53601_profile/PROFILE.db -c Bryobacter_all_anvio_contigs.db -o Bryobacter_all_merged_profile --sample-name Coassembly_Bryobacter_all &> anvi-merge.log &

#Get the clustering results for CONCOCT:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
anvi-export-collection -p Bryobacter_all_merged_profile/PROFILE.db -C CONCOCT -O CONCOCT_binning_results

#Export the splits and coverage information from the contigs database:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-export-splits-and-coverages -p Bryobacter_all_merged_profile/PROFILE.db -o Bryobacter_all_merged_profile/ -c Bryobacter_all_anvio_contigs.db -O Bryobacter_all_coassembly --splits-mode &> anvi-export-splits-and-coverages.log &

#DJS 8/23/18:

#Bin the splits using metabat:
#wd: Sample_Bryobacter_all_megahit_coassembly/METABAT/
comics -- nohup metabat -i Bryobacter_all_coassembly-SPLITS.fa -o Bin -a Bryobacter_all_coassembly-COVs.txt --minContig 2500 --cvExt -t 16 --saveCls --onlyLabel -v &> metabat.log &

#Parse metabat output into format that can be import as a collection into the anvio profile:
#wd: /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_42896_53600_megahit_coassembly/METABAT/
sh Metabat_to_anvio_parser.sh

#Import the Metabat bins as a collection into the merged anvio profile:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_anvio_contigs.db -p Bryobacter_all_merged_profile/PROFILE.db -C METABAT ../METABAT/Metabat_binning_results.txt &> anvi-import-collection_metabat.log &

#Start tetraESOM binning process:

#First generate esome files from splits in ANVIO database:
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 2500 -max 10000 -DIR ESOM_2.5-10kb &> esomwrapper_2.5-10kb.log &

#Opened lrn file (Tetra_esom_2500.lrn) in the GUI and normalized the data using Robust ZT transform. Saved normalized data as Tetra_esom_2500.rzt.lrn

#Run the training script to generate the esomtrn command, using recommended values for rows and columns provided by esomWrapper tool:
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_2.5-10kb/ 
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_2500.rzt.lrn -cls Tetra_esom_2500.cls -names Tetra_esom_2500.names -rows 582 -cols 1164 -epochs 50 &> esomTrain.log

#Then run the command. Replace no_norm.lrn file in command with rzt.lrn file:
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_2.5-10kb/
nohup esomtrn --permute --out Tetra_esom_2500.rzt.582x1164e50.wts -b Tetra_esom_2500.rzt.582x1164e50.bm --cls Tetra_esom_2500.cls --lrn Tetra_esom_2500.rzt.lrn --algorithm kbatch --rows 582 --columns 1164 -bmc 73 --start-radius 291 --epochs 50 -k 0.15 --bmsearch standard --dist euc &>> esomTrain.log &

#DJS 9/5/2018:

#Opened finished map, drew bins:
esomana -c Tetra_esom_2500.cls -l Tetra_esom_2500.rzt.lrn -n Tetra_esom_2500.names -b Tetra_esom_2500.rzt.582x1164e50.bm -w Tetra_esom_2500.rzt.582x1164e50.wts

#Not many drawable bins and references did not cluster well.

#Redo ESOM with contig size window of 4-10 kb:
#Removing the smaller contigs should reduce the signal:noise ratio for the larger assembly.
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 4000 -max 10000 -DIR ESOM_4-10kb &> esomwrapper_4-10kb.log &

#Opened lrn file (Tetra_esom_4000.lrn) in the GUI and normalized the data using Robust ZT transform. Saved normalized data as Tetra_esom_4000.rzt.lrn

#Run the training script to generate the esomtrn command, using recommended values for rows and columns from esomWrapper tool:
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_4-10kb/
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_4000.rzt.lrn -cls Tetra_esom_4000.cls -names Tetra_esom_4000.names -rows 412 -cols 824 -epochs 50 &> esomTrain.log &

#Run the generated training command in the background, manually change the .no_norm.lrn file called in the command to call the appropriate lrn file:
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_4-10kb/
nohup esomtrn --permute --out Tetra_esom_4000.rzt.412x824e50.wts -b Tetra_esom_4000.rzt.412x824e50.bm --cls Tetra_esom_4000.cls --lrn Tetra_esom_4000.rzt.lrn --algorithm kbatch --rows 412 --columns 824 -bmc 52 --start-radius 206 --epochs 50 -k 0.15 --bmsearch standard --dist euc &>> esomTrain.log &

#DJS 9/6/2018:

#Binned using VizBin on personal computer with a contig size window of 2500-10,000 bp and using 4mer frequencies.
#Saved the bins as a collection in anvio.

#DJS 9/8/2018:

#ESOM with contig size window of 4-10 kb finished, opened to draw bins.
#wd: /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_4-10kb/
#The map looks much better! A lot of well delineated bins, and references clustered well.
#Saved binned class file as: Tetra_esom_4000_binned.cls

#Generate the bin fastas:
#wd: /omics/HABs/Microcystis/Erie2014_coAssembly/Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_4-10kb/
nohup sh batch_getClassFasta.sh &> getClassFasta.log &

#Generate the binning results summary file for ANVIO and DAStool:
#wd: Sample_Bryobacter_all_megahit_coassembly/ESOM/ESOM_4-10kb/
sh ESOM_binning_results_parser.sh

#Import the ESOM bins into anvio:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_anvio_contigs.db -p Bryobacter_all_merged_profile/PROFILE.db -C ESOM ../ESOM/ESOM_4-10kb/ESOM_binning_results.txt &> anvi-import-collection_ESOM.log &

#Ran DASTool to dereplicate and summarize all binning results:
#wd: Sample_Bryobacter_all_megahit_coassembly/DASTOOL/
nohup das-tool -i CONCOCT_binning_results.txt,ESOM_binning_results.txt,Metabat_binning_results.txt,VizBin_binning_results.txt -c Bryobacter_all_coassembly-SPLITS.fa -o DASTool_bins --search_engine blast --threads 16 &> das-tool.log &

#Remove the period from DASTool bin ids so that it is compatible with Anvio:
#wd: Sample_Bryobacter_all_megahit_coassembly/DASTOOL/ 
sed 's/results\./results_/g' DASTool_bins_DASTool_scaffolds2bin.txt > DASTool_binning_results.txt

#Import the DASTool bins as an anvio collection:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_anvio_contigs.db -p Bryobacter_all_merged_profile/PROFILE.db -C DASTool ../DASTOOL/DASTool_binning_results.txt &> anvi-import-collection_DASTool.log

#Started refining DASTool bins:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
anvi-refine -c Bryobacter_all_anvio_contigs.db -p Bryobacter_all_merged_profile/PROFILE.db -C DASTool -b Metabat_binning_results_120 --server-only -P 8000
## the b flage will be replaced with the ID for each bin to be refined

#DJS 9/10/18:

#Continued refining bins as above, saved the new DASTool collection:

#Moved data to /geomicro/data22/smitdere/

#Summarize the refined bins:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup anvi-summarize -p Bryobacter_all_merged_profile/PROFILE.db -c Bryobacter_all_anvio_contigs.db --init-gene-coverages -C DASTool -o DASTool_anvi_summarize &> anvi-summarize.log &

#DJS 9/12/18:

#CheckM on Refined DASTool bins:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/
nohup checkm lineage_wf HiQual_DASTool_refined_bin_fastas/ HiQual_DASTool_refined_checkm/ -f HiQual_DASTool_refined_checkm/Bryobacter_all_coassembly.checkm.bin_stats.tsv --tab_table -x fa -t 15 &> checkm.log &

#Get placement of bins in checkm reference tree:
#wd: Sample_Bryobacter_all_megahit_coassembly/ANVIO/

################################################################
#Normalized coassembly
################################################################

#DJS 9/13/18:

#No Microcystis bins from either coassembly with completeness above 70%. Trying to normalize kmers in read dataset to 20X to reduce Microcystis coverage to optimum for assembly.
#Concatenate the forward reads:
#wd: Erie2014_coAssembly/
cat Sample_42896/Sample_42896_adtrim_clean_qtrim_fwd.derep.fastq Sample_49615/Sample_49615_adtrim_clean_qtrim_fwd.derep.fastq Sample_49622/Sample_49622_adtrim_clean_qtrim_fwd.derep.fastq Sample_49625/Sample_49625_adtrim_clean_qtrim_fwd.derep.fastq Sample_49629/Sample_49629_adtrim_clean_qtrim_fwd.derep.fastq Sample_53599/Sample_53599_adtrim_clean_qtrim_fwd.derep.fastq Sample_53600/Sample_53600_adtrim_clean_qtrim_fwd.derep.fastq Sample_53601/Sample_53601_adtrim_clean_qtrim_fwd.derep.fastq > Bryobacter_all_cat_fwd.fastq

#Concatenate the reverse reads:
#wd: Erie2014_coAssembly/
nohup cat Sample_42896/Sample_42896_adtrim_clean_qtrim_rev.derep.fastq Sample_49615/Sample_49615_adtrim_clean_qtrim_rev.derep.fastq Sample_49622/Sample_49622_adtrim_clean_qtrim_rev.derep.fastq Sample_49625/Sample_49625_adtrim_clean_qtrim_rev.derep.fastq Sample_49629/Sample_49629_adtrim_clean_qtrim_rev.derep.fastq Sample_53599/Sample_53599_adtrim_clean_qtrim_rev.derep.fastq Sample_53600/Sample_53600_adtrim_clean_qtrim_rev.derep.fastq Sample_53601/Sample_53601_adtrim_clean_qtrim_rev.derep.fastq > Bryobacter_all_cat_rev.fastq &

#DJS 9/14/18:

#Normalize the kmer depth to 20X for each kmer in the concactenated read file with BBnorm while removing kmers with depth of 5X or lower:
#wd: Erie2014_coAssembly/
nohup bbnorm in1=Bryobacter_all_cat_fwd.fastq in2=Bryobacter_all_cat_rev.fastq out=Byrobacter_all_cat_normalized.fastq target=20 min=5 &> bbnorm.log &

#DJS 9/16/18:

#Assemble the normalized, concactenated reads using megahit:
#wd: Erie2014_coAssembly/ 
comics -- nohup megahit --12 Byrobacter_all_cat_normalized.fastq --k-min 21 --k-max 141 --k-step 12 -t 16 -o Sample_Bryobacter_all_normalized_coassembly --verbose &> Sample_Bryobacter_all_normalized_coassembly.log &

#DJS 9/17/18:

#Change assembly headers to have simple def lines:
#wd: Sample_Bryobacter_all_normalized_coassembly/
sed -e 's/\s/_/g' final.contigs.fa > final.contigs.fixed.fa
sed -e 's/=/_/g' final.contigs.fixed.fa > final.contigs.fixed2.fa
sed 's/\./_/g' final.contigs.fixed2.fa > final.contigs.fixed3.fa
rm final.contigs.fixed.fa
rm final.contigs.fixed2.fa
mv final.contigs.fixed3.fa final.contigs.fa

#Index the assembly:
#wd: Sample_Bryobacter_all_normalized_coassembly/
comics -- nohup omics mapping -a final.contigs.fa --index-only &> Sample_Bryobacter_all_normalized_assembly_index.log &

#DJS 9/19/18:

#Merge the coverage files:
#wd: Sample_Bryobacter_all_normalized_coassembly/
comics -- nohup merge-coverage -a final.contigs.fa Mapping_42896/final.contigs.genomeCovBed.tsv Mapping_49615/final.contigs.genomeCovBed.tsv Mapping_49622/final.contigs.genomeCovBed.tsv Mapping_49625/final.contigs.genomeCovBed.tsv Mapping_49629/final.contigs.genomeCovBed.tsv Mapping_53599/final.contigs.genomeCovBed.tsv Mapping_53600/final.contigs.genomeCovBed.tsv Mapping_53601/final.contigs.genomeCovBed.tsv -o final.contigs.merged.cov &> merge-coverage.log &

#Generate an Anvio contigs database:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-gen-contigs-database -f ../final.contigs.fa -o Bryobacter_all_normalized_anvio_contigs.db -n Bryobacter_all_normalized_kmer_coassembly --split-length 10000 --kmer-size 4 &> anvi-gen-contigs-database.log &

#Get get gene calls for ribosomal RNAs:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-run-hmms -c Bryobacter_all_normalized_anvio_contigs.db --num-threads 15 &> anvi-run-hmms.log &

#DJS 9/20/18:

#Annotate the anvio gene calls with NCBI COGs:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-run-ncbi-cogs -c Bryobacter_all_normalized_anvio_contigs.db --cog-data-dir /geomicro/data2/COMMON/anviCOGs/ -T 16 --search-with diamond --sensitive &> anvi-run-ncbi-cogs.log &

#DJS 9/21/18:

#Get amino acid sequences of all gene calls in the assembly for annotation in KEGG GhostKoala:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-get-sequences-for-gene-calls -c Bryobacter_all_normalized_anvio_contigs.db -o Bryobacter_all_normalized_anvio_gene_calls.faa --get-aa-sequences &> anvi-get-aa-sequences-for-gene-calls.log &

#Remove blank entries, caused from bug with annotating Ribosomal RNAs:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
sed '/^>/ {N; /\n$/d}' Bryobacter_all_normalized_anvio_gene_calls.faa > Bryobacter_all_normalized_anvio_gene_calls.cleaned.faa

#Add the prefix "genecall" to each of the gene headers for GhostKoala:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
sed 's/>/>genecall_/g' Bryobacter_all_normalized_anvio_gene_calls.cleaned.faa > Bryobacter_all_normalized_anvio_gene_calls.cleaned.edited.faa

#File is too large for GhostKoala, so split the file into two equal parts:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
split -b 173m Bryobacter_all_normalized_anvio_gene_calls.cleaned.faa Bryobacter_all_normalized_anvio_gene_calls.faa.

#One gene was split between the two file parts, manually pasted back together. Then submit to GhostKoala for annotation.

#DJS 9/23/18:

#Concatenated the two annotation files for both functions and taxonomy
#Parse the output into a format suitable for anvio:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
python /geomicro/data22/smitdere/Erie2014_coAssembly/GhostKoalaParser/KEGG-to-anvio --KeggDB /geomicro/data22/smitdere/Erie2014_coAssembly/GhostKoalaParser/samples/KO_Orthology_ko00001.txt -i user_ko.cat.txt -o KeggAnnotations-AnviImportable.txt

#Import the functional annotations into the contigs database:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-import-functions -c Bryobacter_all_normalized_anvio_contigs.db -i KeggAnnotations-AnviImportable.txt &> anvi-import-functions.log &

#Parse the GhostKoala taxonomy output into a format readable in anvio:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
python /geomicro/data22/smitdere/Erie2014_coAssembly/GhostKoalaParser/GhostKOALA-taxonomy-to-anvio user.out.cat.top KeggTaxonomy.txt

#Import the taxonomic annotations into anvio:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-import-taxonomy-for-genes -c Bryobacter_all_normalized_anvio_contigs.db -p default_matrix -i KeggTaxonomy.txt &> anvi-import-taxonomy.log &

#Profile the contigs database with mapping data from each sample:
#wd: Sample_Bryobacter_all_normalized_coassembly/
nohup sh batch_anvio_profile.sh &> batch_anvio_profile.log & 

#DJS 9/26/18:

#Merge the anvio profiles, then bin using CONCOCT:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-merge Sample_42896_profile/PROFILE.db Sample_49615_profile/PROFILE.db Sample_49622_profile/PROFILE.db Sample_49625_profile/PROFILE.db Sample_49629_profile/PROFILE.db Sample_53599_profile/PROFILE.db Sample_53600_profile/PROFILE.db Sample_53601_profile/PROFILE.db -c Bryobacter_all_normalized_anvio_contigs.db -o Bryobacter_all_normalized_merged_profile --sample-name Coassembly_Bryobacter_all_normalized &> anvi-merge.log &

#DJS 9/27/18:

#Get the CONCOCT binning results:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
anvi-export-collection -p Bryobacter_all_normalized_merged_profile/PROFILE.db -C CONCOCT -O CONCOCT_binning_results

#Export the splits and coverages from the merged profile:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-export-splits-and-coverages -p Bryobacter_all_normalized_merged_profile/PROFILE.db -o Bryobacter_all_normalized_merged_profile/ -c Bryobacter_all_normalized_anvio_contigs.db -O Bryobacter_all_normalized_coassembly --splits-mode &> anvi-export-splits-and-coverages.log &

#Bin the splits using metabat:
#wd: Sample_Bryobacter_all_normalized_coassembly/METABAT/
comics -- nohup metabat -i Bryobacter_all_normalized_coassembly-SPLITS.fa -o Bin -a Bryobacter_all_normalized_coassembly-COVs.txt --minContig 2500 --cvExt -t 16 --saveCls --onlyLabel -v &> metabat.log &

#Parse the metabat output into a format suitable for anvio database:
#wd: Sample_Bryobacter_all_normalized_coassembly/METABAT/
sh Metabat_to_anvio_parser.sh

#Import the Metabat bins as a collection into the Anvio merged profile:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_normalized_anvio_contigs.db -p Bryobacter_all_normalized_merged_profile/PROFILE.db -C METABAT ../METABAT/Metabat_binning_results.txt &> anvi-import-collection-metabat.log

#Start ESOM binning
#Generate ESOM files (contig size window of 4000-10,000 bp):
#wd: Sample_Bryobacter_all_normalized_coassembly/ESOM/
nohup perl /geomicro/data1/COMMON/scripts/wrappers/ESOM/esomWrapper.pl -p ESOM_ref/ -e fa -min 4000 -max 10000 -DIR ESOM_4-10kb &> esomwrapper_4-10kb.log &

#Opened lrn file (Tetra_esom_4000.lrn) in the GUI and normalized the data using Robust ZT transform. Saved normalized data as Tetra_esom_4000.rzt.lrn

#Ran training command scrip to initialize the training command for ESOM (used the suggested number of rows and columns):
#wd: Sample_Bryobacter_all_normalized_coassembly/ESOM/
nohup perl /geomicro/data1/COMMON/scripts/sandbox/esomTrain.pl -lrn Tetra_esom_4000.rzt.lrn -cls Tetra_esom_4000.cls -names Tetra_esom_4000.names -rows 380 -cols 760 -epochs 50 &> esomTrain.log &

#Run the training command (manually change .no_norm.lrn to rzt.lrn file generated earlier):
#wd: Sample_Bryobacter_all_normalized_coassembly/ESOM/
nohup esomtrn --permute --out Tetra_esom_4000.rzt.380x760e50.wts -b Tetra_esom_4000.rzt.380x760e50.bm --cls Tetra_esom_4000.cls --lrn Tetra_esom_4000.rzt.lrn --algorithm kbatch --rows 380 --columns 760 -bmc 48 --start-radius 190 --epochs 50 -k 0.15 --bmsearch standard --dist euc &>> esomTrain.log &

#DJS 9/29/18:

#Binned using VizBin on personal computer with a contig size window of 2500-10,000 bp and using 4mer frequencies.
#Save as a collection in the anvio merged profile.
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_normalized_anvio_contigs.db -p Bryobacter_all_normalized_merged_profile/PROFILE.db -C VizBin VizBin_binning_results.txt &> anvi-import-collection-vizbin.log &

#DJS 9/30/18:

#ESOM map finished, binned contigs and saved class file as Sample_Bryobacter_all_normalized_coassembly/ESOM/ESOM_4-10kb/Tetra_esom_4000_binned.cls/

#Get ESOM bin fastas:
#wd: Sample_Bryobacter_all_normalized_coassembly/ESOM/ESOM_4-10kb/
nohup sh batch_getClassFasta.sh &> getClassFasta.log &

#Parse ESOM output to generate ESOM binning results file for anvio and DASTool:
#wd: Sample_Bryobacter_all_normalized_coassembly/ESOM/ESOM_4-10kb/
sh ESOM_binning_results_parser.sh

#Import ESOM bins into anvio:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_normalized_anvio_contigs.db -p Bryobacter_all_normalized_merged_profile/PROFILE.db -C ESOM ../ESOM/ESOM_4-10kb/ESOM_binning_results.txt &> anvi-import-collection-esom.log &

#Dereplicate and summarize bins using DASTool:
#wd: Sample_Bryobacter_all_normalized_coassembly/DASTOOL/
nohup das-tool -i CONCOCT_binning_results.txt,ESOM_binning_results.txt,Metabat_binning_results.txt,VizBin_binning_results.txt -c Bryobacter_all_normalized_coassembly-SPLITS.fa -o DASTool_bins --search_engine blast --threads 16 &> das-tool.log &

#Remove the peroid from dastool bin names to be suitable for anvio:
#wd: Sample_Bryobacter_all_normalized_coassembly/DASTOOL/
sed 's/results\./results_/g' DASTool_bins_DASTool_scaffolds2bin.txt > DASTool_binning_results.txt

#Import the Dastool bins into anvio:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-import-collection -c Bryobacter_all_normalized_anvio_contigs.db -p Bryobacter_all_normalized_merged_profile/PROFILE.db -C DASTool ../DASTOOL/DASTool_binning_results.txt &> anvi-import-collection-dastool.log &

#Start refining DASTool bins:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
anvi-refine -c Bryobacter_all_normalized_anvio_contigs.db -p Bryobacter_all_normalized_merged_profile/PROFILE.db -C DASTool -b CONCOCT_binning_results_015 --server-only -P 8000
#replace -b flag with the appropriate bin id

#DJS 10/7/2018:

#Summarized the refined DASTool bins:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup anvi-summarize -p Bryobacter_all_normalized_merged_profile/PROFILE.db -c Bryobacter_all_normalized_anvio_contigs.db --init-gene-coverages -C DASTool -o DASTool_anvi_summarize &> anvi-summarize.log &

#DJS 10/8/2018:

#Ran checkm on the refined, HiQuality DASTool bins:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup checkm lineage_wf HiQual_DASTool_refined_bin_fastas/ HiQual_DASTool_refined_checkm/ -f HiQual_DASTool_refined_bin_fastas/Bryobacter_all_normalized_coassembly.checkm.bin_stats.tsv --tab_table -x fa -t 15 &> checkm.log &

#Get placement of bins in checkm genome tree:
#wd: Sample_Bryobacter_all_normalized_coassembly/ANVIO/
nohup checkm tree_qa -f HiQual_DASTool_refined_checkm/Bryobacter_all_normalized_coassembly_checkm_tree_qa.txt -o 2 --tab_table HiQual_DASTool_refined_checkm/ &> checkm_tree_qa.log &

##############################################################
################   Transcript/Bin analysis   #################
##############################################################

#DJS 10/22/18:

#Start QC of 100 um metatranscriptomic datasets:
#Modified script to accomodate single-end sequencing and removed the dereplication step.
#wd: Metatranscriptomes/
nohup sh bbmap_metaT_qc.sh &> transcript_bbmap_qc.log &

#DJS 10/27/2018:

#Copied all bins with Contamination below 5% (excluding heterogeneity) and Completeness above 90% to Redundant_Final_Erie_Bins directory to begin dRep analysis.
#This included new bins generated from coassemblies with megahit and older bins generated from IDBA on individual assemblies.

#DJS 10/28/2018:

#Run dRep on Redundant bin dataset:
#wd: Erie2014_coAssembly/
nohup dRep dereplicate dRep/ -con 200 -g Redundant_Final_Erie_Bins/*.fa &> dRep.log &

#DJS 11/7/18:

#Moved transcripts up one directory (to reorganize and only keep metagenomic assemblies in the Erie2014_coAssembly folder):
#wd: ./
mv Metatranscriptomes/ ../

#DJS 11/16/2018:

#Ran Phylosift on dereplicated bin dataset:
#wd: dRep/dereplicated_genomes/
nohup sh /geomicro/data22/smitdere/Erie2014_coAssembly/phylosift_wrapper.sh &> phylosift_all.log &

#DJS 1/22/19:

#Reran dRep on all bins with Contamination below 5% (excluding heterogeneity) and Completeness above 75% in the Redudant Erie Bins folder, set the ANI cutoff to form a cluster at 97% ID.
#Only use ANI_mf clustering, skip MASH preclustering.
#wd: Erie2014_coAssembly/
nohup dRep dereplicate dRep/ -con 200 -comp 74.5 -g Redundant_Final_Erie_Bins/*.fa --SkipMash --S_ani 0.97 &> dRep.log &

#DJS 1/27/18:
#Get placement of the winning dRep genomes in the CheckM tree:
#wd: dRep/data/checkM/
nohup checkm tree_qa -f checkM_outdir/dRep_checkm_tree_qa.txt -o 2 --tab_table checkM_outdir/ &> dRep_checkm_tree_qa.log &

#DJS 1/29/19:
#wd: /omics/HABs/Microcystis/2014_Erie_EMIRGE/ 
#Combined Emirge sequences from all samples into one file:
cat Sample_42895_EMIRGE/Sample_42895_emirge.fasta Sample_42896_EMIRGE/Sample_42896_emirge.fasta Sample_49613_EMIRGE/Sample_49613_emirge.fasta Sample_49614_EMIRGE/Sample_49614_emirge.fasta Sample_49615_EMIRGE/Sample_49615_emirge.fasta Sample_49618_EMIRGE/Sample_49618_emirge.fasta Sample_49621_EMIRGE/Sample_49621_emirge.fasta Sample_49622_EMIRGE/Sample_49622_emirge.fasta Sample_49624_EMIRGE/Sample_49624_emirge.fasta Sample_49628_EMIRGE/Sample_49628_emirge.fasta Sample_49629_EMIRGE/Sample_49629_emirge.fasta Sample_49632_EMIRGE/Sample_49632_emirge.fasta Sample_49633_EMIRGE/Sample_49633_emirge.fasta Sample_49636_EMIRGE/Sample_49636_emirge.fasta Sample_49638_EMIRGE/Sample_49638_emirge.fasta Sample_53599_EMIRGE/Sample_53599_emirge.fasta Sample_53600_EMIRGE/Sample_53600_emirge.fasta Sample_53601_EMIRGE/Sample_53601_emirge.fasta Sample_53602_EMIRGE/Sample_53602_emirge.fasta > Emirge_all.fasta

#Copied the file to the Erie2014_coAssembly directory:

#Removed sequences that had ambiguous bases, were longer than 1700 bp, and shorter than 1000 bp
#wd: ./
mothur "#screen.seqs(fasta=Emirge_all.fasta, maxambig=0, minlength=1000, maxlength=1700)"

#Cluster the emirge sequences at 97% nucleotide identity:
#wd: ./
nohup vsearch --cluster_fast Emirge_all.good.fasta --consout Emirge_all.good.clustered.fasta --id 0.97 -iddef 4 --uc Emirge_all.good.cluster_results.txt &> vsearch_emirge.log &

#DJS 1/31/19:

#Aligned the QC'd emirge sequences to the Silva reference database (v132):
#wd: ./
mothur "#align.seqs(fasta=Emirge_all.good.clustered.fasta, reference=../Silva_reference_alignment_v132/silva.nr_v132.align, processors=10)"

#Check the clustered sequence centroids for chimeric sequences using pintail and chimera check in mothur:
#wd: ./
mothur "#chimera.pintail(fasta=Emirge_all.good.clustered.align, template=../Silva_reference_alignment_v132/silva.nr_v132.align, processors=10)"

#DJS 9-Feb-19

#Classify QC'd and clustered emirge sequences using Silvia reference database as a template. Used Wang kmer method in mothur:
#wd: ./
nohup mothur "#classify.seqs(fasta=Emirge_all.good.clustered.fasta, template=../Silva_reference_alignment_v132/silva.nr_v132.align, taxonomy=../Silva_reference_alignment_v132/silva.nr_v132.tax, method=wang, probs=T, processors=10, cutoff=80)" &> Emirge_classify_silva.log &

#Next classify with Freshwater Database:
#wd: ./
nohup mothur "#classify.seqs(fasta=Emirge_all.good.clustered.fasta, template=../FreshTrain30Apr2018SILVAv128/FreshTrain30Apr2018SILVAv128.fasta, taxonomy=../FreshTrain30Apr2018SILVAv128/FreshTrain30Apr2018SILVAv128.taxonomy, method=wang, probs=T, processors=10, cutoff=80)" &> Emirge_classify_FreshwaterDB.log &

#DJS 10-Feb-19:

#Rerun checkm on updated dereplicated genomes (Strain ID set to 95%):
#wd: dRep/
nohup checkm lineage_wf dereplicated_genomes/ Dereplicated_Genomes_CheckM_Strain_AAI_95/ -f Dereplicated_Genomes_CheckM_Strain_AAI_95/dereplicated_genomes.checkm.bin_stats.tsv --tab_table -x fa -t 15 --aai_strain 0.95 &> checkm_Strain_AAI_95.log &
nohup checkm tree_qa -f Dereplicated_Genomes_CheckM_Strain_AAI_95/dRep_checkm_tree_qa.txt -o 2 --tab_table Dereplicated_Genomes_CheckM_Strain_AAI_95/ &> dRep_checkm_strainAAI_0.95tree_qa.log &

#DJS 27-Feb-19:

#Dereplicated the redundant bins using ANIb in pyani software:
#wd ./
comics -- nohup average_nucleotide_identity.py -i ./Redundant_Final_Erie_Bins/ -o ./ pyani_Final_Erie_Bins -m ANIb --workers 30 &> pyani_final_erie_bins.log &

#DJS 4-Mar-19:

#Parse ANIb output with Jacob's scripts:
#First format the ANIb output into a table of pairwise ANI and alignment coverage comparisons:
#wd: pyani_Final_Erie_Bins/
python ../Pairwise_Dereplication-master/Convert_Table.py ANIb_percentage_identity.tab ANIb_alignment_coverage.tab

#Dereplicate at the following cutoffs: ANI = 97%, Coverage = 60%:
#wd: pyani_Final_Erie_Bins/
comics -- nohup python3 ../Pairwise_Dereplication-master/Select_Unique_Genomes.py pairwise_long.tsv checkm_results.tsv 0.97 0.6 &> ANIb_select_unique_genomes_ANI_0.97_COV_0.6.log &

#DJS 12-Mar-19:

#Rerun dRep so that bins are dereplicated at 97% ANI with an alignment coverage of 60% or higher:
#wd: ./
nohup dRep dereplicate dRep/ -con 200 -comp 75 -g Redundant_Final_Erie_Bins/*.fa --SkipMash --S_ani 0.97 --cov_thresh 0.6 &> dRep.log &

#DJS 15-Mar-19:

#Copied over final dereplicated pyani bins into a dereplicated genomes subdirectory in the pyani directory:
#wd: pyani_Final_Erie_Bins/
for i in `cat unique_genomes_0.97ani_0.6cov.txt`; do cp ~/Erie2014_coAssembly/Redundant_Final_Erie_Bins/${i}.fa dereplicated_genomes/; done

#Ran checkm on this bin dataset with strain heterogeneity set to 95% AAI:
#wd: pyani_Final_Erie_Bins/
nohup checkm lineage_wf dereplicated_genomes/ Dereplicated_Genomes_CheckM_Strain_AAI_95/ -f Dereplicated_Genomes_CheckM_Strain_AAI_95/pyani_dereplicated_genomes.checkm.bin_stats.tsv --tab_table -x fa -t 15 --aai_strain 0.95 &> checkm_Strain_AAI_95.log &
nohup checkm tree_qa -f Dereplicated_Genomes_CheckM_Strain_AAI_95/dRep_checkm_tree_qa.txt -o 2 --tab_table Dereplicated_Genomes_CheckM_Strain_AAI_95/ &> dRep_checkm_strainAAI_0.95tree_qa.log &

#DJS 16-Mar-19:

#concatenated all the pyani dereplicated bins to blast against the emirge sequences:
#wd: pyani_Final_Erie_Bins/
cat dereplicated_genomes/*.fa > dereplicated_genomes_all_contigs.fa

#Then run the blast:
nohup blastn -query dereplicated_genomes_all_contigs.fa -db ../Emirge_all.good.clustered.fasta -out dereplicated_genomes_vs_emirge_seqs.blastn -outfmt "7 std qlen slen qcovs" &> dereplicated_genomes_vs_emirge_seqs.log &

#DJS 9-Nov-2020:
#Ran CheckM on the Redundant Bin Dataset:
comics -- nohup checkm lineage_wf Redundant_Final_Erie_Bins/ Redundant_Final_Erie_Bins/CheckM/ -f Redundant_Final_Erie_Bins/CheckM/ --tab_table -x fa -t 30 &> Redundant_Final_Erie_Bins_checkM.log &
