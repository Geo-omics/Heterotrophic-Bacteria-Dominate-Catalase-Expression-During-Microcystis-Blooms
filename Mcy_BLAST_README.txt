##DJS 18-Feb-21##

#Reran mcy blasts with read datasets run through the updated bbmap QC pipeline.

#Copied over reference databases used originally back in 2018
#These contained mcy genes (entire operon) from Microcystis reference genomes in IMG as well as some from WLE 2014 assemblies
#Confirmed by blasting to NCBI that the genes from the assembly are indeed Microcystis. (hit to Microcystis mcy at 99% ID and mcy from Planktothrix and Nodularia at ~70%)
#Also included 16S V4 amplicon sequences that were abundant from Berry et al. 2017 to use as denominator for % toxic Microcystis counts.

#This command blasts the short reads agains the custom databases and filters the blast hits based on set percent identity cutoffs:
nohup bash batch_metaG_blast.sh &> batch_metaG_blast.log &

#This command tallies the read counts to each reference:
#Run it in the metaG_blasts/ directory
nohup bash get_read_counts.sh &> get_read_counts.log &

##DJS 19-Feb-21##
#We are getting 200% toxic Microcystis cells in one sample (49615), which doesn't make sense because mcyD is a single copy gene in all known Microcystis genomes

#Extracting the reads that mapped to the mcyD genes at the set thresholds and aligning to Teal's large database to see if any reads unambiguously map at higher % similarity to other organisms or functions.
#First I need to get the blast output for each sample with no % ID filtering:
#Did this by running the loop in the metaG_blasts/ directory
for i in *_vs_mcy.blastn; do
        echo "Running postblast for $i"
        perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -e 1e-5 -b $i -o ${i}.no_id_cutoff.postblast
done

##DJS 21-Feb-21

#Then run top5 script so each read is only listed once:
for i in *.no_id_cutoff.postblast; do
        echo "Running top hits for $i"
        perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b $i -o $i
done

#Extract reads that mapped to the mcy genes to run through Teal's pipeline:
nohup bash extract_mcy_reads.sh &> extract_mcy_reads.log &

#DJS 24-Feb-21:
#Align the reads to Teal's database using diamond blastx:
nohup bash Diamond_map_Teal_URDB.sh &> Diamond_map_Teal_URDB.log &

#DJS 25-Feb-21:
#Get gene hits and read count data from the diamond alignment:
nohup bash process_diamond_URDB.sh &> process_diamond_URDB.log &

#Also try aligning the reads to URDB with MUSCATO:
nohup bash MUSCATO_map_Teal_URDB.sh &> MUSCATO_map_Teal_URDB.log &

##DJS 3-Mar-21:
#Extract the read counts and URDB hits from the muscato alignment:
nohup bash process_muscato_URDB.sh &> process_muscato_URDB.log &

##Now I want to compare reads that mapped to the V4 reference in Colleen's BBMap analysis to those that mapped in Blast (note that our BLAST results are now identical):

#I made a new directory to do this comparison in:
/geomicro/data22/smitdere/Rerun_mcy_blasts/Colleen_BBMap_v_URDB

#Make a symbolic link to Colleen's BBMap bam file:
ln -s /geomicro/data21/cyancey/LE_2014_MetaG_100x/mcyGenotypes/Special_42895/v4_NoMM_sorted.bam
ln -s /geomicro/data21/cyancey/LE_2014_MetaG_100x/mcyGenotypes/Special_42895/mcy_MM_sorted.bam

##First get a list of reads that mapped with BBmap:
#Index Colleen's bam file:
samtools index v4_NoMM_sorted.bam 

#Extract all the reads that mapped to the Microcystis V4 as a separate bam file:
samtools view -b v4_NoMM_sorted.bam Microcystis_Berry_2017_Otu00005 > v4_NoMM_Microcystis_OTU.bam

#Get a list of reads that mapped to the Microcystis V4:
samtools view -F 4 v4_NoMM_Microcystis_OTU.bam | cut -f1 > bbmap_v4_readlist.txt

#Back one directory up, get a list of reads that BLAST mapped to the Microcystis V4 sequence:
#wd = /geomicro/data22/smitdere/Rerun_mcy_blasts/
nohup bash extract_16S_reads.sh &> extract 16S reads.log &

## DJS 4-Mar-2021 ## 

#Move back into the Colleen_BBMap_v_URDB directory, then get a list of reads in Sample 42896 that mapped to the Microcystis V4 sequence with bbmap but not in BLAST:
#First concactenate the reads that mapped via BLAST
#Need to append the forward and reverse reads markers to the end of each line in the BLAST mapped read lists first:
#Working directory was in the Rerun mcy blast folder:
grep ">" 42895_fwd_mapped_Microcystis_V4.fasta > 42895_fwd_mapped_Microcystis_V4_readlist.txt
grep ">" 42895_rev_mapped_Microcystis_V4.fasta > 42895_rev_mapped_Microcystis_V4_readlist.txt

#Then concactenate into one file back in the Collen BBmap directory:
cat ../42895_fwd_mapped_Microcystis_V4_readlist.txt ../42895_rev_mapped_Microcystis_V4_readlist.txt > 42895_cat_mapped_Microcystis_V4_readlist.txt

#Remove the fasta sequence header >
sed -i 's/^>//g' 42895_cat_mapped_Microcystis_V4_readlist.txt 

#Sort both the BLAST and BBmap read lists:
sort 42895_cat_mapped_Microcystis_V4_readlist.txt > 42895_cat_mapped_Microcystis_V4_readlist_sorted.txt 
sort bbmap_v4_readlist.txt > bbmap_v4_readlist_sorted.txt

#Get a list of reads that mapped in bbmap but not in BLAST:
comm -23 bbmap_v4_readlist_sorted.txt 42895_cat_mapped_Microcystis_V4_readlist_sorted.txt > bbmap_v4_only_readlist.txt

#Extract the reads that only mapped to the Microcystis V4 in bbmap from the fasta file:
perl /geomicro/data1/COMMON/scripts/SeqTools/extractSeqs.pl -l bbmap_v4_only_readlist.txt -f ../Sample_42895_adtrim_clean_qtrim_fwd.derep.fasta -o 42895_fwd_Microcystis_V4_bbmap_mapped_only.fasta
#Note that all sequences were found with this command, which means that no reverse reads mapped uniquely with bbmap.

#Just out of curiosity, how many reads mapped uniquely with either BBmap or BLAST, and how many reads mapped with both aliners?
#Get the number of reads that mapped to just bbmap:
comm -23 bbmap_v4_readlist_sorted.txt 42895_cat_mapped_Microcystis_V4_readlist_sorted.txt | wc -l    #### 109 mapped only with bbmap

#Get the number of reads that mapped to just BLAST:
comm -13 bbmap_v4_readlist_sorted.txt 42895_cat_mapped_Microcystis_V4_readlist_sorted.txt | wc -l    #### 1 read mapped only with BLAST (a reverse read)

#Get the number of reads that mapped with both methods:
comm -12 bbmap_v4_readlist_sorted.txt 42895_cat_mapped_Microcystis_V4_readlist_sorted.txt | wc -l    #### 10 reads mapped with both methods

#Align the reads that only mapped with bbmap to Silva v138 using BLAST:
nohup blastn -query 42895_fwd_Microcystis_V4_bbmap_mapped_only.fasta -db /geomicro/data22/smitdere/Silva_reference_alignment_v138/silva.nr_v138.align.degapped.fixed_headers -out 42895_bbmap_only_Microcystis_V4_vs_Silva_v138.blastn -outfmt "7 std qlen slen qcovs" &> 42895_bbmap_only_Microcystis_V4_vs_Silva_v138.log &

#Filter output so that it is in tabular format and only includes reads that aligned at 97% ID or higher:
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -p 97 -s 50 -e 1e-5 -b 42895_bbmap_only_Microcystis_V4_vs_Silva_v138.blastn -o 42895_bbmap_only_Microcystis_V4_vs_Silva_v138.blastn.postblast

#Concatenate the fwd and rev reads that mapped to the Microcystis V4 sequence with BLAST:
cat ../42895_fwd_mapped_Microcystis_V4.fasta ../42895_rev_mapped_Microcystis_V4.fasta > 42895_cat_mapped_Microcystis_V4.fasta

#Align the reads that mapped with BLAST to Silva v138:
nohup blastn -query 42895_cat_mapped_Microcystis_V4.fasta -db /geomicro/data22/smitdere/Silva_reference_alignment_v138/silva.nr_v138.align.degapped.fixed_headers -out 42895_BLAST_mapped_Microcystis_V4_vs_Silva_v138.blastn -outfmt "7 std qlen slen qcovs" &> 42895_BLAST_mapped_Microcystis_V4_vs_Silva_v138.log &

#Filter the output so that it is in tabular format and only includes reads that aligned at 97% ID or higher:
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -p 97 -s 50 -e 1e-5 -b 42895_BLAST_mapped_Microcystis_V4_vs_Silva_v138.blastn -o 42895_BLAST_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast

#Get the read counts for each BLAST to Silva v138:
nohup bash get_read_counts.sh &> get_read_counts.log &

### DJS 10-Mar-2021 ###

## Next I want to check if there are any reads that are mapping to Microcystis below 97% ID that are true Microcystis hits:  

#I first refiltered the blast results so that the results are tabulated, but there is no percent identity cutoff. Still included an e-value cutoff:  
for i in *_vs_16S.blastn; do
        echo "Running postblast for $i"
        perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -e 1e-5 -b $i -o ${i}.no_id_cutoff.postblast
done

#Then I ran the top5.pl script, so that each read is only counted once:  
for i in *_vs_16S.blastn.no_id_cutoff.postblast; do
        echo "Getting top hits for $i"
	perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b $i -o $i
done

#Extract the reads that mapped to the 16S V4 of Microcystis, no % ID cutoffs:
nohup bash extract_16S_reads_no_pid.sh &> extract_16S_reads_no_pid.log &

#Then, aligned these reads to Silva reference database:  
nohup blastn -query 42895_fwd_mapped_Microcystis_V4_no_pid.fasta -db /geomicro/data22/smitdere/Silva_reference_alignment_v138/silva.nr_v138.align.degapped.fixed_headers -out 42895_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn -outfmt "7 std qlen slen qcovs" &> 42895_fwd_mapped_Microcystis_V4_vs_Silva_v138.log &
nohup blastn -query 42895_rev_mapped_Microcystis_V4_no_pid.fasta -db /geomicro/data22/smitdere/Silva_reference_alignment_v138/silva.nr_v138.align.degapped.fixed_headers -out 42895_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn -outfmt "7 std qlen slen qcovs" &> 42895_rev_mapped_Microcystis_V4_vs_Silva_v138.log &

nohup blastn -query 49615_fwd_mapped_Microcystis_V4_no_pid.fasta -db /geomicro/data22/smitdere/Silva_reference_alignment_v138/silva.nr_v138.align.degapped.fixed_headers -out 49615_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn -outfmt "7 std qlen slen qcovs" &> 49615_fwd_mapped_Microcystis_V4_vs_Silva_v138.log &
nohup blastn -query 49615_rev_mapped_Microcystis_V4_no_pid.fasta -db /geomicro/data22/smitdere/Silva_reference_alignment_v138/silva.nr_v138.align.degapped.fixed_headers -out 49615_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn -outfmt "7 std qlen slen qcovs" &> 49615_rev_mapped_Microcystis_V4_vs_Silva_v138.log &

#Filter each output file so that it is tabular and only includes hits with significant e-values identity or greater.
#Focus on the 42895 sample (which had the biggest disagreement with BBMap) and the 49615 sample (which had 200% toxic Microcystis):  
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b 42895_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn -o 42895_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b 42895_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn -o 42895_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast

perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b 49615_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn -o 49615_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast
perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -s 50 -e 1e-5 -b 49615_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn -o 49615_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast

### DJS 12-Mar-21 ###

#Get just the best hits to Silva in the 49615 sample:  
perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b 49615_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast -o 49615_fwd_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast.tophits
perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b 49615_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast -o 49615_rev_mapped_Microcystis_V4_vs_Silva_v138.blastn.postblast.tophits

## DJS 15-Mar-21 ##

#Created a new directory to compare the URDB gene entries to IMG gene calls in closed Microcystis genomes:
mkdir URDB_hits_vs_Closed_Microcystis_Genomes

#Blast the URDB genes that had significant hits to reads that mapped to mcy genes in several Microcystis isolates to the gene calls.
#The annotations of the URDB genes could be very sparse or unspecific, so I want to see if they are only including microcystin biosynthesis genes or if they are truly genes from other BGCs.
nohup bash batch_blast_URDB_hits.sh &> batch_blast_URDB_hits.log &

###Make a fragment recruitment plot of the mcy BLAST for sample 49615:  
#Working Directory = metaG_blasts/
#First, need to filter the Blast output so that I am only plotting reads that matched at least 95% nuceotdie ID and 80% of read length:
awk -F '\t' '{if($15 >= 80 && $3 >= 95) {print $0}}' Sample_49615_adtrim_clean_qtrim_fwd_vs_mcy.blastn.no_id_cutoff.postblast > Sample_49615_adtrim_clean_qtrim_fwd_vs_mcy.blastn.postblast.filter_for_frag_plot
awk -F '\t' '{if($15 >= 80 && $3 >= 95) {print $0}}' Sample_49615_adtrim_clean_qtrim_rev_vs_mcy.blastn.no_id_cutoff.postblast > Sample_49615_adtrim_clean_qtrim_rev_vs_mcy.blastn.postblast.filter_for_frag_plot
#Remove the last three columns of data, because Robert's script is expecting only 12 (the default number for a BLAST job, I told BLAST to add extra data columns):
cut --complement -f13,14,15 Sample_49615_adtrim_clean_qtrim_fwd_vs_mcy.blastn.postblast.filter_for_frag_plot > Sample_49615_adtrim_clean_qtrim_fwd_vs_mcy.blastn.postblast.filter_for_frag_plot.cut
cut --complement -f13,14,15 Sample_49615_adtrim_clean_qtrim_rev_vs_mcy.blastn.postblast.filter_for_frag_plot > Sample_49615_adtrim_clean_qtrim_rev_vs_mcy.blastn.postblast.filter_for_frag_plot.cut
#Make the fragment recruitment plot:
comics plot-blast-frag-cov -r ../Microcystis_mcy_cluster.fasta Sample_49615_adtrim_clean_qtrim_fwd_vs_mcy.blastn.postblast.filter_for_frag_plot.cut
comics plot-blast-frag-cov -r ../Microcystis_mcy_cluster.fasta Sample_49615_adtrim_clean_qtrim_rev_vs_mcy.blastn.postblast.filter_for_frag_plot.cut

#These plots are hard to interpret with the analysis that I did, mapping simultaneously with multiple references. I'm going to remap to the mcy gene cluster from 1 or 2 strains then make a frag plot for each strain:

#Map Sample 49615 to the two references from single strains with BLAST:
nohup blastn -query Sample_49615_adtrim_clean_qtrim_fwd.derep.fasta -db Microcystis_DIANCHI905_mcy.fasta -out Sample_49615_adtrim_clean_qtrim_fwd_vs_DIANCHI905_mcy.blastn -outfmt "7 std qlen slen qcovs" &> Sample_49615_adtrim_clean_qtrim_fwd_vs_DIANCHI905_mcy.log &
nohup blastn -query Sample_49615_adtrim_clean_qtrim_rev.derep.fasta -db Microcystis_DIANCHI905_mcy.fasta -out Sample_49615_adtrim_clean_qtrim_rev_vs_DIANCHI905_mcy.blastn -outfmt "7 std qlen slen qcovs" &> Sample_49615_adtrim_clean_qtrim_rev_vs_DIANCHI905_mcy.log &
nohup blastn -query Sample_49615_adtrim_clean_qtrim_rev.derep.fasta -db Microcystis_PCC7806_mcy.fasta -out Sample_49615_adtrim_clean_qtrim_rev_vs_PCC7806_mcy.blastn -outfmt "7 std qlen slen qcovs" &> Sample_49615_adtrim_clean_qtrim_rev_vs_PCC7806_mcy.log &
nohup blastn -query Sample_49615_adtrim_clean_qtrim_fwd.derep.fasta -db Microcystis_PCC7806_mcy.fasta -out Sample_49615_adtrim_clean_qtrim_fwd_vs_PCC7806_mcy.blastn -outfmt "7 std qlen slen qcovs" &> Sample_49615_adtrim_clean_qtrim_fwd_vs_PCC7806_mcy.log &

#Tabulate the BLAST output and remove hits below 95% ID and e-value greater than 1e-5:
for i in *_vs_DIANCHI905_mcy.blastn; do perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -p 95 -e 1e-5 -b $i -o ${i}.postblast; done
for i in *_vs_PCC7806_mcy.blastn; do perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -p 95 -e 1e-5 -b $i -o ${i}.postblast; done

#Get just the top hits for each read (only count each read once):
for i in *_vs_DIANCHI905_mcy.blastn.postblast; do perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b $i -o ${i}.tophits; done
for i in *_vs_PCC7806_mcy.blastn.postblast; do perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b $i -o ${i}.tophits; done

#Filter the BLAST output so that we are only including alignments that cover at a minimum 80% of read length:
awk -F '\t' '{if($15 >= 80 && $3 >= 95) {print $0}}' Sample_49615_adtrim_clean_qtrim_fwd_vs_DIANCHI905_mcy.blastn.postblast.tophits > Sample_49615_adtrim_clean_qtrim_fwd_vs_DIANCHI905_mcy.blastn.postblast.filter_for_frag_plot
awk -F '\t' '{if($15 >= 80 && $3 >= 95) {print $0}}' Sample_49615_adtrim_clean_qtrim_rev_vs_DIANCHI905_mcy.blastn.postblast.tophits > Sample_49615_adtrim_clean_qtrim_rev_vs_DIANCHI905_mcy.blastn.postblast.filter_for_frag_plot
awk -F '\t' '{if($15 >= 80 && $3 >= 95) {print $0}}' Sample_49615_adtrim_clean_qtrim_fwd_vs_PCC7806_mcy.blastn.postblast.tophits > Sample_49615_adtrim_clean_qtrim_fwd_vs_PCC7806_mcy.blastn.postblast.filter_for_frag_plot
awk -F '\t' '{if($15 >= 80 && $3 >= 95) {print $0}}' Sample_49615_adtrim_clean_qtrim_rev_vs_PCC7806_mcy.blastn.postblast.tophits > Sample_49615_adtrim_clean_qtrim_rev_vs_PCC7806_mcy.blastn.postblast.filter_for_frag_plot

#Concatenate the output from mapping fwd and rev reads to each reference:
cat Sample_49615_*_vs_DIANCHI905_mcy.blastn.postblast.filter_for_frag_plot > Sample_49615_cat_DIANCHI905_mcy.blastn
cat Sample_49615_*_vs_PCC7806_mcy.blastn.postblast.filter_for_frag_plot > Sample_49615_cat_PCC7806_mcy.blastn

#Remove the last three columns of data for Robert's script:
cut --complement -f13,14,15 Sample_49615_cat_DIANCHI905_mcy.blastn > Sample_49615_cat_DIANCHI905_mcy.blastn.cut
cut --complement -f13,14,15 Sample_49615_cat_PCC7806_mcy.blastn > Sample_49615_cat_PCC7806_mcy.blastn.cut

#Make the frag recruitment plots:
comics plot-blast-frag-cov -r Microcystis_DIANCHI905_mcy.fasta Sample_49615_cat_DIANCHI905_mcy.blastn.cut
comics plot-blast-frag-cov -r Microcystis_PCC7806_mcy.fasta Sample_49615_cat_PCC7806_mcy.blastn.cut

### DJS 16-Mar-2021 ###

#Map the reads that mapped to the mcy database at 100-95 % identity to nr:
nohup blastn -query 49615_fwd_mapped_mcy_100_95.fasta -db /geomicro/data9/flux/reference-data/blast/nt -out 49615_fwd_mapped_mcy_100_95_vs_nt.blastn -outfmt "7 std qlen slen qcovs" &> 49615_fwd_mapped_mcy_100_95_vs_nt.log &
nohup blastn -query 49615_rev_mapped_mcy_100_95.fasta -db /geomicro/data9/flux/reference-data/blast/nt -out 49615_rev_mapped_mcy_100_95_vs_nt.blastn -outfmt "7 std qlen slen qcovs" &> 49615_rev_mapped_mcy_100_95_vs_nt.log &

#I'm concerned that mapping all the hits to NCBI is going to take a very long time, so I'm also going to downsample by 50% and map to nr:
python3 downsample.py 49615_fwd_mapped_mcy_100_95.fasta --single -f 0.5
python3 downsample.py 49615_rev_mapped_mcy_100_95.fasta --single -f 0.5

nohup blastn -query 49615_fwd_mapped_mcy_100_95.0.5.fasta -db /geomicro/data9/flux/reference-data/blast/nt -out 49615_fwd_mapped_mcy_100_95_DS50_vs_nt.blastn -outfmt "7 std qlen slen qcovs" &> 49615_fwd_mapped_mcy_100_95_DS50_vs_nt.log &
nohup blastn -query 49615_rev_mapped_mcy_100_95.0.5.fasta -db /geomicro/data9/flux/reference-data/blast/nt -out 49615_rev_mapped_mcy_100_95_DS50_vs_nt.blastn -outfmt "7 std qlen slen qcovs" &> 49615_rev_mapped_mcy_100_95_DS50_vs_nt.log &

### DJS 17-Mar-2021 ###

#Filter the output to only include hits with an e-value at least 1e-5. Also tabulates output:
for i in *_vs_nt.blastn; do perl /geomicro/data1/COMMON/scripts/BlastTools/postBlast.pl -e 1e-5 -b $i -o ${i}.postblast; done

#Create an output that only includes the top hits for each read:
for i in *_vs_nt.blastn.postblast; do perl /geomicro/data1/COMMON/scripts/BlastTools/top5.pl -t 1 -b $i -o ${i}.tophits; done

### DJS 18-Mar-2021 ###

#Get the NCBI accession numbers for each hit:
nohup bash get_NCBI_accessions.sh &> get_NCBI_accessions.log &

### DJS 19-Mar-2021 ###

#Concactenate the annotation lists and the blast files:
cat 49615_fwd_mapped_mcy_100_95_vs_nt.blastn.postblast.annotation_list.txt 49615_rev_mapped_mcy_100_95_vs_nt.blastn.postblast.annotation_list.txt | sort -u > 49615_cat_mapped_mcy_100_95_vs_nt.blastn.postblast.annotation_list.txt

#Get the number of reads mapped to each gene:
nohup bash get_read_counts_parallelized.sh &> get_read_counts_parallelized.log &

## DJS 26-Mar-2021 ##
#Decided to use these parameters/references for the final alignments and counts:
#95% identity cutoff for 16S V4, but including all cyanobacteria sequences in the reference file so that the BLAST is more competetitve.
#95% identity cutoff for mcy, including all non-Microcystis references that mapped in other cyanobacteria in NCBI in the custom BLAST database.

#Got a list of annotated BLAST hits in NCBI (excluding Microcystis) and the count of reads that they recruited from the related RMarkdown analysis file:
#Filename is Combined49615_vs_NCBI_counts_no_Microcystis.txt
#Extracted just the accession numbers from this file:
awk -F '\t' 'NR>1 {print $2}' Combined49615_vs_NCBI_counts_no_Microcystis.txt > Accession_numbers_to_download.txt

#get a fasta file of these sequences from NCBI:
nohup bash get_NCBI_Seqs.sh &> get_NCBI_Seqs.log &

#Add a new line before each sequence for legibility:
sed -i 's/>/\n&/g' Non_Microcystis_mcy_hits.fasta 

#Get rid of the first empty line:
tail -n +2 Non_Microcystis_mcy_hits.fasta > Non_Microcystis_mcy_hits.fasta.temp && mv Non_Microcystis_mcy_hits.fasta.temp Non_Microcystis_mcy_hits.fasta

#concactenate the Microcysits mcy fasta from colleen and the the non-Microcystis hits from NCBI:
cat mcyGenes_IMG.fna Non_Microcystis_mcy_hits.fasta > Competitive_mcy_database.fasta

#Create the blast database:
makeblastdb -in Competitive_mcy_database.fasta -dbtype nucl

#Extract out all the cyanobacteria reads from the Silva v138 database for the 16S mapping:
#Extracting the sequences from the database that was aligned to the amplicon sequences and trimmed to only include the V4 region sequences:
#The sed commands remove awkward tabs in the fasta headers and degap the sequences
grep -A 1 "Bacteria;Cyanobacteria;" /geomicro/data22/smitdere/Microcystis_Cultures/silva.nr_v138.v4.align | sed 's/-//g' | sed 's/\t/_/g' > Silva.v4.cyanobacteria.fasta

#There are spaces in the fasta header that will make parsing annoying later, so I'm replacing all the white space in the header names with underscores:
sed -E 's/(.) (.)/\1_\2/g' Silva.v4.cyanobacteria.fasta | sed 's/;_//g' > temp.txt && mv temp.txt Silva.v4.cyanobacteria.fasta

#A few sequences did not align properly and are causing problems for the blast database so I'm removing them:
sed '4893,4894d' Silva.v4.cyanobacteria.fasta > temp.txt && mv temp.txt Silva.v4.cyanobacteria.fasta
sed '6872,6873d' Silva.v4.cyanobacteria.fasta > temp.txt && mv temp.txt Silva.v4.cyanobacteria.fasta
sed '8068,8069d' Silva.v4.cyanobacteria.fasta > temp.txt && mv temp.txt Silva.v4.cyanobacteria.fasta
sed '13719,13720d' Silva.v4.cyanobacteria.fasta > temp.txt && mv temp.txt Silva.v4.cyanobacteria.fasta
sed '14594,14595d' Silva.v4.cyanobacteria.fasta > temp.txt && mv temp.txt Silva.v4.cyanobacteria.fasta

#Create the BLAST database:
makeblastdb -in Silva.v4.cyanobacteria.fasta -dbtype nucl

#Run the blast with the new databases and cutoffs:
nohup bash batch_metaG_blast_updated.sh &> batch_metaG_blast_updated.log &

#Remove all the hits with alignment coverage lower than 80% of read length:
for i in *.postblast; do awk -F "\t" '{if($15 >= 80) {print $0}}' $i > ${i}.aln_len_filter; done

## DJS 2-Apr-21 ##

#Testing Robert's script to make sure counts are as expected:
#Strategy is to create a file that is a subset of the first 50 querys in one of the BLAST files. Chose 16S BLAST for Sample 42895:
#Note, need head to spit out 51 lines because of file header:
awk -F '\t' '{print $1}' Sample_42895_adtrim_clean_qtrim_fwd_vs_16S.blastn.postblast | sort -u | head -n 51 > Sample_42895_16S_first_50_queries.txt

#Extract all the hits for the first 50 queries using this for loop:
for i in `cat Sample_42895_16S_first_50_queries.txt`; do awk -F '\t' '{if($1 == "'$i'") {print $0}}' Sample_42895_adtrim_clean_qtrim_fwd_vs_16S.blastn.postblast >> Sample_42895_test_blast.txt; done

#Run Robert's script on the alignment length filtered files:
/geomicro/data2/heinro/work/blast_hit_counts/blast_hit_counts *.aln_len_filter

#Run the script on the smaller test file to make sure reads are being appropriately counted:
/geomicro/data2/heinro/work/blast_hit_counts/blast_hit_counts Sample_42895_test_blast.txt

#Test again for an mcy blast:
#Get the first 50 queries
awk -F '\t' '{print $1}' Sample_42895_adtrim_clean_qtrim_fwd_vs_mcy.blastn.postblast | sort -u | head -n 51 > Sample_42895_mcy_first_50_queries.txt

#Extract all the hits for the first 50 queries using this for loop:
for i in `cat Sample_42895_mcy_first_50_queries.txt`; do awk -F '\t' '{if($1 == "'$i'") {print $0}}' Sample_42895_adtrim_clean_qtrim_fwd_vs_mcy.blastn.postblast >> Sample_42895_mcy_test_blast.txt; done

#Run the script to get the counts for the small mcy test BLAST file:
/geomicro/data2/heinro/work/blast_hit_counts/blast_hit_counts Sample_42895_mcy_test_blast.txt

## DJS 5-Apr-21 ##

#Rerunning Robert's script with the new flags to get the counts broken down my mcy gene:
/geomicro/data2/heinro/work/blast_hit_counts/blast_hit_counts *.aln_len_filter --mcy-break-down --save-map map.txt

#Rename to prevent overwriting later:
mv out.none.csv Read_couts_mcy_cyano16S.aln_len_filter.none.csv
mv out.all.csv Read_couts_mcy_cyano16S.aln_len_filter.all.csv
mv out.part.csv Read_couts_mcy_cyano16S.aln_len_filter.part.csv
