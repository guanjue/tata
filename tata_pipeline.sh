script_bin1='/Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/'
script_bin2='/Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq/pipseq_procap_scripts/bin/'
script_bin3='/Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/pipseq_procap_scripts/bin/'
script_bin='/Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tata/'

### split X_mis to different files
cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="0_mis") print $0}' > TATA_midpoints_HSRhee_0_mis.gff
cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="1_mis") print $0}' > TATA_midpoints_HSRhee_1_mis.gff
cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="2_mis") print $0}' > TATA_midpoints_HSRhee_2_mis.gff

### first expand then liftover
#cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="0_mis" && $7=="+") print $1,$4-501,$4+499,$9,$6,$7; else if ($6=="0_mis" && $7=="-") print $1,$4-501,$4+499,$9,$6,$7}' > TATA_midpoints_HSRhee_0_mis_1000_saccer2.bed
#cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="1_mis" && $7=="+") print $1,$4-501,$4+499,$9,$6,$7; else if ($6=="1_mis" && $7=="-") print $1,$4-501,$4+499,$9,$6,$7}' > TATA_midpoints_HSRhee_1_mis_1000_saccer2.bed
#cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="2_mis" && $7=="+") print $1,$4-501,$4+499,$9,$6,$7; else if ($6=="2_mis" && $7=="-") print $1,$4-501,$4+499,$9,$6,$7}' > TATA_midpoints_HSRhee_2_mis_1000_saccer2.bed

#./liftOver TATA_midpoints_HSRhee_0_mis_1000_saccer2.bed SacCer2toSacCer3_modified.chain TATA_midpoints_HSRhee_0_mis_1000_SacCer3.bed unMapped
#./liftOver TATA_midpoints_HSRhee_1_mis_1000_saccer2.bed SacCer2toSacCer3_modified.chain TATA_midpoints_HSRhee_1_mis_1000_SacCer3.bed unMapped
#./liftOver TATA_midpoints_HSRhee_2_mis_1000_saccer2.bed SacCer2toSacCer3_modified.chain TATA_midpoints_HSRhee_2_mis_1000_SacCer3.bed unMapped

#cat TATA_midpoints_HSRhee_0_mis_1000_SacCer3.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2==1000) print $0}' > TATA_midpoints_HSRhee_0_mis_1000.bed
#cat TATA_midpoints_HSRhee_1_mis_1000_SacCer3.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2==1000) print $0}' > TATA_midpoints_HSRhee_1_mis_1000.bed
#cat TATA_midpoints_HSRhee_2_mis_1000_SacCer3.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2==1000) print $0}' > TATA_midpoints_HSRhee_2_mis_1000.bed

### first liftover then expand
cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="0_mis" && $7=="+") print $1,$4,$4+1,$9,$6,$7; else if ($6=="0_mis" && $7=="-") print $1,$4,$4+1,$9,$6,$7}' > TATA_midpoints_HSRhee_0_mis_1_saccer2.bed
cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="1_mis" && $7=="+") print $1,$4,$4+1,$9,$6,$7; else if ($6=="1_mis" && $7=="-") print $1,$4,$4+1,$9,$6,$7}' > TATA_midpoints_HSRhee_1_mis_1_saccer2.bed
cat TATA-midpoints-HSRhee.gff | awk -F '\t' -v OFS='\t' '{if ($6=="2_mis" && $7=="+") print $1,$4,$4+1,$9,$6,$7; else if ($6=="2_mis" && $7=="-") print $1,$4,$4+1,$9,$6,$7}' > TATA_midpoints_HSRhee_2_mis_1_saccer2.bed

$script_bin'liftOver' TATA_midpoints_HSRhee_0_mis_1_saccer2.bed SacCer2toSacCer3_modified.chain TATA_midpoints_HSRhee_0_mis_1_SacCer3.bed unMapped
$script_bin'liftOver' TATA_midpoints_HSRhee_1_mis_1_saccer2.bed SacCer2toSacCer3_modified.chain TATA_midpoints_HSRhee_1_mis_1_SacCer3.bed unMapped
$script_bin'liftOver' TATA_midpoints_HSRhee_2_mis_1_saccer2.bed SacCer2toSacCer3_modified.chain TATA_midpoints_HSRhee_2_mis_1_SacCer3.bed unMapped

cat TATA_midpoints_HSRhee_0_mis_1_SacCer3.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2==1) print $1,"mRNA",".",$2,$3,$5,$6,"TATA-containing",$4}' > TATA_midpoints_HSRhee_0_mis_1.gff
cat TATA_midpoints_HSRhee_1_mis_1_SacCer3.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2==1) print $1,"mRNA",".",$2,$3,$5,$6,"TATA-containing",$4}' > TATA_midpoints_HSRhee_1_mis_1.gff
cat TATA_midpoints_HSRhee_2_mis_1_SacCer3.bed | awk -F '\t' -v OFS='\t' '{if ($3-$2==1) print $1,"mRNA",".",$2,$3,$5,$6,"TATA-containing",$4}' > TATA_midpoints_HSRhee_2_mis_1.gff

cat TATA_midpoints_HSRhee_0_mis_1.gff | awk -F '\t' -v OFS='\t' '{if ($7=="+") print $1,$4-501,$4+499,$9,$6,$7; else if ($6=="0_mis" && $7=="-") print $1,$4-501,$4+499,$9,$6,$7}' > TATA_midpoints_HSRhee_0_mis_1000.bed
cat TATA_midpoints_HSRhee_1_mis_1.gff | awk -F '\t' -v OFS='\t' '{if ($7=="+") print $1,$4-501,$4+499,$9,$6,$7; else if ($6=="0_mis" && $7=="-") print $1,$4-501,$4+499,$9,$6,$7}' > TATA_midpoints_HSRhee_1_mis_1000.bed
cat TATA_midpoints_HSRhee_2_mis_1.gff | awk -F '\t' -v OFS='\t' '{if ($7=="+") print $1,$4-501,$4+499,$9,$6,$7; else if ($6=="0_mis" && $7=="-") print $1,$4-501,$4+499,$9,$6,$7}' > TATA_midpoints_HSRhee_2_mis_1000.bed




#rm -r for_heatmap
mkdir for_heatmap
cp *bed for_heatmap

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c TATA_midpoints_HSRhee_0_mis_1000.bed -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c TATA_midpoints_HSRhee_1_mis_1000.bed -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_proseq_bam/SRR3031844_1.qt.adpt.36.fastq.mapped.sort.bam.bai -c TATA_midpoints_HSRhee_2_mis_1000.bed -s 0 -n 1 -e false -r 2 -p true -a 0 -t 3 -w 0 -h true -m false

mv *_anti.tabular for_heatmap
mv *_sense.tabular for_heatmap

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tbp1_bam_file/tbp1_53319_hs0.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tbp1_bam_file/tbp1_53319_hs0.bam.bai -c TATA_midpoints_HSRhee_0_mis_1000.bed -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tbp1_bam_file/tbp1_53319_hs0.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tbp1_bam_file/tbp1_53319_hs0.bam.bai -c TATA_midpoints_HSRhee_1_mis_1000.bed -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tbp1_bam_file/tbp1_53319_hs0.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/tbp1_bam_file/tbp1_53319_hs0.bam.bai -c TATA_midpoints_HSRhee_2_mis_1000.bed -s 6 -n 1 -e false -r 0 -p false -a 1 -t 3 -w 0 -h true -m false
mv *_combined.tabular for_heatmap

java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c TATA_midpoints_HSRhee_0_mis_1000.bed -s 0 -n 1 -e false -r 0 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c TATA_midpoints_HSRhee_1_mis_1000.bed -s 0 -n 1 -e false -r 0 -p true -a 1 -t 3 -w 0 -h true -m false
java -jar /Users/universe/Documents/bin/cegr-tools/build/dist/TagPileup.jar -b /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam -i /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/pip_seq_bam/sua7_hs0_Tfiltered_56422_56428_merged.bam.bai -c TATA_midpoints_HSRhee_2_mis_1000.bed -s 0 -n 1 -e false -r 0 -p true -a 1 -t 3 -w 0 -h true -m false
mv *_combined.tabular for_heatmap


cd for_heatmap
python $script_bin'sort_by_reads.py' -f TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined.tabular -o TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular
python $script_bin'sort_by_reads.py' -f TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined.tabular -o TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular
python $script_bin'sort_by_reads.py' -f TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined.tabular -o TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular

python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense.tabular -m 1 -s TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense_sort.tabular
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense.tabular -m 1 -s TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense_sort.tabular
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense.tabular -m 1 -s TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense_sort.tabular

python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti.tabular -m 1 -s TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti_sort.tabular
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti.tabular -m 1 -s TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti_sort.tabular
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti.tabular -m 1 -s TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti_sort.tabular

python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_0_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined.tabular -m 1 -s TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_0_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_1_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined.tabular -m 1 -s TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_1_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_2_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined.tabular -m 1 -s TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_2_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular


#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense_sort blue4 0.5
#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense_sort blue4 0.5
#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense_sort blue4 0.5

#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti_sort red4 0.5
#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti_sort red4 0.5
#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti_sort red4 0.5

### merge BR
#composite -dissolve 50 -transparent-color white TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense_sort.png TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti_sort.png TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sort.png
#composite -dissolve 50 -transparent-color white TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense_sort.png TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti_sort.png TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sort.png
#composite -dissolve 50 -transparent-color white TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense_sort.png TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti_sort.png TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sort.png

#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort green4 10
#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort green4 10
#Rscript $script_bin3'heatmap.R' TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort green4 10



### generate composite plot & remove the genes that have more than 1000 reads at a single position
Rscript $script_bin'composite_3.R' TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_0_mis_1000
Rscript $script_bin'composite_3.R' TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_1_mis_1000
Rscript $script_bin'composite_3.R' TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_2_mis_1000

### plot the background??? composite plot as well
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/composite_bg.R TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_0_mis_1000
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/composite_bg.R TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_1_mis_1000
#Rscript /Volumes/MAC_Data/data/labs/pugh_lab/hs_vinesh/composite_bg.R TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_sense_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_SRR3031844_1_readc_anti_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_sua7_hs0_Tfiltered_56422_56428_merged_read1_combined_sort.tabular TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular TATA_midpoints_HSRhee_2_mis_1000

### sort the bed files like the heatmap cdt files for the next fasta file generation
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_0_mis_1000.bed -m 4 -s TATA_midpoints_HSRhee_0_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_0_mis_1000_sort.bed
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_1_mis_1000.bed -m 4 -s TATA_midpoints_HSRhee_1_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_1_mis_1000_sort.bed
python $script_bin'vlookup.py' -t TATA_midpoints_HSRhee_2_mis_1000.bed -m 4 -s TATA_midpoints_HSRhee_2_mis_1000_tbp1_53319_hs0_read1_combined_sort.tabular -n 1 -o TATA_midpoints_HSRhee_2_mis_1000_sort.bed

### generate fasta files for different windows
#rm -r fasta
mkdir fasta
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_0_mis_1000_sort.bed -fo TATA_midpoints_HSRhee_0_mis_1000_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_1_mis_1000_sort.bed -fo TATA_midpoints_HSRhee_1_mis_1000_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_2_mis_1000_sort.bed -fo TATA_midpoints_HSRhee_2_mis_1000_sort.fa -s 
mv *fa fasta


cat TATA_midpoints_HSRhee_0_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-50,($2+$3)/2+50,$4,$5,$6}' > TATA_midpoints_HSRhee_0_mis_100_sort.bed
cat TATA_midpoints_HSRhee_1_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-50,($2+$3)/2+50,$4,$5,$6}' > TATA_midpoints_HSRhee_1_mis_100_sort.bed
cat TATA_midpoints_HSRhee_2_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-50,($2+$3)/2+50,$4,$5,$6}' > TATA_midpoints_HSRhee_2_mis_100_sort.bed

bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_0_mis_100_sort.bed -fo TATA_midpoints_HSRhee_0_mis_100_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_1_mis_100_sort.bed -fo TATA_midpoints_HSRhee_1_mis_100_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_2_mis_100_sort.bed -fo TATA_midpoints_HSRhee_2_mis_100_sort.fa -s 
mv *fa fasta

cat TATA_midpoints_HSRhee_0_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-25,($2+$3)/2+25,$4,$5,$6}' > TATA_midpoints_HSRhee_0_mis_50_sort.bed
cat TATA_midpoints_HSRhee_1_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-25,($2+$3)/2+25,$4,$5,$6}' > TATA_midpoints_HSRhee_1_mis_50_sort.bed
cat TATA_midpoints_HSRhee_2_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-25,($2+$3)/2+25,$4,$5,$6}' > TATA_midpoints_HSRhee_2_mis_50_sort.bed

bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_0_mis_50_sort.bed -fo TATA_midpoints_HSRhee_0_mis_50_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_1_mis_50_sort.bed -fo TATA_midpoints_HSRhee_1_mis_50_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_2_mis_50_sort.bed -fo TATA_midpoints_HSRhee_2_mis_50_sort.fa -s 
mv *fa fasta

cat TATA_midpoints_HSRhee_0_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-27,($2+$3)/2+27,$4,$5,$6}' > TATA_midpoints_HSRhee_0_mis_54_sort.bed
cat TATA_midpoints_HSRhee_1_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-27,($2+$3)/2+27,$4,$5,$6}' > TATA_midpoints_HSRhee_1_mis_54_sort.bed
cat TATA_midpoints_HSRhee_2_mis_1000_sort.bed | awk -F '\t' -v OFS='\t' '{print $1,($2+$3)/2-27,($2+$3)/2+27,$4,$5,$6}' > TATA_midpoints_HSRhee_2_mis_54_sort.bed

bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_0_mis_54_sort.bed -fo TATA_midpoints_HSRhee_0_mis_54_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_1_mis_54_sort.bed -fo TATA_midpoints_HSRhee_1_mis_54_sort.fa -s 
bedtools getfasta -fi /Volumes/MAC_Data/data/labs/pugh_lab/master_ref/procap_pipseq_procap_pip_sort_250/sacCer3.fa -bed TATA_midpoints_HSRhee_2_mis_54_sort.bed -fo TATA_midpoints_HSRhee_2_mis_54_sort.fa -s 
mv *fa fasta

mkdir fasta/dnashape
cd fasta/dnashape
### works only after generate the DNA shape cdt files
Rscript $script_bin'composite_1.R' TATA_midpoints_HSRhee_0_mis_54_sort_MGW.cdt TATA_midpoints_HSRhee_1_mis_54_sort_MGW.cdt TATA_midpoints_HSRhee_2_mis_54_sort_MGW.cdt
Rscript $script_bin'composite_1.R' TATA_midpoints_HSRhee_0_mis_54_sort_PTwist.cdt TATA_midpoints_HSRhee_1_mis_54_sort_PTwist.cdt TATA_midpoints_HSRhee_2_mis_54_sort_PTwist.cdt
Rscript $script_bin'composite_1.R' TATA_midpoints_HSRhee_0_mis_54_sort_Roll.cdt TATA_midpoints_HSRhee_1_mis_54_sort_Roll.cdt TATA_midpoints_HSRhee_2_mis_54_sort_Roll.cdt
Rscript $script_bin'composite_1.R' TATA_midpoints_HSRhee_0_mis_54_sort_HTwist.cdt TATA_midpoints_HSRhee_1_mis_54_sort_HTwist.cdt TATA_midpoints_HSRhee_2_mis_54_sort_HTwist.cdt
cd -

cd ..




