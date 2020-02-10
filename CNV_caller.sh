##Files located in
  ~/GWAS_project/fastp_trimmed_A6

##Files needed go through scripts
##Unique
for Filename in *trimmedaligned.bam; do name=`ls $Filename | cut -d"." -f1`;
~/programs/sambamba-0.6.8 view -t 8 -h -f bam $Filename -o $name"_unique.bam" ;
done
##dup_removed

for Filename in *.bam; do name=`ls $Filename | cut -d"." -f1` ;
~/programs/sambamba-0.6.8  markdup -r -t 8 $Filename $name".DUPremoved.bam";
done
## add_replace

for Filename in *unique.DUPremoved.bam; do name=`ls $Filename | cut -d"." -f1`;
java -jar ~/programs/picard.jar AddOrReplaceReadGroups I=$Filename
O=$name"_fixed.bam" RGLB=libl RGPL=illumina RGPU=until RGSM=20 ; done
## sort
for filename in *fixed.bam; do samtools sort $filename -o $filename'.sorted.bam';
done

## Change sorted to sorted1
for i in *.sorted.bam;
 do name=`echo ${i} | cut -d'.' -f1,2`;
 mv ${i} $name.sorted1.bam; done

## Append the header for bam files
for i in *.bam.sorted1.bam;
 do name=`echo ${i} | cut -d'.' -f1,2`;
 code=`echo ${i} | cut -d'_' -f1,2`;
 samtools view -H ${i} > header.sam;
 sed "s/SM:20/SM:$code/g" header.sam > header_corrected.sam;
 samtools reheader header_corrected.sam ${i%.*}.bam > $name.sorted.bam; done

rm *.sorted1*


## experiment with re-indexing sucessful
for i in *.sorted.bam; do
samtools index -b ${i}; done}

##CNV calling in GATK

##Generate consectutive bins of 1000bp from the reference
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar PreprocessIntervals
-R revisedAssemblyUnmasked.fa
--bin-length 1000
--padding 250
-O preprocessed_intervals.interval_list

##Annotate intervals with GC content
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar AnnotateIntervals
-R revisedAssemblyUnmasked.fa
-L preprocessed_intervals.interval_list
--interval-merging-rule OVERLAPPING_ONLY
-O annotated_intervals.tsv

## Filter intervals produced by preprocessed intervals
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar FilterIntervals
-L preprocessed_intervals.interval_list
--annotated-intervals annotated_intervals.tsv
--interval-merging-rule OVERLAPPING_ONLY
-O filtered_intervals.interval_list

##  gatk CollectReadCounts
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar CollectReadCounts
-I 5270-S2_S2_trimmedaligned_unique_fixed.bam.sorted1.bam
-L filtered_intervals.interval_list
--interval-merging-rule OVERLAPPING_ONLY
-O 5270-S2_S2.counts.hdf5

for i in *.sorted.bam;
 do code=`echo ${i} | cut -d'_' -f1,2`;
 ../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar CollectReadCounts
 -I ${i}
 -L filtered_intervals.interval_list
 --interval-merging-rule OVERLAPPING_ONLY
 -O $code.counts.hdf5; done
 
 #Need to activate mini conda environment
## Determine Germline CONTIG PLOIDY, can use an input list or alternatively add every filename in individually 
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar DetermineGermlineContigPloidy --input 5270-S1_S1.counts.hdf5  
--contig-ploidy-priors ploidy_priors.tsv --output ../Little_Campbell_Stream_CNV --output-prefix normal_cohort

## GermlineCNVCaller, again you can use an interval list
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar GermlineCNVCaller --run-mode COHORT -L filtered_intervals.interval_list 
--interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ../Little_Campbell_Stream_CNV/normal_cohort-calls 
--input 5270-S1_S1.counts.hdf5 --output ../Little_Campbell_Stream_CNV --output-prefix normal_cohort_run

## PostProcessGermlineCNVCalls
../../programs/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar PostprocessGermlineCNVCalls --contig-ploidy-calls normal_cohort-calls 
--calls-shard-path normal_cohort_run-calls --model-shard-path normal_cohort_run-model 
--output-genotyped-intervals sample_x_genotyped_intervals.vcf --output-genotyped-segments sample_x_genotyped_segments.vcf

