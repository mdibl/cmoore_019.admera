### 3-REAP-paired reads


#####################################################
### USER-DEFINED CONFIGURATION — EDIT THESE! ###
#####################################################

script_path=/compbio/scratch/cnobrega/3_REAP/scripts
fastq_path=/compbio/data/ClaireMoore/cmoore_019.admera
work_path=/compbio/scratch/cnobrega/3_REAP/cmoore_019.admera/output
sample_file=/compbio/scratch/cnobrega/3_REAP/cmoore_019.admera/sample.txt

nrThreads=12 # cores

star_index_path=/compbio/scratch/cnobrega/3_REAP/ref/star_index_hg38 # for Step 3, STAR index. # If you have no star index, please build one first.
refPAS_file=/compbio/scratch/cnobrega/3_REAP/ref/human.PAS.hg38.txt # for Step 6, PAS annotation file from PolyA_DBv3.2
chromsizes=/compbio/scratch/cnobrega/3_REAP/ref/hg38.chrom.sizes # for Step 7, Chromosome sizes file

echo "Starting 3-REAP pipeline..."
#=============================================================
### step 1. UMI extraction
#=============================================================

#echo "Step 1: UMI extraction"
#
#cd $fastq_path
#while IFS=$'\t' read -r sample fastq_1 fastq_2; do
#	outpath=${work_path}/s1_fastq_umi/${sample}
#	mkdir -p ${outpath}
#	umi_tools extract --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN \
#	-I $fastq_1 \
#	--read2-in=$fastq_2 \
#	-L ${outpath}/extraction_log \
#	-S ${outpath}/${sample}_R1.fastq.gz \
#	--read2-out=${outpath}/${sample}_R2.fastq.gz
#done < $sample_file

#=============================================================
### step 2. Trimming adapters and polyA/polyT sequences
#=============================================================
echo "Step 2: Adapter and polyA/polyT trimming"

cd $fastq_path
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
	echo "Processing $sample read 1..."
	echo "Path: $fastq_1"
	outpath=${work_path}/s2_fastq_cutadapt/${sample}
	mkdir -p ${outpath}
	cutadapt --quiet -j 24 -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 $fastq_1 | \
	cutadapt --quiet -j 24 -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
	cutadapt --quiet -j 24 -m 20 -O 3 -a "r1polyA=A{18}" - | \
	cutadapt --quiet -j 24 -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" -o ${outpath}/${sample}_R1.fastq.gz --discard-trimmed -
done < $sample_file

while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    echo "Processing $sample read 2..."
    echo "Path: $fastq_2"
    outpath=${work_path}/s2_fastq_cutadapt/${sample}
    mkdir -p ${outpath}
    cutadapt --quiet -j 24 -m 18 -O 3 --nextseq-trim=10 -a "rd2Adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=3;max_error_rate=0.100000" $fastq_2 | \
    cutadapt --quiet -j 24 -m 18 -O 20 -g "rd2Adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=20" -o ${outpath}/${sample}_R2_trimAdapt.fastq.gz --discard-trimmed -
    ### Reverse reads contain 5'Ts, remove remaining 5'Ts before the alignment
    python ${script_path}/trim_N-5T_fq.py --rawfastq ${outpath}/${sample}_R2_trimAdapt.fastq.gz --out_dir ${outpath} --random_NT_len 0
done < $sample_file

#=============================================================
### step 3. Genome alignment using paired reads
#=============================================================
# Since the fastq files were preprocessed separately, they need to be repaired to re-establish read pairing before alignment. repair.sh is a tool of bbmap.
# If you have no star index, please build one first.
echo "Step 3: Repair read pairs and align with STAR"
# pair reads repair
inpath=${work_path}/s2_fastq_cutadapt
outpath=${work_path}/s2_fastq_cutadapt_pair; mkdir -p ${outpath}
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
repair.sh overwrite=t in1=${inpath}/${sample}/${sample}_R1.fastq.gz in2=${inpath}/${sample}/${sample}_R2_trimAdapt.5Ttrimmed.fastq.gz \
           out1=${outpath}/${sample}/${sample}_R1.fastq.gz out2=${outpath}/${sample}/${sample}_R2.fastq.gz outs=${outpath}/${sample}/${sample}_unpaired.fastq.gz \
           repair
done < $sample_file

echo "Aligning with STAR..."
inpath=${work_path}/s2_fastq_cutadapt_pair
cd ${work_path}
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
	outpath=${work_path}/s3_star_align/${sample}; mkdir -p ${outpath}
	STAR --runThreadN ${nrThreads} --readFilesCommand zcat --genomeDir ${star_index_path} \
	--readFilesIn ${inpath}/${sample}/${sample}_R1.fastq.gz ${inpath}/${sample}/${sample}_R2.fastq.gz --outReadsUnmapped None --outFilterType BySJout --limitBAMsortRAM 4000000000 \
	--outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/ \
	--outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2
	
	# delete STAR temp dir
	rm -rf ${outpath}/_STARtmp
	
	# delete intermediate FASTQ files after STAR alignment
	rm -f ${work_path}/s2_fastq_cutadapt/${sample}/*.fastq.gz
	rm -f ${work_path}/s2_fastq_cutadapt_pair/${sample}/*unpaired.fastq.gz
done < $sample_file

## Statistical alignment rate
echo "Compiling STAR alignment statistics..."
cd ${work_path}/s3_star_align
output_file="stats_STAR_summary.csv"
echo "Sample,Total Reads,Uniquely Mapped Reads (%),Mapped to Multiple Loci (%),Unmapped: Too Many Mismatches (%),Unmapped: Too Short (%),Unmapped: Other (%)" > $output_file
for logfile in */*Log.final.out; do
    #sample=$(basename $(dirname $logfile))
    sample=${logfile%%Log.final.out}
    total_reads=$(grep "Number of input reads" $logfile | awk '{print $NF}')
    uniquely_mapped=$(grep "Uniquely mapped reads %" $logfile | awk '{print $NF}')
    multiple_mapped=$(grep "% of reads mapped to multiple loci" $logfile | awk '{print $NF}')
    unmapped_mismatches=$(grep "% of reads unmapped: too many mismatches" $logfile | awk '{print $NF}')
    unmapped_short=$(grep "% of reads unmapped: too short" $logfile | awk '{print $NF}')
    unmapped_other=$(grep "% of reads unmapped: other" $logfile | awk '{print $NF}')
    echo "$sample,$total_reads,$uniquely_mapped,$multiple_mapped,$unmapped_mismatches,$unmapped_short,$unmapped_other" >> $output_file
done

#=============================================================
### step 4. Removing PCR duplication
#=============================================================
echo "Step 4: Deduplication"
inpath=${work_path}/s3_star_align

while IFS=$'\t' read -r sample fastq_1 fastq_2; do
	outpath=${work_path}/s4_bam_dedup/${sample}; mkdir -p ${outpath}
	
	# ensure BAM index exists
	if [ ! -f "$inpath/$sample/Aligned.sortedByCoord.out.bam.bai" ] || [ "$inpath/$sample/Aligned.sortedByCoord.out.bam" -nt "$inpath/$sample/Aligned.sortedByCoord.out.bam.bai" ]; then
    	samtools index "${inpath}/${sample}/Aligned.sortedByCoord.out.bam"
	fi
	
	# create a symbolic link since not using UMIs
	ln -s ${inpath}/${sample}/Aligned.sortedByCoord.out.bam ${outpath}/${sample}.dedup.bam
	
	# umi_tools dedup --paired -I ${inpath}/${sample}/Aligned.sortedByCoord.out.bam -S ${outpath}/${sample}.dedup.bam \
	# --multimapping-detection-method=NH --output-stats=${outpath}/deduplicated.txt --log=${outpath}/deduplication.log
done < $sample_file

#=============================================================
### step 5. Bam file to Bed
#=============================================================
# Extract R2 reads that are successfully paired and aligned to the reference genome from the BAM file, then convert them to BED format.
echo "Step 5: Convert BAM to BED"
outdir=${work_path}/s5_R2_bed; mkdir -p $outdir
cd $outdir
bam_path=${work_path}/s4_bam_dedup
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    samtools view -f 128 -F 4 -b ${bam_path}/${sample}/${sample}.dedup.bam > ${bam_path}/${sample}/${sample}.dedup_R2.bam
    samtools index ${bam_path}/${sample}/${sample}.dedup_R2.bam
    bedtools bamtobed -cigar -i ${bam_path}/${sample}/${sample}.dedup_R2.bam > ${sample}.dedup_R2.bed
    sort -k 1,1 ${sample}.dedup_R2.bed > ${sample}.dedup_R2.sorted.bed
    gzip ${sample}.dedup_R2.sorted.bed
done < $sample_file

wc -l *.sorted.bed >> stats_deduped_reads.txt

# samtools view: -f 128, select read2 (second in pair); -F 4, Exclude mates that are unmapped
#=============================================================
### step 6. Defining LAPs (last aligned positions) and matching PASs
#=============================================================
# -misM 2: Allow up to 2 nt soft clipping at the ends of aligned reads.
# -dist 24: Assign LAPs to PASs in PolyA_DB within a ±24-nt window.
echo "Step 6: LAP/PAS analysis"
indir=${work_path}/s5_R2_bed
outdir=${work_path}/s6_LAP; mkdir -p $outdir

while IFS=$'\t' read -r sample fastq_1 fastq_2; do
	Rscript ${script_path}/genome_LAPandPAS_define.R -bedLAP <(zcat $indir/${sample}.dedup_R2.sorted.bed.gz) -out $outdir/${sample} -refPAS $refPAS_file -misM 2 -dist 24
done < $sample_file

# combine the PAS quantification table
mkdir -p ${work_path}/s6_LAP/PAS_quant; cd ${work_path}/s6_LAP/PAS_quant
mv ${work_path}/s6_LAP/*cluster.all.reads.csv ./
Rscript ${script_path}/combine_all_sample_PAS_count_tables.R -csv ./ -out ./cluster.all.reads.csv --pattern "_cluster\\.all\\.reads\\.csv$"
Rscript ${script_path}/combine_all_sample_PAS_count_tables.R -csv ${outdir} -out ./stats_PAS.csv --pattern ".*_stats\\.csv$"

#=============================================================
### step 7. Generating bigwig files for PAS usage visualization
#=============================================================
echo "Step 7: Generate bigWig tracks"
indir=${work_path}/s6_LAP

for type in PASS_bw_LAP24 PASS_bw_LAP24_positon PASS_bw_LAP24_PAS mapped_read
do
outdir=${work_path}/s7_bigwig/${type}; mkdir -p $outdir; cd $outdir
echo "Working on ${type}..."

while IFS=$'\t' read -r sample fastq_1 fastq_2; do
  echo "Working on $sample..."
  ## Count total read number
  totalReadNum=`wc -l $indir/${sample}_${type}.bed | sed s/[[:blank:]].*//`
  echo "for file ${sample}_PASS_bw.bed, TotalReadNum=$totalReadNum"

  zcat $indir/${sample}_${type}.bed | sort -k 1,1 > $indir/${sample}_${type}.bed.sorted

  ## it is strand-specific 
  ## Generate bedgraph file
  echo "Generate bedgraph files for + and - strands..."
  genomeCoverageBed -bg -split -i $indir/${sample}_${type}.bed.sorted -strand '+' -g $chromsizes > $sample.plus.bedgraph
  genomeCoverageBed -bg -split -i $indir/${sample}_${type}.bed.sorted -strand '-' -g $chromsizes > $sample.minus.bedgraph
  
  gzip $indir/${sample}_${type}.bed.sorted
  
  ## Normalize bedgraph counts
  echo "Normalize bedgraph counts..."
  ${script_path}/norm_bedgraph.pl -t $totalReadNum -i "$sample.plus.bedgraph"
  ${script_path}/norm_bedgraph.pl -t $totalReadNum -i "$sample.minus.bedgraph"

    ## give minus strand negative value
  awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2, $3, -$4}' $sample.minus.bedgraph.normalized > $sample.minus.bedgraph.normalized1 

  ## Convert to bigwig file
  echo "Convert to bigwig file..."
  bedGraphToBigWig $sample.plus.bedgraph.normalized  $chromsizes plus.$sample.$type.bw
  bedGraphToBigWig $sample.minus.bedgraph.normalized  $chromsizes minus.$sample.$type.bw
  bedGraphToBigWig $sample.minus.bedgraph.normalized1  $chromsizes mi.$sample.$type.bw

done < $sample_file
chmod 775 *.bw #change file permission
rm -f *.bedgraph
rm -f *.bedgraph.normalized
rm -f *.bedgraph.normalized1
done

echo "Pipeline complete!"
