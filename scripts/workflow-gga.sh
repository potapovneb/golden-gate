#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 8

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

### update shell $PATH to include required tools
module load smrtlink-5.0.0
module load samtools-1.3.1
module load ncbi-blast-2.3.0+

### switch to working directory
mkdir -p "$rundir"
cd "$rundir"

# ### looking up sequencing data
# echo ""
# find "$collectionPathUri"/Analysis_Results -type f -name "*.bax.h5" | sort -u > input.fofn
# echo "Task 1 completed at $(date)"

# ### convert
# echo ""
# bax2bam --fofn=input.fofn -o movie
# echo "Task 2 completed at $(date)"

# ### build ccs
# echo ""
# ccs --reportFile=subreads_ccs.csv --logFile=subreads_ccs.log --numThreads=8 --minPasses=1 movie.subreads.bam subreads_ccs.bam
# echo "Task 3 completed at $(date)"

# ### cleanup
# echo ""
# rm -f movie.scraps.bam*
# rm -f movie.subreads.bam*
# echo "Task 4 completed at $(date)"

# ### save consensus reads in fasta format
# echo ""
# "$root"/bin/bam2fq.pl -fasta subreads_ccs.bam subreads_ccs.fasta
# "$root"/bin/bam2csv.pl subreads_ccs.bam subreads_ccs.sts.csv
# echo "Task 5 completed at $(date)"

# ### create blast database
# echo ""
# mkdir db
# makeblastdb -dbtype nucl -in subreads_ccs.fasta -input_type fasta -title goldengate -out db/goldengate -logfile db/goldengate.log
# echo "Task 5 completed at $(date)"

# ### map fragments to assemblies
# echo ""
# BLASTDB="$rundir"/db
# outfmt="qseqid qlen sseqid slen qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand qcovs qcovhsp qcovus"
# blastn \
#     -query "$inserts" \
#     -db goldengate \
#     -out blast_result.tmp \
#     -outfmt "10 $outfmt" \
#     -num_alignments 1000000 \
#     -num_threads 8
# echo $outfmt | tr ' ' , > blast_result.csv
# cat blast_result.tmp >> blast_result.csv
# echo "Task 6 completed at $(date)"

### process golden gate assemb,lies
echo ""
mkdir -p "$rundir"/results
cd "$rundir"/results
"$root"/bin/assemble.pl \
    --symmetric \
    --nobreak 3,4,5 \
    --keep 3,4,5,6,7,8 \
    --pro-max 0 \
    --epi-max 0 \
    --discard \
    ../{blast_result.csv,subreads_ccs.fasta,subreads_ccs.sts.csv} "$inserts"
    echo "Task 7 completed at $(date)"

echo ""
echo "Finished on $(date)"
