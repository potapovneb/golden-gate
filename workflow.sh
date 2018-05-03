#!/bin/bash

samples=$1
command=$2

PROJDIR=$PWD

while read line
do
    if [[ ! $line =~ "SampleID" ]] && [[ ! $line =~ "#" ]]
    then
	### sample details
	sampleId=`echo $line | cut -d, -f1`
	collectionPathUri=`echo $line | cut -d, -f5`
	
	### preview run info
	echo ""
	echo "sampleId=$sampleId"
	echo "collectionPathUri=$collectionPathUri"

	### output directory
	rundir=`printf "%s/samples/%05i" $PROJDIR $sampleId`
	mkdir -p "$rundir"

	# ### run GGA script
	# qsub \
	#     -v root="$PROJDIR",rundir="$rundir",inserts="$PROJDIR/references/inserts.fasta",collectionPathUri="$collectionPathUri" \
	#     -N "gga$sampleId" \
	#     -o "$rundir"/workflow.log \
	#     -j yes \
	#     "$PROJDIR"/scripts/workflow-gga.sh

    	### run GGA script
	qsub \
	    -v root="$PROJDIR",rundir="$rundir",inserts="$PROJDIR/references/inserts.fasta",collectionPathUri="$collectionPathUri" \
	    -N "gga$sampleId" \
	    -o /dev/null \
	    -j yes \
	    "$PROJDIR"/scripts/workflow-gga.sh
    fi
done < $samples
