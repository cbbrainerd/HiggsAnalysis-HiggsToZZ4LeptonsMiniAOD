#!/bin/bash


if [ "$1" == "" ] || [ "$2" == "" ] || [ "$4" == "" ] || [ "$4" == "" ]; then
    echo "Please provide arguments to the script: site configuration, data type and MC type"
    echo "Usage bash loopcheck_bkg_<finalstate>.sh <arg1> <arg2> <arg3> <arg4>"
    exit      
fi



echo "$1 configuration";
echo "$2 data"
echo "$3 simulation"
echo "$4 site"

SCERN="CERN";
SFNAL="FNAL";
SDESY="DESY";
SBARI="BARI";

###### Signal

n=0;
m=0;

mkdir -p jobs$4

echo "Reading sig_input_$3_AN.txt file"

cp sig_input_$3_AN.txt sig_input.txt
nlines=`wc -l sig_input.txt | awk '{print $1}'`;

mkdir -p SigCards$4$3

while [ $n -lt ${nlines} ]; do
    (( n = n + 1 ))
    (( m = ${nlines} - n ))
    echo $n $m
    rm -f SigCards$4$3/sig_input_${n}.txt
    cat sig_input.txt | head -1 > SigCards$4$3/sig_input_${n}.txt
    samplename=`cat SigCards$4$3/sig_input_${n}.txt | awk '{print $1}'`
    cat sig_input.txt | tail -n $m >  sig_input_tmp.txt
    mv  sig_input_tmp.txt sig_input.txt
    echo $samplename
    rm -f jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh	

    if [ $1 = ${SCERN} ]; then
     cat submit_HZZ4LeptonsAnalysis_CERN.sh | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?HZZ4LeptonsAnalysis?RunReferenceAnalysis?g" | sed "s?jobdir?jobs$4_25ns?g" | sed "s?histodir?histos$4_25ns?g" | sed "s?output?output_${samplename}?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?SigCards$4$3/sig_input_${n}.txt?g" | sed "s?_log?_${samplename}_$4.log?g" > jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh
    elif  [ $1 = ${SFNAL} ]; then 
     mkdir -p /eos/uscms/store/user/`whoami`/80X/jobs$4_25ns
     mkdir -p /eos/uscms/store/user/`whoami`/80X/histos$4_25ns

     cat submit_HZZ4LeptonsAnalysis_FNAL.sh | sed "s?site?$1?g" | sed "s?CMSSW_SEARCH_PATH_DIR?${CMSSW_SEARCH_PATH}?g" | sed "s?CMSSW_BASE_DIR?${CMSSW_BASE}?g" | sed "s?path?$PATH?g"  | sed "s?lib:?$LD_LIBRARY_PATH:?g" | sed "s?4mu?$4?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?HZZ4LeptonsAnalysis?RunReferenceAnalysis?g" | sed "s?jobdir?jobs$4_25ns?g" | sed "s?histodir?histos$4_25ns?g" | sed "s?output?output_${samplename}?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?sig_input_${n}.txt?g" | sed "s?_log?_${samplename}_$4.log?g" > jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh
     cat condor_template.cfg  | sed "s?4mu?$4?g" | sed "s?submit_HZZ4LeptonsAnalysis_BARI?submit_ReferenceAnalysis_sig_${samplename}_$4?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?SigCards$4$3/sig_input_${n}.txt?g" | sed "s?mail?`whoami`?g" > jobs$4/condor_ReferenceAnalysis_sig_${samplename}_$4.cfg

    elif  [ $1 = ${SDESY} ]; then
     cat submit_HZZ4LeptonsAnalysis_DESY.sh | sed "s?site?$1?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?HZZ4LeptonsAnalysis?RunReferenceAnalysis?g" | sed "s?jobdir?jobs$4_25ns?g" | sed "s?histodir?histos$4_25ns?g" | sed "s?output?output_${samplename}?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?SigCards$4$3/sig_input_${n}.txt?g" | sed "s?_log?_${samplename}_$4.log?g" > jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh
    elif  [ $1 = ${SBARI} ]; then
	cat submit_HZZ4LeptonsAnalysis_BARI.sh | sed "s?site?$1?g" | sed "s?CMSSW_BASE_DIR?${CMSSW_BASE}?g" | sed "s?path?$PATH?g"  | sed "s?lib:?$LD_LIBRARY_PATH:?g" | sed "s?4mu?$4?g" | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?HZZ4LeptonsAnalysis?RunReferenceAnalysis?g" | sed "s?jobdir?jobs$4_25ns?g" | sed "s?histodir?histos$4_25ns?g" | sed "s?output?output_${samplename}?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?sig_input_${n}.txt?g" | sed "s?_log?_${samplename}_$4.log?g" > jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh
	cat condor_template.cfg  | sed "s?4mu?$4?g" | sed "s?submit_HZZ4LeptonsAnalysis_BARI?submit_ReferenceAnalysis_sig_${samplename}_$4?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?SigCards$4$3/sig_input_${n}.txt?g" | sed "s?mail?`whoami`?g" > jobs$4/condor_ReferenceAnalysis_sig_${samplename}_$4.cfg
    else 
     cat submit_HZZ4LeptonsAnalysis.sh | sed "s?mc?$3?g" |sed "s?year?$2?g" | sed "s?HZZ4LeptonsAnalysis?RunReferenceAnalysis?g" | sed "s?jobdir?jobs$4_25ns?g" | sed "s?histodir?histos$4_25ns?g" | sed "s?output?output_${samplename}?g" | sed "s?RunReferenceAnalysis?RunReference$4_signal?g" | sed "s?sig_input_h150.txt?SigCards$4$3/sig_input_${n}.txt?g" | sed "s?_log?_${samplename}_$4.log?g" > jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh
    fi

    chmod u+xr  jobs$4/submit_ReferenceAnalysis_sig_${samplename}_$4.sh

    cd jobs$4
    
    if [ $1 = ${SCERN} ]; then
       echo "Submitting jobs via LSF at CERN"
       bsub -q 8nh  submit_ReferenceAnalysis_sig_${samplename}_$4.sh
    elif  [ $1 = ${SFNAL} ]; then
       echo "Submitting jobs via CONDOR at FNAL"
      condor_submit  condor_ReferenceAnalysis_sig_${samplename}_$4.cfg
    elif  [ $1 = ${SDESY} ]; then
       echo "Submitting jobs via SGE"
       qsub submit_ReferenceAnalysis_sig_${samplename}_$4.sh   
    elif  [ $1 = ${SBARI} ]; then
       echo "Submitting jobs via CONDOR at BARI"
       condor_submit -name ettore  condor_ReferenceAnalysis_sig_${samplename}_$4.cfg
    else
       echo "Submitting jobs via PBS"
       qsub -q local submit_ReferenceAnalysis_sig_${samplename}_$4.sh
    fi 
    cd ..
done 

