#!/bin/bash


mkdir -p /lustre/cms/store/user/reham/MonoHiggs_2017/MonoHiggs_SR/94X/jobs4e_25ns
mkdir -p /lustre/cms/store/user/reham/MonoHiggs_2017/MonoHiggs_SR/94X/histos4e_25ns

echo "Running HtoZZto4Leptons Analysis with executables RunRunReference4e_bkg_bkg"
source /cvmfs/cms.cern.ch/cmsset_default.sh

export LD_LIBRARY_PATH=/lustre/home/reham/TEST_ELE/CMSSW_9_4_10/biglib/slc6_amd64_gcc630:/lustre/home/reham/TEST_ELE/CMSSW_9_4_10/lib/slc6_amd64_gcc630:/lustre/home/reham/TEST_ELE/CMSSW_9_4_10/external/slc6_amd64_gcc630/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_10/biglib/slc6_amd64_gcc630:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_10/lib/slc6_amd64_gcc630:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_10/external/slc6_amd64_gcc630/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/llvm/4.0.1/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/lib:$LD_LIBRARY_PATH
export PATH=/cvmfs/cms.cern.ch/share/overrides/bin:/lustre/home/reham/TEST_ELE/CMSSW_9_4_10/bin/slc6_amd64_gcc630:/lustre/home/reham/TEST_ELE/CMSSW_9_4_10/external/slc6_amd64_gcc630/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_10/bin/slc6_amd64_gcc630:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_10/external/slc6_amd64_gcc630/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/llvm/4.0.1/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/bin:/cvmfs/cms.cern.ch/common:/cvmfs/cms.cern.ch/bin:/usr/lib64/qt-3.3/bin:/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lpp/mmfs/bin:/opt/puppetlabs/bin:$PATH


export CMSSW_BASE=/lustre/home/reham/TEST_ELE/CMSSW_9_4_10
melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc493/
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR`;
    cd ${workdir};
else 
    workdir=`echo $PWD`;
    cd ${workdir};
fi

savedir=`echo /lustre/cms/store/user/reham/MonoHiggs_2017/MonoHiggs_SR/94X/histos4e_25ns`

echo "Working dir is $workdir"
#echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

echo "Compiling the macros"
bash compilereference.sh 4e


./RunReference4e_bkg ./sig_input_h150.txt 1 ./bkg_input_1.txt 1 ./data_input.txt 1 BARI NO Fall17 >& ${workdir}/RunReference4e_bkg_bkg_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_4e.log
cp -f ${workdir}/RunReference4e_bkg_bkg_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_4e.log /lustre/cms/store/user/reham/MonoHiggs_2017/MonoHiggs_SR/94X/jobs4e_25ns/RunReference4e_bkg_bkg_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_4e.log

mv ${workdir}/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root    ${savedir}/.
mv ${workdir}/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_bnn.txt ${savedir}/.
mv ${workdir}/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_bnn.root ${savedir}/.
mv ${workdir}/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_txt.txt ${savedir}/.
mv ${workdir}/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_txt_vbf.txt ${savedir}/.

if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
 rm -f $_CONDOR_SCRATCH_DIR/*
fi
