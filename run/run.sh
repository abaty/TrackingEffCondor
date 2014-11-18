if [ $# -ne 1 ]
then
  echo "Usage: ./run.sh <condor_iteration>"
  exit 1
fi

echo | awk -v i=$1 '{print "./run.exe "i" "}' 
echo | awk -v i=$1 '{print "./run.exe "i" "}' | bash

echo | awk -v tag=$4 -v user=$USER '{print "mv eff_corr*.root /net/hisrv0001/home/"user"/factorized_corrections/CMSSW_5_3_12_patch3/src/condor_trk_corr/final_hists/"}' | bash

echo | awk -v tag=$4 -v user=$USER '{print "mv eff_pt*.root /net/hisrv0001/home/"user"/factorized_corrections/CMSSW_5_3_12_patch3/src/condor_trk_corr/final_hists/"}'
echo | awk -v tag=$4 -v user=$USER '{print "mv eff_pt*.root /net/hisrv0001/home/"user"/factorized_corrections/CMSSW_5_3_12_patch3/src/condor_trk_corr/final_hists/"}' | bash
rm *.root
echo "job done successfully"
