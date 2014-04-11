if [ $# -ne 1 ]
then
  echo "Usage: ./run.sh <condor_iteration>"
  exit 1
fi

echo | awk -v i=$1 '{print "./run.exe "i" "}' 
echo | awk -v i=$1 '{print "./run.exe "i" "}' | bash
echo | awk -v tag=$4 -v user=$USER '{print "mkdir -p /net/hisrv0001/home/"user"/scratch_proxy/HIRun2013/merged"}'
echo | awk -v tag=$4 -v user=$USER '{print "mkdir -p /net/hisrv0001/home/"user"/scratch_proxy/HIRun2013/merged"}' | bash
echo | awk -v tag=$4 -v user=$USER '{print "mkdir -p /net/hisrv0001/home/"user"/scratch_proxy/HIRun2013/umerged"}'
echo | awk -v tag=$4 -v user=$USER '{print "mkdir -p /net/hisrv0001/home/"user"/scratch_proxy/HIRun2013/unmerged"}' | bash
echo | awk -v tag=$4 -v user=$USER '{print "mv eff_pt*.root /net/hisrv0001/home/"user"/tracking_eff_corrections/CMSSW_5_3_12_patch3/src/condor_trk_corr/final_hists/"}'
echo | awk -v tag=$4 -v user=$USER '{print "mv eff_pt*.root /net/hisrv0001/home/"user"/tracking_eff_corrections/CMSSW_5_3_12_patch3/src/condor_trk_corr/final_hists/"}' | bash
rm *
echo "job done successfully"
