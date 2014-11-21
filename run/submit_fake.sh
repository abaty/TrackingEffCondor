if [ $# -ne 0 ]
then
  echo "Usage: ./psort.sh <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax>"
  exit 1
fi

now="trkfake_corr_$(date +"%Y_%m_%d__%H_%M_%S")"
njobs=29

mkdir $now
cp centrality_weights.root $now
cp run_fake.sh $now
cp FFJEC_correction_PF_akVs3Calo_pt2.root $now
cat run_fake.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@user_flag@$USER@g" |  sed "s@arglist@ @g" | sed "s@transfer_filelist@run_fake.exe@g" | sed "s@njobs@$njobs@g" > $now/run_fake.condor

NAME="run_fake.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"
cp run_fake.exe $now
rm run_fake.exe
echo
cat $now/run_fake.condor
echo 
echo condor_submit $now/run_fake.condor
echo
# condor_submit $now/run.condor

