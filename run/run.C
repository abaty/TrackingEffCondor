#include "../derivation/plot_efficiency_cent.C"
#include "../ntupler/track_ntupler_cent.C"

//i is the iteration being run by the condor job
void run(int i = 0){

  float ptmin[14]=   {0.5,0.5,0.5,0.5,0.5,1 ,1 ,1 ,1 ,1  ,3 ,3 ,3  ,8};
  float ptmax[14]=   {1  ,1  ,1  ,1  ,1  ,3 ,3 ,3 ,3 ,3  ,8 ,8 ,8  ,300};
  float centmin[14]= {0  ,10 ,20 ,30 ,50 ,0 ,10,20,30,50 ,0 ,10,20 ,0};
  float centmax[14]= {10 ,20 ,30 ,50 ,100,10,20,30,50,100,10,20,100,100}; 
  
  int nevents =     200;
  int ncent_step=   4;
  int naccept_step= 4;
  int npt_step=     4;
  int nrmin_step=   3;

  int icent_step=0;
  int iaccept_step=0;
  int ipt_step=0;
  int irmin_step=0;
  int istep=0;
  		
  while(icent_step<=ncent_step && iaccept_step<=naccept_step && ipt_step<=npt_step && irmin_step<=nrmin_step){
    std::cout << "cent:"<< icent_step << "  accept: " << iaccept_step << "  pt: " << ipt_step << "  rmin: " << irmin_step << std::endl;
    track_ntupler_cent(icent_step,iaccept_step,ipt_step,irmin_step,ptmin[i],ptmax[i],centmin[i],centmax[i],nevents);
    plot_efficiency_cent(icent_step,iaccept_step,ipt_step,irmin_step,ptmin[i],ptmax[i],centmin[i],centmax[i],nevents);
    if(istep%4==0) icent_step++;
    if(istep%4==1) iaccept_step++;
    if(istep%4==2) ipt_step++;
    if(istep%4==3) irmin_step++;
    istep++;
  }
  if((istep-1)%4==0) icent_step--;
  if((istep-1)%4==1) iaccept_step--;
  if((istep-1)%4==2) ipt_step--;
  if((istep-1)%4==3) irmin_step--;
  plot_efficiency_cent(icent_step,iaccept_step,ipt_step,irmin_step,ptmin[i],ptmax[i],centmin[i],centmax[i],nevents,1);
}

int main(int argc, char *argv[]){
  if(argc != 2){
    std::cout << "Usage: runcorr <condor_iter>" << std::endl;
    return 1;
  }
  run(std::atoi(argv[1]));
  return 0;
}


