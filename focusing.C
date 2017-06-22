{

#define CHANNEL CH4
#define DEBUG 1
  gROOT->ProcessLine("gErrorIgnoreLevel=2001");
  char* file = (char *)  "Run6_HPK80D_KU_RED.rtct";
	PSTCT meas(file,0,2);
  meas.PrintInfo();
  enum{
    CH1=0,
    CH2=1,
    CH3=2,
    CH4=3
  } channel;

  int i=0, j=0;
  int x0=0, z0=0;
  int *dx;
  dx = (int *) malloc(sizeof(int)*meas.Nz);
  bool foundStart=false, foundEnd=false;
  double startSignal = -1*0.1*meas.GetHA(CHANNEL, meas.Nx*3/4, 0, meas.Nz/2, 0, 0)->GetMinimum(); 
  double fullSignal = -1*0.9*meas.GetHA(CHANNEL, meas.Nx*3/4, 0, meas.Nz/2, 0, 0)->GetMinimum(); 

  for(i=0;i<meas.Nz;i++){
    foundStart = false;
    foundEnd = false;
    for(j=0;j<meas.Nx;j++){
      TH1F * t1;
      t1 = meas.GetHA(CHANNEL, j, 0, i, 0, 0);
      t1->GetXaxis()->SetRange(300, 400);
      if(!foundStart){
        foundStart = startSignal<=-1*t1->GetMinimum();
        x0 = j;
      }
      if(!foundEnd){
        foundEnd = fullSignal<=-1*t1->GetMinimum();
        dx[i] = j-x0;
      }
    }
  }

  int min = 1000, indx = 0;
  for(i=0;i<meas.Nz;i++){
    if(dx[i]<min){
      min = dx[i];
      indx = i;
    }
  }
  meas.GetHA(CHANNEL, 200, 0, 20, 0, 0)->Draw("COLZ");
  cout << "Focal point at z="<<meas.z0+meas.dz*indx<<" with spot size="<<min<<endl;
}