{
  /****************************************************************************
   * File: focusing.C
   * Version: 1.0
   * Author: Ric Rodriguez
   * Algorithm Credit: Federico Siviero
   * Function: Calculates the focal point optical coordinate and beam spot size
   * Example usage
   * ~bash$: root
   * root[0]: gSystem->Load("TCTAnalyse.sl");
   * root[1]: .x focusing.C;
   ****************************************************************************/


	#define FILEPATH "data.rtct"

  gROOT->ProcessLine("gErrorIgnoreLevel=2001");
  char * file = (char * ) FILEPATH;
  PSTCT *meas;
  meas = new PSTCT(file, 0,2);
  meas->PrintInfo();
  //meas->CorrectBaseLine();
  int SCANAXIS;
  if(meas->Nx>1){
    SCANAXIS=0;
  }else{
    SCANAXIS=1;
  }
  int CHANNEL;
  for(int i=0;i<4;i++){
    if(meas->WFOnOff[i]){
      CHANNEL = i;
    }
  }

  enum {
    CH1 = 0,
    CH2 = 1,
    CH3 = 2,
    CH4 = 3
  }
  channel;
  enum {
    X=0,
    Y=1
  }axis;

  // Set up input data buffers
  int i = 0, j = 0;
  int axis0 = 0, z0 = 0;
  int * dAxis;
  int start=0, end=0, binSize=0;
  dAxis = new int[meas->Nz];
  std::vector<int> distance;

  bool foundStart = false, foundEnd = false;

  // Judiciously choose a waveform which the laser hits entirely as a reference
  // In this case we assume that the scan is from on the metal to off the metal
  // We also assume that our focusing axis is Z
  double startSignal = 0, fullSignal = 0;

  switch (SCANAXIS) {
    TH1F *t1;
  case X:
    startSignal =0.1 * meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz / 2, 0, 0)->GetMaximum();
    fullSignal = 0.9 * meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz / 2, 0, 0)->GetMaximum();
    t1 =  meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz / 2, 0, 0);
    binSize =meas->NP/t1->GetXaxis()->GetXmax();
    if(t1->GetMaximumBin()-binSize*5<0){
      start =0;
    }else{
      start = t1->GetMaximumBin()-binSize*5;
    }
    if(t1->GetMaximumBin()+binSize*5>meas->NP){
      end = meas->NP;
    }else{
      end = t1->GetMaximumBin()+binSize*5;
    }

    break;
  case Y:
    startSignal = 0.1 * meas->GetHA(CHANNEL, 0, meas->Ny * 3 / 4, meas->Nz / 2, 0, 0)->GetMaximum();
    fullSignal = 0.9 * meas->GetHA(CHANNEL, 0, meas->Ny * 3 / 4, meas->Nz / 2, 0, 0)->GetMaximum();
    t1 =  meas->GetHA(CHANNEL, 0, meas->Ny * 3 / 4, meas->Nz / 2, 0, 0);
    binSize =meas->NP/t1->GetXaxis()->GetXmax();
    if(t1->GetMaximumBin()-binSize*5<0){
      start =0;
    }else{
      start = t1->GetMaximumBin()-binSize*5;
    }
    if(t1->GetMaximumBin()+binSize*5>meas->NP){
      end = meas->NP;
    }else{
      end = t1->GetMaximumBin()+binSize*5;
    }

    break;
  default:
    startSignal = 15;
    fullSignal = 120;
  }

  for (i = 0; i < meas->Nz; i++) {
    foundStart = false;
    foundEnd = false;

    switch (SCANAXIS) {
    case X:
      for (j = 0; j < meas->Ny; j++) {
        TH1F * t1;
        t1 = meas->GetHA(CHANNEL, j, 0, i, 0, 0);
        t1->GetXaxis()->SetRange(start, end);
        if (!foundStart) {
          foundStart = startSignal <= t1->GetMaximum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= t1->GetMaximum();
          dAxis[i] = j - axis0;
          distance.push_back(j-axis0);
        }
      }
      break;

    case Y:
      for (j = 0; j < meas->Ny; j++) {
        TH1F * t1;
        t1 = meas->GetHA(CHANNEL, 0, j, i, 0, 0);

        //t1->GetXaxis()->SetRange(start, end);
        if (!foundStart) {
          foundStart = startSignal <= t1->GetMaximum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= t1->GetMaximum();
          dAxis[i] = j - axis0;
          distance.push_back(j-axis0);
        }

      }
      break;

    default:
      for (j = 0; j < meas->Ny; j++) {
        TH1F * t1;
        t1 = meas->GetHA(CHANNEL, 0, j, i, 0, 0);
        t1->GetXaxis()->SetRange(start, end);
        if (!foundStart) {
          foundStart = startSignal <= t1->GetMaximum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= t1->GetMaximum();
          dAxis[i] = j - axis0;
          distance.push_back(j-axis0);
        }
      }
      break;
    }
  }

  // Iterate through each Z-coordinates and find the shortest distance from 10% to 90% signal
  int min = 1000, indx = 0;
  for (i = 0; i < meas->Nz; i++) {
    if (dAxis[i] < min) {
      min = dAxis[i];
      indx = i;
    }
  }
  double * vals, *bins;
  bins = new double[meas->Ny];
  vals = new double[meas->Ny+1];
  for(int i = 0;i<meas->Ny;i++){
    TH1F *t;
    t = meas->GetHA(CHANNEL, 0, i, indx, 0);
    vals[i] = t->GetMaximum();
    bins[i] = meas->dy*i+meas->y0;
  }
  TCanvas *c1 = new TCanvas("c1", "Amplitude vs Y-Position");
  c1->SetCanvasSize(1000, 500);
  c1->SetWindowSize(1050, 550);
  c1->SetGrid();
  TGraph * gr = new TGraph(meas->Ny, bins, vals);
  gr->Draw("AC");
  gr->SetLineColor(4);
  gr->SetTitle("Amplitude vs Y-Position");
  gr->GetYaxis()->SetTitle("Amplitude [mV]");
  gr->GetXaxis()->SetTitle("Y-Axis Position [um]");

  // Print result
  cout << "Focal point at z=" <<indx*meas->dz+meas->z0 << "um with spot size=" << (min+1)*meas->dy << endl;
  delete [] dAxis;
}
