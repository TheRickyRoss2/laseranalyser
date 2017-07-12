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


  // User defined path to waveform file
	#define FILEPATH "Run28.rtct"

  gROOT->ProcessLine("gErrorIgnoreLevel=2001;");
  char * file = (char * ) FILEPATH;
  PSTCT *meas;

  meas = new PSTCT(file, 0,2);

  // Choose the scanning axis with higher resolution
  int SCANAXIS;

  if(meas->Nx>meas->Ny){
    SCANAXIS=0;
  }else{
    SCANAXIS=1;
  }

  // Dynamically select the channel with the highest resolution
  int CHANNEL;
  TH2F * waves;
  double refMax=0, refMin=0;

  for(int i=0;i<4;i++){

    if(meas->WFOnOff[i]){

      waves = meas->Draw(i, 0, 2, 0, 0, 0, 0, 2000);
      double wavemin = waves->GetMinimum();
      waves = meas->Draw(i, 0, 1, 0, 0, 0, 0, 2000);
      double wavemax = waves->GetMaximum();

      if(wavemin<refMax || wavemax>refMax){

        CHANNEL = i;
        refMin = wavemin;
        refMax = wavemax;
      }
    }
  }

  enum {
    CH1 = 0,
    CH2 = 1,
    CH3 = 2,
    CH4 = 3
  }channel;

  enum {
    X=0,
    Y=1
  }axis;

  // Set up input data buffers
  int i = 0, j = 0;
  int axis0 = 0, axis1=0, z0 = 0;
  int * dAxis;
  int start=0, end=0, binSize=0;
  dAxis = new int[meas->Nz];

  bool foundStart = false, foundEnd = false, metalToPad=true;

  // Judiciously choose a waveform which the laser hits entirely as a reference
  // We assume that our focusing axis is Z
  // In this section we also determine the direction of the pulse
  // In this section we also determine the time interval of the pulse for charge collection
  // Additionally we determine whether direction the scan is going
  double startSignal = 0, fullSignal = 0;
  bool positivePulse = true;

  switch (SCANAXIS) {
    TH1F *t1;

  case X:
    // Choose a reference signal around 3/4 of the way down
    // Standard metal->pad scan with positive pulse shape
    startSignal =0.1 * meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz / 2, 0, 0)->GetMaximum();
    fullSignal = 0.9 * meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz / 2, 0, 0)->GetMaximum();

    // If a low signal is found try looking for a negative pulse
    // metal->pad scan with negative pulse shape
    if(fullSignal<50){
      startSignal = -0.1*meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz/2, 0, 0)->GetMinimum();
      fullSignal = -0.9*meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz/2, 0, 0)->GetMinimum();
      positivePulse = false;
    }

    // If still no signal is found then check for a signal 3/4 of the way up to the pad
    // pad->metal scan with positive pulse shape
    if(fullSignal<50){
      startSignal = 0.1*meas->GetHA(CHANNEL, meas->Nx*1/4, 0, meas->Nz/2, 0, 0)->GetMaximum();
      fullSignal = 0.9*meas->GetHA(CHANNEL, meas->Nx*1/4, 0, meas->Nz/2, 0, 0)->GetMaximum();
      metalToPad = false;
      positivePulse = true;
    }

    // Determine the reference waveform based on the type of scan determined previously
    if(metalToPad){
      t1 =  meas->GetHA(CHANNEL, meas->Nx*3/4, 0, meas->Nz / 2, 0, 0);
    }else{
      t1 =  meas->GetHA(CHANNEL, meas->Nx*1/4, 0, meas->Nz / 2, 0, 0);
    }

    // Determine the step value of the bins
    binSize =meas->NP/t1->GetXaxis()->GetXmax();

    // Cut before beginning of pulse
    if(t1->GetMinimumBin()-binSize*5<0){
      start =0;
    }else{
      start = t1->GetMinimumBin()-binSize*5;
    }

    // Cut after end of pulse
    if(t1->GetMaximumBin()+binSize*5>meas->NP){
      end = meas->NP;
    }else{
      end = t1->GetMaximumBin()+binSize*5;
    }

    break;

  case Y:
    startSignal = 0.1 * meas->GetHA(CHANNEL, 0, meas->Ny * 3 / 4, meas->Nz / 2, 0, 0)->GetMaximum();
    fullSignal = 0.9 * meas->GetHA(CHANNEL, 0, meas->Ny * 3 / 4, meas->Nz / 2, 0, 0)->GetMaximum();
    if(fullSignal<50){
      startSignal = -0.1*meas->GetHA(CHANNEL, 0, meas->Ny*3/4, meas->Nz/2, 0, 0)->GetMinimum();
      fullSignal = -0.9*meas->GetHA(CHANNEL, 0, meas->Ny*3/4, meas->Nz/2, 0, 0)->GetMinimum();
      positivePulse = false;
    }
    if(fullSignal<50){
      startSignal = 0.1*meas->GetHA(CHANNEL, 0, meas->Ny*1/4, meas->Nz/2, 0, 0)->GetMaximum();
      fullSignal = 0.9*meas->GetHA(CHANNEL, 0, meas->Ny*1/4, meas->Nz/2, 0, 0)->GetMaximum();
      metalToPad = false;
      positivePulse = true;
    }

    if(metalToPad){
      t1 =  meas->GetHA(CHANNEL, 0, meas->Ny * 3 / 4, meas->Nz / 2, 0, 0);
    }else{
      t1 =  meas->GetHA(CHANNEL, 0, meas->Ny * 1 / 4, meas->Nz / 2, 0, 0);
    }
    binSize =meas->NP/t1->GetXaxis()->GetXmax();
    if(t1->GetMinimumBin()-binSize*5<0){
      start =0;
    }else{
      start = t1->GetMinimumBin()-binSize*5;
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
      for (j = 0; j < meas->Nx; j++) {
        TH1F * t1;
        t1 = meas->GetHA(CHANNEL, j, 0, i, 0, 0);
        t1->GetXaxis()->SetRange(start, end);
        if (!foundStart) {
          if(positivePulse){
            if(metalToPad){
              foundStart = startSignal <= t1->GetMaximum();
            }else{
              foundStart = startSignal >= t1->GetMaximum();
            }
          }else{
            foundStart = startSignal <= -1*t1->GetMinimum();
          }
          axis0 = j;
        }

        if (!foundEnd) {
          if(positivePulse){
            if(metalToPad){
              foundEnd = fullSignal <= t1->GetMaximum();
            }else{
              foundEnd = fullSignal >= t1->GetMaximum();
            }
          }else{
            foundEnd = fullSignal <= -1*t1->GetMinimum();
          }
          axis1 = j;
        }

      }
      if(metalToPad){
        dAxis[i] = axis1-axis0;
      }else{
        dAxis[i] = axis0-axis1;
      }
      break;

    case Y:
      for (j = 0; j < meas->Ny; j++) {
        TH1F * t1;
        t1 = meas->GetHA(CHANNEL, 0, j, i, 0, 0);
        t1->GetXaxis()->SetRange(start, end);

        if (!foundStart) {
          if(positivePulse){
            if(metalToPad){
              foundStart = startSignal <= t1->GetMaximum();
            }else{
              foundStart = startSignal >= t1->GetMaximum();
            }
          }else{
            foundStart = startSignal <= -1*t1->GetMinimum();
          }
          axis0 = j;
        }

        if (!foundEnd) {
          if(positivePulse){
            if(metalToPad){
              foundEnd = fullSignal <= t1->GetMaximum();
            }else{
              foundEnd = fullSignal >= t1->GetMaximum();
            }
          }else{
            foundEnd = fullSignal <= -1*t1->GetMinimum();
          }
          axis1 = j;
        }

        if(metalToPad){
          dAxis[i] = axis1-axis0;
        }else{
          dAxis[i] = axis0-axis1;
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
  int arrSize=0;
  if(SCANAXIS==X){
    arrSize = meas->Nx;
  }else{
    arrSize = meas->Ny;
  }
  double * vals, *bins;
  bins = new double[arrSize];
  vals = new double[arrSize];
  for(int i = 0;i<arrSize;i++){
    TH1F *t;
    if(SCANAXIS==X){
      t = meas->GetHA(CHANNEL, i, 0, indx, 0);
      bins[i] = meas->dx*i+meas->x0;
    }else{
      t = meas->GetHA(CHANNEL, 0, i, indx, 0);
      bins[i] = meas->dy*i+meas->y0;
    }
    if(positivePulse){
      vals[i] = t->GetMaximum();
    }else{
      vals[i] = -1*t->GetMinimum();
    }
  }
  TCanvas *c1 = new TCanvas("c1", "Amplitude vs Y-Position");
  c1->SetCanvasSize(1000, 500);
  c1->SetWindowSize(1050, 550);
  c1->SetGrid();
  TGraph * gr = new TGraph(arrSize, bins, vals);
  gr->Draw("AC");
  gr->SetLineColor(2);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(2);
  gr->SetMarkerSize(0.25);
  if(SCANAXIS==X){
    gr->SetTitle("Amplitude vs X-Position");
    gr->GetXaxis()->SetTitle("X-Axis Position [um]");
  }else{
    gr->SetTitle("Amplitude vs Y-Position");
    gr->GetXaxis()->SetTitle("Y-Axis Position [um]");
  }
  gr->GetYaxis()->SetTitle("Amplitude [mV]");

  // Print result
  if(SCANAXIS==X){
    cout << "Focal point at z=" <<indx*meas->dz+meas->z0 << "um with spot size=" << (min)*meas->dx << endl;
  }else{
    cout << "Focal point at z=" <<indx*meas->dz+meas->z0 << "um with spot size=" << (min)*meas->dy << endl;
  }
  delete [] dAxis;
}
