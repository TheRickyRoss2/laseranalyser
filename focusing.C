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

  // User defined variables
  #define CHANNEL CH2 // Select which channel the pad is linked to
	#define FILEPATH "../../fredData.rtct"
	#define SCANAXIS Y
  #define tStart 13 // Approx start time of pulse
  #define tEnd 19 // Approx end time of pulse
  // Aux macros
  # define time(X)(402 * (X)/20) // DRS4 time->bin conversion

  // Load in waveform file
  gROOT->ProcessLine("gErrorIgnoreLevel=2001");
  char * file = (char * ) FILEPATH;
  PSTCT meas(file,0, 2);
  meas.PrintInfo();

  enum {
    X = 0,
      Y = 1
  }

  axis;

  enum {
    CH1 = 0,
      CH2 = 1,
      CH3 = 2,
      CH4 = 3
  }
  channel;

  // Set up input data buffers
  int i = 0, j = 0;
  int axis0 = 0, z0 = 0;
  int * dAxis;
  dAxis = (int * ) malloc(sizeof(int) * meas.Nz);
  bool foundStart = false, foundEnd = false;

  // Judiciously choose a waveform which the laser hits entirely as a reference
  // In this case we assume that the scan is from on the metal to off the metal
  // We also assume that our focusing axis is Z
  double startSignal = 0, fullSignal = 0;

  switch (SCANAXIS) {
  case X:
    startSignal = 0.1 * meas.GetHA(CHANNEL, meas.Nx * 3 / 4, 0, meas.Nz / 2, 0, 0)->GetMaximum();
    fullSignal = 0.9 * meas.GetHA(CHANNEL, meas.Nx * 3 / 4, 0, meas.Nz / 2, 0, 0)->GetMaximum();
    break;
  case Y:
    startSignal = 17;//0.1 * meas.GetHA(CHANNEL, 0, meas.Ny * 3 / 4, meas.Nz / 2, 0, 0)->GetMaximum();
    fullSignal = 190;//0.9 * meas.GetHA(CHANNEL, 0, meas.Ny * 3 / 4, meas.Nz / 2, 0, 0)->GetMaximum();
    cout <<"ss"<<startSignal<<"fs"<<fullSignal<<endl;
    break;
  default:
    startSignal = 15;
    fullSignal = 120;
  }

  for (i = 0; i < meas.Nz; i++) {
    cout<<"Z="<<i<<endl;
    foundStart = false;
    foundEnd = false;

    switch (SCANAXIS) {
    case X:
      for (j = 0; j < meas.Nx; j++) {
        TH1F * t1;
        t1 = meas.GetHA(CHANNEL, j, 0, i, 0, 0);
        t1->GetXaxis()->SetRange(time(tStart), time(tEnd));
        if (!foundStart) {
          foundStart = startSignal >= t1->GetMaximum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal >= t1->GetMaximum();
          dAxis[i] = j - axis0;
        }
      }
      break;

    case Y:
      for (j = 0; j < meas.Ny; j++) {
        TH1F * t1;
        t1 = meas.GetHA(CHANNEL, 0, j, i, 0, 0);
        t1->GetXaxis()->SetRange(time(tStart), time(tEnd));
        cout <<"Max of current graph:"<<t1->GetMaximum()<<endl;
        if (!foundStart) {
          foundStart = startSignal <= t1->GetMaximum();
          axis0 = j;
          if(foundStart)cout<<"Start "<<j<<endl;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= t1->GetMaximum();
          dAxis[i] = j - axis0;
          if(foundStart)cout<<"End "<<dAxis[i]<<endl;
        }
      }
      break;

    default:
      for (j = 0; j < meas.Nx; j++) {
        TH1F * t1;
        t1 = meas.GetHA(CHANNEL, j, 0, i, 0, 0);
        t1->GetXaxis()->SetRange(time(tStart), time(tEnd));
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
    cout <<endl;
  }

  // Iterate through each Z-coordinates and find the shortest distance from 10% to 90% signal
  int min = 1000, indx = 0;
  for (i = 0; i < meas.Nz; i++) {
    cout << dAxis[i]<<endl;
    if (dAxis[i] < min) {
      min = dAxis[i];
      indx = i;
    }
  }

  // Print result
  cout << "Focal point at z=" <<indx*meas.dz+meas.z0 << " with spot size=" << min*meas.dy << endl;
  free(dAxis);
}
