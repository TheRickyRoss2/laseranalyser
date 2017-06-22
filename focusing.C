{
  /****************************************************************************
   * File: focusing.C
   * Version: 1.0
   * Author: Ric Rodriguez
   * Function: Calculates the focal point optical coordinate and beam spot size
   * Example usage
   * ~bash$: root
   * root[0]: gSystem->Load("TCTAnalyse.sl");
   * root[1]: .x focusing.C;
   ****************************************************************************/

  // User defined variables
  #define CHANNEL CH4 // Select which channel the pad is linked to
  #define FILEPATH "Run6_HPK80D_KU.rtct" // Path to PSTCT waveform file
  #define SCANAXIS X // Scanning axis for waveforms
  #define tStart 60 // Approx start time of pulse
  #define tEnd 80 // Approx end time of pulse

  // Aux macros
  # define time(X)(1024 * (X) / 200) // DRS4 time->bin conversion

  // Load in waveform file
  gROOT - > ProcessLine("gErrorIgnoreLevel=2001");
  char * file = (char * ) FILEPATH;
  PSTCT meas(file, 0, 2);
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
    startSignal = -1 * 0.1 * meas.GetHA(CHANNEL, meas.Nx * 3 / 4, 0, meas.Nz / 2, 0, 0) - > GetMinimum();
    fullSignal = -1 * 0.9 * meas.GetHA(CHANNEL, meas.Nx * 3 / 4, 0, meas.Nz / 2, 0, 0) - > GetMinimum();
    break;
  case Y:
    startSignal = -1 * 0.1 * meas.GetHA(CHANNEL, 0, meas.Ny * 3 / 4, meas.Nz / 2, 0, 0) - > GetMinimum();
    fullSignal = -1 * 0.9 * meas.GetHA(CHANNEL, 0, meas.Ny * 3 / 4, meas.Nz / 2, 0, 0) - > GetMinimum();
    break;
  default:
    startSignal = 15;
    fullSignal = 120;
  }

  for (i = 0; i < meas.Nz; i++) {
    foundStart = false;
    foundEnd = false;

    switch (SCANAXIS) {
    case X:
      for (j = 0; j < meas.Nx; j++) {
        TH1F * t1;
        t1 = meas.GetHA(CHANNEL, j, 0, i, 0, 0);
        t1 - > GetXaxis() - > SetRange(time(tStart), time(tEnd));
        if (!foundStart) {
          foundStart = startSignal <= -1 * t1 - > GetMinimum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= -1 * t1 - > GetMinimum();
          dAxis[i] = j - axis0;
        }
      }
      break;

    case Y:
      for (j = 0; j < meas.Ny; j++) {
        TH1F * t1;
        t1 = meas.GetHA(CHANNEL, 0, j, i, 0, 0);
        t1 - > GetXaxis() - > SetRange(time(tStart), time(tEnd));
        if (!foundStart) {
          foundStart = startSignal <= -1 * t1 - > GetMinimum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= -1 * t1 - > GetMinimum();
          dAxis[i] = j - axis0;
        }
      }
      break;

    default:
      for (j = 0; j < meas.Nx; j++) {
        TH1F * t1;
        t1 = meas.GetHA(CHANNEL, j, 0, i, 0, 0);
        t1 - > GetXaxis() - > SetRange(time(tStart), time(tEnd));
        if (!foundStart) {
          foundStart = startSignal <= -1 * t1 - > GetMinimum();
          axis0 = j;
        }
        if (!foundEnd) {
          foundEnd = fullSignal <= -1 * t1 - > GetMinimum();
          dAxis[i] = j - axis0;
        }
      }
      break;
    }
  }

  // Iterate through each Z-coordinates and find the shortest distance from 10% to 90% signal
  int min = 1000, indx = 0;
  for (i = 0; i < meas.Nz; i++) {
    if (dAxis[i] < min) {
      min = dAxis[i];
      indx = i;
    }
  }

  // Print result
  cout << "Focal point at z=" << meas.z0 + meas.dz * indx << " with spot size=" << min << endl;
  free(dAxis);
}