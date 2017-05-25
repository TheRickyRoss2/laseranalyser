{

	gSystem->Load("~/Downloads/TCTAnalyse.V2.0/TCTAnalyse.sl");
	
  char * file = (char *) "./example.rtct";
	
  gROOT->ProcessLine(".x LoadPSTCT.C");
	stct.PrintInfo();
  stct.Draw(1,0,2,0,0,0,65,85)->Draw("SURF2Z");
}
