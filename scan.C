{
	Int_t k=0;
	Int_t numO, numS;	
	MeasureWF *wf[100];
	TGraph *cc[100];
	Float_t width[100], pos[100];
	Float_t optical_axis_co[100];
	Float_t Os;
	
	stct.CorrectBaseLine();
	stct.PrintInfo();	
	Int_t optic_axis=2;
	Int_t scanning_axis=0;
	Float_t LowLim=0;
	Float_t HiLim=100;
	Float_t FWHM=10;
	Float_t Level = 0;
	switch(optic_axis)
	{
		case 0:numO=stct.Nx-1;Os=stct.dx;break;
		case 1:numO=stct.Ny-1;Os=stct.dy;break;
		case 2:numO=stct.Nz-1;Os=stct.dz;break;
	}	
	switch(scanning_axis)
	{
		case 0:numS=stct.Nx-1;break;
		case 1:numS=stct.Ny-1;break;
		case 2:numS=stct.Nz-1;break;
	}
	for(int j=0;j<numO;j++){
		switch(optic_axis)
		{
		case 0: wf[j]=stct.Projection(1,scanning_axis,j,0,0,0,0,numS); break;
		case 1: wf[j]=stct.Projection(1,scanning_axis,0,j,0,0,0,numS); break;
		case 2: wf[j]=stct.Projection(1,scanning_axis,0,0,j,0,0,numS); break;
		}
	 	wf[j]->SetTemperatures(stct.T);
       wf[j]->ScaleHisto(-0.001); //use any value to set the scale
       wf[j]->CorrectBaseLine(1); //correct base line - not needed
       cc[j]=wf[j]->CCE(0,80);   //integrate the charge in time window 0-80 ns
       wf[j]->CCE(0,80)->Draw();
     }

  //--------------------------------------------------------------------
  //define fit function
  
  //  TF1 *ff=new TF1("ff","TMath::Erf((x-[0])/[1])*[2]+[3]",LowLim,HiLim);
  TF1 *ff=new TF1("ff","TMath::Erf((x-[0])/[1])*[2]+[3]",LowLim,HiLim);
  ff->SetParameter(0,Level); 
  ff->SetParameter(1,FWHM); 
  ff->SetParameter(2,8); 
  ff->SetParameter(3,4);

  //-----------------------------------------------------
  //fitting

  for(int j=0;j<=numO;j++)
   {
     //cc[j]->Fit("ff","R");
     width[j]=ff->GetParameter(1)*2.35/TMath::Sqrt(2);
     pos[j]=ff->GetParameter(0);
     optical_axis_co[j]=j*Os;
   }

  //-----------------------------------------------------
  //plot graphs at different positions along optical axis
/*
  for(int j=1;j<=numO;j++)
    {
      cc[j]->SetLineColor(j%8+1);
      if(j==1)
	{
	  cc[j]->Draw("AL"); 
	  cc[j]->GetHistogram()->GetXaxis()->SetTitle("scanning distance [#mum]");
	}
      else cc[j]->Draw("L");
    }

  //-----------------------------------------------------
  //draw the gaussian beam profile 

 TCanvas c2;
 TGraph *FWHMg=new TGraph(numO+1,optical_axis_co,width);
 FWHMg->SetMarkerStyle(21); 
 FWHMg->Draw("APL");
 FWHMg->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
 FWHMg->GetHistogram()->GetYaxis()->SetTitle("FWHM [#mum]");
TCanvas c3;
 TGraph *POSg=new TGraph(numO+1,optical_axis_co,pos);
 POSg->SetMarkerStyle(21);
 POSg->Draw("APL");
 POSg->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
 POSg->GetHistogram()->GetYaxis()->SetTitle("position of the edge [#mum]");
 POSg->Fit("pol1");
*/
}
