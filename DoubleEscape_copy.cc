#include "DoubleEscape.hh"

int main(int argc, char* argv[])
{
  if (argc < 3 || argc > 4) {std::cout << "Usage: DoubleEscape 'filename.[root/dat]' output.root <WaveformFile.root>" << std::endl; return 0;}

  fname = argv[1];
  oname = argv[2];

  WriteWaveforms = false;
  if (argc == 4) {WriteWaveforms = true;}

  // get input filetype
  std::string filename = std::string(fname);
  int id = filename.find_last_of(".") + 1;
  std::string ext = filename.substr(id);

  TChain *t = new TChain("tree","t");

  if (ext.compare("dat") == 0) {
     std::string fnameTMP;
     ifstream ifs(fname);
     while (getline(ifs, fnameTMP)) {
        std::cout << "Adding " << fnameTMP << " to the chain" << std::endl;
        t->Add(fnameTMP.c_str());
     }
  }
  else {t->Add(fname);}

  ED = 0;

  // Set EXOEventData branch address
  t->SetBranchAddress("EventBranch",&ED);

  int nentries = t->GetEntries();
  std::cout << "Number of events: " << nentries << std::endl;

  // create 2D graphs
  gr = new TGraph2D(3);
  grPosEE = new TGraph(2*nentries);
  grPosDEP = new TGraph(nentries);

  grPosEE->SetMarkerStyle(2);
  grPosEE->SetMarkerColor(kBlack);
  grPosDEP->SetMarkerStyle(2);
  grPosDEP->SetMarkerColor(kRed);

  // initialize histograms
  hCSC = new TH1F("hCSC","Scintillation",200,0,10000);

  char hEPCLName[50];
  char hEPCLTitle[50];
  char hESUMName[50];
  char hESUMTitle[50];
  char hWSName[50];
  char hWSTitle[50];
  char hISName[50];
  char hISTitle[50];
  char hEPCLMeanTitle[50];
  char hEPCLSigmaNTitle[50];
  char hEPCLSigmaPTitle[50];

  hEPCL[0] = new TH1F("hEPCL","rec. energy (not purity corrected)",150,0,3000);
  hESUM[0] = new TH1F("hESUM","rec. sum energy (not purity corrected)",200,0,4000);
  hWS[0] = new TH1F("hWS","Weighted sum #alpha=0.1 (not purity corrected)",300,0,6000);
  hIS[0] = new TH2F("hIS","Ion vs. Scint (not purity corrected)",200,0,4000,200,0,10000);

  for (int i = 1; i < 12; i++) {
     sprintf(hEPCLName,"hEPCL%i",i);
     sprintf(hEPCLTitle,"rec. energy (%i #mus)",150+i*20);

     sprintf(hESUMName,"hESUM%i",i);
     sprintf(hESUMTitle,"rec. sum energy (%i #mus)",150+i*20);

     sprintf(hWSName,"hWS%i",i);
     sprintf(hWSTitle,"Weighted sum #alpha=0.1 (%i #mus)",150+i*20);

     sprintf(hISName,"hIS%i",i);
     sprintf(hISTitle,"Ion vs. Scint (%i #mus)",150+i*20);

     hEPCL[i] = new TH1F(hEPCLName,hEPCLTitle,150,0,3000);
     hESUM[i] = new TH1F(hESUMName,hESUMTitle,200,0,4000);
     hWS[i] = new TH1F(hWSName,hWSTitle,300,0,6000);
     hIS[i] = new TH2F(hISName,hISTitle,200,0,4000,200,0,10000);
  }

  hCHI2 = new TH1F("hCHI2","chi2 distribution",250,0,2);
  hPHI = new TH1F("hPHI","phi distribution",100,-1,1);
  hDist = new TH1F("hDist","distance to 511 energy deposit",50,0,400);

  hNUW = new TH1I("hNUW","Number of U wire signals",4,0,4);

  hEPCLTrue = new TH1F("hEPCLTrue","True energy",150,0,3000);

  // set the output tree
  oTree = new TTree("tree","output tree");
  oTree->Branch("fChargeCluster",&OUT_CC);
  oTree->Branch("fScintillationCluster",&OUT_SC);
  oTree->Branch("fUWireWaveform1",&OUT_U1);
  oTree->Branch("fUWireWaveform2",&OUT_U2);
  oTree->Branch("fUWireWaveform3",&OUT_U3);
  oTree->Branch("fUWireWaveform4",&OUT_U4);

  nPointsDEP = 0;
  nPointsEE = 0;

  // load first entra to get trigger time
  t->GetEntry(0);
  double purFitP0 = 293.8;
  double purFitP1 = -1.68;
  double purTime = double(ED->fEventHeader.fTriggerSeconds - 1304146800.0) / 3600.0 / 24.0;
  elifeMean = purFitP1*purTime + purFitP0;
  elifeSigmaN = elifeMean - 2.0;
  elifeSigmaP = elifeMean + 2.0;

  std::cout << "Electron lifetime: " << elifeMean << "(" << elifeSigmaN << "," << elifeSigmaP << ")" << std::endl;

  sprintf(hEPCLMeanTitle,"rec. energy (%i ns)", int(elifeMean));
  sprintf(hEPCLSigmaNTitle,"rec. energy (%i ns)", int(elifeSigmaN));
  sprintf(hEPCLSigmaPTitle,"rec. energy (%i ns)", int(elifeSigmaP));

  hEPCLMean = new TH1F("hEPCLMean",hEPCLMeanTitle,150,0,3000);
  hEPCLSigmaN = new TH1F("hEPCLSigmaN",hEPCLSigmaNTitle,150,0,3000);
  hEPCLSigmaP = new TH1F("hEPCLSigmaP",hEPCLSigmaPTitle,150,0,3000);

  // loop over all entries in tree
  for (evtID = 0; evtID < nentries; evtID++) {
     t->GetEntry(evtID);

     if (evtID % 1000 == 0) {std::cout << evtID << " events processed" << std::endl;}

     DoAnalysis();
  }

  TFile *oFile = new TFile(oname,"RECREATE");

  hCHI2->Write();
  hPHI->Write();
  hCSC->Write();
  hNUW->Write();
  hDist->Write();
  hEPCLTrue->Write();

  TCanvas *c1 = new TCanvas();
  grPosEE->Draw("AP");
  grPosDEP->Draw("Psame");
  c1->Write();

  hEPCLMean->Write();
  hEPCLSigmaN->Write();
  hEPCLSigmaP->Write();

  for (int i = 0; i < 12; i++) {
     hEPCL[i]->Write();
     hESUM[i]->Write();
     hWS[i]->Write();
     hIS[i]->Write();
  }

  oFile->Close();

  if (WriteWaveforms) {
     char *oName2 = argv[3];

     TFile *oFile2 = new TFile(oName2,"RECREATE");
     oTree->Write();
     oFile2->Close();
  }

  return 1;
}

void DoAnalysis()
{
  int nsc = ED->GetNumScintillationClusters();

  // find scintillatino clusters with 3 associated charge clusters
  for (int scID = 0; scID < nsc; scID++) {
     EXOScintillationCluster *scint_cluster = ED->GetScintillationCluster(scID);

     int ncl = scint_cluster->GetNumChargeClusters();

     // select 3 site events
     if (ncl != 3) {continue;}
     if (scint_cluster->fTime > 1900000) {continue;}
  
     // get positions and energies of the three charge clusters
     double cX[3];
     double cY[3];
     double cZ[3];
     double eccl[3];
     double dtcl[3];

     bool GoodEvent = true;
     for (int clID = 0; clID < 3; clID++) {
        EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(clID);
  
        // all three coordinates must be reconstructed
        if (!(charge_cluster->fIs3DCluster)) {GoodEvent = false; break;}
  
        cX[clID] = charge_cluster->fX;
        cY[clID] = charge_cluster->fY;
        cZ[clID] = charge_cluster->fZ;
        eccl[clID] = charge_cluster->fCorrectedEnergy;
        dtcl[clID] = charge_cluster->fDriftTime;
 
        // apply fiducial cut
        double R = TMath::Sqrt(cX[clID]*cX[clID] + cY[clID]*cY[clID]);
        if (R > 163 || cZ[clID] < -180 || cZ[clID] > 180) {GoodEvent = false; break;}
     }
  
     if (!GoodEvent) {continue;}

     // set data points
     for (int k = 0; k < 3; k++) {
        gr->SetPoint(k,cX[k],cY[k],cZ[k]);
     }

     double n1;
     double n2;
     double chi2Fit = FitDataPoints();
     double phi = GetPhi(eccl,cX,cY,cZ,&n1,&n2);

     if (phi > -0.99) {continue;}

     hCHI2->Fill(chi2Fit);

     if (chi2Fit > 0.5) {continue;}

     hDist->Fill(n1);
     hDist->Fill(n2);

     int maxID;
     double maxCluster = max(eccl, 3, &maxID);

     grPosDEP->SetPoint(nPointsDEP,cX[maxID],cY[maxID]);
     for (int i = 0; i < 3; i++) {
        if (i == maxID) {continue;}
        grPosEE->SetPoint(nPointsEE,cX[i],cY[i]);
        nPointsEE++;
     }
     nPointsDEP++;

     double csc = scint_cluster->fCountsSumOnAPDPlaneOne + scint_cluster->fCountsSumOnAPDPlaneTwo;

     hCSC->Fill(csc);

     double alpha = 0.1;
     double w_correction = 1.0;
     double p1 = 1.0;
     double p2 = 0.0;

     hESUM[0]->Fill(eccl[0]+eccl[1]+eccl[3]);
     hIS[0]->Fill(eccl[0]+eccl[1]+eccl[3],csc);
     hWS[0]->Fill((1.0-alpha)*(eccl[0]+eccl[1]+eccl[3]) + alpha*csc);

     for (int i = 0; i < 3; i++) {
        hEPCLMean->Fill(eccl[i]*TMath::Exp(dtcl[i]/elifeMean/1000.0));
        hEPCLSigmaN->Fill(eccl[i]*TMath::Exp(dtcl[i]/elifeSigmaN/1000.0));
        hEPCLSigmaP->Fill(eccl[i]*TMath::Exp(dtcl[i]/elifeSigmaP/1000.0));

        hEPCLTrue->Fill(eccl[i]*TMath::Exp(dtcl[i]/elifeMean/1000.0) * 0.916 + 85);

        hEPCL[0]->Fill(eccl[i]);

        for (int k = 1; k < 12; k++) {
           double elife = 150000 + (k-1)*20000;
           double epcl = eccl[i]*TMath::Exp(dtcl[i]/elife) * p1 + p2;

           if (i == 0) {
              double ESum = (eccl[0]*TMath::Exp(dtcl[0]/elife) + eccl[1]*TMath::Exp(dtcl[1]/elife) + eccl[2]*TMath::Exp(dtcl[2]/elife)) * w_correction;
              hESUM[k]->Fill(ESum);
              hIS[k]->Fill(ESum,csc);
              hWS[k]->Fill((1.0-alpha)*ESum + alpha*csc);
           }

           hEPCL[k]->Fill(epcl);
        }
     }

     // get the double escape cluster
     //int maxID;
     //double maxCluster = max(eccl, 3, &maxID);
 
     OUT_CC = scint_cluster->GetChargeClusterAt(maxID);
     OUT_SC = scint_cluster;

     if (OUT_CC->fDetectorHalf > 1) {continue;}

     int nUWires = OUT_CC->GetNumUWireSignals();
     hNUW->Fill(nUWires);

     // fill puritiy corrected energy
     double dt = OUT_CC->fDriftTime;
     double ep = (OUT_CC->fCorrectedEnergy) * TMath::Exp(dt/245000);
     OUT_CC->fPurityCorrectedEnergy = ep;

     if (WriteWaveforms) {
        if (nUWires > 2) {continue;}
        if (nUWires == 1) {
           int ch = OUT_CC->GetUWireSignalChannelAt(0);

           // dont take events at the boundary
           if (ch == 0 || ch == 37 || ch == 76 || ch == 113) {continue;}

           OUT_U1 = ED->GetWaveformData()->GetWaveformWithChannel(ch-1);
           OUT_U2 = ED->GetWaveformData()->GetWaveformWithChannel(ch);
           OUT_U3 = ED->GetWaveformData()->GetWaveformWithChannel(ch+1);
        }
        if (nUWires == 2) {
           int ch1TMP = OUT_CC->GetUWireSignalChannelAt(0);
           int ch2TMP = OUT_CC->GetUWireSignalChannelAt(1);

           int ch1 = 0;
           int ch2 = 0;
           if (ch1TMP < ch2TMP) {ch1 = ch1TMP; ch2 = ch2TMP;}
           else {ch1 = ch2TMP; ch2 = ch1TMP;}

           if (ch1 == 0 || ch2 == 37 || ch1 == 76 || ch2 == 113) {continue;}

           // make sure the two channels are adjecent
           if (ch2 - ch1 != 1) {continue;}

           OUT_U1 = ED->GetWaveformData()->GetWaveformWithChannel(ch1-1);
           OUT_U2 = ED->GetWaveformData()->GetWaveformWithChannel(ch1);
           OUT_U3 = ED->GetWaveformData()->GetWaveformWithChannel(ch2);
           OUT_U4 = ED->GetWaveformData()->GetWaveformWithChannel(ch2+1);
        }

        oTree->Fill();
     }
  }

  return;
}

double FitDataPoints()
{
  TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
  min->SetObjectFit(gr);
  min->SetFCN(SumDistance2);

  double pStart[4] = {1,1,1,1};
  min->SetParameter(0,"x0",pStart[0],0.01,0,0);
  min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  min->SetParameter(2,"y0",pStart[2],0.01,0,0);
  min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

  Double_t arglist[10];
  arglist[0] = -1;
  min->ExecuteCommand("SET PRINT", arglist, 1);

  arglist[0] = 0;
  min->ExecuteCommand("SET NOW", arglist,0);

  arglist[0] = 1000; // number of function calls
  arglist[1] = 0.001; // tolerance
  min->ExecuteCommand("MIGRAD",arglist,2);

  /*int nvpar,nparx;
  double amin,edm,errdef;
  min->GetStats(amin,edm,errdef,nvpar,nparx);
  min->PrintResults(1,amin);*/
  /*gr->Draw("p0");
  gr->GetXaxis()->SetRangeUser(-200,200);
  gr->GetYaxis()->SetRangeUser(-200,200);
  gr->GetZaxis()->SetRangeUser(-200,200);*/

  // get fit parameters
  double parFit[4];
  for (int i = 0; i < 4; i++) {parFit[i] = min->GetParameter(i);}

  // draw the fitted line
  l = new TPolyLine3D(2);
  double t1 = -200;
  double t2 = 200;
  double x1,y1,z1,x2,y2,z2;
  line(t1,parFit,x1,y1,z1);
  line(t2,parFit,x2,y2,z2);

  l->SetPoint(0,x1,y1,z1);
  l->SetPoint(1,x2,y2,z2);
  l->SetLineColor(kRed);

  return GetChi2(parFit);
}

double GetPhi(double *eccl, double *cX, double *cY, double *cZ, double *n1, double *n2)
{
  int maxID;
  double maxCluster = max(eccl, 3, &maxID);

  XYZVector *x1 = 0;
  XYZVector *x2 = 0;
  if (maxID == 0) {
     x1 = new XYZVector(cX[1]-cX[0],cY[1]-cY[0],cZ[1]-cZ[0]);
     x2 = new XYZVector(cX[2]-cX[0],cY[2]-cY[0],cZ[2]-cZ[0]);
  }

  if (maxID == 1) {
     x1 = new XYZVector(cX[0]-cX[1],cY[0]-cY[1],cZ[0]-cZ[1]);
     x2 = new XYZVector(cX[2]-cX[1],cY[2]-cY[1],cZ[2]-cZ[1]);
  }

  if (maxID == 2) {
     x1 = new XYZVector(cX[0]-cX[2],cY[0]-cY[2],cZ[0]-cZ[2]);
     x2 = new XYZVector(cX[1]-cX[2],cY[1]-cY[2],cZ[1]-cZ[2]);
  }

  double n1TMP;
  double n2TMP;

  n1TMP = TMath::Sqrt(x1->Dot(*x1));
  n2TMP = TMath::Sqrt(x2->Dot(*x2));

  *n1 = n1TMP;
  *n2 = n2TMP;

  XYZVector u1 = x1->Unit();
  XYZVector u2 = x2->Unit();

  double CosT = u1.Dot(u2);

  return CosT;
}

// define the parameteric line equation
void line(double t, double *p, double &x, double &y, double &z)
{
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

// calculate distance line-point
double distance2(double x,double y,double z, double *p)
{
   // distance line point is D= | (xp-x0) cross  ux |
   // where ux is direction of line and x0 is a point in the line (like t = 0)
   XYZVector xp(x,y,z);
   XYZVector x0(p[0], p[2], 0. );
   XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
   XYZVector u = (x1-x0).Unit();
   double d2 = ((xp-x0).Cross(u)) .Mag2();

   return d2;
}

// function to be minimized
void SumDistance2(int &, double *, double &sum, double *par, int )
{
   // the TGraph must be a global variable
   TGraph2D *gr_local = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double *x = gr_local->GetX();
   double *y = gr_local->GetY();
   double *z = gr_local->GetZ();
   int npoints = gr_local->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      sum += d;
   }
}

double GetChi2(double *par)
{
  // the TGraph must be a global variable
  TGraph2D *gr_local = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  assert(gr != 0);
  double *x = gr_local->GetX();
  double *y = gr_local->GetY();
  double *z = gr_local->GetZ();
  int npoints = gr_local->GetN();
  double chi2 = 0.0;
  for (int i  = 0; i < npoints; ++i) {
     double d = distance2(x[i],y[i],z[i],par);
     chi2 += d / 64;

  }

  if (chi2 == 0.0) {return -999.9;}

  return chi2;
}

double max(double *data, int size, int *ID)
{
  double max = -1000000.0;
  int idTMP = 0;

  for (int i = 0; i < size; i++) {
     if (data[i] > max) {max = data[i]; idTMP = i;}
  }

  *ID = idTMP;

  return max;
}

