#include "DoubleEscape.hh"

int main(int argc, char* argv[])
{
  std::cout << argc << std::endl;
  if (argc < 3 || argc > 4) {std::cout << "Usage: DoubleEscape 'filename.[root/dat]' output.root [PrintEvents?]" << std::endl; return 0;}

  fname = argv[1];
  oname = argv[2];

  if (argc == 3) {PrintEvents = false;}
  else {PrintEvents = bool(atoi(argv[3]));}

  WriteWaveforms = false;

  // get input filetype
  std::string filename = std::string(fname);
  int id = filename.find_last_of(".") + 1;
  std::string ext = filename.substr(id);

  t = new TChain("tree","t");

  if (ext.compare("dat") == 0) {
     std::string fnameTMP;
     ifstream ifs(fname);
     while (getline(ifs, fnameTMP)) {
        std::cout << "Adding " << fnameTMP << " to the chain" << std::endl;
        t->Add(fnameTMP.c_str());
     }
  }
  else {t->Add(fname);}

  // initialize
  init();

  int nentries = t->GetEntries();
  std::cout << "Number of events: " << nentries << std::endl;

  // loop over all entries in tree
  for (evtID = 0; evtID < nentries; evtID++) {
     t->GetEntry(evtID);

     if (evtID % 1000 == 0 && !PrintEvents) {std::cout << evtID << " events processed" << "  Purity: " << elife << std::endl;}

     // check if this is a noise event
     //if (ED->fEventHeader.fTaggedAsNoise == 1) {continue;}
     //if (ED->fEventHeader.fIndividualTriggerRequest == 0 && ED->fEventHeader.fSumTriggerRequest == 0) {continue;}

     DoAnalysis();
  }

  WriteResults();

  return 1;
}

void init() {
  // Set EXOEventData branch address
  ED = 0;
  t->SetBranchAddress("EventBranch",&ED);

  // create 2D graph
  gr = new TGraph2D(3);

  hCHI2 = new TH1F("hCHI2","chi2 distribution",250,0,2);
  hPHI = new TH1F("hPHI","phi distribution",100,-1,1);

  // set 511 gammas output tree
  oTreeLP = new TTree("treeLP","511 gammas");
  oTreeLP->Branch("fecclLP",&fecclLP,"fecclLP/d");
  oTreeLP->Branch("fepclLP",&fepclLP,"fepclLP/D");
  oTreeLP->Branch("fetclLP",&fetclLP,"fetclLP/D");
  oTreeLP->Branch("fecclLP_gc",&fecclLP_gc,"fecclLP_gc/d");
  oTreeLP->Branch("fepclLP_gc",&fepclLP_gc,"fepclLP_gc/D");
  oTreeLP->Branch("fetclLP_gc",&fetclLP_gc,"fetclLP_gc/D");
  oTreeLP->Branch("fdtclLP",&fdtclLP,"fdtclLP/D");
  oTreeLP->Branch("fXLP",&fXLP,"fXLP/D");
  oTreeLP->Branch("fYLP",&fYLP,"fYLP/D");
  oTreeLP->Branch("fZLP",&fZLP,"fZLP/D");
  oTreeLP->Branch("fDist",&fDist,"fDist/D");

  // set the double escape peak output tree
  oTreeUP = new TTree("treeUP","Double escape peak");
  oTreeUP->Branch("fnUWires",&fnUWires,"fnUWires/I");
  oTreeUP->Branch("fecclUP",&fecclUP,"fecclUP/D");
  oTreeUP->Branch("fepclUP",&fepclUP,"fepclUP/D");
  oTreeUP->Branch("fetclUP",&fetclUP,"fetclUP/D");
  oTreeUP->Branch("fecclUP_gc",&fecclUP_gc,"fecclUP_gc/D");
  oTreeUP->Branch("fepclUP_gc",&fepclUP_gc,"fepclUP_gc/D");
  oTreeUP->Branch("fetclUP_gc",&fetclUP_gc,"fetclUP_gc/D");
  oTreeUP->Branch("fdtclUP",&fdtclUP,"fdtclUP/D");
  oTreeUP->Branch("fXUP",&fXUP,"fXUP/D");
  oTreeUP->Branch("fYUP",&fYUP,"fYUP/D");
  oTreeUP->Branch("fZUP",&fZUP,"fZUP/D");
  oTreeUP->Branch("fcsc",&fcsc,"fcsc/D");
  oTreeUP->Branch("fesum",&fesum,"fesum/D");

  // set the waveform output tree
  oTreeWF = new TTree("treeWF","output tree");
  oTreeWF->Branch("fChargeCluster",&OUT_CC);
  oTreeWF->Branch("fScintillationCluster",&OUT_SC);
  oTreeWF->Branch("fUWireWaveform1",&OUT_U1);
  oTreeWF->Branch("fUWireWaveform2",&OUT_U2);
  oTreeWF->Branch("fUWireWaveform3",&OUT_U3);
  oTreeWF->Branch("fUWireWaveform4",&OUT_U4);

  // load first entry to get trigger time
  /*t->GetEntry(0);
  double purFitP0 = 293.8;
  double purFitP1 = -1.68;
  double purTime = double(ED->fEventHeader.fTriggerSeconds - 1304146800.0) / 3600.0 / 24.0;
  elife = purFitP1*purTime + purFitP0;

  std::cout << "Electron lifetime: " << elife << std::endl;*/

  return;
}

void WriteResults() {
  TFile *oFile = new TFile(oname,"RECREATE");

  hCHI2->Write();
  hPHI->Write();

  oTreeLP->Write();
  oTreeUP->Write();
  oTreeWF->Write();

  oFile->Close();

  return;
}

void DoAnalysis()
{
  // get the current purity
  //double purFitP0 = 293.8;
  //double purFitP1 = -1.68;
  //double purFitP0 = 286;
  //double purFitP1 = -1.427;
  double purFitP0;
  double purFitP1;
  double purFitP2;
  double purFitP3;
  double purFitP4;

  double purTime = double(ED->fEventHeader.fTriggerSeconds - 1304146800.0) / 3600.0 / 24.;;

  if (purTime < 58) {
     purFitP0 = -284.596;
     purFitP1 = 53.6978;
     purFitP2 = -1.88664;
     purFitP3 = 0.0269101;
     purFitP4 = -0.000133772;
  }
  if (purTime >= 58 && purTime < 81.6) {
     purFitP0 = 14068.5;
     purFitP1 = -908.011;
     purFitP2 = 21.8864;
     purFitP3 = -0.230994;
     purFitP4 = 0.00090631;
  }
  if (purTime >= 81.6 && purTime < 94.0) {
     purFitP0 = -9011.55;
     purFitP1 = 115.417;
     purFitP2 = 0.0;
     purFitP3 = 0.0;
     purFitP4 = 0.0;
  }
  if (purTime >= 94.0 && purTime < 102.5) {
     purFitP0 = 2000.0;
     purFitP1 = 0.0;
     purFitP2 = 0.0;
     purFitP3 = 0.0;
     purFitP4 = 0.0;
  }
  if (purTime >= 102.5 && purTime < 113.0) {
     purFitP0 = -1208000.0;
     purFitP1 = 34380.0;
     purFitP2 = -325.9;
     purFitP3 = 1.03;
     purFitP4 = 0.0;
  }
  if (purTime >= 113.0 && purTime < 129.6) {
     purFitP0 = -48740.0;
     purFitP1 = 805.0;
     purFitP2 = -3.259;
     purFitP3 = 0.0;
     purFitP4 = 0.0;
  }
  if (purTime >= 129.6 && purTime < 142.0) {
     purFitP0 = -29510.0;
     purFitP1 = 230.1;
     purFitP2 = 0.0;
     purFitP3 = 0.0;
     purFitP4 = 0.0;
  }
  if (purTime >= 142.0) {
     purFitP0 = 3300.0;
     purFitP1 = 0.0;
     purFitP2 = 0.0;
     purFitP3 = 0.0;
     purFitP4 = 0.0;
  }

  elife = purFitP4*purTime*purTime*purTime*purTime + purFitP3*purTime*purTime*purTime + purFitP2*purTime*purTime + purFitP1*purTime + purFitP0;

  int nsc = ED->GetNumScintillationClusters();

  // find scintillation clusters with 3 associated charge clusters
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
     double eccl_gc[3];
     double dtcl[3];
     int nUWires[3];

     bool GoodEvent = true;
     for (int clID = 0; clID < 3; clID++) {
        EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(clID);
  
        // all three coordinates must be reconstructed
        //if (!(charge_cluster->fIs3DCluster)) {GoodEvent = false; break;}
  
        cX[clID] = charge_cluster->fX;
        cY[clID] = charge_cluster->fY;
        cZ[clID] = charge_cluster->fZ;

        if (cX[clID] == -999 || cY[clID] == -999) {GoodEvent = false; break;}

        eccl[clID] = charge_cluster->fPurityCorrectedEnergy;
        //eccl[clID] = charge_cluster->fRawEnergy;
        dtcl[clID] = charge_cluster->fDriftTime;
        nUWires[clID] = charge_cluster->GetNumUWireSignals();
 
        int nU = charge_cluster->GetNumUWireSignals();
        eccl_gc[clID] = 0.0;
        for (int UID = 0; UID < nU; UID++) {
           EXOUWireSignal *u_signal = charge_cluster->GetUWireSignalAt(UID);

           int UCHID = u_signal->fChannel;
           double U_corr = 1.0;
           if (UCHID >= 0 && UCHID < 114) {U_corr = GainCorrection[UCHID];}
           double eU = u_signal->fCorrectedEnergy * U_corr;
           eccl_gc[clID] += eU;
        }

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

     hPHI->Fill(phi);

     // cut on phi angle
     if (phi > -0.99) {continue;}

     hCHI2->Fill(chi2Fit);

     // cut on chi2
     if (chi2Fit > 0.5) {continue;}

     // get the double escape charge cluster
     int maxID;
     double maxCluster = max(eccl, 3, &maxID);

     OUT_CC = scint_cluster->GetChargeClusterAt(maxID);
     OUT_SC = scint_cluster;

     fnUWires = OUT_CC->GetNumUWireSignals();
     fecclUP = eccl[maxID];
     fepclUP = eccl[maxID] * TMath::Exp(dtcl[maxID]/elife/1000.0);
     fetclUP = fepclUP * 0.916 + 85;
     fdtclUP = dtcl[maxID];
     fecclUP_gc = eccl_gc[maxID];
     fepclUP_gc = eccl_gc[maxID] * TMath::Exp(dtcl[maxID]/elife/1000.0);
     fetclUP_gc = fepclUP_gc * 0.916 + 85;
     fXUP = cX[maxID];
     fYUP = cY[maxID];
     fZUP = cZ[maxID];
     fesum = (eccl[0] * TMath::Exp(dtcl[0]/elife/1000.0) + eccl[1] * TMath::Exp(dtcl[1]/elife/1000.0) + eccl[2] * TMath::Exp(dtcl[2]/elife/1000.0)) * 0.916 + 85;
     fcsc = scint_cluster->fRawEnergy;

     //if (nUWires[maxID] == 1) {oTreeUP->Fill();}
     oTreeUP->Fill();
    
     if (PrintEvents) {std::cout << ED->fEventNumber << "\t" << scID << "\t" << maxID << std::endl;}

     // get 511 gammas
     bool DistFirst = true;
     for (int i = 0; i < 3; i++) {
        if (i == maxID) {continue;}

        fecclLP = eccl[i];
        fepclLP = eccl[i] * TMath::Exp(dtcl[i]/elife/1000.0);
        fetclLP = fepclLP * 0.916 + 85;
        fecclLP_gc = eccl_gc[i];
        fepclLP_gc = eccl_gc[i] * TMath::Exp(dtcl[i]/elife/1000.0);
        fetclLP_gc = fepclLP_gc * 0.916 + 85;
        fdtclLP = dtcl[i];
        fXLP = cX[i];
        fYLP = cY[i];
        fZLP = cZ[i];
        if (DistFirst) {fDist = n1; DistFirst = false;}
        else {fDist = n2;}

        //if (nUWires[i] == 1) {oTreeLP->Fill();}
        oTreeLP->Fill();
     }

     if (OUT_CC->fDetectorHalf > 1) {continue;}

     // check if the waveforms are in the file
     EXOWaveformData *wf_data = ED->GetWaveformData();
     wf_data->Decompress();

     if (wf_data->fNumSamples == 0) {continue;}

     // fill puritiy corrected energy
     double dt = OUT_CC->fDriftTime;
     double ep = (OUT_CC->fCorrectedEnergy) * TMath::Exp(dt/elife/1000.0);
     OUT_CC->fPurityCorrectedEnergy = ep;

     if (fnUWires > 2) {continue;}
     if (fnUWires == 1) {
        int ch = OUT_CC->GetUWireSignalChannelAt(0);

        // dont take events at the boundary
        if (ch == 0 || ch == 37 || ch == 76 || ch == 113) {continue;}

        OUT_U1 = ED->GetWaveformData()->GetWaveformWithChannel(ch-1);
        OUT_U2 = ED->GetWaveformData()->GetWaveformWithChannel(ch);
        OUT_U3 = ED->GetWaveformData()->GetWaveformWithChannel(ch+1);
     }
     if (fnUWires == 2) {
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

     oTreeWF->Fill();
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

