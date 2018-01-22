#include "MinimizePlot.hh"

int main(int argc, char* argv[])
{
  if (argc != 3) {std::cout << "Usage: MinimizePlot 'filename.[root/dat]' output.root" << std::endl; return 0;}

  fname = argv[1];
  oname = argv[2];

  t = new TChain("tree","t");

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

  nentries = t->GetEntries();

  // create 2D graph
  gr = new TGraph2D(3);

  hEPCL = new TH1F("hEPCL","rec. energy",150,0,3000);

  ArraySize = 0;

  // loop over all entries in tree and create histogram
  for (evtID = 0; evtID < nentries; evtID++) {
     t->GetEntry(evtID);

     DoAnalysis();
  }

  int nrStep = 200;
  double par1Limits[2] = {0.84,0.92};
  double par2Limits[2] = {60,170};

  hMinimize = new TH2F("hMinimize","Minimization Plot",nrStep,par2Limits[0],par2Limits[1],nrStep,par1Limits[0],par1Limits[1]);


  double val = 0.0;
  for (int i = 0; i < nrStep; i++) {
     double par1 = par1Limits[0] + i * (par1Limits[1] - par1Limits[0]) / nrStep;
     for (int k = 0; k < nrStep; k++) {
        double par2 = par2Limits[0] + k * (par2Limits[1] - par2Limits[0]) / nrStep;

        val = PeakDiff(par1, par2);
        hMinimize->SetBinContent(k+1,i+1,val);
     }
  }

  TFile *oFile = new TFile(oname,"RECREATE");
  hEPCL->Write();
  hMinimize->Write();
  oFile->Close();

  return 1;
}

double PeakDiff(double p1, double p2)
{
  hEPCL->Clear();
  hEPCL = new TH1F("hEPCL","rec. energy",150,0,3000);

  //std::cout << "creating histogram with parameters: p1 = " << p1 << ", p2 = " << p2 << " ... ";

  double elife = 250000;
  for (int i = 0; i < ArraySize; i++) {
     double ecorr = ECCL[i] * TMath::Exp(DTCL[i]/elife) * p1 + p2;

     hEPCL->Fill(ecorr);
  }

  //std::cout << " done." << std::endl;

  double refP1_1 = 300.0;
  double refP2_1 = 600.0;
  double refP1_2 = 1500.0;
  double refP2_2 = 1900.0;

  TF1 *fit1 = new TF1("fit1","gaus",refP1_1*p1+p2,refP2_1*p1+p2);
  TF1 *fit2 = new TF1("fit2","gaus",refP1_2*p1+p2,refP2_2*p1+p2);

  //std::cout << "fitting histogram ... ";

  hEPCL->Fit("fit1","Rq");
  hEPCL->Fit("fit2","Rq");

  //std::cout << "done." << std::endl;

  double x1 = fit1->GetParameter(1);
  double x2 = fit2->GetParameter(1);
  double s1 = fit1->GetParameter(2);
  double s2 = fit2->GetParameter(2);

  double trueX1 = 511;
  double trueX2 = 1592;

  double val = (x1 - trueX1)*(x1 - trueX1)/(s1*s1) + (x2 - trueX2)*(x2 - trueX2)/(s2*s2);

  //std::cout << val << std::endl;

  if (val > 10000) {return 0.0;}

  return val;
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
 
        // apply position cut
        if (cX[clID] < -180 || cX[clID] > 180 || cY[clID] < -180 || cY[clID] > 180 || cZ[clID] < -180 || cZ[clID] > 180) {GoodEvent = false; break;}
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
     if (chi2Fit > 0.5) {continue;}

     /*int maxID;
     double maxCharge = max(eccl, 3, &maxID);

     if (maxCharge < 900.0) {continue;}
     for (int i = 0; i < 3; i++) {
        if (i == maxID) {continue;}
        if (eccl[i] > 900.0) {GoodEvent = false; break;}
     }

     if (!GoodEvent) {continue;}*/

     for (int i = 0; i < 3; i++) {
        ECCL[ArraySize] = eccl[i];
        DTCL[ArraySize] = dtcl[i];

        ArraySize++;
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
  min->ExecuteCommand("SET NOW", arglist, 0);

  arglist[0] = 1000; // number of function calls
  arglist[1] = 0.001; // tolerance
  min->ExecuteCommand("MIGRAD",arglist,2);

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

