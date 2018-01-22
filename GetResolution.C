double fitFunction(double *x, double *par)
{
  double val = errf(x,par) + gauss(x,&par[1]);

  return val;
}

double errf(double *x, double *par)
{
  return par[0]*TMath::Erfc((x[0]-par[2]) / (TMath::Sqrt(2)*par[3]));
  //return par[0] * (1 - 1 / (1 + TMath::Exp(-par[1]*(x[0] - par[2]))));
}

double gauss(double *x, double *par)
{
  return par[0]*TMath::Gaus(x[0],par[1],par[2]);
}

void GetResolution()
{
  //TFile *f = new TFile("AllTh_1573_1940.root","READ");
  //TFile *f = new TFile("2421.root","READ");
  //TFile *f = new TFile("2424_2448_new.root","READ");
  //TTree *tLP = (TTree*)f->Get("treeLP");
  //TTree *tUP = (TTree*)f->Get("treeUP");

  TChain *tLP = new TChain("treeLP");
  TChain *tUP = new TChain("treeUP");

  /*tLP->Add("2424_new.root");
  tLP->Add("2426_new.root");
  tLP->Add("2431_new.root");
  tLP->Add("2432_new.root");
  tLP->Add("2433_new.root");
  tLP->Add("2434_new.root");
  tLP->Add("2447_new.root");
  tLP->Add("2448_new.root");

  tUP->Add("2424_new.root");
  tUP->Add("2426_new.root");
  tUP->Add("2431_new.root");
  tUP->Add("2432_new.root");
  tUP->Add("2433_new.root");
  tUP->Add("2434_new.root");
  tUP->Add("2447_new.root");
  tUP->Add("2448_new.root");*/
  
  tLP->Add("../analysis/200412/*.root");
  tUP->Add("../analysis/200412/*.root");

  TH1F *hLP = new TH1F("hLP","511 peak",100,0,1500);
  TH1F *hUP = new TH1F("hUP","Double escape peak",150,0,3000);

  //tLP->Draw("0.9154*fepclLP+82.9>>hLP");
  //tUP->Draw("0.9154*fepclUP+82.9>>hUP");
  //tLP->Draw("0.9278*fepclLP+59.81>>hLP");
  //tUP->Draw("0.9278*fepclUP+59.81>>hUP");
  //tLP->Draw("fepclLP>>hLP");
  //tUP->Draw("fepclUP>>hUP");
  
  // Reprocessed data (01/30/12)
  //tLP->Draw("fecclLP*0.9466 + 74.5>>hLP");
  //tUP->Draw("(fecclUP*0.9466 + 74.5)*0.973>>hUP");
  
  // Reprocessed data (02/07/12)
  //tLP->Draw("fecclLP*0.9414 + 69.35>>hLP");
  //tUP->Draw("(fecclUP*0.9414 + 69.35)*0.9675>>hUP");
  //tLP->Draw("fecclLP*0.9334 + 68.57>>hLP");
  //tUP->Draw("fecclUP*0.9334 + 68.57>>hUP");
  
  // Reprocessed data (04/20/12)
  tLP->Draw("fecclLP*0.924 + 72.7>>hLP");
  tUP->Draw("fecclUP*0.924 + 72.7>>hUP");

  TF1 *fit1 = new TF1("fit1",fitFunction,200,650,4);
  //TF1 *fit2 = new TF1("fit2","gaus",1500,1900);
  TF1 *fit2 = new TF1("fit2","gaus",1480,1780);
  TF1 *FitGauss = new TF1("FitGauss",gauss,0,1000,3);
  TF1 *FitErrf = new TF1("FitErrf",errf,0,1000,4);

  FitGauss->SetLineColor(kGreen);
  FitErrf->SetLineColor(kBlue);
  //fit1->SetLineColor(kRed);
  //fit2->SetLineColor(kRed);
  //fit1->SetLineWidth(1);
  //fit2->SetLineWidth(1);

  fit1->SetParNames("A","B","x0","#sigma");
  fit1->SetParameters(4,10,511,60);

  fit1->SetParLimits(0,0,10000);
  fit1->SetParLimits(1,0,10000);
  fit1->SetParLimits(2,400,600);
  fit1->SetParLimits(3,10,200);

  hLP->Fit("fit1","r");
  hUP->Fit("fit2","r");

  double par1[4];
  double par2[3];
  double *par1Err;
  double *par2Err;

  fit1->GetParameters(par1);
  fit2->GetParameters(par2);
  par1Err = fit1->GetParErrors();
  par2Err = fit2->GetParErrors();

  FitGauss->SetParameters(&par1[1]);
  FitErrf->SetParameters(par1);

  double resolution1 = par1[3]/par1[2];
  double resolution2 = par2[2]/par2[1];
  double resolution1_err = TMath::Sqrt((par1[3]/par1[2]**2)**2 * par1Err[2]**2 + (1.0/par1[2])**2 * par1Err[3]**2);
  double resolution2_err = TMath::Sqrt((par2[2]/par2[1]**2)**2 * par2Err[1]**2 + (1.0/par2[1])**2 * par2Err[2]**2);

  //cout << "Resolution: Peak1 = " << par1[3] << " (" << par1[3]/par1[2]*100 << "%)  Peak2 = " << par2[2] << " (" << par2[2]/par2[1]*100 << "%)" << endl;
  cout << "Resolution: Peak1 = " << par1[3] << " (" << par1[3]/par1[2] << " +- " << resolution1_err << ")  Peak2 = " << par2[2]/par2[1] << " (" << par2[2]/par2[1] << " +- " << resolution2_err << ")" << endl;

  TCanvas *c1 = new TCanvas();
  hLP->Draw();
  fit1->Draw("same");
  FitErrf->Draw("same");
  FitGauss->Draw("same");

  TCanvas *c2 = new TCanvas();
  hUP->Draw();
  fit2->Draw("same");

  return;
}
