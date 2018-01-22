double fitFunction(double *x, double *par)
{
  double val = errf(x,par) + gauss(x,&par[1]) + par[4];

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

void GetPeakValue()
{
  const int nFilesC = 12;
  const int nFilesA = 5;
  const int N =  nFilesC + nFilesA;

  string fnamesC[nFilesC] = {"1573","1595","1605","1619","1638","1647","1918","1926","1930","1931","1932","1937"};
  string fnamesA[nFilesA] = {"1587","1599","1631","1642","1929"};

  TGraphErrors *gr1C = new TGraphErrors(nFilesC);
  TGraphErrors *gr1A = new TGraphErrors(nFilesA);
  TGraphErrors *gr2C = new TGraphErrors(nFilesC);
  TGraphErrors *gr2A = new TGraphErrors(nFilesA);

  gr1C->SetTitle("511 peak");
  gr1A->SetTitle("511 peak");
  gr2C->SetTitle("double escape peak");
  gr2A->SetTitle("double escape peak");

  gr1C->GetXaxis()->SetTitle("run number");
  gr1C->GetYaxis()->SetTitle("peak position");
  gr1A->GetXaxis()->SetTitle("run number");
  gr1A->GetYaxis()->SetTitle("peak position");

  gr2C->GetXaxis()->SetTitle("run number");
  gr2C->GetYaxis()->SetTitle("peak position");
  gr2A->GetXaxis()->SetTitle("run number");
  gr2A->GetYaxis()->SetTitle("peak position");

  gr1C->SetMarkerStyle(20);
  gr1C->SetMarkerColor(kBlack);

  gr1A->SetMarkerStyle(20);
  gr1A->SetMarkerColor(kRed);

  gr2C->SetMarkerStyle(20);
  gr2C->SetMarkerColor(kBlack);

  gr2A->SetMarkerStyle(20);
  gr2A->SetMarkerColor(kRed);

  double *Peak1Mean = new double[N];
  double *Peak2Mean = new double[N];

  // create histograms
  //TH1F *hLP = new TH1F("hLP","511 peak",150,0,3000);
  //TH1F *hUP = new TH1F("hUP","Double escape peak",150,0,3000);

  // define fit functions
  TF1 *fitLP = new TF1("fitLP","gaus",350,600);
  TF1 *fitUP = new TF1("fitUP","gaus",1500,1900);

/*  fitLP->SetParNames("A","B","x0","#sigma");
  fitLP->SetParameters(4,10,511,60);

  fitLP->SetParLimits(0,0,10000);
  fitLP->SetParLimits(1,0,10000);
  fitLP->SetParLimits(2,400,600);
  fitLP->SetParLimits(3,20,200);
*/
  cout << "Run\tPeak1\tPeak2" << endl;

  // loop over cathode runs
  int count = 0;
  for (int i = 0; i < nFilesC; i++) {
     string name = fnamesC[i] + ".root";
     TFile *f = new TFile(name.c_str(),"READ");

     TTree *tLP = (TTree*)f->Get("treeLP");
     TTree *tUP = (TTree*)f->Get("treeUP");

     // create histograms
     TH1F *hLP = new TH1F("hLP","511 peak",150,0,3000);
     TH1F *hUP = new TH1F("hUP","Double escape peak",150,0,3000);

     // fill histograms
     tLP->Draw("fepclLP>>hLP");
     tUP->Draw("fepclUP>>hUP");

     // fit the peaks
     hLP->Fit("fitLP","rq");
     hUP->Fit("fitUP","rq");

     cout << fnamesC[i] << "\t" << fitLP->GetParameter(1) << " +- " << fitLP->GetParError(1) << "\t" << fitUP->GetParameter(1) << " +- " << fitUP->GetParError(1) << endl;

     // set data points
     gr1C->SetPoint(i,atoi(fnamesC[i].c_str()),fitLP->GetParameter(1));
     gr2C->SetPoint(i,atoi(fnamesC[i].c_str()),fitUP->GetParameter(1));

     // set errors
     gr1C->SetPointError(i,0,fitLP->GetParError(1));
     gr2C->SetPointError(i,0,fitUP->GetParError(1));

     Peak1Mean[count] = fitLP->GetParameter(1);
     Peak2Mean[count] = fitUP->GetParameter(1);

     count++;
  }

  // loop over anode runs
  for (int i = 0; i < nFilesA; i++) {
     string name = fnamesA[i] + ".root";
     TFile *f = new TFile(name.c_str(),"READ");

     TTree *tLP = (TTree*)f->Get("treeLP");
     TTree *tUP = (TTree*)f->Get("treeUP");

     // create histograms
     TH1F *hLP = new TH1F("hLP","511 peak",75,0,3000);
     TH1F *hUP = new TH1F("hUP","Double escape peak",75,0,3000);

     // fill histograms
     tLP->Draw("fepclLP>>hLP");
     tUP->Draw("fepclUP>>hUP");

     // fit the peaks
     hLP->Fit("fitLP","rqn");
     hUP->Fit("fitUP","rqn");

     cout << fnamesA[i] << "\t" << fitLP->GetParameter(1) << " +- " << fitLP->GetParError(1) << "\t" << fitUP->GetParameter(1) << " +- " << fitUP->GetParError(1) << endl;

     // set data points
     gr1A->SetPoint(i,atoi(fnamesA[i].c_str()),fitLP->GetParameter(1));
     gr2A->SetPoint(i,atoi(fnamesA[i].c_str()),fitUP->GetParameter(1));

     // set errors
     gr1A->SetPointError(i,0,fitLP->GetParError(1));
     gr2A->SetPointError(i,0,fitUP->GetParError(1));

     Peak1Mean[count] = fitLP->GetParameter(1);
     Peak2Mean[count] = fitUP->GetParameter(1);

     count++;
  }

  // calculate mean and rms
  double mean1 = 0;
  double mean2 = 0;
  double rms21 = 0;
  double rms22 = 0;
  for (int i = 0; i < N; i++) {
     mean1 += Peak1Mean[i] / N;
     mean2 += Peak2Mean[i] / N;
  }

  for (int i = 0; i < N; i++) {
     rms21 += (Peak1Mean[i] - mean1) * (Peak1Mean[i] - mean1);
     rms22 += (Peak2Mean[i] - mean2) * (Peak2Mean[i] - mean2);
  }

  double rms1 = TMath::Sqrt(rms21/(N-1));
  double rms2 = TMath::Sqrt(rms22/(N-1));

  cout << mean1 << " +- " << rms1 << "  " << mean2 << " +-" << rms2 << endl;

  // create lines that inidcate mean and rms
  TLine *lMean1 = new TLine(1567,mean1,1940,mean1);
  TLine *lSigN1 = new TLine(1567,mean1-rms1,1940,mean1-rms1);
  TLine *lSigP1 = new TLine(1567,mean1+rms1,1949,mean1+rms1);

  TLine *lMean2 = new TLine(1567,mean2,1940,mean2);
  TLine *lSigN2 = new TLine(1567,mean2-rms2,1940,mean2-rms2);
  TLine *lSigP2 = new TLine(1567,mean2+rms2,1940,mean2+rms2);

  lMean1->SetLineColor(16);
  lMean2->SetLineColor(16);

  lSigN1->SetLineColor(16);
  lSigP1->SetLineColor(16);
  lSigN2->SetLineColor(16);
  lSigP2->SetLineColor(16);

  lSigN1->SetLineStyle(2);
  lSigP1->SetLineStyle(2);
  lSigN2->SetLineStyle(2);
  lSigP2->SetLineStyle(2);

  // text for result output
  char result1[50];
  char result2[50];

  sprintf(result1,"%.2f +- %.2f",mean1,rms1);
  sprintf(result2,"%.2f +- %.2f",mean2,rms2);

  TText *resultText1 = new TText(1610,450,result1);
  TText *resultText2 = new TText(1610,1700,result2);

  TCanvas *c1 = new TCanvas();
  gr1C->Draw("AP");
  gr1A->Draw("Psame");
  lMean1->Draw("same");
  lSigN1->Draw("same");
  lSigP1->Draw("same");
  resultText1->Draw("same");

  TCanvas *c2 = new TCanvas();
  gr2C->Draw("AP");
  gr2A->Draw("Psame");
  lMean2->Draw("same");
  lSigN2->Draw("same");
  lSigP2->Draw("same");
  resultText2->Draw("same");

  return;
}
