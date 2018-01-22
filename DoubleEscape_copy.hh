#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TH1.h>

#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOChargeCluster.hh"
#include "EXOUtilities/EXOScintillationCluster.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXO3DView.hh"

using namespace ROOT::Math;

void DoAnalysis();
double FitDataPoints();
double GetPhi(double *eccl, double *cX, double *cY, double *cZ, double *n1, double *n2);
void line(double t, double *p, double &x, double &y, double &z);
double distance2(double x,double y,double z, double *p);
void SumDistance2(int &, double *, double &sum, double *par, int );
double GetChi2(double *par);
double max(double *data, int size, int *ID);

char *fname;
char *oname;

bool WriteWaveforms;
TTree *oTree;
EXOChargeCluster *OUT_CC;
EXOScintillationCluster *OUT_SC;
const EXOWaveform *OUT_U1;
const EXOWaveform *OUT_U2;
const EXOWaveform *OUT_U3;
const EXOWaveform *OUT_U4;

int nPointsDEP;
int nPointsEE;
TGraph *grPosEE;
TGraph *grPosDEP;

int evtID;
EXOEventData *ED;
TGraph2D *gr;
TPolyLine3D *l;

double elifeMean;
double elifeSigmaN;
double elifeSigmaP;

TH1F *hCHI2;
TH1F *hPHI;
TH1F *hDist;
TH1F *hEPCLMean;
TH1F *hEPCLSigmaN;
TH1F *hEPCLSigmaP;
TH1F *hEPCLTrue;
TH1F *hEPCL[12];
TH1F *hESUM[12];
TH1F *hCSC;
TH1F *hWS[12];
TH2F *hIS[12];
TH1I *hNUW;
