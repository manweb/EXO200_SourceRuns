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
#include <TFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TGraph2D.h>
#include <TH1.h>

#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOChargeCluster.hh"
#include "EXOUtilities/EXOScintillationCluster.hh"
#include "EXOUtilities/EXO3DView.hh"

using namespace ROOT::Math;

void init();
void PeakDiff(int &, double *, double &fval, double *p, int );
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

TChain *t;
int nentries;

int evtID;
EXOEventData *ED;
TGraph2D *gr;
TPolyLine3D *l;
TH1F *hEPCL;
int ArraySize;
double ECCL[100000];
double DTCL[100000];
