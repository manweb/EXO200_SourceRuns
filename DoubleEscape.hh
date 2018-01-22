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

void init();
void WriteResults();
void DoAnalysis();
double FitDataPoints();
double GetPhi(double *eccl, double *cX, double *cY, double *cZ, double *n1, double *n2);
void line(double t, double *p, double &x, double &y, double &z);
double distance2(double x,double y,double z, double *p);
void SumDistance2(int &, double *, double &sum, double *par, int );
double GetChi2(double *par);
double max(double *data, int size, int *ID);

double GainCorrection[114] ={1.0063,0.99631,0.95355,0.98478,1.0008,0.97925,0.97718,0.97029,0.95922,0.9997,1.0057,1.003,0.99045,1.0057,1.0056,1.0054,0.94392,1.0269,1.012,1.0023,0.99201,1.0132,1.003,1.0087,1.0061,1.0144,0.98746,0.99365,0.99455,0.9931,1.0019,0.97818,0.99212,0.97474,0.89061,0.99739,0.99313,0.90054,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0075,1.0083,1.0017,1.0243,1.0395,1.0274,1.0304,1.0064,0.99619,1.027,0.99832,1.0113,0.99871,1.0122,1.0006,0.98523,0.99921,1.0162,1.0269,0.99464,0.97529,0.99449,0.99837,0.98738,1.0005,0.98692,1.0051,0.99625,0.98746,1.028,1.0159,1.0279,1.0015,1.0108,1.0035,1.009,0.99671,0.99012};

char *fname;
char *oname;

bool PrintEvents;

TChain *t;

TTree *oTreeLP;
TTree *oTreeUP;
TTree *oTreeWF;

bool WriteWaveforms;
EXOChargeCluster *OUT_CC;
EXOScintillationCluster *OUT_SC;
const EXOWaveform *OUT_U1;
const EXOWaveform *OUT_U2;
const EXOWaveform *OUT_U3;
const EXOWaveform *OUT_U4;

int evtID;
EXOEventData *ED;
TGraph2D *gr;

double elife;

double fecclLP;
double fepclLP;
double fetclLP;
double fecclLP_gc;
double fepclLP_gc;
double fetclLP_gc;
double fdtclLP;
double fXLP;
double fYLP;
double fZLP;
double fDist;

int fnUWires;
double fecclUP;
double fepclUP;
double fetclUP;
double fecclUP_gc;
double fepclUP_gc;
double fetclUP_gc;
double fdtclUP;
double fXUP;
double fYUP;
double fZUP;
double fcsc;
double fesum;

TH1F *hCHI2;
TH1F *hPHI;
