void PlotWaveform()
{

TFile * f = new TFile("DoubleEscapePeak_1619_re.root");
TTree *t = (TTree*)f->Get("tree");

EXOWaveform *wf = 0;
t->SetBranchAddress("fUWireWaveform1",&wf);

t->GetEntry(0);
wf->Decompress();
int *data = wf->GetData();

TH1I *h = new TH1I("h","histo",2048,0,2048);

for (int i = 0; i < 2047; i++) {h->SetBinContent(i+1,data[i]);}

h->Draw();

return;
}
