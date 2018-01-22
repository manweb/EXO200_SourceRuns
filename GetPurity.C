void GetPurity()
{
  TTree *t = new TTree;
  t->ReadFile("Manuel_Th_elife.dat");

  int runStart;
  int runStop;
  double elife;
  double elifeError;

  t->SetBranchAddress("runStart",&runStart);
  t->SetBranchAddress("runStop",&runStop);
  t->SetBranchAddress("Elife",&elife);
  t->SetBranchAddress("ElifeError",&elifeError);

  int Entries = t->GetEntries();

  double p0 = 286.0;
  double p1 = -1.427;

  cout << "time\tpurity" << endl;
  cout << "*****************************" << endl;
  for (int i=0; i<Entries; i++) {
     t->GetEntry(i);

     double time = ((runStart-1304146800)+(runStop-runStart)/2)/3600/24;
     double purity = p1*time + p0;

     cout << time << "\t" << purity << endl;
  }
}

