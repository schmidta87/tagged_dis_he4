#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"

using namespace std;

const int nArgs=4;

int main(int argc, char** argv)
{
  if (argc != nArgs)
    {
      cerr << "Wrong number of arguments.\n"
      << "Call the program with:\n"
      << "    tagged_dis_gen /path/to/output/file [N events] [E_beam (GeV)]\n\n";
      return -1;
    }

  const int Nevents=atoi(argv[2]);
  const double Ebeam=atof(argv[3]);

  TFile * outfile = new TFile(argv[1],"RECREATE");
  TTree * outtree = new TTree("gen","Generated Events");

  double theta_e, phi_e, mom_e, theta_n, phi_n, mom_n, weight;
  outtree->Branch("theta_e",&theta_e,"theta_e/D");
  outtree->Branch("phi_e",&phi_e,"phi_e/D");
  outtree->Branch("mom_e",&mom_e,"mom_e/D");
  outtree->Branch("theta_n",&theta_n,"theta_n/D");
  outtree->Branch("phi_n",&phi_n,"phi_n/D");
  outtree->Branch("mom_n",&mom_n,"mom_n/D");
  outtree->Branch("weight",&weight,"weight/D");

  for (int event=0 ; event < Nevents ; event++)
    {
      theta_e=1.;
      phi_e=2.;
      mom_e=3.;
      theta_n=4.;
      phi_n=5.;
      mom_n=6.;
      weight=1.;

      outtree->Fill();
    }

  outtree->Write();
  outfile->Close();
  
  cout << "\n\tAll done generating " << Nevents << " events.\n\n";
  return 0;
}
