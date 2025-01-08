#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH2D.h"
#include "TSystem.h"

double ssth(double Enu, double dmsq_vac, double ssth12_vac, double ssth13)
{
  double rhoY = 0.090; // rho density of production in kg/cm^3, Y = N(e)/N(n+p)
  double A = 1.53e-4 * rhoY * Enu;
  A *= (1 - ssth13);
  double th12 = asin(sqrt(ssth12_vac));
  double s2th_vac = sin(2 * th12);
  double ret = 0.5 * (1 + (A - dmsq_vac * sqrt(1 - s2th_vac * s2th_vac)) / sqrt(pow(dmsq_vac * sqrt(1 - s2th_vac * s2th_vac) - A, 2.) + pow(dmsq_vac * s2th_vac, 2.)));
  return ret;
}

// int main(double ssth12, double ssth13, double dmsq12)
int main(int argc, char *argv[])
{
  double ssth12 = 0;
  double ssth13 = 0;
  double dmsq12 = 0;

  if (argc >= 1)
    ssth12 = (double)atof(argv[1]);
  if (argc >= 2)
    ssth13 = (double)atof(argv[2]);
  if (argc >= 3)
    dmsq12 = (double)atof(argv[3]);

  //gSystem->Load("libThreeProb_3.10.a");
  gSystem->Load("install/lib/libProb3plusplus.so");
  TH1::AddDirectory(0);

  int NBinsEnergy = 300;
  int NBinsNadir = 1001;

  double EnergyEdge[NBinsEnergy + 1];
  for (int i = 0; i <= NBinsEnergy; i++)
    EnergyEdge[i] = 0.03 * double(i) / NBinsEnergy;

  double NadirEdge[NBinsNadir + 1];
  for (int i = 0; i <= NBinsNadir - 1; i++)
    NadirEdge[i] = -1 + double(i) / (NBinsNadir - 1);
  NadirEdge[NBinsNadir] = +1;

  bool kSquared = true; // are we using sin^2(x) variables?

  double dm2 = dmsq12;
  double Theta13 = ssth13; // bf = 0.0214
  double Theta12 = ssth12;
  double Theta23 = 0.5;
  double DM2 = 2.5e-3;
  double delta = 0;
  int kNuBar = 1;

  std::cout << "Using:" << "\n"
            << "\tdm2:\t " << dm2 << "\n"
            << "\tTheta13: " << Theta13 << "\n"
            << "\tTheta12: " << Theta12 << "\n"
            << "\tTheta23: " << Theta23 << "\n"
            << "\tDM2:\t " << DM2 << "\n"
            << "\tdelta:   " << delta << "\n"
            << "\tknubar:  " << kNuBar << "\n";

  NeutrinoPropagator *myNu;
  BargerPropagator *bNu;

  bNu = new BargerPropagator();
  // flavor is default, but produced as mass states in sun
  bNu->UseMassEigenstates(true);

  // Octant for Theta23 in sin2(2x) mode
  bNu->SetDefaultOctant(23, 2);
  // use the standard barger
  myNu = bNu;

  // Save the name of the file in a string with dm2 and sin12 as floats in scientific notation
  TString FilePath = "./oscilograms/";
  TString FileName = "osc_probability_dm2_" + TString::Format("%.3e", dm2) + "_sin13_" + TString::Format("%.3e", ssth13) + "_sin12_" + TString::Format("%.3e", ssth12);

  TH2D *hsurv = new TH2D("hsurv", "P_{#nu_{e}#rightarrow#nu_{e}}",
                         NBinsEnergy, EnergyEdge, NBinsNadir, NadirEdge);
  TH2D *hsurv_1 = new TH2D("hsurv_1", "P_{#nu_{1}#rightarrow#nu_{e}}",
                           NBinsEnergy, EnergyEdge, NBinsNadir, NadirEdge);
  TH2D *hsurv_2 = new TH2D("hsurv_2", "P_{#nu_{2}#rightarrow#nu_{e}}",
                           NBinsEnergy, EnergyEdge, NBinsNadir, NadirEdge);
  TH2D *hsurv_3 = new TH2D("hsurv_3", "P_{#nu_{2}#rightarrow#nu_{e}}",
                           NBinsEnergy, EnergyEdge, NBinsNadir, NadirEdge);
  myNu->SetMNS(Theta12, Theta13, Theta23, dm2, DM2, delta,
      1.0, kSquared, kNuBar);
  for (int i = 1; i <= hsurv->GetNbinsX(); i++)
  {
    double E = hsurv->GetXaxis()->GetBinCenter(i);
    for (int j = 1; j <= hsurv->GetNbinsY(); j++)
    {
      double n = hsurv->GetYaxis()->GetBinCenter(j);

      //myNu->SetMNS(Theta12, Theta13, Theta23, dm2, DM2, delta,
                   //E, kSquared, kNuBar);
      myNu->SetEnergy( E ); 

      myNu->DefinePath(n, 1.47e8);
      myNu->propagate(1 * kNuBar);

      double surv_1 = myNu->GetProb(1, 1);
      double surv_2 = myNu->GetProb(2, 1);
      double surv_3 = myNu->GetProb(3, 1);

      // double f_3 = pow(asin(sqrt(Theta13)),2.);
      //double f_2 = ssth(E * 1e3, dm2, Theta12, Theta13) * (1 - Theta13);
      double f_2 = myNu->GetSinSqTheta12Sun();
      double f_3 = Theta13;
      double surv = (1 - f_2 - f_3) * surv_1 + f_2 * surv_2 + f_3 * surv_3;

      //   if (n>0) std::cout << E << "  " << 1-f_2-f_3 << "  "<< f_2 <<"   " << f_3 << "  " << surv << std::endl;

      hsurv_1->Fill(E, n, surv_1);
      hsurv_2->Fill(E, n, surv_2);
      hsurv_3->Fill(E, n, surv_3);

      hsurv->Fill(E, n, surv);
    }
  }
  // gStyle->SetNdivisions(506,"XYZ");
  // hsurv->Draw("colz");
  // hsurv->ProjectionY("a",150,150)->Draw("hist");
  // TFile *out = new TFile("out.root","recreate");
  // hsurv  ->Write("hsurv_e");
  // hsurv_1->Write("hsurv_1");
  // hsurv_2->Write("hsurv_2");
  // hsurv_3->Write("hsurv_3");

  TFile *f = new TFile(FilePath + FileName + ".root", "RECREATE");
  hsurv->Write();
  f->Close();
  // Make a colorful print to the terminal including the name of the file
  std::cout << "\033[1;32m" << "Saved to: " << FilePath << "\nFile: " << FileName << ".root" << "\033[0m" << "\n\n";
  return 0;
}

// Write a reference list for colors in the terminal with full format
// Color  |   Code
// Black  |   0;30
// Blue   |   0;34
// Green  |   0;32
// Cyan   |   0;36
// Red    |   0;31
// Purple |   0;35
// Brown  |   0;33
