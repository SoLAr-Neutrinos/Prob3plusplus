#include <iostream>
#include <cmath>
#include <vector>
#include "BargerPropagator.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

// Parameters for neutrino propagation
double Theta23 = 0.5, DM2 = 2.5e-3, delta = 0;
bool kSquared = true;
int kNuBar = 1;

// Function to interpolate the value from a TH1F
double interpolateValueTH1F(double n, TH1F *hist)
{
    int bin = hist->GetXaxis()->FindBin(n);

    // Handle edges: clamp to valid bin range
    if (bin <= 0)
        return hist->GetBinContent(1);
    if (bin >= hist->GetNbinsX())
        return hist->GetBinContent(hist->GetNbinsX());

    // Linear interpolation
    double x1 = hist->GetXaxis()->GetBinCenter(bin);
    double y1 = hist->GetBinContent(bin);
    double x2 = hist->GetXaxis()->GetBinCenter(bin + 1);
    double y2 = hist->GetBinContent(bin + 1);

    if (x2 == x1)
        return y1; // Avoid division by zero

    return y1 + (y2 - y1) * (n - x1) / (x2 - x1);
}

double ssth(double Enu, double dmsq_vac, double ssth12_vac, double ssth13)
{
    double rhoY = 0.090; // density of production in kg/cm^3
    double A = 1.53e-4 * rhoY * Enu;
    A *= (1 - ssth13);
    double th12 = asin(sqrt(ssth12_vac));
    double s2th_vac = sin(2 * th12);
    double ret = 0.5 * (1 + (A - dmsq_vac * sqrt(1 - s2th_vac * s2th_vac)) /
                                sqrt(pow(dmsq_vac * sqrt(1 - s2th_vac * s2th_vac) - A, 2.) +
                                     pow(dmsq_vac * s2th_vac, 2.)));
    return ret;
}

int main(int argc, char **argv)
{
    // Check for correct number of arguments
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <test_histogram.root> <histogram name>" << std::endl;
        return 1;
    }

    gSystem->Load("libThreeProb_3.10.a");

    // Open test histogram file
    const char *filename = argv[1];
    const char *histname = argv[2];
    TFile *inputFile = TFile::Open(filename);
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return 1;
    }

    TH2D *testHist = (TH2D *)inputFile->Get(histname);
    if (!testHist)
    {
        std::cerr << "Error: Histogram '" << histname << "' not found in file " << filename << std::endl;
        inputFile->Close();
        return 1;
    }

    // Retrieve the nadir scatter plot (TGraph2D)
    TFile *nadirFile = TFile::Open("nadir.root");
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error: Cannot open nadir.root" << std::endl;
        return 1;
    }

    TH1F *nadirHist = dynamic_cast<TH1F *>(nadirFile->Get("nadir"));
    if (!nadirHist)
    {
        std::cerr << "Error: Nadir scatter plot not found in file." << std::endl;
        return 1;
    }

    BargerPropagator *myNu = new BargerPropagator();
    myNu->UseMassEigenstates(true);

    // Chi-squared calculation function
    auto ChiSquared = [&](const double *coeff)
    {
        double Theta12 = coeff[0]; // Parameters to vary
        double Theta13 = coeff[1];
        double dm2 = coeff[2];
        double chi2 = 0.0;

        for (int i = 1; i <= testHist->GetNbinsX(); i++)
        {
            double E = testHist->GetXaxis()->GetBinCenter(i);
            for (int j = 1; j <= testHist->GetNbinsY(); j++)
            {
                double n = testHist->GetYaxis()->GetBinCenter(j);
                double observed = testHist->GetBinContent(i, j);

                myNu->SetMNS(Theta12, Theta13, Theta23, dm2, DM2, delta, E, kSquared, kNuBar);
                myNu->DefinePath(n, 25.00);
                myNu->propagate(kNuBar);

                double surv_1 = myNu->GetProb(1, 1);
                double surv_2 = myNu->GetProb(2, 1);
                double surv_3 = myNu->GetProb(3, 1);

                double f_2 = ssth(E * 1e3, DM2, Theta12, Theta13) * (1 - Theta13);
                double f_3 = Theta13;

                double nadirConv = interpolateValueTH1F(n, nadirHist);

                double expected = nadirConv * (1 - f_2 - f_3) * surv_1 + f_2 * surv_2 + f_3 * surv_3;

                if (observed > 0)
                {
                    chi2 += pow(observed - expected, 2) / observed;
                }
            }
        }

        return chi2;
    };

    // Define Functor
    ROOT::Math::Functor Func(ChiSquared, 2);

    // Set up minimizer
    std::unique_ptr<ROOT::Math::Minimizer> Minimizer(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

    Minimizer->SetMaxFunctionCalls(5000);
    Minimizer->SetTolerance(0.001);
    Minimizer->SetFunction(Func);

    // Initial values and limits
    Minimizer->SetVariable(0, "Theta12", 0, 0.005);
    Minimizer->SetVariable(1, "Theta13", 0, 0.005);
    Minimizer->SetVariable(2, "m2", 3.5e-5, 1e-6);
    Minimizer->SetVariableLimits(1, 1, 9.5e-5);

    // Perform minimization
    Minimizer->Minimize();

    // Retrieve results
    double theta12_opt = Minimizer->X()[0];
    double theta13_opt = Minimizer->X()[1];
    double m2_opt = Minimizer->X()[2];
    double minChi2 = Minimizer->MinValue();

    if (Minimizer->Status())
        std::cout << "Fit failed" << std::endl;

    std::cout << "Optimized Theta12: " << theta12_opt << std::endl;
    std::cout << "Optimized Theta13: " << theta13_opt << std::endl;
    std::cout << "Optimized m2: " << m2_opt << std::endl;
    std::cout << "Minimum Chi-Squared: " << minChi2 << std::endl;

    // Clean up
    delete myNu;
    inputFile->Close();
    return 0;
}
