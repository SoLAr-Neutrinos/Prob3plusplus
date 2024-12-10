#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include <vector>
#include "TFile.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

bool kSquared = true; // are we using sin^2(x) variables?
double Theta23 = 0.5;
double DM2 = 2.5e-3;
double delta = 0;
int kNuBar = 1;

double ssth12 = 0.303;
double ssth13 = 0.021;
double dmsq12 = 7.4e-5;
double exposure = 6e5;
TString reaction = "";
TString process = "";
TString qvar = "Ecol";

int NBinsEnergy = 0;
int NBinsNadir = 0;
double *EnergyEdge = nullptr;
double *NadirEdge = nullptr;

TH2F *testHist = nullptr;

TH2F *smear = nullptr;
TGraph *xsec = nullptr;
TH1F *nadirHist = nullptr;
TGraph *spectrum = nullptr;

// Function to get the base file name without extension
TString getBaseFileName(const TString &filename)
{
    // Find the last dot in the file name
    Ssiz_t lastDot = filename.Last('.');
    // Find the last slash (both Unix and Windows separators)
    Ssiz_t lastSlash = filename.Last('/') > filename.Last('\\') ? filename.Last('/') : filename.Last('\\');

    // If there's no dot or the dot is before the last slash, there's no extension
    if (lastDot == kNPOS || lastDot < lastSlash)
    {
        return filename(lastSlash + 1, filename.Length() - lastSlash - 1);
    }

    // Extract the base file name without extension
    return filename(lastSlash + 1, lastDot - lastSlash - 1);
}

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

TH2F *applySmearingMatrix(TH2F *inputHistogram, TH2F *smearingMatrix)
{
    // Validate inputs
    if (!inputHistogram || !smearingMatrix)
    {
        printf("Error: Input histogram or smearing matrix is null.\n");
        return nullptr;
    }

    // Check compatibility between input histogram and smearing matrix
    if (inputHistogram->GetNbinsX() != smearingMatrix->GetNbinsX())
    {
        printf("Error: Input histogram x-axis bins do not match smearing matrix rows.\n");
        return nullptr;
    }

    // Prepare the output histogram H'(x', y)
    TString outputName = TString(inputHistogram->GetName()) + "_smeared";
    TH2F *smearedHistogram = new TH2F(outputName, "Smeared Histogram",
                                      smearingMatrix->GetNbinsY(), // Smeared energy bins (x')
                                      smearingMatrix->GetYaxis()->GetXmin(),
                                      smearingMatrix->GetYaxis()->GetXmax(),
                                      inputHistogram->GetNbinsY(), // Same y-axis as input
                                      inputHistogram->GetYaxis()->GetXmin(),
                                      inputHistogram->GetYaxis()->GetXmax());

    // Perform the matrix multiplication for each y-bin
    for (int yBin = 1; yBin <= inputHistogram->GetNbinsY(); ++yBin)
    {
        // Extract the input slice (true energy) for the current y-bin
        TH1D *sliceTrue = inputHistogram->ProjectionX("sliceTrue", yBin, yBin);

        // Create a temporary slice for the smeared (output) x-axis
        TH1D sliceSmeared("sliceSmeared", "Smeared X-axis",
                          smearingMatrix->GetNbinsY(),
                          smearingMatrix->GetYaxis()->GetXmin(),
                          smearingMatrix->GetYaxis()->GetXmax());

        // Perform the matrix multiplication: sum over true energy bins
        for (int xPrimeBin = 1; xPrimeBin <= smearingMatrix->GetNbinsY(); ++xPrimeBin)
        {
            double smearedValue = 0.0;
            for (int xBin = 1; xBin <= smearingMatrix->GetNbinsX(); ++xBin)
            {
                double smearingFactor = smearingMatrix->GetBinContent(xBin, xPrimeBin); // Note the transposed access
                double inputValue = sliceTrue->GetBinContent(xBin);
                smearedValue += smearingFactor * inputValue;
            }
            sliceSmeared.SetBinContent(xPrimeBin, smearedValue);
        }

        // Fill the smeared values back into the output histogram
        for (int xPrimeBin = 1; xPrimeBin <= sliceSmeared.GetNbinsX(); ++xPrimeBin)
        {
            double value = sliceSmeared.GetBinContent(xPrimeBin);
            smearedHistogram->SetBinContent(xPrimeBin, yBin, value);
        }

        // Cleanup
        delete sliceTrue;
    }

    smearedHistogram->Scale(inputHistogram->Integral() / smearedHistogram->Integral());

    return smearedHistogram;
}

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

void saveResults(const TString &fileName, const std::vector<TH2F *> &histograms)
{
    TFile *f = new TFile(fileName + ".root", "RECREATE");
    for (auto hist : histograms)
    {
        hist->Write();
    }
    f->Close();
    std::cout << "\033[1;32m" << "Saved to: " << fileName << ".root" << "\033[0m" << "\n\n";
}

std::vector<TH2F *> calculateSurvivalHistograms(double Theta12, double Theta13, double dm2)
{
    // Assuming testHist and modelHist are global or accessible within this scope
    extern int NBinsEnergy;
    extern int NBinsNadir;
    extern double *EnergyEdge;
    extern double *NadirEdge;
    extern TH1F *nadirHist;
    extern TGraph *spectrum;
    extern TGraph *xsec;
    extern TH2F *smear;
    extern double exposure;

    TH2F *hsurv = new TH2F("hsurv", "P_{#nu_{e}#rightarrow#nu_{e}}",
                           NBinsEnergy, EnergyEdge, NBinsNadir, NadirEdge);
    TH2F *hsurvSolar = new TH2F("hsurvSolar", "Solar P_{#nu_{e}#rightarrow#nu_{e}}",
                                NBinsEnergy, EnergyEdge, NBinsNadir, NadirEdge);

    BargerPropagator *myNu = new BargerPropagator();
    // flavor is default, but produced as mass states in sun
    myNu->UseMassEigenstates(true);
    // Octant for Theta23 in sin2(2x) mode
    myNu->SetDefaultOctant(23, 2);

    for (int i = 1; i <= hsurv->GetNbinsX(); i++)
    {
        double E = hsurv->GetXaxis()->GetBinCenter(i);
        double E_GeV = E * 1e3; // Convert to GeV
        for (int j = 1; j <= hsurv->GetNbinsY(); j++)
        {
            double n = hsurv->GetYaxis()->GetBinCenter(j);

            myNu->SetMNS(Theta12, Theta13, Theta23, dm2, DM2, delta,
                         E_GeV, kSquared, kNuBar);

            myNu->DefinePath(n, 25.00);
            myNu->propagate(1 * kNuBar);

            double surv_1 = myNu->GetProb(1, 1);
            double surv_2 = myNu->GetProb(2, 1);
            double surv_3 = myNu->GetProb(3, 1);

            double f_2 = ssth(E, dm2, Theta12, Theta13) * (1 - Theta13);
            double f_3 = Theta13;
            double surv = (1 - f_2 - f_3) * surv_1 + f_2 * surv_2 + f_3 * surv_3;
            hsurv->Fill(E, n, surv);

            double nadirEval = interpolateValueTH1F(n, nadirHist);

            double modelSurv = surv * nadirEval;

            if (spectrum != nullptr)
            {
                modelSurv *= spectrum->Eval(E);
            }
            if (xsec != nullptr)
            {
                modelSurv *= xsec->Eval(E);
            }
            hsurvSolar->Fill(E, n, modelSurv);
        }
    }

    hsurvSolar->Scale(exposure / hsurvSolar->Integral());

    TH2F *smeared = applySmearingMatrix(hsurvSolar, smear);

    return {hsurv, hsurvSolar, smeared};
}

double ChiSquared(const double *params)
{
    double theta12 = params[0];
    double theta13 = params[1];
    double m2 = params[2];

    // std::cout << "Theta12: " << theta12 << std::endl;
    // std::cout << "Theta13: " << theta13 << std::endl;
    // std::cout << "m2: " << m2 << std::endl;

    extern TH2F *testHist;

    std::vector<TH2F *> histograms = calculateSurvivalHistograms(theta12, theta13, m2);
    TH2F *modelHist = histograms[2];

    double chi2 = 0.0;
    for (int i = 1; i <= testHist->GetNbinsX(); i++)
    {
        for (int j = 1; j <= testHist->GetNbinsY(); j++)
        {
            double observed = testHist->GetBinContent(i, j);
            double expected = modelHist->GetBinContent(i, j);

            if (observed > 0)
            {
                chi2 += pow(observed - expected, 2) / expected;
            }
        }
    }
    return chi2;
}

std::vector<double> minimize()
{
    // Define Functor
    ROOT::Math::Functor Func(ChiSquared, 3);

    // Set up minimizer
    std::unique_ptr<ROOT::Math::Minimizer> Minimizer(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize"));

    Minimizer->SetMaxFunctionCalls(1e6);
    Minimizer->SetTolerance(0.001);
    Minimizer->SetFunction(Func);

    // Initial values and limits
    Minimizer->SetVariable(0, "Theta12", 0, 0.005);
    Minimizer->SetVariable(1, "Theta13", 0, 0.005);
    Minimizer->SetVariable(2, "m2", 3.5e-5, 1e-6);
    Minimizer->SetVariableLimits(1.0, 1.0, 9.5e-5);

    // Perform minimization
    Minimizer->Minimize();

    // Retrieve results
    double theta12_opt = Minimizer->X()[0];
    double theta13_opt = Minimizer->X()[1];
    double m2_opt = Minimizer->X()[2];
    double minChi2 = Minimizer->MinValue();

    Minimizer->PrintResults();

    // std::cout << "\nOptimized Theta12: " << theta12_opt << std::endl;
    // std::cout << "Optimized Theta13: " << theta13_opt << std::endl;
    // std::cout << "Optimized m2: " << m2_opt << std::endl;
    // std::cout << "Minimum Chi-Squared: " << minChi2 << std::endl;

    return {theta12_opt, theta13_opt, m2_opt, minChi2};
}

int main(int argc, char *argv[])
{
    if (argc > 2)
        reaction = argv[2];
    if (argc > 3)
        process = argv[3];
    if (argc > 4)
        ssth12 = (double)atof(argv[4]);
    if (argc > 5)
        ssth13 = (double)atof(argv[5]);
    if (argc > 6)
        dmsq12 = (double)atof(argv[6]);
    if (argc > 7)
        qvar = argv[7];
    // if (argc > 8)
    //     exposure = (double)atof(argv[8]);

    // Open test histogram file
    const char *inputName = argv[1];
    TFile *inputFile = TFile::Open(inputName);
    if (!inputFile || inputFile->IsZombie())
    {
        std::cerr << "Error: Cannot open file " << inputName << std::endl;
        return 1;
    }

    testHist = (TH2F *)inputFile->Get(qvar);
    if (!testHist)
    {
        std::cerr << "Error: Histogram '" << qvar << "' not found in file " << inputName << std::endl;
        inputFile->Close();
        return 1;
    }
    exposure = testHist->Integral();
    // testHist->Scale(exposure / testHist->Integral());

    NBinsEnergy = testHist->GetNbinsX();
    NBinsNadir = testHist->GetNbinsY();

    std::cout << "Nadir bins: " << NBinsNadir << std::endl;
    std::cout << "Energy bins: " << NBinsEnergy << std::endl;

    // Allocate memory for the arrays
    EnergyEdge = new double[NBinsEnergy + 1];
    NadirEdge = new double[NBinsNadir + 1];

    // Get the bin edges for the X axis (Energy)
    for (int i = 0; i <= NBinsEnergy; ++i)
    {
        EnergyEdge[i] = testHist->GetXaxis()->GetBinLowEdge(i + 1);
        // std::cout << EnergyEdge[i] << std::endl;
    }

    // Get the bin edges for the Y axis (Nadir)
    for (int i = 0; i <= NBinsNadir; ++i)
    {
        NadirEdge[i] = testHist->GetYaxis()->GetBinLowEdge(i + 1);
        // std::cout << NadirEdge[i] << std::endl;
    }

    std::cout << "Energy Edge: " << EnergyEdge[0] << " - " << EnergyEdge[NBinsEnergy] << std::endl;
    std::cout << "Nadir Edge: " << NadirEdge[0] << " - " << NadirEdge[NBinsNadir] << std::endl;

    // ----------------------------------------------------------

    // Retrieve the nadir probability distribution (TH1F)
    TFile *nadirFile = TFile::Open("nadir.root");
    if (!nadirFile || nadirFile->IsZombie())
    {
        std::cerr << "Error: Cannot open nadir.root" << std::endl;
        return 1;
    }
    nadirHist = dynamic_cast<TH1F *>(nadirFile->Get("nadir"));

    if (!nadirHist)
    {
        std::cerr << "Error: Nadir histogram not found in file." << std::endl;
        return 1;
    }
    nadirHist->Scale(1.0 / nadirHist->Integral());
    // ----------------------------------------------------------

    // Retrieve the solar neutrino spectra (TGraph)
    TFile *spectraFile = nullptr;
    if (process != "")
    {
        spectraFile = TFile::Open("solar_spectra.root");
        if (!spectraFile || spectraFile->IsZombie())
        {
            std::cerr << "Error: Cannot open solar_spectra.root" << std::endl;
            return 1;
        }
        spectrum = dynamic_cast<TGraph *>(spectraFile->Get(process));
        if (!spectrum)
        {
            std::cerr << "Error: Spectrum Graph not found in file." << std::endl;
            return 1;
        }
        spectrum->Scale(1.0 / spectrum->Integral());
    }
    // ----------------------------------------------------------

    // Retrieve the LAr cross-section (TGraph)
    TFile *xsecFile = nullptr;
    if (reaction != "")
    {
        xsecFile = TFile::Open("xsec.root");
        if (!xsecFile || xsecFile->IsZombie())
        {
            std::cerr << "Error: Cannot open xsec.root" << std::endl;
            return 1;
        }
        xsec = dynamic_cast<TGraph *>(xsecFile->Get(reaction));
        if (!xsec)
        {
            std::cerr << "Error: Cross-section Graph not found in file." << std::endl;
            return 1;
        }
        xsec->Scale(1.0 / xsec->Integral());
    }
    // ----------------------------------------------------------

    // Retrieve the smearing matrix (TH2F)
    TFile *smearFile = TFile::Open("smear.root");
    if (!smearFile || smearFile->IsZombie())
    {
        std::cerr << "Error: Cannot open smear.root" << std::endl;
        return 1;
    }
    smear = dynamic_cast<TH2F *>(smearFile->Get("E->" + qvar));
    if (!smear)
    {
        std::cerr << "Error: " << "E->" + qvar << " Smearing matrix not found in file." << std::endl;
        return 1;
    }
    smear->Scale(1.0 / smear->Integral());
    // ----------------------------------------------------------

    gSystem->Load("libThreeProb_3.10.a");
    TH1::AddDirectory(0);

    std::cout << "\nStarting minimization..." << std::endl;

    std::vector<double> optimum = minimize();
    ssth12 = optimum[0];
    ssth13 = optimum[1];
    dmsq12 = optimum[2];

    std::cout << "\nUsing:" << "\n"
              << "\tdm2:\t " << dmsq12 << "\n"
              << "\tTheta13: " << ssth13 << "\n"
              << "\tTheta12: " << ssth12 << "\n"
              << "\tTheta23: " << Theta23 << "\n"
              << "\tDM2:\t " << DM2 << "\n"
              << "\tdelta:   " << delta << "\n"
              << "\tknubar:  " << kNuBar << "\n";

    std::vector<TH2F *> histograms = calculateSurvivalHistograms(ssth12, ssth13, dmsq12);

    // Save the name of the file in a string with dm2 and sin12 as floats in scientific notation
    TString FileName = "./models/" + getBaseFileName(inputName) + "_" + qvar + "_dm2_" + TString::Format("%.3e", dmsq12) + "_sin13_" + TString::Format("%.3e", ssth13) + "_sin12_" + TString::Format("%.3e", ssth12);
    if (process != "")
    {
        FileName += "_" + process;
    }
    if (reaction != "")
    {
        FileName += "_" + reaction;
    }

    saveResults(FileName, histograms);

    if (smearFile != nullptr)
        smearFile->Close();
    if (spectraFile != nullptr)
        spectraFile->Close();

    nadirFile->Close();
    inputFile->Close();

    // Free the allocated memory
    delete[] EnergyEdge;
    delete[] NadirEdge;

    return 0;
}