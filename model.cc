#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include <vector>
#include <set>
#include "TFile.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"
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
TString qvar = "Ecal";
TString inputName = "";

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
            for (int xBin = 1; xBin <= sliceTrue->GetNbinsX(); ++xBin)
            {
                double inputBinCenter = sliceTrue->GetXaxis()->GetBinCenter(xBin);
                int smearingBin = smearingMatrix->GetXaxis()->FindBin(inputBinCenter);
                if (smearingBin < 1 || smearingBin > smearingMatrix->GetNbinsX())
                    continue; // Skip if the bin is out of bounds
                if (smearingMatrix->GetYaxis()->GetBinCenter(xPrimeBin) < 0)
                    continue;

                double smearingFactor = smearingMatrix->GetBinContent(smearingBin, xPrimeBin); // Note the transposed access
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
    TH2F *hsolar = new TH2F("hsolar", "Solar P_{#nu_{e}#rightarrow#nu_{e}}",
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
            hsolar->Fill(E, n, modelSurv);
        }
    }

    // hsolar->Scale(exposure / hsolar->Integral());

    TH2F *smeared = applySmearingMatrix(hsolar, smear);

    smeared->Scale(exposure / smeared->Integral());

    return {hsurv, hsolar, smeared};
}

TH2F *histChi2(TH2F *hist1, TH2F *hist2)
{
    // Determine the common binning range and step size
    double xMin = std::min(hist1->GetXaxis()->GetXmin(), hist2->GetXaxis()->GetXmin());
    double xMax = std::max(hist1->GetXaxis()->GetXmax(), hist2->GetXaxis()->GetXmax());
    double yMin = std::min(hist1->GetYaxis()->GetXmin(), hist2->GetYaxis()->GetXmin());
    double yMax = std::max(hist1->GetYaxis()->GetXmax(), hist2->GetYaxis()->GetXmax());

    int nBinsX = max(hist1->GetNbinsX(), hist2->GetNbinsX());
    int nBinsY = max(hist1->GetNbinsY(), hist2->GetNbinsY());

    double xStep = (xMax - xMin) / nBinsX;
    double yStep = (yMax - yMin) / nBinsY;

    // Create a new histogram with the common binning
    TH2F *result = new TH2F("chi2", "Chi2", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    // Subtract the bin contents of the two histograms
    for (int i = 1; i <= nBinsX; ++i)
    {
        for (int j = 1; j <= nBinsY; ++j)
        {
            double xCenter = xMin + (i - 0.5) * xStep;
            double yCenter = yMin + (j - 0.5) * yStep;

            int bin1 = hist1->FindBin(xCenter, yCenter);
            int bin2 = hist2->FindBin(xCenter, yCenter);

            double content1 = (hist1->GetXaxis()->GetBinLowEdge(i) >= hist1->GetXaxis()->GetXmin() &&
                               hist1->GetXaxis()->GetBinUpEdge(i) <= hist1->GetXaxis()->GetXmax() &&
                               hist1->GetYaxis()->GetBinLowEdge(j) >= hist1->GetYaxis()->GetXmin() &&
                               hist1->GetYaxis()->GetBinUpEdge(j) <= hist1->GetYaxis()->GetXmax())
                                  ? hist1->GetBinContent(bin1)
                                  : 0.0;

            double content2 = (hist2->GetXaxis()->GetBinLowEdge(i) >= hist2->GetXaxis()->GetXmin() &&
                               hist2->GetXaxis()->GetBinUpEdge(i) <= hist2->GetXaxis()->GetXmax() &&
                               hist2->GetYaxis()->GetBinLowEdge(j) >= hist2->GetYaxis()->GetXmin() &&
                               hist2->GetYaxis()->GetBinUpEdge(j) <= hist2->GetYaxis()->GetXmax())
                                  ? hist2->GetBinContent(bin2)
                                  : 0.0;
            if (content2 != 0)
            {
                result->SetBinContent(i, j, pow((content1 - content2), 2) / content2);
            }
            else
            {
                // Handle the case where content2 is zero, e.g., set the bin content to zero
                result->SetBinContent(i, j, 0);
            }
        }
    }

    return result;
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

    TH2F *resultHist = histChi2(testHist, modelHist);
    double chi2 = resultHist->Integral();
    // std::cout << "theta12: " << theta12 << " theta13: " << theta13 << " m2: " << m2 << " Chi2: " << chi2 << std::endl;
    return chi2;
}

TH2F *histLogLikelihood(TH2F *hist1, TH2F *hist2)
{
    // Determine the common binning range and step size
    double xMin = std::min(hist1->GetXaxis()->GetXmin(), hist2->GetXaxis()->GetXmin());
    double xMax = std::max(hist1->GetXaxis()->GetXmax(), hist2->GetXaxis()->GetXmax());
    double yMin = std::min(hist1->GetYaxis()->GetXmin(), hist2->GetYaxis()->GetXmin());
    double yMax = std::max(hist1->GetYaxis()->GetXmax(), hist2->GetYaxis()->GetXmax());

    int nBinsX = max(hist1->GetNbinsX(), hist2->GetNbinsX());
    int nBinsY = max(hist1->GetNbinsY(), hist2->GetNbinsY());

    double xStep = (xMax - xMin) / nBinsX;
    double yStep = (yMax - yMin) / nBinsY;

    // Create a new histogram with the common binning
    TH2F *result = new TH2F("log-likelihood", "-ln(L)", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    // Subtract the bin contents of the two histograms
    for (int i = 1; i <= nBinsX; ++i)
    {
        for (int j = 1; j <= nBinsY; ++j)
        {
            double xCenter = xMin + (i - 0.5) * xStep;
            double yCenter = yMin + (j - 0.5) * yStep;

            int bin1 = hist1->FindBin(xCenter, yCenter);
            int bin2 = hist2->FindBin(xCenter, yCenter);

            double observed = (hist1->GetXaxis()->GetBinLowEdge(i) >= hist1->GetXaxis()->GetXmin() &&
                               hist1->GetXaxis()->GetBinUpEdge(i) <= hist1->GetXaxis()->GetXmax() &&
                               hist1->GetYaxis()->GetBinLowEdge(j) >= hist1->GetYaxis()->GetXmin() &&
                               hist1->GetYaxis()->GetBinUpEdge(j) <= hist1->GetYaxis()->GetXmax())
                                  ? hist1->GetBinContent(bin1)
                                  : 0.0;

            double expected = (hist2->GetXaxis()->GetBinLowEdge(i) >= hist2->GetXaxis()->GetXmin() &&
                               hist2->GetXaxis()->GetBinUpEdge(i) <= hist2->GetXaxis()->GetXmax() &&
                               hist2->GetYaxis()->GetBinLowEdge(j) >= hist2->GetYaxis()->GetXmin() &&
                               hist2->GetYaxis()->GetBinUpEdge(j) <= hist2->GetYaxis()->GetXmax())
                                  ? hist2->GetBinContent(bin2)
                                  : 0.0;

            double logL = 0.0;
            if (expected > 0)
            {
                logL = -observed * log(expected) + expected; // Poisson likelihood
            }
            else if (observed > 0)
            {
                logL = std::numeric_limits<double>::infinity(); // Handle invalid cases
            }
            // std::cout << likelihood << std::endl;
            result->SetBinContent(i, j, logL);
        }
    }

    return result;
}

double logLikelihood(const double *params)
{
    double theta12 = (params[0]);
    double theta13 = (params[1]);
    double m2 = (params[2]);

    // std::cout << "Theta12: " << theta12 << std::endl;
    // std::cout << "Theta13: " << theta13 << std::endl;
    // std::cout << "m2: " << m2 << std::endl;

    extern TH2F *testHist;

    std::vector<TH2F *> histograms = calculateSurvivalHistograms(theta12, theta13, m2);
    TH2F *modelHist = histograms[2];

    TH2F *resultHist = histLogLikelihood(testHist, modelHist);
    double l = resultHist->Integral();
    // std::cout << "theta12: " << theta12 << " theta13: " << theta13 << " m2: " << m2 << " L: " << l << std::endl;

    return l;
}

std::vector<std::vector<double>> minimize()
{

    TString fileName = "models/" + getBaseFileName(inputName) + "_" + qvar + "_minimizer_plots.root";
    TFile *rootFile = new TFile(fileName, "RECREATE");

    // Helper lambda for minimization
    auto performMinimization = [](bool fixParam0, bool fixParam1, bool fixParam2)
    {
        ROOT::Math::Functor Func(logLikelihood, 3);
        std::unique_ptr<ROOT::Math::Minimizer> Minimizer(
            ROOT::Math::Factory::CreateMinimizer("Minuit", "Simplex"));

        Minimizer->SetMaxFunctionCalls(1e6);
        Minimizer->SetMaxIterations(1e6);
        Minimizer->SetTolerance(1e-5);
        Minimizer->SetPrecision(1e-12);
        Minimizer->SetPrintLevel(1);
        Minimizer->SetStrategy(2);
        Minimizer->SetFunction(Func);

        // Set up fixed or variable parameters
        if (fixParam0)
        {
            Minimizer->SetFixedVariable(0, "Theta12", ssth12);
        }
        else
        {
            Minimizer->SetVariable(0, "Theta12", (ssth12), 1e-3);
            Minimizer->SetVariableLimits(0, (0.001), (1));
        }

        if (fixParam1)
        {
            Minimizer->SetFixedVariable(1, "Theta13", ssth13);
        }
        else
        {
            Minimizer->SetVariable(1, "Theta13", (ssth13), 1e-3);
            Minimizer->SetVariableLimits(1, (0.001), (1));
        }

        if (fixParam2)
        {
            Minimizer->SetFixedVariable(2, "m2", dmsq12);
        }
        else
        {
            Minimizer->SetVariable(2, "m2", (dmsq12), 1e-6);
            Minimizer->SetVariableLimits(2, (3.5e-5), (9.5e-5));
        }

        Minimizer->Minimize();

        // Retrieve optimized parameters
        double theta12_opt = Minimizer->X()[0];
        double theta13_opt = Minimizer->X()[1];
        double m2_opt = Minimizer->X()[2];
        double maxLikelihood = -Minimizer->MinValue();

        // Print optimized parameters
        std::cout << "\nOptimized Parameters:\n"
                  << "Theta12: " << theta12_opt << "\n"
                  << "Theta13: " << theta13_opt << "\n"
                  << "m2: " << m2_opt << "\n"
                  << "Max Likelihood: " << maxLikelihood << "\n\n";

        // Create unique labels based on fixed parameters
        std::string fixedLabel = "";
        if (fixParam0)
            fixedLabel += "_Theta12_fixed";
        if (fixParam1)
            fixedLabel += "_Theta13_fixed";
        if (fixParam2)
            fixedLabel += "_m2_fixed";

        // Covariance and correlation matrices
        std::string covMatrixName = "CovarianceMatrix" + fixedLabel;
        std::string corrMatrixName = "CorrelationMatrix" + fixedLabel;
        TH2F *covMatrixHist = new TH2F(covMatrixName.c_str(), (covMatrixName + ";Param1;Param2").c_str(), 3, 0, 3, 3, 0, 3);
        TH2F *corrMatrixHist = new TH2F(corrMatrixName.c_str(), (corrMatrixName + ";Param1;Param2").c_str(), 3, 0, 3, 3, 0, 3);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                covMatrixHist->SetBinContent(i + 1, j + 1, Minimizer->CovMatrix(i, j));
                corrMatrixHist->SetBinContent(i + 1, j + 1, Minimizer->Correlation(i, j));
            }
        }

        covMatrixHist->Write();
        corrMatrixHist->Write();

        // Contour plots
        auto generateContour = [&](int var1, int var2, const char *var1Name, const char *var2Name)
        {
            unsigned int npoints = 40;
            std::vector<double> xi(npoints), xj(npoints);
            if (Minimizer->Contour(var1, var2, npoints, xi.data(), xj.data()))
            {
                std::string contourName = "Contour_" + std::string(var1Name) + "_" + std::string(var2Name) + fixedLabel;
                std::string contourTitle = "Contour: " + std::string(var1Name) + " vs " + std::string(var2Name) + ";" + std::string(var1Name) + ";" + std::string(var2Name);
                auto graph = new TGraph(npoints, xi.data(), xj.data());
                graph->SetName(contourName.c_str());
                graph->SetTitle(contourTitle.c_str());
                graph->Write();
            }
        };

        generateContour(0, 1, "Theta12", "Theta13");
        generateContour(0, 2, "Theta12", "m2");
        generateContour(1, 2, "Theta13", "m2");

        return std::vector<double>{theta12_opt, theta13_opt, m2_opt, maxLikelihood};
    };

    // Perform four minimizations
    std::vector<std::vector<double>> results;
    results.push_back(performMinimization(false, false, false));
    results.push_back(performMinimization(true, false, false));
    results.push_back(performMinimization(false, true, false));
    results.push_back(performMinimization(false, false, true));

    rootFile->Close();

    return results;
}

int main(int argc, char *argv[])
{
    if (argc > 2)
        qvar = argv[2];
    if (argc > 3)
        reaction = argv[3];
    if (argc > 4)
        process = argv[4];
    if (argc > 5)
        ssth12 = (double)atof(argv[5]);
    if (argc > 6)
        ssth13 = (double)atof(argv[6]);
    if (argc > 7)
        dmsq12 = (double)atof(argv[7]);
    // if (argc > 8)
    //     exposure = (double)atof(argv[8]);

    // Open test histogram file
    inputName = argv[1];
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

    // std::cout << "Nadir bins: " << NBinsNadir << std::endl;
    // std::cout << "Energy bins: " << NBinsEnergy << std::endl;

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

    // std::cout << "Energy Edge: " << EnergyEdge[0] << " - " << EnergyEdge[NBinsEnergy] << std::endl;
    // std::cout << "Nadir Edge: " << NadirEdge[0] << " - " << NadirEdge[NBinsNadir] << std::endl;

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
    // smear->Scale(1.0 / smear->Integral());
    // ----------------------------------------------------------

    gSystem->Load("libThreeProb_3.10.a");
    TH1::AddDirectory(0);

    std::cout << "\nStarting minimization..." << std::endl;

    // Obtain all optima
    std::vector<std::vector<double>> allOptima = minimize();

    // Set to track processed optima
    std::set<std::string> processedOptima;

    // Loop over all optima
    for (const auto &optimum : allOptima)
    {
        // Generate a unique key for the current optimum
        std::string key = std::to_string(optimum[0]) + "_" + std::to_string(optimum[1]) + "_" + std::to_string(optimum[2]);

        // Check if this optimum has already been processed
        if (processedOptima.find(key) != processedOptima.end())
        {
            continue; // Skip if already processed
        }

        // Mark this optimum as processed
        processedOptima.insert(key);

        // Extract parameters
        double ssth12 = optimum[0];
        double ssth13 = optimum[1];
        double dmsq12 = optimum[2];
        double likelihood = optimum[3];

        // Log the parameters
        std::cout << "\nUsing:" << "\n"
                  << "\tdm2:\t " << dmsq12 << "\n"
                  << "\tTheta13: " << ssth13 << "\n"
                  << "\tTheta12: " << ssth12 << "\n"
                  << "\tLog-likelihood: " << likelihood << "\n\n";

        // Perform calculations based on the optimum
        std::vector<TH2F *> histograms = calculateSurvivalHistograms(ssth12, ssth13, dmsq12);

        TH2F *chi2Hist = histChi2(testHist, histograms[2]);
        TH2F *likelihoodHist = histLogLikelihood(testHist, histograms[2]);

        likelihoodHist->Scale(-1); // Convert to ln(L)

        if (chi2Hist)
        {
            histograms.push_back(chi2Hist);
        }
        if (likelihoodHist)
        {
            histograms.push_back(likelihoodHist);
        }

        // Save the name of the file with parameters in scientific notation
        TString FileName = "./models/" + getBaseFileName(inputName) + "_" + qvar + "_dm2_" + TString::Format("%.3e", dmsq12) + "_sin13_" + TString::Format("%.3e", ssth13) + "_sin12_" + TString::Format("%.3e", ssth12);
        if (process != "")
        {
            FileName += "_" + process;
        }
        if (reaction != "")
        {
            FileName += "_" + reaction;
        }

        // Save the results
        saveResults(FileName, histograms);
    }

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