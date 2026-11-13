#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TList.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TPaveText.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

struct PeakFitResult {
    double peakPosition = 0;
    double fitMin = 0, fitMax = 0;
    double initialSigma = 0;
    double rangeMean = 0, rangeSigma = 0; // Pre-fit parameters: the final range is rangeMean +/- 2*rangeSigma.
    double amplitude = std::numeric_limits<double>::quiet_NaN();
    double mean = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double amplitudeError = 0;
    double meanError = 0;        // Mean uncertainty from TMinuit/HESSE, distinct from the peak width sigma.
    double sigmaError = 0;       // Fit uncertainty of the peak width sigma.
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int status = -1;             // MIGRAD: 0 means convergence; -1 means the fit has not been performed.
    int covarianceStatus = 0;    // TMinuit: 3 means a full, accurate covariance matrix.
};

// Pass the current histogram and bin range through the fixed TMinuit FCN interface.
namespace {
TH1 *fitHistogram = nullptr;
int fitFirstBin = 0, fitLastBin = 0;

void GaussianChi2(int &, double *, double &chi2, double *parameters, int) {
    chi2 = 0;
    for (int bin = fitFirstBin; bin <= fitLastBin; ++bin) {
        const double x = fitHistogram->GetBinCenter(bin);
        const double y = fitHistogram->GetBinContent(bin);
        const double error = fitHistogram->GetBinError(bin);
        const double z = (x - parameters[1]) / parameters[2];
        const double model = parameters[0] * std::exp(-0.5 * z * z);
        const double variance = error > 0 ? error * error : 1.0;
        chi2 += (y - model) * (y - model) / variance;
    }
}
}

// Find peaks, estimate mu and sigma in the half-height range, then fit within mu +/- 2*sigma.
// Return only converged fits with reliable errors, sorted by increasing fitted mean.
std::vector<PeakFitResult> FindAndFitPeaks(TH1 &histogram) {
    const int maxPeaks = histogram.GetNbinsX();
    const double sigma = 2.0;       // Peak-search width parameter in bins.
    const double searchThreshold = 1.0e-6; // Loose candidate-search threshold; do not impose a 5% peak-height cut.
    const double minimumPeakCounts = 100; // Require the raw bin count at each candidate position to be strictly greater than 100.
    std::vector<double> positions;
    std::vector<PeakFitResult> results;

    if (histogram.GetDimension() == 1 && histogram.GetMaximum() > 0) {
        TSpectrum spectrum(maxPeaks);
        const int count = spectrum.Search(&histogram, sigma, "goff", searchThreshold);
        for (int i = 0; i < count; ++i) {
            const double position = spectrum.GetPositionX()[i];
            if (histogram.GetBinContent(histogram.FindBin(position)) > minimumPeakCounts) {
                positions.push_back(position);
            }
        }
        std::sort(positions.begin(), positions.end());
    }

    for (double position : positions) {
        PeakFitResult result;
        result.peakPosition = position;
        int peakBin = histogram.FindBin(position);
        // TSpectrum candidates may lie a few bins away from the raw local maximum.
        while (peakBin > 1 && histogram.GetBinContent(peakBin - 1) > histogram.GetBinContent(peakBin)) --peakBin;
        while (peakBin < histogram.GetNbinsX() && histogram.GetBinContent(peakBin + 1) > histogram.GetBinContent(peakBin)) ++peakBin;
        const double height = histogram.GetBinContent(peakBin);
        const double halfHeight = 0.5 * height;
        int left = peakBin, right = peakBin;
        while (left > 1 && histogram.GetBinContent(left) > halfHeight) --left;
        while (right < histogram.GetNbinsX() && histogram.GetBinContent(right) > halfHeight) ++right;
        if (histogram.GetBinContent(left) > halfHeight || histogram.GetBinContent(right) > halfHeight) {
            // Skip boundary peaks without both half-height crossings.
            continue;
        }

        // Interpolate linearly between adjacent bin centers to locate half-height crossings.
        const auto crossing = [&](int a, int b) {
            const double ya = histogram.GetBinContent(a), yb = histogram.GetBinContent(b);
            return histogram.GetBinCenter(a) + (halfHeight - ya) / (yb - ya)
                * (histogram.GetBinCenter(b) - histogram.GetBinCenter(a));
        };
        result.fitMin = crossing(left, left + 1);
        result.fitMax = crossing(right - 1, right);
        fitFirstBin = left + 1;
        fitLastBin = right - 1;
        result.ndf = fitLastBin - fitFirstBin + 1 - 3;
        if (result.ndf <= 0) {
            continue;
        }

        // Use the count-weighted sample standard deviation in the fit range as the initial sigma.
        // Each bin count is the number of samples at that bin center.
        double count = 0, sumX = 0, sumX2 = 0;
        for (int bin = fitFirstBin; bin <= fitLastBin; ++bin) {
            const double x = histogram.GetBinCenter(bin), n = histogram.GetBinContent(bin);
            count += n;
            sumX += n * x;
            sumX2 += n * x * x;
        }
        if (count <= 1) {
            continue;
        }
        result.initialSigma = std::sqrt(std::max(0.0, (sumX2 - sumX * sumX / count) / (count - 1)));
        const double binWidth = histogram.GetBinWidth(peakBin);
        result.initialSigma = std::max(result.initialSigma, binWidth * 0.1);

        result.amplitude = height;
        result.mean = position;
        result.sigma = result.initialSigma;
        bool fitSucceeded = false;
        for (int pass = 0; pass < 2; ++pass) {
            if (pass == 1) {
                result.rangeMean = result.mean;
                result.rangeSigma = result.sigma;
                result.fitMin = std::max(histogram.GetXaxis()->GetXmin(), result.rangeMean - 2 * result.rangeSigma);
                result.fitMax = std::min(histogram.GetXaxis()->GetXmax(), result.rangeMean + 2 * result.rangeSigma);
                fitFirstBin = std::max(1, histogram.FindFixBin(result.fitMin));
                fitLastBin = std::min(histogram.GetNbinsX(), histogram.FindFixBin(result.fitMax));
                if (histogram.GetBinCenter(fitFirstBin) < result.fitMin) ++fitFirstBin;
                if (histogram.GetBinCenter(fitLastBin) > result.fitMax) --fitLastBin;
                result.ndf = fitLastBin - fitFirstBin + 1 - 3;
                if (result.ndf <= 0) break;
            }
            fitHistogram = &histogram;
            TMinuit minuit(3);
            minuit.SetPrintLevel(-1);
            minuit.SetFCN(GaussianChi2);
            minuit.SetErrorDef(1.0); // chi-square: Delta chi2 = 1.
            // Equal lower/upper limits (0, 0) leave each parameter unbounded in TMinuit.
            minuit.DefineParameter(0, "amplitude", result.amplitude, std::max(1.0, height * 0.01), 0, 0);
            minuit.DefineParameter(1, "mean", result.mean, binWidth * 0.1, 0, 0);
            minuit.DefineParameter(2, "sigma", result.sigma, binWidth * 0.1, 0, 0);
            double arguments[2] = {5000, 0.01};
            minuit.mnexcm("MIGRAD", arguments, 2, result.status);
            int hesseStatus = 0;
            minuit.mnexcm("HESSE", arguments, 1, hesseStatus);
            minuit.GetParameter(0, result.amplitude, result.amplitudeError);
            minuit.GetParameter(1, result.mean, result.meanError);
            minuit.GetParameter(2, result.sigma, result.sigmaError);
            double edm = 0, errorDefinition = 0;
            int variableParameters = 0, totalParameters = 0;
            minuit.mnstat(result.chi2, edm, errorDefinition, variableParameters, totalParameters, result.covarianceStatus);
            if (result.status != 0 || hesseStatus != 0 || result.covarianceStatus != 3
                || !std::isfinite(result.amplitude) || !std::isfinite(result.mean)
                || !std::isfinite(result.sigma) || !std::isfinite(result.chi2)
                || !std::isfinite(result.meanError) || !std::isfinite(result.sigmaError)
                || result.sigma <= 0 || result.meanError <= 0 || result.sigmaError <= 0) break;
            if (pass == 1) fitSucceeded = true;
        }
        if (!fitSucceeded) continue;
        results.push_back(result);

        // Attach successful fit curves to the histogram for plotting and ROOT output.
        if (result.status == 0) {
            const std::string name = std::string(histogram.GetName()) + "_gaussian_" + std::to_string(results.size());
            if (auto *previous = histogram.GetListOfFunctions()->FindObject(name.c_str())) {
                histogram.GetListOfFunctions()->Remove(previous);
                delete previous;
            }
            auto *curve = new TF1(name.c_str(), "gaus", result.fitMin, result.fitMax);
            curve->SetParameters(result.amplitude, result.mean, result.sigma);
            curve->SetParError(0, result.amplitudeError);
            curve->SetParError(1, result.meanError);
            curve->SetParError(2, result.sigmaError);
            curve->SetChisquare(result.chi2);
            curve->SetNDF(result.ndf);
            curve->SetLineColor(kRed);
            histogram.GetListOfFunctions()->Add(curve);
        }
    }

    std::sort(results.begin(), results.end(), [](const PeakFitResult &a, const PeakFitResult &b) {
        return a.mean < b.mean;
    });
    std::cout << "\n" << histogram.GetName() << ": " << results.size()
              << " successfully fitted peaks (" << positions.size() << " candidates)\n";
    for (size_t i = 0; i < results.size(); ++i) {
        const auto &r = results[i];
        std::cout << "  Peak " << i + 1 << ": candidate=" << r.peakPosition
                  << ", range=[" << r.fitMin << ", " << r.fitMax << "]"
                  << " (pre-fit mean=" << r.rangeMean << ", sigma=" << r.rangeSigma << ")"
                  << ", initial sigma=" << r.initialSigma << '\n'
                  << "    A=" << r.amplitude << " +/- " << r.amplitudeError
                  << ", mean=" << r.mean << " +/- " << r.meanError
                  << ", sigma=" << r.sigma << " +/- " << r.sigmaError
                  << ", chi2/ndf=" << r.chi2 << "/" << r.ndf
                  << ", MIGRAD status=" << r.status << ", covariance status=" << r.covarianceStatus << '\n';
    }
    return results;
}

// Plot the returned peak-fit results with indices starting at zero.
// Indices enumerate retained peaks; missing or duplicate peaks prevent an absolute p.e. interpretation.
void PlotPeakResults(const std::vector<PeakFitResult> &peaks, const std::string &prefix,
                     const std::string &gain) {
    if (peaks.empty()) {
        std::cout << gain << ": no successful fits, skipping peak-index plots\n";
        return;
    }
    TGraphErrors meanGraph, sigmaGraph;
    meanGraph.SetName(("g_" + gain + "_mean_vs_index").c_str());
    sigmaGraph.SetName(("g_" + gain + "_sigma_vs_index").c_str());
    meanGraph.SetTitle((gain + ": fitted mean;" + (gain == "HG" ? "nPEs" : "Peak index (starting at 0)")
                        + ";Mean [ADC counts]").c_str());
    sigmaGraph.SetTitle((gain + ": fitted sigma;Peak index (starting at 0);Sigma [ADC counts]").c_str());
    const size_t meanPointCount = gain == "HG" ? std::min(size_t(6), peaks.size()) : peaks.size();
    for (size_t i = 0; i < peaks.size(); ++i) {
        if (i < meanPointCount) {
            meanGraph.SetPoint(i, i, peaks[i].mean);
            // nPEs is a discrete index with no x error; y errors are the Gaussian mean uncertainties.
            meanGraph.SetPointError(i, 0, peaks[i].meanError);
            std::cout << gain << " mean graph point: index=" << i
                      << ", mean=" << peaks[i].mean << ", y_error=" << peaks[i].meanError << '\n';
        }
        sigmaGraph.SetPoint(i, i, peaks[i].sigma);
        sigmaGraph.SetPointError(i, 0, peaks[i].sigmaError);
    }
    TCanvas canvas("peak_canvas", "Peak fit results", 1100, 750);
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.13);
    meanGraph.SetMarkerStyle(20);
    meanGraph.Draw("AP");
    // For HG, plot the first six peaks (nPEs=0..5) and exclude the first from the linear fit.
    TF1 linearFit("hg_mean_linear", "[0]+[1]*x", 1, meanPointCount - 1);
    TPaveText fitText(0.15, 0.66, 0.61, 0.89, "NDC");
    if (gain == "HG" && meanPointCount >= 4) {
        // TGraphErrors::Fit uses the y uncertainties in the weighted chi-square fit by default.
        // Do not use equal-weight options W/WW: chi2 = sum((mean_i - line_i)^2 / meanError_i^2).
        const int fitStatus = meanGraph.Fit(&linearFit, "QR");
        double weightedChi2 = 0;
        for (size_t i = 1; i < meanPointCount; ++i) {
            const double residual = (peaks[i].mean - linearFit.Eval(i)) / peaks[i].meanError;
            weightedChi2 += residual * residual;
        }
        std::cout << "HG mean linear fit (excluding nPEs=0): status=" << fitStatus
                  << ", intercept=" << linearFit.GetParameter(0) << " +/- " << linearFit.GetParError(0)
                  << ", slope=" << linearFit.GetParameter(1) << " +/- " << linearFit.GetParError(1)
                  << ", chi2/ndf=" << linearFit.GetChisquare() << "/" << linearFit.GetNDF() << '\n';
        std::cout << "HG weighted chi2 cross-check (using meanError): " << weightedChi2 << '\n';
        fitText.SetFillStyle(0);
        fitText.SetBorderSize(0);
        fitText.SetTextAlign(12);
        fitText.SetTextFont(42);
        fitText.SetTextSize(0.028);
        fitText.AddText("Mean = intercept + slope #times nPEs");
        fitText.AddText("Weighted fit, nPEs = 1..5; errors: #delta#mu");
        fitText.AddText(Form("Intercept = %.3f #pm %.3f ADC counts", linearFit.GetParameter(0), linearFit.GetParError(0)));
        fitText.AddText(Form("Slope = %.3f #pm %.3f ADC counts/p.e.", linearFit.GetParameter(1), linearFit.GetParError(1)));
        fitText.AddText(Form("#chi^{2}/ndf = %.2f / %d", linearFit.GetChisquare(), linearFit.GetNDF()));
        if (fitStatus != 0) fitText.AddText(Form("Fit failed: status %d", fitStatus));
        fitText.Draw();
    }
    canvas.SaveAs((prefix + "_PHA_" + gain + "_mean_vs_index.png").c_str());
    canvas.SaveAs((prefix + "_PHA_" + gain + "_mean_vs_index.pdf").c_str());
    meanGraph.Write();
    canvas.Clear();
    sigmaGraph.SetMarkerStyle(20);
    sigmaGraph.Draw("AP");
    canvas.SaveAs((prefix + "_PHA_" + gain + "_sigma_vs_index.png").c_str());
    canvas.SaveAs((prefix + "_PHA_" + gain + "_sigma_vs_index.pdf").c_str());
    sigmaGraph.Write();
}

int main(int argc, char **argv) {
    // Input and output paths can also be supplied on the command line.
    const std::string inputPath = argc > 1 ? argv[1] : "data/Run94_list.csv";
    const std::string outputDir = argc > 2 ? argv[2] : "plots";
    std::ifstream input(inputPath);
    if (!input) {
        std::cerr << "Cannot open: " << inputPath << '\n';
        return 1;
    }
    std::filesystem::create_directories(outputDir);
    const std::string prefix = outputDir + "/" + std::filesystem::path(inputPath).stem().string();
    gROOT->SetBatch(true);
    gStyle->SetOptStat(1110);

    // 1. Define three histograms: time in microseconds and PHA in ADC counts.
    //    Keep all 8192 ADC channels, with one channel per bin.
    TH1D hTime("h_TStamp_us", "Trigger time;TStamp [#mus];Entries", 300, 0, 6.0e8);
    TH1D hLG("h_PHA_LG", "Low gain charge spectrum;PHA LG [ADC counts];Entries", 8192, -0.5, 8191.5);
    TH1D hHG("h_PHA_HG", "High gain charge spectrum;PHA HG [ADC counts];Entries", 8192, -0.5, 8191.5);

    // 2. Read the CSV line by line using the fixed column order:
    // TStamp_us,Trg_Id,Board_Id,Num_Chs,ChannelMask,CH_Id,DataType,PHA_LG,PHA_HG
    // Read ChannelMask and DataType as strings without analyzing them.
    std::string line;
    unsigned long long entries = 0, lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        if (line.rfind("//", 0) == 0 || line.rfind("TStamp_us,", 0) == 0) continue;

        for (char &character : line) {
            if (character == ',') character = ' ';
        }
        std::istringstream row(line);
        double time;
        unsigned long long trigger;
        int board, numChannels, channel, phaLG, phaHG;
        std::string channelMask, dataType;
        if (!(row >> time >> trigger >> board >> numChannels >> channelMask
                  >> channel >> dataType >> phaLG >> phaHG)) {
            std::cerr << "Invalid CSV row at line " << lineNumber << '\n';
            return 1;
        }

        // 3. Fill the histograms with the values from this row.
        hTime.Fill(time);
        hLG.Fill(phaLG);
        hHG.Fill(phaHG);
        ++entries;
    }
    if (entries == 0) {
        std::cerr << "No data rows in: " << inputPath << '\n';
        return 1;
    }

    // Find and fit peaks in the full HG and LG spectra before zooming the plot axes.
    const std::vector<PeakFitResult> peaksLG = FindAndFitPeaks(hLG);
    const std::vector<PeakFitResult> peaksHG = FindAndFitPeaks(hHG);

    // 4. Save each plot as PNG and PDF, and save the complete histograms in a ROOT file.
    //    Zoom only the displayed PHA range; all 8192 bins remain available.
    hLG.GetXaxis()->SetRangeUser(-0.5, 499.5);
    hHG.GetXaxis()->SetRangeUser(-0.5, 999.5);
    TH1D *histograms[3] = {&hTime, &hLG, &hHG};
    const char *names[3] = {"TStamp_us", "PHA_LG", "PHA_HG"};
    TFile output((prefix + "_histograms.root").c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "Cannot create ROOT output\n";
        return 1;
    }
    TCanvas canvas("canvas", "SiPM test", 1100, 750);
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.13);
    for (int i = 0; i < 3; ++i) {
        histograms[i]->SetLineWidth(2);
        histograms[i]->Draw(i >= 1 ? "" : "HIST");
        canvas.SaveAs((prefix + "_" + names[i] + ".png").c_str());
        canvas.SaveAs((prefix + "_" + names[i] + ".pdf").c_str());
        // Use the full range in the ROOT file and the zoomed range in the images.
        histograms[i]->GetXaxis()->SetRange(0, 0);
        histograms[i]->Write();
        std::cout << names[i] << ": entries=" << histograms[i]->GetEntries()
                  << ", underflow=" << histograms[i]->GetBinContent(0)
                  << ", overflow=" << histograms[i]->GetBinContent(histograms[i]->GetNbinsX() + 1)
                  << '\n';
    }
    // 5. Plot mean and sigma versus index with uncertainties from the successful fits.
    PlotPeakResults(peaksLG, prefix, "LG");
    PlotPeakResults(peaksHG, prefix, "HG");
    output.Close();
    std::cout << "Processed " << entries << " rows. Output: " << outputDir << '\n';
    return 0;
}
