#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TList.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TPaveText.h>

#include <algorithm>
#include <cmath>
#include <deque>
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
    double rangeMean = 0, rangeSigma = 0; // Pre-fit parameters: the final range is rangeMean +/- 1.5*rangeSigma.
    double amplitude = std::numeric_limits<double>::quiet_NaN();
    double mean = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double amplitudeError = 0;
    double meanError = 0;        // Mean uncertainty from TMinuit/HESSE, distinct from the peak width sigma.
    double sigmaError = 0;       // Fit uncertainty of the peak width sigma.
    double meanSigmaCovariance = 0; // Covariance used to propagate sigma/mean uncertainty.
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int status = -1;             // MIGRAD: 0 means convergence; -1 means the fit has not been performed.
    int covarianceStatus = 0;    // TMinuit: 3 means a full, accurate covariance matrix.
};

// Pass the current histogram and bin range through the fixed TMinuit FCN interface.
namespace {
TH1 *fitHistogram = nullptr;
int fitFirstBin = 0, fitLastBin = 0;
// int fitPeakCount = 0; // Global fitting is temporarily disabled.

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

// Sum of fitPeakCount Gaussian peaks plus B*exp(-x/tau).
// Global fitting is temporarily disabled.
// void MultiGaussianChi2(int &, double *, double &chi2, double *parameters, int) {
//     chi2 = 0;
//     for (int bin = fitFirstBin; bin <= fitLastBin; ++bin) {
//         const double x = fitHistogram->GetBinCenter(bin);
//         const double backgroundAmplitude = parameters[3 * fitPeakCount];
//         const double decayLength = parameters[3 * fitPeakCount + 1];
//         if (decayLength == 0) { chi2 = 1e100; return; }
//         double model = backgroundAmplitude * std::exp(-x / decayLength);
//         if (!std::isfinite(model)) { chi2 = 1e100; return; }
//         for (int peak = 0; peak < fitPeakCount; ++peak) {
//             if (parameters[3 * peak + 2] == 0) { chi2 = 1e100; return; }
//             const double z = (x - parameters[3 * peak + 1]) / parameters[3 * peak + 2];
//             model += parameters[3 * peak] * std::exp(-0.5 * z * z);
//         }
//         const double error = fitHistogram->GetBinError(bin);
//         const double residual = fitHistogram->GetBinContent(bin) - model;
//         chi2 += residual * residual / (error > 0 ? error * error : 1.0);
//     }
// }
}

// Find peaks, estimate mu and sigma in the half-height range, then fit within mu +/- 1.5*sigma.
// Return only converged fits with reliable errors, sorted by increasing fitted mean.
std::vector<PeakFitResult> FindAndFitPeaks(TH1 &histogram, std::vector<double> *candidatePositions = nullptr) {
    const int maxPeaks = histogram.GetNbinsX();
    const double sigma = 8.0;       // Peak-search width parameter in bins.
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
    if (candidatePositions) *candidatePositions = positions;
    std::cout << histogram.GetName() << " candidate positions (>100 raw bin counts):";
    for (double position : positions) std::cout << " " << position;
    std::cout << '\n';

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
                result.fitMin = std::max(histogram.GetXaxis()->GetXmin(), result.rangeMean - 1.5 * result.rangeSigma);
                result.fitMax = std::min(histogram.GetXaxis()->GetXmax(), result.rangeMean + 1.5 * result.rangeSigma);
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
            double covariance[9] = {};
            minuit.mnemat(covariance, 3);
            result.meanSigmaCovariance = covariance[1 * 3 + 2];
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
            curve->SetBit(TF1::kNotDraw); // Retain individual fits without drawing them.
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

// Keep the independent fits above; use them only as seeds for a simultaneous fit.
// Stop the simultaneous fit at the last seed mean plus 1.5 seed sigmas.
// Global fitting is temporarily disabled.
// std::vector<PeakFitResult> FitMultipleGaussians(
//     TH1 &histogram, const std::vector<PeakFitResult> &individualFits) {
//     // Use every successful independent fit; no manually selected peak count.
//     const auto &seeds = individualFits;
//     const int peakCount = static_cast<int>(seeds.size());
//     if (seeds.empty()) {
//         std::cout << histogram.GetName() << ": no independent fits available for the global fit\n";
//         return {};
//     }
//
//     const double fitMin = histogram.GetXaxis()->GetXmin();
//     const double fitMax = std::min(histogram.GetXaxis()->GetXmax(),
//                                    seeds.back().mean + 1.5 * seeds.back().sigma);
//     fitHistogram = &histogram;
//     fitPeakCount = peakCount;
//     fitFirstBin = 1;
//     fitLastBin = std::min(histogram.GetNbinsX(), histogram.FindFixBin(fitMax));
//     if (histogram.GetBinCenter(fitLastBin) > fitMax) --fitLastBin;
//     const int parameterCount = 3 * peakCount + 2;
//     const int ndf = fitLastBin - fitFirstBin + 1 - parameterCount;
//     if (ndf <= 0) return {};
//     TMinuit minuit(parameterCount);
//     minuit.SetPrintLevel(-1);
//     minuit.SetFCN(MultiGaussianChi2);
//     minuit.SetErrorDef(1.0);
//     for (int peak = 0; peak < peakCount; ++peak) {
//         const auto &seed = seeds[peak];
//         std::cout << histogram.GetName() << " global seed " << peak
//                   << ": A=" << seed.amplitude << ", mean=" << seed.mean << ", sigma=" << seed.sigma << '\n';
//         minuit.DefineParameter(3 * peak, Form("A%d", peak), seed.amplitude,
//                                std::max(1.0, seed.amplitude * 0.01), 0, 0);
//         minuit.DefineParameter(3 * peak + 1, Form("mean%d", peak), seed.mean, 0.1, 0, 0);
//         minuit.DefineParameter(3 * peak + 2, Form("sigma%d", peak), seed.sigma, 0.1, 0, 0);
//     }
//     const double backgroundInitial = std::max(1.0, histogram.GetMaximum() * 0.1);
//     const double decayInitial = (fitMax - fitMin) / 2;
//     minuit.DefineParameter(3 * peakCount, "background_B", backgroundInitial,
//                            std::max(1.0, backgroundInitial * 0.01), 0, 0);
//     minuit.DefineParameter(3 * peakCount + 1, "decay_tau", decayInitial, decayInitial * 0.01, 0, 0);
//     double arguments[2] = {60000, 0.01};
//     int status = 0, hesseStatus = 0;
//     minuit.mnexcm("MIGRAD", arguments, 2, status);
//     minuit.mnexcm("HESSE", arguments, 1, hesseStatus);
//     double chi2 = 0, edm = 0, errorDefinition = 0;
//     int nVariable = 0, nTotal = 0, covarianceStatus = 0;
//     minuit.mnstat(chi2, edm, errorDefinition, nVariable, nTotal, covarianceStatus);
//     std::cout << histogram.GetName() << " global " << peakCount << "-Gaussian + exponential fit: status=" << status
//               << ", HESSE=" << hesseStatus << ", covariance=" << covarianceStatus
//               << ", chi2/ndf=" << chi2 << "/" << ndf
//               << ", range=[" << fitMin << ", " << fitMax << "]\n";
//     if (!std::isfinite(chi2)) return {};
//     const bool fitSucceeded = status == 0 && hesseStatus == 0 && covarianceStatus == 3;
//     if (!fitSucceeded)
//         std::cout << "WARNING: global fit did not converge reliably; draw the current model only,"
//                      " and do not return it as a successful fit\n";
//
//     const std::string name = std::string(histogram.GetName()) + "_multi_gaussians";
//     std::string formula;
//     for (int peak = 0; peak < peakCount; ++peak) {
//         if (peak > 0) formula += "+";
//         formula += "gaus(" + std::to_string(3 * peak) + ")";
//     }
//     formula += "+[" + std::to_string(3 * peakCount) + "]*exp(-x/["
//              + std::to_string(3 * peakCount + 1) + "])";
//     auto *total = new TF1(name.c_str(), formula.c_str(), fitMin, fitMax);
//     total->SetNpx(8192);
//     for (int parameter = 0; parameter < parameterCount; ++parameter) {
//         double value = 0, error = 0;
//         minuit.GetParameter(parameter, value, error);
//         if (!std::isfinite(value)) { delete total; return {}; }
//         total->SetParameter(parameter, value);
//         total->SetParError(parameter, error);
//     }
//     total->SetChisquare(chi2);
//     total->SetNDF(ndf);
//     total->SetLineColor(kBlue + 1);
//     total->SetLineWidth(3);
//     total->SetTitle(Form("%d Gaussians + exponential%s", peakCount, fitSucceeded ? "" : " (not converged)"));
//     histogram.GetListOfFunctions()->Add(total);
//
//     std::vector<double> covariance(parameterCount * parameterCount);
//     minuit.mnemat(covariance.data(), parameterCount);
//
//     // Decompose the simultaneous model, not the earlier independent fits.
//     // Use peakCount signal curves plus one exponential-background curve.
//     for (int component = 0; component <= peakCount; ++component) {
//         const std::string componentName = name + "_component_" + std::to_string(component);
//         TF1 *curve = nullptr;
//         if (component < peakCount) {
//             curve = new TF1(componentName.c_str(), "gaus", fitMin, fitMax);
//             curve->SetParameters(total->GetParameter(3 * component),
//                                  total->GetParameter(3 * component + 1),
//                                  total->GetParameter(3 * component + 2));
//             curve->SetTitle(Form("Signal %d", component));
//             curve->SetLineStyle(2);
//         } else {
//             curve = new TF1(componentName.c_str(), "[0]*exp(-x/[1])", fitMin, fitMax);
//             curve->SetParameters(total->GetParameter(3 * peakCount), total->GetParameter(3 * peakCount + 1));
//             curve->SetTitle("Exponential background");
//             curve->SetLineStyle(7);
//         }
//         Float_t red = 0, green = 0, blue = 0;
//         TColor::HSV2RGB(360.0f * component / (peakCount + 1), 0.85f, 0.75f, red, green, blue);
//         curve->SetLineColor(TColor::GetColor(red, green, blue));
//         curve->SetLineWidth(2);
//         curve->SetNpx(4096);
//         histogram.GetListOfFunctions()->Add(curve);
//     }
//
//     std::vector<PeakFitResult> results;
//     for (int peak = 0; peak < peakCount; ++peak) {
//         PeakFitResult result;
//         result.peakPosition = seeds[peak].peakPosition;
//         result.amplitude = total->GetParameter(3 * peak);
//         result.amplitudeError = total->GetParError(3 * peak);
//         result.mean = total->GetParameter(3 * peak + 1);
//         result.meanError = total->GetParError(3 * peak + 1);
//         result.sigma = std::abs(total->GetParameter(3 * peak + 2));
//         result.sigmaError = total->GetParError(3 * peak + 2);
//         result.meanSigmaCovariance = covariance[(3 * peak + 1) * parameterCount + 3 * peak + 2]
//                                   * std::copysign(1.0, total->GetParameter(3 * peak + 2));
//         result.fitMin = total->GetXmin();
//         result.fitMax = total->GetXmax();
//         result.chi2 = chi2;
//         result.ndf = total->GetNDF();
//         result.status = status;
//         result.covarianceStatus = covarianceStatus;
//         results.push_back(result);
//         std::cout << "  Global peak " << peak << ": A=" << result.amplitude << " +/- " << result.amplitudeError
//                   << ", mean=" << result.mean << " +/- " << result.meanError
//                   << ", sigma=" << result.sigma << " +/- " << result.sigmaError << '\n';
//     }
//     std::cout << "  Background B=" << total->GetParameter(3 * peakCount)
//               << " +/- " << total->GetParError(3 * peakCount)
//               << ", tau=" << total->GetParameter(3 * peakCount + 1)
//               << " +/- " << total->GetParError(3 * peakCount + 1) << " ADC counts\n";
//     if (total->GetParameter(3 * peakCount) <= 0 || total->GetParameter(3 * peakCount + 1) <= 0)
//         std::cout << "WARNING: unconstrained background parameters do not describe a positive decaying background\n";
//     if (!fitSucceeded) return {};
//     return results;
// }

// Plot the returned peak-fit results with indices starting at zero.
// Indices enumerate retained peaks; missing or duplicate peaks prevent an absolute p.e. interpretation.
void PlotPeakResults(const std::vector<PeakFitResult> &peaks, const std::string &prefix,
                     const std::string &gain) {
    if (peaks.empty()) {
        std::cout << gain << ": no successful fits, skipping peak-index plots\n";
        return;
    }
    TGraphErrors meanGraph, resolutionGraph;
    meanGraph.SetName(("g_" + gain + "_mean_vs_index").c_str());
    resolutionGraph.SetName(("g_" + gain + "_resolution_vs_nPEs").c_str());
    meanGraph.SetTitle((gain + ": fitted mean;" + (gain == "HG" ? "nPEs" : "Peak index (starting at 0)")
                        + ";Mean [ADC counts]").c_str());
    resolutionGraph.SetTitle((gain + ": relative resolution;nPEs;Relative resolution #sigma/#mu [%]").c_str());
    const size_t meanPointCount = peaks.size();
    for (size_t i = 0; i < peaks.size(); ++i) {
        if (i < meanPointCount) {
            meanGraph.SetPoint(i, i, peaks[i].mean);
            // nPEs is a discrete index with no x error; y errors are the Gaussian mean uncertainties.
            meanGraph.SetPointError(i, 0, peaks[i].meanError);
            std::cout << gain << " mean graph point: index=" << i
                      << ", mean=" << peaks[i].mean << ", y_error=" << peaks[i].meanError << '\n';
        }
        // Resolution uses the raw fitted mean, without pedestal subtraction.
        // Propagate both parameter errors and their covariance for R=100*sigma/mean.
        if (peaks[i].mean > 0) {
            const double resolution = 100 * peaks[i].sigma / peaks[i].mean;
            const double derivativeSigma = 100 / peaks[i].mean;
            const double derivativeMean = -100 * peaks[i].sigma / (peaks[i].mean * peaks[i].mean);
            const double variance = std::pow(derivativeSigma * peaks[i].sigmaError, 2)
                                  + std::pow(derivativeMean * peaks[i].meanError, 2)
                                  + 2 * derivativeSigma * derivativeMean * peaks[i].meanSigmaCovariance;
            const double resolutionError = std::sqrt(std::max(0.0, variance));
            const int point = resolutionGraph.GetN();
            resolutionGraph.SetPoint(point, i, resolution);
            resolutionGraph.SetPointError(point, 0, resolutionError);
            std::cout << gain << " resolution point: nPEs=" << i << ", sigma/mean="
                      << resolution << " +/- " << resolutionError << " percent\n";
        }
    }
    TCanvas canvas("peak_canvas", "Peak fit results", 1100, 750);
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.13);
    meanGraph.SetMarkerStyle(20);
    meanGraph.Draw("AP");
    // Plot all successful peaks and exclude the first from the HG linear fit.
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
        fitText.AddText(Form("Weighted fit, nPEs = 1..%zu; errors: #delta#mu", meanPointCount - 1));
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
    if (resolutionGraph.GetN() > 0) {
        resolutionGraph.SetMarkerStyle(20);
        resolutionGraph.Draw("AP");
        canvas.SaveAs((prefix + "_PHA_" + gain + "_resolution_vs_nPEs.png").c_str());
        canvas.SaveAs((prefix + "_PHA_" + gain + "_resolution_vs_nPEs.pdf").c_str());
        resolutionGraph.Write();
    }
    std::filesystem::remove(prefix + "_PHA_" + gain + "_sigma_vs_index.png");
    std::filesystem::remove(prefix + "_PHA_" + gain + "_sigma_vs_index.pdf");
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

    // 1. Define three histograms: time in seconds and PHA in ADC counts.
    //    Keep all 8192 ADC channels, with one channel per bin.
    TH1D hTime("h_TStamp_us", "Trigger time;TStamp [s];Entries / 2 s", 300, 0, 600);
    // Resolve recorded bursts and gaps: 10-second window, 0.01-second bins.
    TH1D hTimeFirst("h_TStamp_first10s", "Recorded events: first 10 s;TStamp [s];Trigger rates [kHz]",
                    1000, 0, 10);
    // Retain only the latest 10 seconds of timestamps, not the entire CSV.
    std::deque<double> recentTimestamps;
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
        const double timeSeconds = time * 1e-6; // CSV timestamps are in microseconds.
        hTime.Fill(timeSeconds);
        if (timeSeconds >= 0 && timeSeconds < 10) hTimeFirst.Fill(timeSeconds);
        recentTimestamps.push_back(timeSeconds);
        while (recentTimestamps.front() < timeSeconds - 10) recentTimestamps.pop_front();
        hLG.Fill(phaLG);
        hHG.Fill(phaHG);
        ++entries;
    }
    if (entries == 0) {
        std::cerr << "No data rows in: " << inputPath << '\n';
        return 1;
    }

    // Anchor the final window to the last recorded timestamp, not the preset run duration.
    const double lastWindowEnd = std::max(10.0, recentTimestamps.back());
    const double lastWindowStart = lastWindowEnd - 10;
    TH1D hTimeLast("h_TStamp_last10s", "Recorded events: last 10 s;TStamp [s];Trigger rates [kHz]",
                   1000, lastWindowStart, lastWindowEnd);
    for (double timestamp : recentTimestamps) {
        // ROOT excludes the upper axis edge; include the final event in the last bin.
        hTimeLast.Fill(std::min(timestamp, std::nextafter(lastWindowEnd, lastWindowStart)));
    }

    // Find and fit peaks in the full HG and LG spectra before zooming the plot axes.
    // const std::vector<PeakFitResult> peaksLG = FindAndFitPeaks(hLG);
    const std::vector<PeakFitResult> peaksHG = FindAndFitPeaks(hHG);

    // Simultaneous fits supplement, rather than replace, the individual fits.
    // const auto globalPeaksLG = FitMultipleGaussians(hLG, peaksLG);
    // const auto globalPeaksHG = FitMultipleGaussians(hHG, peaksHG);

    // 4. Save charge spectra with linear/logarithmic y axes and preserve all fit functions.
    //    Zoom only the displayed PHA range; all 8192 bins remain available.
    hLG.GetXaxis()->SetRangeUser(-0.5, 499.5);
    hHG.GetXaxis()->SetRangeUser(-0.5, 999.5);
    TH1D *chargeHistograms[2] = {&hLG, &hHG};
    const char *chargeNames[2] = {"PHA_LG", "PHA_HG"};
    TFile output((prefix + "_histograms.root").c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "Cannot create ROOT output\n";
        return 1;
    }
    TCanvas canvas("canvas", "SiPM test", 1100, 750);
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.13);
    // Keep the overview as raw counts per 2-second bin. Convert only the
    // first/last 10-second plots to recorded rates in kHz.
    hTimeFirst.Scale(1.0 / (hTimeFirst.GetBinWidth(1) * 1000));
    hTimeLast.Scale(1.0 / (hTimeLast.GetBinWidth(1) * 1000));

    // Draw and save the timestamp histogram separately. No fit or log axis.
    canvas.SetLogy(false);
    hTime.SetLineWidth(2);
    hTime.SetLineColor(kBlack);
    hTime.Draw("HIST");
    canvas.SaveAs((prefix + "_TStamp_us.png").c_str());
    canvas.SaveAs((prefix + "_TStamp_us.pdf").c_str());
    hTime.Write();
    std::cout << "TStamp_us: entries=" << hTime.GetEntries()
              << ", underflow=" << hTime.GetBinContent(0)
              << ", overflow=" << hTime.GetBinContent(hTime.GetNbinsX() + 1) << '\n';

    // Draw the first 10 seconds with 0.01-second bins.
    canvas.Clear();
    canvas.SetLogy(false);
    hTimeFirst.SetLineColor(kBlack);
    hTimeFirst.SetLineWidth(2);
    hTimeFirst.SetMinimum(0);
    hTimeFirst.Draw("HIST");
    canvas.SaveAs((prefix + "_TStamp_first10s.png").c_str());
    canvas.SaveAs((prefix + "_TStamp_first10s.pdf").c_str());
    hTimeFirst.Write();

    // Draw the final 10 seconds using the same bin width and linear y axis.
    canvas.Clear();
    hTimeLast.SetLineColor(kBlack);
    hTimeLast.SetLineWidth(2);
    hTimeLast.SetMinimum(0);
    hTimeLast.Draw("HIST");
    canvas.SaveAs((prefix + "_TStamp_last10s.png").c_str());
    canvas.SaveAs((prefix + "_TStamp_last10s.pdf").c_str());
    hTimeLast.Write();
    std::cout << "Final time window: [" << lastWindowStart << ", " << lastWindowEnd
              << "] s, entries=" << hTimeLast.GetEntries() << '\n';

    // Draw the charge histograms below, independently of the timestamp plots.
    TLegend legend(0.60, 0.46, 0.89, 0.84);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.027);
    // Each charge histogram produces independent-fit images on linear and log axes.
    // The global-fit view is temporarily disabled.
    for (int chargeIndex = 0; chargeIndex < 2; ++chargeIndex) {
        TH1D &histogram = *chargeHistograms[chargeIndex];
        const std::string histogramLabel = chargeNames[chargeIndex];
        histogram.SetLineWidth(2);
        histogram.SetLineColor(kBlack);

        // Select the fit overlay: 0 = independent peaks, 1 = global model.
        for (int fitViewIndex = 0; fitViewIndex < 1; ++fitViewIndex) {

            // Select the vertical axis: 0 = linear, 1 = logarithmic.
            for (int yAxisScaleIndex = 0; yAxisScaleIndex < 2; ++yAxisScaleIndex) {
                const bool useLogYAxis = yAxisScaleIndex == 1;

                // 1. Clear the previous plot, configure the axis, and draw data.
                // HIST suppresses automatic drawing of attached fit functions;
                // the selected curves are added explicitly below using SAME.
                canvas.Clear();
                canvas.SetLogy(useLogYAxis);
                histogram.SetMinimum(useLogYAxis ? 0.5 : 0);
                histogram.Draw("HIST");
                legend.Clear();

                legend.AddEntry(&histogram, "Data", "l");
                if (fitViewIndex == 0) {
                    // 2a. Overlay the individual Gaussian fits. Their saved
                    // names start at 1 and form a consecutive sequence.
                    for (int independentPeakNumber = 1; ; ++independentPeakNumber) {
                        const std::string fitName = std::string(histogram.GetName())
                            + "_gaussian_" + std::to_string(independentPeakNumber);
                        auto *independentCurve = histogram.GetFunction(fitName.c_str());
                        if (!independentCurve) break;
                        independentCurve->ResetBit(TF1::kNotDraw);
                        independentCurve->Draw("SAME");
                        if (independentPeakNumber == 1)
                            legend.AddEntry(independentCurve, "Independent Gaussian fits", "l");
                    }

                }
                // else {
                //     // 2b. Overlay the global model's signal and background
                //     // components, then draw the total fit above them.
                //     const std::string totalFitName = std::string(histogram.GetName()) + "_multi_gaussians";
                //     auto *totalFit = histogram.GetFunction(totalFitName.c_str());
                //     if (totalFit) {
                //         // Each Gaussian has three parameters; the exponential
                //         // background has two. Component N is the background.
                //         const int numberOfSignalPeaks = (totalFit->GetNpar() - 2) / 3;
                //         for (int componentIndex = 0; componentIndex <= numberOfSignalPeaks; ++componentIndex) {
                //             const std::string componentName = totalFitName
                //                 + "_component_" + std::to_string(componentIndex);
                //             if (auto *componentCurve = histogram.GetFunction(componentName.c_str())) {
                //                 componentCurve->Draw("SAME");
                //                 legend.AddEntry(componentCurve, componentCurve->GetTitle(), "l");
                //             }
                //         }
                //         totalFit->Draw("SAME");
                //         legend.AddEntry(totalFit, totalFit->GetTitle(), "l");
                //     } else {
                //         // LG fitting is disabled, so its global view shows
                //         // only the charge data and this explanatory entry.
                //         legend.AddEntry(static_cast<TObject *>(nullptr), "Global fit unavailable", "");
                //     }
                // }
                legend.Draw();

                // 3. Save the selected fit/scale combination as PNG and PDF.
                // Example: Run94_list_PHA_HG_individual_log.png.
                std::string imageSuffix;
                imageSuffix = fitViewIndex == 0 ? "_individual" : "_global";
                imageSuffix += useLogYAxis ? "_log" : "_linear";
                const std::string imagePath = prefix + "_" + histogramLabel + imageSuffix;
                canvas.SaveAs((imagePath + ".png").c_str());
                canvas.SaveAs((imagePath + ".pdf").c_str());
            }
        }

        // Keep independent fits in the ROOT file, but disable their automatic
        // display so they do not get mixed with the global-fit components.
        for (int independentPeakNumber = 1; ; ++independentPeakNumber) {
            const std::string fitName = std::string(histogram.GetName())
                + "_gaussian_" + std::to_string(independentPeakNumber);
            auto *independentCurve = histogram.GetFunction(fitName.c_str());
            if (!independentCurve) break;
            independentCurve->SetBit(TF1::kNotDraw);
        }

        // 5. Save the complete histogram and attached functions to the ROOT file.
        // Reset the display zoom so every stored ADC bin remains visible.
        histogram.SetMinimum(0);
        histogram.GetXaxis()->SetRange(0, 0);
        histogram.Write();
        std::cout << histogramLabel << ": entries=" << histogram.GetEntries()
                  << ", underflow=" << histogram.GetBinContent(0)
                  << ", overflow=" << histogram.GetBinContent(histogram.GetNbinsX() + 1)
                  << '\n';
    }
    // 5. Plot mean and sigma versus index with uncertainties from the successful fits.
    // PlotPeakResults(peaksLG, prefix, "LG");
    PlotPeakResults(peaksHG, prefix, "HG");
    output.Close();
    std::cout << "Processed " << entries << " rows. Output: " << outputDir << '\n';
    return 0;
}
