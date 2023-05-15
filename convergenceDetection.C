/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Nikola Majksner, Wouter Remmerie, AirShaper
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "convergenceDetection.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "cartesianCS.H"
#include <vector>
#include <algorithm>
#include <map>

using namespace Foam::functionObjects::fieldValues;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(convergenceDetection, 0);
        addToRunTimeSelectionTable(functionObject, convergenceDetection, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::convergenceDetection::convergenceDetection(
    const word &name,
    const Time &runTime,
    const dictionary &dict)
    : forces(name, runTime, dict, false),
      windowNormalization_(0.5),
      windowGradient_(0.5),
      windowGradientAveraging_(0.99),
      windowEvaluation_(0.333),
      windowEvaluationAveraging_(0.333),
      thresholdConvergence_(0.075),
      thresholdAveraging_(0.025),
      iterationMaxConvergence_(4000),
      iterationMaxAveraging_(2000),
      convergenceStabilityFactor_(2),
      iterationMinConvergence_(200),
      iterationMinAveraging_(200),
      convergenceFound_(false),
      simulationFinished_(false),
      forcedConvergence_(false),
      currentIteration_(0),
      averagingStartedAt_(0),
      patches_(wordRes()),
      rho_("rhoInf"),
      rhoInf_(1.225),
      CofR_(Zero),
      normalizedDataMeanConverged_(),
      porousObjectsNames_(),
      rotatingObjectsNames_(),
      porousNames_(),
      rotatingNames_(),
      porousObjects_(),
      rotatingObjects_(),
      coordSysPtrC_(nullptr),
      convergenceData_(),
      gradData_(),
      gradDataAveraging_(),
      convergenceNormalized_(),
      dataNormalized_(),
      convergenceNormalizedAveraging_(),
      totalForceFilePtr_(nullptr),
      gradientFilePtr_(nullptr),
      gradientAveragedFilePtr_(nullptr)

{
    read(dict);

    setCoordinateSystem(dict);

    if (dict.found("porous"))
    {
        porousObjectsNames_ = dict.subDict("porous").toc();
        for (const word &name : porousObjectsNames_)
        {
            porousObjects_.push_back(new surfaceFieldValue(name, runTime, dict.subDict("porous").subDict(name)));
            const dictionary &porousSubDict = dict.subDict("porous").subDict(name);
            word porousName = porousSubDict.lookupOrDefault<word>("name", name);
            porousNames_.push_back(porousName);
        }
    }

    if (dict.found("rotation_forces"))
    {
        rotatingObjectsNames_ = dict.subDict("rotation_forces").toc();
        for (const word &name : rotatingObjectsNames_)
        {
            rotatingObjects_.push_back(new forces(name, runTime, dict.subDict("rotation_forces").subDict(name)));
            coordSysPtrC_.clear();
            coordSysPtrC_.reset(new coordSystem::cartesian(dict.subDict("rotation_forces").subDict(name)));
            const dictionary &rotatingSubDict = dict.subDict("rotation_forces").subDict(name);
            word rotatingName = rotatingSubDict.lookupOrDefault<word>("name", name);
            rotatingNames_.push_back(rotatingName);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::convergenceDetection::read(const dictionary &dict)
{
    forces::read(dict);

    dict.readIfPresent("windowNormalization", windowNormalization_);
    dict.readIfPresent("windowGradient", windowGradient_);
    dict.readIfPresent("windowGradientAveraging", windowGradientAveraging_);
    dict.readIfPresent("windowEvaluation", windowEvaluation_);
    dict.readIfPresent("windowEvaluationAveraging", windowEvaluationAveraging_);
    dict.readIfPresent("thresholdConvergence", thresholdConvergence_);
    dict.readIfPresent("thresholdAveraging", thresholdAveraging_);
    dict.readIfPresent("iterationMaxConvergence", iterationMaxConvergence_);
    dict.readIfPresent("iterationMaxAveraging", iterationMaxAveraging_);
    dict.readIfPresent("convergenceStabilityFactor_", convergenceStabilityFactor_);
    dict.readIfPresent("iterationMinConvergence", iterationMinConvergence_);
    dict.readIfPresent("iterationMinAveraging", iterationMinAveraging_);

    if (!exists(time().globalPath() + "/system/averaging"))
    {
        FatalError
            << "Averaging file does not exists, please create averaging file in system/averaging and include it in controlDict (#include \"averaging\")" << nl
            << "https://www.openfoam.com/documentation/guides/latest/doc/guide-fos-field-fieldAverage.html" << nl
            << abort(FatalError);

        return false;
    }

    return true;
}

bool Foam::functionObjects::convergenceDetection::execute()
{
    forces::execute();

    if (!porousObjects_.empty())
    {
        unsigned porousCount = 0;
        for (surfaceFieldValue *sfv : porousObjects_)
        {
            sfv->write();
            scalar area = sfv->getResult<scalar>("areaNormalIntegrate(" + porousNames_[porousCount] + ",U)");

            convergenceData_[porousNames_[porousCount]].push_back(area);
            porousCount++;
        }
    }

    if (!rotatingObjects_.empty())
    {
        unsigned roCount = 0;
        for (forces *f : rotatingObjects_)
        {
            f->calcForcesMoments();
            const vector forces = f->forceEff();

            const double thrust = coordSysPtrC_->localVector(forces)[0];

            convergenceData_[rotatingNames_[roCount]].push_back(thrust);
            roCount++;
        }
    }

    const vector sumForce = forces::forceEff();

    const double totalForce = sqrt(pow(sumForce[0], 2) + pow(sumForce[1], 2) + pow(sumForce[2], 2));

    convergenceData_["forces"].push_back(totalForce);

    currentIteration_ = time().value();

    // wait for a few iterations before starting gradient calculation
    if (currentIteration_ >= 5)
    {

        std::map<std::string, double> caclculatedGradient = calculateGradient();
        for (const auto &cd : caclculatedGradient)
        {
            gradData_[cd.first].push_back(cd.second);
        }
    }

    // need few iterations before starting calculation
    if (currentIteration_ >= iterationMinConvergence_)
    {
        if (!convergenceFound_)
        {
            convergenceEvaluation();
        }
        // Start to check at iterationMinAveraging to avoid calculate each iteration
        if (currentIteration_ >= iterationMinAveraging_ &&
            !simulationFinished_ &&
            convergenceFound_)
        {
            std::map<std::string, double> cmg(convergenceMaxGradient());

            bool shouldBeKickedOut = std::any_of(cmg.begin(), cmg.end(), [this](const auto &cmg)
                                                 { return cmg.second >= thresholdConvergence_; });

            if ((shouldBeKickedOut && !forcedConvergence_) &&
                (currentIteration_ < iterationMaxConvergence_))
            {
                stopAveraging();
            }
            std::map<std::string, double> caclculatedGradientAveraging = calculateGradientAveraging();
            for (const auto &c : caclculatedGradientAveraging)
            {
                gradDataAveraging_[c.first].push_back(c.second);
            }
            isExploded();
            isFinished();
        }
    }

    return true;
}

bool Foam::functionObjects::convergenceDetection::end()
{
    return true;
}

bool Foam::functionObjects::convergenceDetection::write()
{
    if (writeToFile())
    {

        createDataFile();
        writeDataFile(totalForceFilePtr_(), convergenceData_["forces"].back());
        if (currentIteration_ >= 5)
        {
            writeDataFile(gradientFilePtr_(), gradData_["forces"].back());
        }
        if (convergenceFound_)
        {
            writeDataFile(gradientAveragedFilePtr_(), gradDataAveraging_["forces"].back());
        }
        forces::write();
    }
    return true;
}

// * * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::convergenceDetection::createDataFile()
{
    if (!totalForceFilePtr_.valid())
    {
        totalForceFilePtr_ = createFile("totalForce");
        writeDataFileHeader("totalForce", "total (x,y,z)", totalForceFilePtr_());
    }

    if (!gradientFilePtr_.valid())
    {
        gradientFilePtr_ = createFile("gradient");
        writeDataFileHeader("Gradient", "Gradient", gradientFilePtr_());
    }

    if (!gradientAveragedFilePtr_.valid())
    {
        gradientAveragedFilePtr_ = createFile("gradientAveraged");
        writeDataFileHeader("GradientAveraged", "GradientAveraged", gradientAveragedFilePtr_());
    }
}

void Foam::functionObjects::convergenceDetection::writeDataFileHeader(
    const word &header,
    const word &tabbedName,
    OFstream &os) const
{
    writeHeader(os, header);
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, tabbedName);

    os << endl;
}

void Foam::functionObjects::convergenceDetection::writeDataFile(
    OFstream &os, double &value) const
{
    writeCurrentTime(os);

    writeValue(os, value);

    os << endl;
}

bool Foam::functionObjects::convergenceDetection::reachedMaxIterations()
{
    return currentIteration_ >= (iterationMaxConvergence_ + iterationMaxAveraging_);
}

bool Foam::functionObjects::convergenceDetection::minIterationsForAveraging()
{
    return currentIteration_ >= static_cast<int>(1.1 * averagingStartedAt_) &&
           currentIteration_ >= static_cast<int>(iterationMinAveraging_ + averagingStartedAt_);
}

bool Foam::functionObjects::convergenceDetection::isAveraged()
{
    std::map<std::string, double> amg = averagingMaxGradient();

    bool allDataAveragedEnough = std::all_of(amg.begin(), amg.end(), [this](const auto &amg)
                                             { return amg.second < thresholdAveraging_; });

    bool averagingMaxGradientGtZero = std::all_of(amg.begin(), amg.end(), [](const auto &amg)
                                                  { return amg.second > 0.0; });
    return allDataAveragedEnough &&
           averagingMaxGradientGtZero &&
           currentIteration_ > static_cast<int>(1.25 * averagingStartedAt_) &&
           !simulationFinished_;
}

void Foam::functionObjects::convergenceDetection::isFinished()
{
    if ((isAveraged() ||
         reachedMaxIterations()) &&
        minIterationsForAveraging())
    {
        std::map<std::string, double> amg = averagingMaxGradient();
        Info << "###########" << endl;
        Info << "Simulation should stop!!!!" << endl;
        for (auto const &a : amg)
        {
            Info << "Gradient of " << a.first << ": " << a.second << endl;
        }
        Info << "Condition for averaging: " << thresholdAveraging_ << endl;
        Info << "###########" << endl;

        time().stopAt(Time::saWriteNow);
        toggleAveraging(false);
        simulationFinished_ = true;
    }
}

void Foam::functionObjects::convergenceDetection::stopAveraging()
{
    convergenceFound_ = false;
    toggleAveraging(false);
    gradDataAveraging_.clear();
    if (exists(time().globalPath() + "/averaging"))
    {
        rm(time().globalPath() + "/averaging");
    }
    Info << "## Getting out of averaging ###" << endl;
}

void Foam::functionObjects::convergenceDetection::isExploded()
{
    bool anyNormalizedMeanDataGtZero = std::any_of(normalizedDataMeanConverged_.begin(), normalizedDataMeanConverged_.end(), [](const auto &ndmc)
                                                   { return ndmc.second > 0; });

    if (anyNormalizedMeanDataGtZero && convergenceFound_)
    {
        double test1 = normalizedDataMeanConverged_["forces"] / convergenceStabilityFactor_;
        double test2 = normalizedDataMeanConverged_["forces"] * convergenceStabilityFactor_;
        double lastIterationForces = convergenceData_["forces"].back();

        if (lastIterationForces < test1 or lastIterationForces > test2)
        {
            if (Pstream::master())
            {
                FatalErrorInFunction
                    << "Forces Exploded, mean force value converged: " << lastIterationForces << nl
                    << exit(FatalError);
            }
        }
    }
}

void Foam::functionObjects::convergenceDetection::convergenceEvaluation()
{
    std::map<std::string, double> cmg(convergenceMaxGradient());

    bool allDataConverged = std::all_of(cmg.begin(), cmg.end(), [this](const auto &cmg)
                                        { return cmg.second < thresholdConvergence_; });

    if ((allDataConverged &&
         currentIteration_ >= iterationMinConvergence_ &&
         !convergenceFound_) ||
        (currentIteration_ >= iterationMaxConvergence_ &&
         !convergenceFound_))
    {
        convergenceFound_ = true;
        if (currentIteration_ >= iterationMaxConvergence_)
        {
            forcedConvergence_ = true;
        }
        word forcedConvergence = forcedConvergence_ ? "Forced " : "";
        Info << "#####################################" << endl;
        Info << forcedConvergence << "Convergence Found" << endl;
        Info << "#####################################" << endl;
        for (auto const &c : cmg)
        {
            Info << "Gradient of " << c.first << ": " << c.second << endl;
        }
        Info << "Condition for convergence: " << thresholdConvergence_ << endl;
        for (auto const &c : cmg)
        {
            Info << "Condition for " << c.first << "is : " << (c.second < thresholdConvergence_) << endl;
        }

        averagingStartedAt_ = currentIteration_;

        int windowDataNormalization = static_cast<int>(currentIteration_ * windowNormalization_);

        std::map<std::string, std::vector<double>> normalizedData;

        for (auto const &d : convergenceData_)
        {
            normalizedData[d.first].insert(normalizedData[d.first].end(), d.second.end() - windowDataNormalization, d.second.end());
        }

        normalizedDataMeanConverged_ = meanValue(normalizedData);

        OFstream os(time().globalPath() + "/averaging");
        os << averagingStartedAt_ << endl;

        toggleAveraging(true);
    }
}

void Foam::functionObjects::convergenceDetection::toggleAveraging(bool toggle)
{
    IOdictionary averaging(IOobject(
        "averaging",
        time().caseSystem(),
        obr_,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE));
    averaging.subDict("fieldAverage1").set("timeStart", toggle ? currentIteration_ : 0);
    averaging.subDict("fieldAverage1").set("enabled", toggle);
    averaging.regIOobject::write();
}

std::map<std::string, std::vector<double>> Foam::functionObjects::convergenceDetection::divideData(std::map<std::string, double> dataMean)
{
    std::map<std::string, std::vector<double>> dividedData;

    for (const auto &cd : convergenceData_)
    {
        for (auto itr : cd.second)
            dividedData[cd.first].push_back(itr / dataMean[cd.first]);
    }

    return dividedData;
}

std::map<std::string, double> Foam::functionObjects::convergenceDetection::meanValue(std::map<std::string, std::vector<double>> normalizedData)
{
    std::map<std::string, double> meanData;

    for (const auto &nd : normalizedData)
    {
        double sum = 0.0;
        for (auto itr : nd.second)
            sum += itr;
        meanData[nd.first] = sum / nd.second.size();
    }

    return meanData; // for transient change this
}

std::map<std::string, double> Foam::functionObjects::convergenceDetection::calculateGradientAveraging()
{
    int averagingSize = currentIteration_ - averagingStartedAt_;
    std::map<std::string, double> gradient;

    // allow some data to come in
    if (averagingSize >= 5)
    {
        double maxIterationNormalized = static_cast<double>(currentIteration_) / static_cast<double>(averagingStartedAt_);

        convergenceNormalizedAveraging_ = divideData(normalizedDataMeanConverged_);

        double start = 1 + (maxIterationNormalized / averagingSize);
        double end = maxIterationNormalized + maxIterationNormalized / averagingSize;
        double step = (maxIterationNormalized - 1) / averagingSize;

        std::vector<double> xAxisNormalized = arange(start, end, step);

        std::map<std::string, std::vector<double>> windowDataForGradient;

        int windowDataNormalization = static_cast<int>(xAxisNormalized.size() * windowGradientAveraging_);

        for (auto const &c : convergenceNormalizedAveraging_)
        {
            windowDataForGradient[c.first].insert(windowDataForGradient[c.first].end(), c.second.end() - windowDataNormalization, c.second.end());
        }

        std::vector<double> windowXAxisForGradient(xAxisNormalized.end() - windowDataNormalization, xAxisNormalized.end());

        if (windowDataForGradient["forces"].size() < 5)
        {
            for (const auto &d : windowDataForGradient)
            {
                gradient[d.first] = 0.0;
            }
            return gradient;
        }

        for (const auto &d : windowDataForGradient)
        {
            gradient[d.first] = polyfit(windowXAxisForGradient, d.second, 1, {0});
        }

        return gradient;
    }
    for (const auto &n : normalizedDataMeanConverged_)
    {
        gradient[n.first] = 0.0;
    }
    return gradient;
}

std::map<std::string, double> Foam::functionObjects::convergenceDetection::calculateGradient()
{
    int windowDataNormalization = static_cast<int>(currentIteration_ * windowNormalization_);
    int windowDataGradient = static_cast<int>(currentIteration_ * windowGradient_);

    std::map<std::string, std::vector<double>> normalizedData;

    for (auto const &d : convergenceData_)
    {
        normalizedData[d.first].insert(normalizedData[d.first].end(), d.second.end() - windowDataNormalization, d.second.end());
    }

    std::map<std::string, double> dataMean = meanValue(normalizedData);

    dataNormalized_ = divideData(dataMean);

    // take this into account for transient
    double start = 1.0f / currentIteration_;
    double end = 1.0f + 1.0f / currentIteration_;
    double step = 1.0f / currentIteration_;

    std::vector<double> xAxisNormalized = arange(start, end, step);

    std::map<std::string, std::vector<double>> windowDataForGradient;

    for (auto const &d : dataNormalized_)
    {
        windowDataForGradient[d.first].insert(windowDataForGradient[d.first].end(), d.second.end() - windowDataGradient, d.second.end());
    }

    std::vector<double> windowXAxisForGradient(xAxisNormalized.end() - windowDataGradient, xAxisNormalized.end());

    std::map<std::string, double> gradient;

    for (const auto &data : windowDataForGradient)
    {
        gradient[data.first] = polyfit(windowXAxisForGradient, data.second, 1, {0});
    }

    return gradient;
}

// DRY
std::map<std::string, double> Foam::functionObjects::convergenceDetection::convergenceMaxGradient()
{
    int evaluationWindow = static_cast<int>(windowEvaluation_ * gradData_["forces"].size());
    std::map<std::string, std::vector<double>> dataForEvaluation;
    std::map<std::string, double> maxGradient;
    std::map<std::string, double> defaultGradient;

    for (const auto &gd : gradData_)
    {
        dataForEvaluation[gd.first].insert(dataForEvaluation[gd.first].end(), gd.second.end() - evaluationWindow, gd.second.end());
    }

    if (dataForEvaluation["forces"].size() > 0)
    {
        for (const auto &afe : dataForEvaluation)
        {
            double maxValue = *std::max_element(afe.second.begin(), afe.second.end());
            double minValue = *std::min_element(afe.second.begin(), afe.second.end());

            maxGradient[afe.first] = max(maxValue, fabs(minValue));
        }

        return maxGradient;
    }
    for (const auto &gd : gradData_)
    {
        defaultGradient[gd.first] = 1.0;
    }

    return defaultGradient;
}

// DRY
std::map<std::string, double> Foam::functionObjects::convergenceDetection::averagingMaxGradient()
{
    int evaluationWindow = static_cast<int>(
        windowEvaluationAveraging_ * gradDataAveraging_["forces"].size());

    std::map<std::string, std::vector<double>> dataForEvaluation;
    std::map<std::string, double> maxGradient;
    std::map<std::string, double> defaultGradient;

    for (const auto &gd : gradDataAveraging_)
    {
        dataForEvaluation[gd.first].insert(dataForEvaluation[gd.first].end(), gd.second.end() - evaluationWindow, gd.second.end());
    }

    if (dataForEvaluation["forces"].size() > 0)
    {
        for (const auto &afe : dataForEvaluation)
        {
            double maxValue = *std::max_element(afe.second.begin(), afe.second.end());
            double minValue = *std::min_element(afe.second.begin(), afe.second.end());

            maxGradient[afe.first] = max(maxValue, fabs(minValue));
        }

        return maxGradient;
    }
    for (const auto &gd : gradDataAveraging_)
    {
        defaultGradient[gd.first] = 1.0;
    }

    return defaultGradient;
}

std::vector<double> Foam::functionObjects::convergenceDetection::arange(double start, double stop, double step)
{
    std::vector<double> values;
    for (double value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}
// https://gist.github.com/chrisengelsma/108f7ab0a746323beaaf7d6634cf4add
double Foam::functionObjects::convergenceDetection::polyfit(
    const std::vector<double> x,
    const std::vector<double> y,
    const int order,
    std::vector<double> coeffs)
{
    // The size of xValues and yValues should be same
    if (x.size() != y.size())
    {
        throw std::runtime_error("The size of x & y arrays are different");
        return {};
    }
    // The size of xValues and yValues cannot be 0, should not happen
    if (x.size() == 0 || y.size() == 0)
    {
        throw std::runtime_error("The size of x or y arrays is 0");
        return {};
    }

    size_t N = x.size();
    int n = order;
    int np1 = n + 1;
    int np2 = n + 2;
    int tnp1 = 2 * n + 1;
    double tmp;

    // X = vector that stores values of sigma(xi^2n)
    std::vector<double> X(tnp1);
    for (int i = 0; i < tnp1; ++i)
    {
        X[i] = 0;
        for (size_t j = 0; j < N; ++j)
            X[i] += static_cast<double>(pow(x[j], i));
    }

    // a = vector to store final coefficients.
    std::vector<double> a(np1);

    // B = normal augmented matrix that stores the equations.
    std::vector<std::vector<double>> B(np1, std::vector<double>(np2, 0));

    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            B[i][j] = X[i + j];

    // Y = vector to store values of sigma(xi^n * yi)
    std::vector<double> Y(np1);
    for (int i = 0; i < np1; ++i)
    {
        Y[i] = static_cast<double>(0);
        for (size_t j = 0; j < N; ++j)
        {
            Y[i] += static_cast<double>(pow(x[j], i) * y[j]);
        }
    }

    // Load values of Y as last column of B
    for (int i = 0; i <= n; ++i)
        B[i][np1] = Y[i];

    n += 1;
    int nm1 = n - 1;

    // Pivotisation of the B matrix.
    for (int i = 0; i < n; ++i)
        for (int k = i + 1; k < n; ++k)
            if (B[i][i] < B[k][i])
                for (int j = 0; j <= n; ++j)
                {
                    tmp = B[i][j];
                    B[i][j] = B[k][j];
                    B[k][j] = tmp;
                }

    // Performs the Gaussian elimination.
    // (1) Make all elements below the pivot equals to zero
    //     or eliminate the variable.
    for (int i = 0; i < nm1; ++i)
        for (int k = i + 1; k < n; ++k)
        {
            double t = B[k][i] / B[i][i];
            for (int j = 0; j <= n; ++j)
                B[k][j] -= t * B[i][j]; // (1)
        }

    // Back substitution.
    // (1) Set the variable as the rhs of last equation
    // (2) Subtract all lhs values except the target coefficient.
    // (3) Divide rhs by coefficient of variable being calculated.
    for (int i = nm1; i >= 0; --i)
    {
        a[i] = B[i][n]; // (1)
        for (int j = 0; j < n; ++j)
            if (j != i)
                a[i] -= B[i][j] * a[j]; // (2)
        a[i] /= B[i][i];                // (3)
    }

    coeffs.resize(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        coeffs[i] = a[i];

    return static_cast<double>(coeffs[1]);
}

// ************************************************************************* //