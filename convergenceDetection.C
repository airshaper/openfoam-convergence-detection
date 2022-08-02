/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
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
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include <vector>
#include <algorithm>

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
      forceStabilityFactor_(2),
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
      normalizedForcesMeanConverged_(0),
      forcesData_(0),
      gradArray_(0),
      gradArrayAveraging_(),
      forcesNormalized_(),
      forcesNormalizedAveraging_(),
      totalForceFilePtr_(nullptr),
      gradientFilePtr_(nullptr),
      gradientAveragedFilePtr_(nullptr)

{
    read(dict);
    setCoordinateSystem(dict);
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
    dict.readIfPresent("forceStabilityFactor", forceStabilityFactor_);
    dict.readIfPresent("iterationMinConvergence", iterationMinConvergence_);
    dict.readIfPresent("iterationMinAveraging", iterationMinAveraging_);

    /*
    if (!exists(time().caseSystem() + "/averaging"))
    {
        if (Pstream::master())
        {
            FatalError
                << "Averaging file does not exists, please create averaging file in system/averaging and include it in controlDict (#include \"averaging\")" << nl
                << "https://www.openfoam.com/documentation/guides/latest/doc/guide-fos-field-fieldAverage.html" << nl
                << abort(FatalError);

            return false;
        }
    }
    */

    return true;
}

bool Foam::functionObjects::convergenceDetection::execute()
{

    forces::execute();

    const vector sumForce = forces::forceEff();

    const double totalForce = sqrt(pow(sumForce[0], 2) + pow(sumForce[1], 2) + pow(sumForce[2], 2));

    forcesData_.push_back(totalForce);

    currentIteration_ = time().value();

    // wait for a few iterations before starting gradient calculation
    if (currentIteration_ >= 5)
    {
        double caclculatedGradient = calculateGradient();
        gradArray_.push_back(caclculatedGradient);
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
            if (convergenceMaxGradient() >= thresholdConvergence_ && !forcedConvergence_)
            {
                stopAveraging();
            }
            double caclculatedGradientAveraging = calculateGradientAveraging();
            gradArrayAveraging_.push_back(caclculatedGradientAveraging);
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
        writeDataFile(totalForceFilePtr_(), forcesData_.back());
        if (currentIteration_ >= 5)
        {
            writeDataFile(gradientFilePtr_(), gradArray_.back());
        }
        if (convergenceFound_)
        {
            writeDataFile(gradientAveragedFilePtr_(), gradArrayAveraging_.back());
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
    return averagingMaxGradient() < thresholdAveraging_ &&
           averagingMaxGradient() > 0.0 &&
           currentIteration_ > static_cast<int>(1.25 * averagingStartedAt_) &&
           !simulationFinished_;
}

void Foam::functionObjects::convergenceDetection::isFinished()
{

    if (isAveraged() ||
        reachedMaxIterations())
    {
        if (minIterationsForAveraging())
        {
            Info << "###########" << endl;
            Info << "Simulation should stop!!!!" << endl;
            Info << "Gradient: " << averagingMaxGradient() << endl;
            Info << "Condition for averaging: " << thresholdAveraging_ << endl;
            Info << "###########" << endl;

            time().stopAt(Time::saWriteNow);
            toggleAveraging(false);
            simulationFinished_ = true;
        }
    }
}

void Foam::functionObjects::convergenceDetection::stopAveraging()
{
    convergenceFound_ = false;
    toggleAveraging(false);
    gradArrayAveraging_.clear();
    if (exists(time().globalPath() + "/averaging"))
    {
        rm(time().globalPath() + "/averaging");
    }
    Info << "## Getting out of averaging ###" << endl;
}

void Foam::functionObjects::convergenceDetection::isExploded()
{
    if (normalizedForcesMeanConverged_ > 0 && convergenceFound_)
    {
        double test1 = normalizedForcesMeanConverged_ / forceStabilityFactor_;
        double test2 = normalizedForcesMeanConverged_ * forceStabilityFactor_;
        double lastIterationForces = forcesData_.back();

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
    if ((convergenceMaxGradient() < thresholdConvergence_ &&
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
        Info << "Gradient: " << convergenceMaxGradient() << endl;
        Info << "Condition for convergence: " << thresholdConvergence_ << endl;
        Info << "Condition is: " << (convergenceMaxGradient() < thresholdConvergence_) << endl;

        averagingStartedAt_ = currentIteration_;

        int windowForcesNormalization = static_cast<int>(currentIteration_ * windowNormalization_);

        std::vector<double> normalizedForces(forcesData_.end() - windowForcesNormalization, forcesData_.end());

        normalizedForcesMeanConverged_ = meanValue(normalizedForces);

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

std::vector<double> Foam::functionObjects::convergenceDetection::divideForces(double d)
{
    std::vector<double> dividedForces = {0};
    for (auto itr : forcesData_)
        dividedForces.push_back(itr / d);

    return dividedForces;
}

double Foam::functionObjects::convergenceDetection::meanValue(std::vector<double> normalizedForces)
{
    double sum = 0.0;
    for (auto itr : normalizedForces)
        sum += itr;

    return sum / normalizedForces.size(); // for transient change this
}

double Foam::functionObjects::convergenceDetection::calculateGradientAveraging()
{
    int averagingSize = currentIteration_ - averagingStartedAt_;

    // allow some data to come in
    if (averagingSize >= 5)
    {
        double maxIterationNormalized = static_cast<double>(currentIteration_) / static_cast<double>(averagingStartedAt_);

        forcesNormalizedAveraging_ = divideForces(normalizedForcesMeanConverged_);

        double start = 1 + (maxIterationNormalized / averagingSize);
        double end = maxIterationNormalized + maxIterationNormalized / averagingSize;
        double step = (maxIterationNormalized - 1) / averagingSize;

        std::vector<double> xAxisNormalized = arange(start, end, step);

        int windowForcesNormalization = static_cast<int>(xAxisNormalized.size() * windowGradientAveraging_);

        std::vector<double> windowForcesForGradient(forcesNormalizedAveraging_.end() - windowForcesNormalization, forcesNormalizedAveraging_.end());

        std::vector<double> windowXAxisForGradient(xAxisNormalized.end() - windowForcesNormalization, xAxisNormalized.end());

        if (windowXAxisForGradient.size() < 5)
        {
            return 0.0;
        }

        std::vector<double> gradient = polyfit(windowXAxisForGradient, windowForcesForGradient, 1, {0});

        return static_cast<double>(gradient[1]);
    }
    return 0.0;
}

double Foam::functionObjects::convergenceDetection::calculateGradient()
{
    int windowForcesNormalization = static_cast<int>(currentIteration_ * windowNormalization_);
    int windowForcesGradient = static_cast<int>(currentIteration_ * windowGradient_);

    // for transient get last 50% of time - not of iterations
    std::vector<double> normalizedForces(forcesData_.end() - windowForcesNormalization, forcesData_.end());

    double forcesMean = meanValue(normalizedForces);

    forcesNormalized_ = divideForces(forcesMean);

    // take this into account for transient
    double start = 1.0f / currentIteration_;
    double end = 1.0f + 1.0f / currentIteration_;
    double step = 1.0f / currentIteration_;

    // Info << "Start: " << start << " End: " << end << " Step: " << step << endl;

    std::vector<double> xAxisNormalized = arange(start, end, step);

    std::vector<double> windowForcesForGradient(forcesNormalized_.end() - windowForcesGradient, forcesNormalized_.end());

    std::vector<double> windowXAxisForGradient(xAxisNormalized.end() - windowForcesGradient, xAxisNormalized.end());

    std::vector<double> gradient = polyfit(windowXAxisForGradient, windowForcesForGradient, 1, {0});

    return static_cast<double>(gradient[1]);
}

// DRY
double Foam::functionObjects::convergenceDetection::convergenceMaxGradient()
{
    int evaluationWindow = static_cast<int>(windowEvaluation_ * gradArray_.size());
    std::vector<double> arrayForEvaluation(gradArray_.end() - evaluationWindow, gradArray_.end());

    if (arrayForEvaluation.size() > 0)
    {
        double maxValue = *std::max_element(arrayForEvaluation.begin(), arrayForEvaluation.end());
        double minValue = *std::min_element(arrayForEvaluation.begin(), arrayForEvaluation.end());
        return max(maxValue, fabs(minValue));
    }
    return 1.0;
}

// DRY
double Foam::functionObjects::convergenceDetection::averagingMaxGradient()
{
    int evaluationWindow = static_cast<int>(
        windowEvaluationAveraging_ * gradArrayAveraging_.size());

    std::vector<double> arrayForEvaluation(gradArrayAveraging_.end() - evaluationWindow, gradArrayAveraging_.end());

    if (arrayForEvaluation.size() > 0)
    {
        double maxValue = *std::max_element(arrayForEvaluation.begin(), arrayForEvaluation.end());
        double minValue = *std::min_element(arrayForEvaluation.begin(), arrayForEvaluation.end());
        return max(maxValue, fabs(minValue));
    }
    return 1.0;
}

std::vector<double> Foam::functionObjects::convergenceDetection::arange(double start, double stop, double step)
{
    std::vector<double> values;
    for (double value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}
// https://gist.github.com/chrisengelsma/108f7ab0a746323beaaf7d6634cf4add
std::vector<double> Foam::functionObjects::convergenceDetection::polyfit(
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

    return coeffs;
}

// ************************************************************************* //