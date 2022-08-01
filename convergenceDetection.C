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
      windowPolynom_(0.5),
      windowPolynomAveraging_(0.99),
      windowEvaluation_(0.333),
      windowEvaluationAveraging_(0.333),
      conditionConvergence_(0.075),
      conditionAveraging_(0.025),
      maxStepConvergence_(4000),
      maxStepAveraging_(2000),
      forceStabilityFactor_(2),
      convergenceMinIter_(200),
      averagingMinIter_(200),
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
      polyVector_(0),
      polyVectorAveraging_(),
      forcesNormalized_(),
      forcesNormalizedAveraging_(),
      totalForceFilePtr_(nullptr),
      polynomGradFilePtr_(nullptr),
      polynomGradAveragedFilePtr_(nullptr)

{
    read(dict);
    setCoordinateSystem(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::convergenceDetection::read(const dictionary &dict)
{

    forces::read(dict);

    dict.readIfPresent("windowNormalization", windowNormalization_);
    dict.readIfPresent("windowPolynom", windowPolynom_);
    dict.readIfPresent("windowPolynomAveraging", windowPolynomAveraging_);
    dict.readIfPresent("windowEvaluation", windowEvaluation_);
    dict.readIfPresent("windowEvaluationAveraging", windowEvaluationAveraging_);
    dict.readIfPresent("conditionConvergence", conditionConvergence_);
    dict.readIfPresent("conditionAveraging", conditionAveraging_);
    dict.readIfPresent("maxStepConvergence", maxStepConvergence_);
    dict.readIfPresent("maxStepAveraging", maxStepAveraging_);
    dict.readIfPresent("forceStabilityFactor", forceStabilityFactor_);
    dict.readIfPresent("convergenceMinIter", convergenceMinIter_);
    dict.readIfPresent("averagingMinIter", averagingMinIter_);

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

    // wait for a few iterations before starting polynom grad calculation
    if (currentIteration_ >= 5)
    {
        double caclculatedPolynomGrad = calculatePolynomGrad();
        polyVector_.push_back(caclculatedPolynomGrad);
    }

    // need few iterations before starting calculation
    if (currentIteration_ >= convergenceMinIter_)
    {
        if (!convergenceFound_)
        {
            checkConvergence();
        }
        // Start to check at averagingMinIter to avoid calculate each iteration
        if (currentIteration_ >= averagingMinIter_ &&
            !simulationFinished_ &&
            convergenceFound_)
        {
            if (checkCriteriaForConvergence() >= conditionConvergence_ && !forcedConvergence_)
            {
                stopAveraging();
            }
            checkIfForcesExploded();
            checkIfFinished();
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
            writeDataFile(polynomGradFilePtr_(), polyVector_.back());
        }
        if (convergenceFound_)
        {
            writeDataFile(polynomGradAveragedFilePtr_(), polyVectorAveraging_.back());
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

    if (!polynomGradFilePtr_.valid())
    {
        polynomGradFilePtr_ = createFile("polynomGrad");
        writeDataFileHeader("polynomGrad", "polynomGrad", polynomGradFilePtr_());
    }

    if (!polynomGradAveragedFilePtr_.valid())
    {
        polynomGradAveragedFilePtr_ = createFile("polynomGradAveraged");
        writeDataFileHeader("polynomGradAveraged", "polynomGradAveraged", polynomGradAveragedFilePtr_());
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
    return currentIteration_ >= (maxStepConvergence_ + maxStepAveraging_);
}

bool Foam::functionObjects::convergenceDetection::minIterationsForAveraging()
{
    return currentIteration_ >= static_cast<int>(1.1 * averagingStartedAt_) &&
           currentIteration_ >= static_cast<int>(averagingMinIter_ + averagingStartedAt_);
}

bool Foam::functionObjects::convergenceDetection::checkAveragingCriteria()
{
    return checkCriteriaForAveraging() < conditionAveraging_ &&
           checkCriteriaForAveraging() > 0.0 &&
           currentIteration_ > static_cast<int>(1.25 * averagingStartedAt_) &&
           !simulationFinished_;
}

void Foam::functionObjects::convergenceDetection::checkIfFinished()
{
    double caclculatedPolynomGradAveraging = calculatePolynomGradAveraging();

    polyVectorAveraging_.push_back(caclculatedPolynomGradAveraging);

    if (checkAveragingCriteria() ||
        reachedMaxIterations())
    {
        if (minIterationsForAveraging())
        {
            Info << "###########" << endl;
            Info << "Simulation should stop!!!!" << endl;
            Info << "Polynom Grad: " << checkCriteriaForAveraging() << endl;
            Info << "Condition for averaging: " << conditionAveraging_ << endl;
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
    polyVectorAveraging_.clear();
    if (exists(time().globalPath() + "/averaging"))
    {
        rm(time().globalPath() + "/averaging");
    }
    Info << "## Getting out of averaging ###" << endl;
}

void Foam::functionObjects::convergenceDetection::checkIfForcesExploded()
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

void Foam::functionObjects::convergenceDetection::checkConvergence()
{
    if ((checkCriteriaForConvergence() < conditionConvergence_ &&
         currentIteration_ >= convergenceMinIter_ &&
         !convergenceFound_) ||
        (currentIteration_ >= maxStepConvergence_ &&
         !convergenceFound_))
    {
        convergenceFound_ = true;
        if (currentIteration_ >= maxStepConvergence_)
        {
            forcedConvergence_ = true;
        }
        word forcedConvergence = forcedConvergence_ ? "Forced " : "";
        Info << "#####################################" << endl;
        Info << forcedConvergence << "Convergence Found" << endl;
        Info << "#####################################" << endl;
        Info << "Polynom Grad: " << checkCriteriaForConvergence() << endl;
        Info << "Condition for convergence: " << conditionConvergence_ << endl;
        Info << "Condition is: " << (checkCriteriaForConvergence() < conditionConvergence_) << endl;

        averagingStartedAt_ = currentIteration_;

        int normalizationForcesWindow = static_cast<int>(currentIteration_ * windowNormalization_);

        std::vector<double> normalizedForces(forcesData_.end() - normalizationForcesWindow, forcesData_.end());

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

double Foam::functionObjects::convergenceDetection::calculatePolynomGradAveraging()
{
    int averagingSize = currentIteration_ - averagingStartedAt_;

    // allow some data to come in
    if (averagingSize >= 5)
    {
        double xMax = static_cast<double>(currentIteration_) / static_cast<double>(averagingStartedAt_);

        forcesNormalizedAveraging_ = divideForces(normalizedForcesMeanConverged_);

        double start = 1 + (xMax / averagingSize);
        double end = xMax + xMax / averagingSize;
        double step = (xMax - 1) / averagingSize;

        std::vector<double> xAxisGeneral = arange(start, end, step);

        int windowForcesPolynom = static_cast<int>(xAxisGeneral.size() * windowPolynomAveraging_);

        std::vector<double> polynomForces(forcesNormalizedAveraging_.end() - windowForcesPolynom, forcesNormalizedAveraging_.end());

        std::vector<double> xAxisPolynom(xAxisGeneral.end() - windowForcesPolynom, xAxisGeneral.end());

        if (xAxisPolynom.size() < 5)
        {
            return 0.0;
        }

        std::vector<double> polynom = polyfit(xAxisPolynom, polynomForces, 1, {0});

        return static_cast<double>(polynom[1]);
    }
    return 0.0;
}

double Foam::functionObjects::convergenceDetection::calculatePolynomGrad()
{
    int normalizationForcesWindow = static_cast<int>(currentIteration_ * windowNormalization_);
    int windowForcesPolynom = static_cast<int>(currentIteration_ * windowPolynom_);

    // for transient get last 50% of time - not of iterations
    std::vector<double> normalizedForces(forcesData_.end() - normalizationForcesWindow, forcesData_.end());

    double forcesMean = meanValue(normalizedForces);

    forcesNormalized_ = divideForces(forcesMean);

    // take this into account for transient
    double start = 1.0f / currentIteration_;
    double end = 1.0f + 1.0f / currentIteration_;
    double step = 1.0f / currentIteration_;

    // Info << "Start: " << start << " End: " << end << " Step: " << step << endl;

    std::vector<double> xAxisGeneral = arange(start, end, step);

    std::vector<double> polynomForces(forcesNormalized_.end() - windowForcesPolynom, forcesNormalized_.end());

    std::vector<double> xAxisPolynom(xAxisGeneral.end() - windowForcesPolynom, xAxisGeneral.end());

    // Info << "xAxisPolynom: " << xAxisPolynom << endl;
    // Info << "polynomForces: " << polynomForces << endl;

    std::vector<double> polynom = polyfit(xAxisPolynom, polynomForces, 1, {0});

    return static_cast<double>(polynom[1]);
}

// DRY
double Foam::functionObjects::convergenceDetection::checkCriteriaForConvergence()
{
    int evaluationWindow = static_cast<int>(windowEvaluation_ * polyVector_.size());
    std::vector<double> vectorForEvaluation(polyVector_.end() - evaluationWindow, polyVector_.end());

    if (vectorForEvaluation.size() > 0)
    {
        double maxValue = *std::max_element(vectorForEvaluation.begin(), vectorForEvaluation.end());
        double minValue = *std::min_element(vectorForEvaluation.begin(), vectorForEvaluation.end());
        return max(maxValue, fabs(minValue));
    }
    return 1.0;
}

// DRY
double Foam::functionObjects::convergenceDetection::checkCriteriaForAveraging()
{
    int evaluationWindow = static_cast<int>(
        windowEvaluationAveraging_ * polyVectorAveraging_.size());

    std::vector<double> vectorForEvaluation(polyVectorAveraging_.end() - evaluationWindow, polyVectorAveraging_.end());

    if (vectorForEvaluation.size() > 0)
    {
        double maxValue = *std::max_element(vectorForEvaluation.begin(), vectorForEvaluation.end());
        double minValue = *std::min_element(vectorForEvaluation.begin(), vectorForEvaluation.end());
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