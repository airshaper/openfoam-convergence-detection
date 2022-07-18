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
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"
#include <vector>
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    static std::vector<double> arange(double start, double stop, double step)
    {
        std::vector<double> values;
        for (double value = start; value < stop; value += step)
            values.push_back(value);
        return values;
    }
    // https://gist.github.com/chrisengelsma/108f7ab0a746323beaaf7d6634cf4add
    static std::vector<double> polyfit(
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
            for (int j = 0; j < N; ++j)
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
            for (int j = 0; j < N; ++j)
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
      windowNormalization_(dict.getOrDefault<scalar>("windowNormalization", 0.5)),
      windowPolynom_(dict.getOrDefault<scalar>("windowPolynom", 0.5)),
      windowPolynomAveraging_(dict.getOrDefault<scalar>("windowPolynomAveraging", 0.99)),
      windowEvaluation_(dict.getOrDefault<scalar>("windowEvaluation", 0.333)),
      windowEvaluationAveraging_(dict.getOrDefault<scalar>("windowEvaluationAveraging", 0.333)),
      conditionConvergence_(dict.getOrDefault<scalar>("conditionConvergence", 0.075)),
      conditionAveraging_(dict.getOrDefault<scalar>("conditionAveraging", 0.025)),
      maxStepConvergence_(dict.getOrDefault<scalar>("maxStepConvergence", 2500)),
      maxStepAveraging_(dict.getOrDefault<scalar>("maxStepAveraging", 2500)),
      forceStabilityFactor_(dict.getOrDefault<scalar>("forceStabilityFactor", 2)),
      convergenceMinIter_(dict.getOrDefault<scalar>("convergenceMinIter", 200)),
      averagingMinIter_(dict.getOrDefault<scalar>("averagingMinIter", 200)),
      convergenceFound_(false),
      averagingStartedAt_(0),
      normalizedForcesMeanConverged_(0),
      forces_data_(0),
      polyVector_(0),
      polyVectorAveraging_(),
      forcesDict_(dictionary(IFstream("system/forces")()))
//   averagingDict_(dictionary(IFstream("system/averaging")())),

{
    read(dict);
}

// * * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::convergenceDetection::read(const dictionary &dict)
{
    // dict.readEntry("windowNormalization", windowNormalization_);

    // Info << "Window: " << windowNormalization_ << endl;
    //  dict.readEntry("labelData", labelData_);
    //  dict.readIfPresent("wordData", wordData_);
    //  dict.readEntry("scalarData", scalarData_);

    return true;
}

bool Foam::functionObjects::convergenceDetection::execute()
{
    forces f("forces", mesh_.time(), forcesDict_.subDict("forces"));

    f.calcForcesMoments();

    const vector sumForce = f.forceEff();

    const double totalForce = sqrt(pow(sumForce[0], 2) + pow(sumForce[1], 2) + pow(sumForce[2], 2));

    forces_data_.push_back(totalForce);

    // need few iterations before starting calculation
    if (mesh_.time().timeToUserTime(mesh_.time().value()) >= 5)
    {
        double caclculatedPolynomGrad = calculatePolynomGrad();

        polyVector_.push_back(caclculatedPolynomGrad);

        if (checkCriteriaForConvergence() < conditionConvergence_ && forces_data_.size() > convergenceMinIter_ && !convergenceFound_)
        {
            convergenceFound_ = true;
            Info << "###########" << endl;
            Info << "Convergence Found" << endl;
            Info << "###########" << endl;

            averagingStartedAt_ = forces_data_.size();

            turnOnAveraging();
        }

        int normalizationForcesWindow = static_cast<int>(forces_data_.size() * windowNormalization_);

        std::vector<double> normalizedForces(forces_data_.end() - normalizationForcesWindow, forces_data_.end());

        normalizedForcesMeanConverged_ = meanValue(normalizedForces);

        if (normalizedForcesMeanConverged_ > 0 && convergenceFound_)
        {

            double test1 = normalizedForcesMeanConverged_ / forceStabilityFactor_;
            double test2 = normalizedForcesMeanConverged_ * forceStabilityFactor_;
            double lastIterationForces = forces_data_.back();

            if (lastIterationForces < test1 or lastIterationForces > test2)
            {
                FatalErrorInFunction
                    << "Forces Exploded, mean force converged" << lastIterationForces << nl
                    << abort(FatalError);
            }

            double caclculatedPolynomGradAveraging = calculatePolynomGradAveraging();

            polyVectorAveraging_.push_back(caclculatedPolynomGradAveraging);
        }

        if (convergenceFound_)
        {
            Info << "checkCriteriaForAveraging" << endl;
            if (checkCriteriaForAveraging() < conditionAveraging_ && checkCriteriaForAveraging() > 0.0 and forces_data_.size() > static_cast<int>(1.25 * averagingStartedAt_))
            {
                Info << "###########" << endl;
                Info << "Simulation should stop!!!!" << endl;
                Info << "###########" << endl;
            }
        }
    }

    return true;
}

void Foam::functionObjects::convergenceDetection::turnOnAveraging()
{
    IOdictionary averaging(IOobject(
        "averaging",
        mesh_.time().caseSystem(),
        mesh_,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE));
    averaging.subDict("fieldAverage1").set("timeStart", forces_data_.size());
    averaging.subDict("fieldAverage1").set("enabled", true);
    averaging.regIOobject::write();
}

std::vector<double> Foam::functionObjects::convergenceDetection::divideForcesByScalar(double d)
{
    std::vector<double> dividedForces = {0};
    for (auto itr : forces_data_)
        dividedForces.push_back(itr / d);

    return dividedForces;
}

double Foam::functionObjects::convergenceDetection::meanValue(std::vector<double> normalizedForces)
{
    double sum = 0.0;
    for (auto itr : normalizedForces)
        sum += itr;

    return sum / normalizedForces.size();
}

double Foam::functionObjects::convergenceDetection::calculatePolynomGradAveraging()
{
    double averagingSize = forces_data_.size() / averagingStartedAt_;

    // allow some data to come in
    if (averagingSize > 20)
    {
        double xMax = forces_data_.size() / averagingStartedAt_;
        std::vector<double> forcesNormalized = divideForcesByScalar(normalizedForcesMeanConverged_);

        std::vector<double>
            xAxisGeneral = arange(
                1 + (xMax / averagingSize),
                (xMax + xMax / averagingSize),
                (xMax / averagingSize));

        int windowForcesPolynom = xAxisGeneral.size() * windowPolynomAveraging_;

        std::vector<double> normalizedForces(forcesNormalized.end() - windowForcesPolynom, forcesNormalized.end());

        std::vector<double> xAxisPolynom(xAxisGeneral.size() - windowForcesPolynom, xAxisGeneral.size());

        if (xAxisPolynom.size() < 5)
        {
            return 0.0;
        }

        std::vector<double> polynom = polyfit(xAxisPolynom, normalizedForces, 1, {0});

        return static_cast<double>(fabs(polynom[1]));
    }
    return 0.0;
}

double
Foam::functionObjects::convergenceDetection::calculatePolynomGrad()
{
    int normalizationForcesWindow = static_cast<int>(forces_data_.size() * windowNormalization_);
    int windowForcesPolynom = static_cast<int>(forces_data_.size() * windowPolynom_);

    std::vector<double> normalizedForces(forces_data_.end() - normalizationForcesWindow, forces_data_.end());

    double forcesMean = meanValue(normalizedForces);

    std::vector<double> forcesNormalized = divideForcesByScalar(forcesMean);

    double start = 1.0f / forces_data_.size();
    double end = 1.0f + 1.0f / forces_data_.size();
    double step = 1.0f / forces_data_.size();

    // Info << "Start: " << start << " End: " << end << " Step: " << step << endl;

    std::vector<double> xAxisGeneral = arange(start, end, step);

    std::vector<double> polynomForces(forcesNormalized.end() - windowForcesPolynom, forcesNormalized.end());

    std::vector<double> xAxisPolynom(xAxisGeneral.end() - windowForcesPolynom, xAxisGeneral.end());

    // Info << "xAxisPolynom: " << xAxisPolynom << endl;
    // Info << "polynomForces: " << polynomForces << endl;

    std::vector<double> polynom = polyfit(xAxisPolynom, polynomForces, 1, {0});

    return static_cast<double>(fabs(polynom[1]));
}

// DRY
double Foam::functionObjects::convergenceDetection::checkCriteriaForConvergence()
{
    int evaluationWindow = static_cast<int>(windowEvaluation_ * polyVector_.size());
    std::vector<double> vectorForEvaluation(polyVector_.end() - evaluationWindow, polyVector_.end());

    if (vectorForEvaluation.size() > 0)
    {
        return *std::max_element(vectorForEvaluation.begin(), vectorForEvaluation.end());
    }
    return 0.0;
}

// DRY
double Foam::functionObjects::convergenceDetection::checkCriteriaForAveraging()
{
    int evaluationWindow = static_cast<int>(
        windowEvaluationAveraging_ * polyVectorAveraging_.size());

    std::vector<double> vectorForEvaluation(polyVectorAveraging_.size() - evaluationWindow, polyVectorAveraging_.size());

    if (vectorForEvaluation.size() > 0)
    {
        return *std::max_element(vectorForEvaluation.begin(), vectorForEvaluation.end());
    }
    return 0.0;
}

bool Foam::functionObjects::convergenceDetection::end()
{
    return true;
}

bool Foam::functionObjects::convergenceDetection::write()
{
    return true;
}

// ************************************************************************* //