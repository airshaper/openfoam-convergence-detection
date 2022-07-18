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

std::vector<float> linspace(float start, float end, float step)
{
    vector<float> range;
    float delta = (end-start)/float(step-1);
    for(int i=0; i<step; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}

static std::vector<double> divideVectorByScalar(std::vector<double> &forces_, double d) 
{
    std::vector<double> dividedForces = {0};
    for (auto itr : forces_) {
         dividedForces.push_back(itr / d);
    }

    return dividedForces;
}

static double meanValue(std::vector<double> forces_)

{
    double sum = 0.0;
    for (auto itr : forces_) {
         sum += itr;
    }

    return sum / forces_.size();
}

static double calculate_polynom_grad
(
    std::vector<double> forces_, 
    int normalizationForcesWindow_, 
    double windowPolynom_
)
{
    std::vector<double> normalizedForces(forces_.end() - normalizationForcesWindow_, forces_.end());

    double forcesMean = meanValue(normalizedForces);

    std::vector<double> forcesNormalized = divideVectorByScalar(forces_, forcesMean);

    Info << "Forces size: " << forces_.size() << endl;

    float start = 1.0/forces_.size();
    float end = 2.0/forces_.size();
    float step = 1.0/forces_.size();

    Info << "Start: " << start << " End: " << end << " Step: " << step << endl;

    std::vector<float> v = linspace(start, end, step);
    // auto L = xAxisGeneral.length();

    // std::vector<double> normForces = normalizedForces;

    // 

    Info << "Axis: " << v << endl;

    // std::vector<float> v(100);

    //std::iota(v.begin(), v.end(), spacer(1.132))
    

    // Info << v.begin() << v.end() << endl;

    // std::vector<double> forcesNormalized(forces_.end() - windowPolynom_, forces_.end());

    // Info << normalizedForces << endl;
    
    return 0.0;
}
namespace functionObjects
{
    defineTypeNameAndDebug(convergenceDetection, 0);
    addToRunTimeSelectionTable(functionObject, convergenceDetection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::convergenceDetection::convergenceDetection
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces(name, runTime, dict, false),
    windowNormalization_(dict.getOrDefault<scalar>("windowNormalization", 0.5)),
    windowPolynom_(dict.getOrDefault<scalar>("windowPolynom", 0.5)),
    windowPolynomAveraging_(dict.getOrDefault<scalar>("windowPolynomAveraging", 0.99)),
    windowEvaluation_(dict.getOrDefault<scalar>("windowEvaluation", 0.333)),
    windowEvaluationAveraging_(dict.getOrDefault<scalar>("windowEvaluationAveraging", 0.333)),
    conditionConvergence_(dict.getOrDefault<scalar>("conditionConvergence", 0.075)),
    conditionAveraging_(dict.getOrDefault<scalar>("conditionAveraging", 0.025)),
    maxStepConvergence_(dict.getOrDefault<scalar>("maxStepConvergence", 2500)),
    maxStepAveraging_(dict.getOrDefault<scalar>("maxStepAveraging", 2500)),
    forceStabilityDactor_(dict.getOrDefault<scalar>("forceStabilityDactor", 2)),
    convergenceMinIter_(dict.getOrDefault<scalar>("convergenceMinIter", 200)),
    averagingMinIter_(dict.getOrDefault<scalar>("averagingMinIter", 200)),
    forces_data_(0),
    polyVector_(0),
    forcesDict_(dictionary(IFstream("system/forces")()))
{
    read(dict);}

// * * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::convergenceDetection::read(const dictionary& dict)
{
    //dict.readEntry("windowNormalization", windowNormalization_);

    //Info << "Window: " << windowNormalization_ << endl;
    // dict.readEntry("labelData", labelData_);
    // dict.readIfPresent("wordData", wordData_);
    // dict.readEntry("scalarData", scalarData_);

    return true;
}


bool Foam::functionObjects::convergenceDetection::execute()
{
    forces f("forces", mesh_.time(), forcesDict_.subDict("forces"));

    f.calcForcesMoments();

    const vector sumForce = f.forceEff();

    const scalar totalForce = sqrt(pow(sumForce[0], 2) + pow(sumForce[1], 2) + pow(sumForce[2], 2));

    forces_data_.push_back(totalForce);

    int normalizationForcesWindow = static_cast<int>(forces_data_.size() * windowNormalization_);
    int windowForcesPolynom = int(forces_data_.size() * windowPolynom_);


    if (mesh_.time().timeToUserTime(mesh_.time().value()) > 1) {
        double caclculated = calculate_polynom_grad(forces_data_, normalizationForcesWindow, windowForcesPolynom);
        // Info << "caclculated: " << caclculated << endl;
    }

    
    return true;
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
