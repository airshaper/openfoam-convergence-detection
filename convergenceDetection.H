/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

Class
    Foam::functionObjects::convergenceDetection

Group
    grpFieldFunctionObjects

Description
    <minimal description of the function object>

    <equation>
    \f[
        x = x_{ref}^x + \rho \omega
    \f]

    <variable-explanation table>
    where
    \vartable
      \rho        | <explanation> [units, e.g. kg/m3]
      \omega      | \f$ \nabla \cdot \vec U \f$
      ...         | ...
    \endvartable

    <inline equation>
    where \f$ x_k \f$ is ...

    <input-output table>
    \table
      Operand      | Type                 | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
               --> |$FOAM_CASE/\<time\>/\<inpField\>s
      output file  | dat                  <!--
              -->  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field | volScalarField       | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    convergenceDetection1
    {
        // Mandatory entries (unmodifiable)
        type           convergenceDetection;
        libs           (convergenceDetectionFunctionObject);

        // Mandatory entries (runtime modifiable)
        ...

        // Mandatory (inherited) entries (unmodifiable)
        ...

        // Mandatory (inherited) entries (runtime unmodifiable)
        ...

        // Optional entries (unmodifiable)
        ...

        // Optional entries (runtime modifiable)
        boolData       <bool>;
        labelData      <label>;
        wordData       <word>;
        scalarData     <scalar>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                      | Type | Req'd | Dflt
      type       | Type name: convergenceDetection        | word |  yes  | -
      libs       | Library name: convergenceDetectionFunctionObject  <!--
             --> | word  | yes | -
      boolData   | <explanation>                    | bool |  yes  | -
      labelData  | <explanation>                    | label |  yes | -
      wordData   | <explanation>                    | word |  yes  | -
      scalarData | <explanation>                    | scalar |  no | 1.0
      wordListData  | <explanation>                 | wordList | yes | -
    \endtable

    Options for the \c ENTRY entry:
    \verbatim
        <option1>
        <option2> | <explanation>
      ...
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link fieldsExpression.H \endlink
     - \link writeFile.H \endlink
     ...

    <if \c postProcess is applicable>
    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func convergenceDetection
    \endverbatim

    <if \c postProcess is not applicable>
    Usage by the \c postProcess utility is not available.

Note
    - <note1>
    - <note2>
    ...

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::convergenceDetection
    ...

SourceFiles
    convergenceDetection.C
    convergenceDetectionTEMPLATES.C
    ...

\*---------------------------------------------------------------------------*/

#ifndef convergenceDetection_H
#define convergenceDetection_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {

        /*---------------------------------------------------------------------------*\
                           Class convergenceDetection Declaration
        \*---------------------------------------------------------------------------*/

        class convergenceDetection
            : public forces
        {

            scalar windowNormalization_;
            scalar windowPolynomAveraging_;
            scalar windowPolynom_;
            scalar windowEvaluation_;
            scalar windowEvaluationAveraging_;
            scalar conditionConvergence_;
            scalar conditionAveraging_;
            scalar maxStepConvergence_;
            scalar maxStepAveraging_;
            scalar forceStabilityFactor_;
            scalar convergenceMinIter_;
            scalar averagingMinIter_;
            bool convergenceFound_;
            bool simulationFinished_;
            int currentIteration_;
            int averagingStartedAt_;
            scalar normalizedForcesMeanConverged_;
            std::vector<double>
                forces_data_;
            std::vector<double> polyVector_;
            std::vector<double> polyVectorAveraging_;
            dictionary forcesDict_;
            // dictionary averagingDict_;

        public:
            //- Runtime type information
            TypeName("convergenceDetection");

            // Constructors

            //- Construct from Time and dictionary
            convergenceDetection(
                const word &name,
                const Time &runTime,
                const dictionary &dict);

            //- No copy construct
            convergenceDetection(const convergenceDetection &) = delete;

            //- No copy assignment
            void operator=(const convergenceDetection &) = delete;

            //- Destructor
            virtual ~convergenceDetection() = default;

            // Member Functions

            //- Read the convergenceDetection data
            virtual bool read(const dictionary &dict);

            //- Execute, currently does nothing
            virtual bool execute();

            //- Execute at the final time-loop, currently does nothing
            virtual bool end();

            //- Write the convergenceDetection
            virtual bool
            write();

        protected:
            // virtual double meanValue(std::vector<double>& forces_);
            // virtual double calculate_polynom_grad(std::vector<double> forces_);

            virtual std::vector<double> divideForcesByScalar(double d);

            virtual double meanValue(std::vector<double> normalizedForces);

            virtual double calculatePolynomGrad();
            virtual double calculatePolynomGradAveraging();

            virtual double checkCriteriaForConvergence();
            virtual double checkCriteriaForAveraging();

            virtual void
            turnOnAveraging();
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "convergenceDetectionTemplates.C"
// #endif

#endif

// ************************************************************************* //
