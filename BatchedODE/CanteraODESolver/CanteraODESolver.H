/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::CanteraODESolver

Description
    Abstract base-class batch-solution of ODE systems

SourceFiles
    CanteraODESolver.C

\*---------------------------------------------------------------------------*/

#ifndef CanteraODESolver_H
#define CanteraODESolver_H

#include "ODESystem.H"
#include "UniformField.H"
#include "typeInfo.H"
#include "autoPtr.H"
// cantera includes
#include "cantera/IdealGasMix.h"
#include "cantera/zerodim.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ODESolver Declaration
\*---------------------------------------------------------------------------*/

class CanteraODESolver
{

protected:

    // Protected data

        //- Reference to ODESystem
        const ODESystem& odes_;

        //- Maximum size of the ODESystem
        const label maxN_;

        //- Size of the ODESystem (adjustable)
        mutable label n_;

        //- Absolute convergence tolerance per step
        scalarField absTol_;

        //- Relative convergence tolerance per step
        scalarField relTol_;

        //- The maximum number of sub-steps allowed for the integration step
        label maxSteps_;

        //- The gas object
        Cantera::IdealGasMix gas_;

        //- The reactor
        Cantera::IdealGasConstPressureReactor reac_;

        //- The reactor network
        Cantera::ReactorNet net_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        CanteraODESolver(const CanteraODESolver&);

        //- Disallow default bitwise assignment
        void operator=(const CanteraODESolver&);


public:

    friend class ODESystem;

    //- Runtime type information
    TypeName("CanteraODESolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            CanteraODESolver,
            dictionary,
            (const ODESystem& ode, const dictionary& dict),
            (ode, dict)
        );


    // Constructors

        //- Construct for given ODESystem
        CanteraODESolver(const ODESystem& ode, const dictionary& dict);

        //- Construct for given ODESystem specifying tolerances
        CanteraODESolver
        (
            const ODESystem& ode,
            const scalarField& absTol,
            const scalarField& relTol
        );


    // Selectors

        //- Select null constructed
        static autoPtr<CanteraODESolver> New
        (
            const ODESystem& ode,
            const dictionary& dict
        );


    //- Destructor
    virtual ~CanteraODESolver()
    {}


    // Member Functions

        //- Return the number of equations to solve
        inline label nEqns() const;

        //- Return access to the absolute tolerance field
        inline scalarField& absTol();

        //- Return access to the relative tolerance field
        inline scalarField& relTol();

        //- Resize the ODE solver
        virtual bool resize() = 0;

        template<class Type>
        static inline void resizeField(UList<Type>& f, const label n);

        template<class Type>
        inline void resizeField(UList<Type>& f) const;

        //- Solve num problems up to deltaT
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CanteraODESolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
