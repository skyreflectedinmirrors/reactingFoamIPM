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
    Foam::BatchedODESolver

Description
    Abstract base-class batch-solution of ODE systems

SourceFiles
    BatchedODESolver.C

\*---------------------------------------------------------------------------*/

#ifndef BatchedODESolver_H
#define BatchedODESolver_H

#include "BatchedODESystem.H"
#include "UniformField.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ODESolver Declaration
\*---------------------------------------------------------------------------*/

class BatchedODESolver
{

protected:

    // Protected data

        //- Reference to ODESystem
        const BatchedODESystem& odes_;

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


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        BatchedODESolver(const BatchedODESolver&);

        //- Disallow default bitwise assignment
        void operator=(const BatchedODESolver&);


public:

    friend class BatchedODESystem;

    //- Runtime type information
    TypeName("BatchedODESolver");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            BatchedODESolver,
            dictionary,
            (const BatchedODESystem& ode, const dictionary& dict),
            (ode, dict)
        );


    // Constructors

        //- Construct for given ODESystem
        BatchedODESolver(const BatchedODESystem& ode, const dictionary& dict);

        //- Construct for given ODESystem specifying tolerances
        BatchedODESolver
        (
            const BatchedODESystem& ode,
            const scalarField& absTol,
            const scalarField& relTol
        );


    // Selectors

        //- Select null constructed
        static autoPtr<BatchedODESolver> New
        (
            const BatchedODESystem& ode,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BatchedODESolver()
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
        virtual void integrate
        (
            label num,
            const scalarField& deltaT,
            scalarField& phi,
            const scalarField& p,
            scalarField& dtChem
        ) const = 0;

        //- Solve num problems up to deltaT
        virtual void integrate
        (
            label num,
            const UniformField<scalar>& deltaT,
            scalarField& phi,
            const scalarField& p,
            scalarField& dtChem
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BatchedODESolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
