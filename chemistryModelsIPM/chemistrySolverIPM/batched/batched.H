/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::batched

Description
    An batch-evaluated ODE solver for chemistry

SourceFiles
    batched.C

\*---------------------------------------------------------------------------*/

#ifndef batched_H
#define batched_H

#include "BatchedODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class batched Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class batched
:
    public ChemistryModel
{
    // Private data

        dictionary coeffsDict_;

        mutable autoPtr<BatchedODESolver> odeSolver_;



public:

    //- Runtime type information
    TypeName("batchedODE");


    // Constructors

        //- Construct from thermo
        batched(typename ChemistryModel::reactionThermo& thermo);


    //- Destructor
    virtual ~batched();


    // Member Functions

        //- Interface to the Batched ODE integrator
        virtual void integrate
        (
            label num,
            const scalarField& deltaT,
            scalarField& phi,
            const scalarField& p,
            scalarField& dtChem
        ) const;

        //- Interface to the Batched ODE integrator
        virtual void integrate
        (
            label num,
            const UniformField<scalar>& deltaT,
            scalarField& phi,
            const scalarField& p,
            scalarField& dtChem
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "batched.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
