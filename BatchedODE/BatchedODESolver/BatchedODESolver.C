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

\*---------------------------------------------------------------------------*/

#include "BatchedODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BatchedODESolver, 0);
    defineRunTimeSelectionTable(BatchedODESolver, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BatchedODESolver::BatchedODESolver(const BatchedODESystem& ode, const dictionary& dict)
:
    odes_(ode),
    maxN_(ode.nEqns()),
    n_(ode.nEqns()),
    absTol_(n_, dict.lookupOrDefault<scalar>("absTol", small)),
    relTol_(n_, dict.lookupOrDefault<scalar>("relTol", 1e-4)),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 10000))
{}


Foam::BatchedODESolver::BatchedODESolver
(
    const BatchedODESystem& ode,
    const scalarField& absTol,
    const scalarField& relTol
)
:
    odes_(ode),
    maxN_(ode.nEqns()),
    n_(ode.nEqns()),
    absTol_(absTol),
    relTol_(relTol),
    maxSteps_(10000)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::BatchedODESolver::resize()
{
    if (odes_.nEqns() != n_)
    {
        if (odes_.nEqns() > maxN_)
        {
            FatalErrorInFunction
                << "Specified number of equations " << odes_.nEqns()
                << " greater than maximum " << maxN_
                << abort(FatalError);
        }

        n_ = odes_.nEqns();

        resizeField(absTol_);
        resizeField(relTol_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
