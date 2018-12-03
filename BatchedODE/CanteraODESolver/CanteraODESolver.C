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

#include "CanteraODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CanteraODESolver, 0);
    defineRunTimeSelectionTable(CanteraODESolver, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CanteraODESolver::CanteraODESolver(const ODESystem& ode, const dictionary& dict)
:
    odes_(ode),
    maxN_(ode.nEqns()),
    n_(ode.nEqns()),
    absTol_(n_, dict.lookupOrDefault<scalar>("absTol", small)),
    relTol_(n_, dict.lookupOrDefault<scalar>("relTol", 1e-4)),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 10000)),
    gas_(word(dict.lookup("ct_mechanism")))
{
    // add the gas to the reactor
    reac_.insert(gas_);

    // and initialize the reactor network
    net_.setTolerances(absTol_[0], relTol_[0]);
    net_.addReactor(reac_);
    net_.integrator().setMaxSteps(maxSteps_);
    net_.reinitialize();
}


Foam::CanteraODESolver::CanteraODESolver
(
    const ODESystem& ode,
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

bool Foam::CanteraODESolver::resize()
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


//- Solve num problems up to deltaT
void Foam::CanteraODESolver::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT
)
{
    // reset state
    this->gas_.setState_TPX(T, p, &c[0]);
    // reset reactor & net
    this->reac_.syncState();
    this->net_.setInitialTime(0);
    this->net_.reinitialize();
    // and solver
    this->net_.advance(deltaT);
    // and copy back
    T = this->reac_.temperature();
    p = this->reac_.pressure();
    this->reac_.contents().getConcentrations(&c[0]);

}


// ************************************************************************* //
