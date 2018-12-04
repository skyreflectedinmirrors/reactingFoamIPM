/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "CanteraSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CanteraSolver, 0);
    addToRunTimeSelectionTable(CanteraODESolver, CanteraSolver, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CanteraSolver::CanteraSolver(const ODESystem& ode, const dictionary& dict)
:
    CanteraODESolver(ode, dict),
    gas_(string(dict.lookup("ct_mechanism")))
{
    // add the gas to the reactor
    reac_.insert(gas_);

    // and initialize the reactor network
    net_.setTolerances(relTol_[0], absTol_[0]);
    net_.setMaxErrTestFails(5);
    net_.addReactor(reac_);
    net_.integrator().setMaxSteps(maxSteps_);
    net_.reinitialize();
}

//- Solve num problems up to deltaT
void Foam::CanteraSolver::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
)
{
    // reset state
    this->gas_.setConcentrations(&c[0]);
    this->gas_.setState_TP(T, p);
    // reset reactor & net
    this->reac_.syncState();
    this->net_.setInitialTime(0);
    this->net_.setMaxTimeStep(deltaT);
    this->net_.reinitialize();
    // and solver
    this->net_.advance(this->net_.time() + deltaT);
    // and copy back
    T = this->reac_.temperature();
    p = this->reac_.pressure();
    this->reac_.contents().getConcentrations(&c[0]);
    // advance a single extra step to get an estimation of the current-stepsize
    subDeltaT = this->net_.step();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::CanteraSolver::resize()
{
    return false;
}


// ************************************************************************* //
