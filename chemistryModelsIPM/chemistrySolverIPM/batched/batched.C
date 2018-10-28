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

\*---------------------------------------------------------------------------*/

#include "batched.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::batched<ChemistryModel>::batched(typename ChemistryModel::reactionThermo& thermo)
:
    ChemistryModel(thermo),
    coeffsDict_(this->subDict("odeCoeffs")),
    odeSolver_(BatchedODESolver::New(*this, coeffsDict_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::batched<ChemistryModel>::~batched()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::batched<ChemistryModel>::integrate
(
    label num,
    const scalarField& deltaT,
    scalarField& phi,
    const scalarField& p
) const
{
    // Reset the size of the batched system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        NotImplemented;
    }

    odeSolver_->integrate(num, deltaT, phi, p);
}

template<class ChemistryModel>
void Foam::batched<ChemistryModel>::integrate
(
    label num,
    const UniformField<scalar>& deltaT,
    scalarField& phi,
    const scalarField& p
) const
{
    // Reset the size of the batched system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        NotImplemented;
    }

    odeSolver_->integrate(num, deltaT, phi, p);
}


// ************************************************************************* //
