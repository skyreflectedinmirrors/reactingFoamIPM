/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "edcIPM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Simply call base EDC constructor

template<class ReactionThermo>
Foam::combustionModels::edcIPM<ReactionThermo>::edcIPM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    EDC<ReactionThermo>(modelType, thermo, turb, combustionProperties)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::edcIPM<ReactionThermo>::~edcIPM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::edcIPM<ReactionThermo>::correct()
{
    MPI_PControl(1, "EDC_correction");
    EDC<ReactionThermo>::correct();
    MPI_PControl(-1, "EDC_correction");
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::edcIPM<ReactionThermo>::R(volScalarField& Y) const
{
    MPI_PControl(1, "EDC_fuel_consumption");
    Foam::tmp<Foam::fvScalarMatrix>& R = EDC<ReactionThermo>::R(Y);
    MPI_PControl(-1, "EDC_fuel_consumption");
    return R;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::edcIPM<ReactionThermo>::Qdot() const
{
    MPI_PControl(1, "EDC_heat_release_rate");
    Foam::tmp<Foam::volScalarField>& tQdot = EDC<ReactionThermo>::Qdot();
    MPI_PControl(-1, "EDC_heat_release_rate");
}

// ************************************************************************* //
