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

#include "pasrIPM.H"
// include MPI header for PControl
#include "mpi.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::pasrIPM<ReactionThermo>::pasrIPM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    PaSR<ReactionThermo>(modelType, thermo, turb, combustionProperties)
{
    Info << "Using IPM-profiled version of PaSR." << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::pasrIPM<ReactionThermo>::~pasrIPM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::pasrIPM<ReactionThermo>::correct()
{
    MPI_Pcontrol(1, "PaSR_correction");
    PaSR<ReactionThermo>::correct();
    MPI_Pcontrol(-1, "PaSR_correction");
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::pasrIPM<ReactionThermo>::R(volScalarField& Y) const
{
    MPI_Pcontrol(1, "PaSR_fuel_consumption");
    const Foam::tmp<Foam::fvScalarMatrix>& R = PaSR<ReactionThermo>::R(Y);
    MPI_Pcontrol(-1, "PaSR_fuel_consumption");
    return R;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::pasrIPM<ReactionThermo>::Qdot() const
{
    MPI_Pcontrol(1, "PaSR_heat_release_rate");
    const Foam::tmp<Foam::volScalarField>& tQdot = PaSR<ReactionThermo>::Qdot();
    MPI_Pcontrol(-1, "PaSR_heat_release_rate");
    return tQdot;
}

// ************************************************************************* //
