/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "laminarIPM.H"
#include "mpi.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminarIPM<ReactionThermo>::laminarIPM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties)
{
    Info << "Using IPM-profiled version of laminar." << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminarIPM<ReactionThermo>::~laminarIPM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::laminarIPM<ReactionThermo>::correct()
{
    MPI_Pcontrol(1, "laminar_correction");
    laminar<ReactionThermo>::correct();
    MPI_Pcontrol(-1, "laminar_correction");
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::laminarIPM<ReactionThermo>::R(volScalarField& Y) const
{
    MPI_Pcontrol(1, "laminar_fuel_consumption");
    const Foam::tmp<Foam::fvScalarMatrix>& R = laminar<ReactionThermo>::R(Y);
    MPI_Pcontrol(-1, "laminar_fuel_consumption");
    return R;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminarIPM<ReactionThermo>::Qdot() const
{
    MPI_Pcontrol(1, "laminar_heat_release_rate");
    const Foam::tmp<Foam::volScalarField>& tQdot = laminar<ReactionThermo>::Qdot();
    MPI_Pcontrol(-1, "laminar_heat_release_rate");
    return tQdot;
}

// ************************************************************************* //
