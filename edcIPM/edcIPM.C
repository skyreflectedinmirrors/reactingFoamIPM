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

// ipm_control forward declare
// there's no way we can get to any of these function w/o IPM being initialized
extern "C" {
    int ipm_control(const int ctl, char *cmd, void *data);
}

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
    ipm_control(1, const_cast<char*>("EDC_correction"), 0);
    EDC<ReactionThermo>::correct();
    ipm_control(-1, const_cast<char*>("EDC_correction"), 0);
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::edcIPM<ReactionThermo>::R(volScalarField& Y) const
{
    ipm_control(1, const_cast<char*>("EDC_fuel_consumption"), 0);
    const Foam::tmp<Foam::fvScalarMatrix>& R = EDC<ReactionThermo>::R(Y);
    ipm_control(-1, const_cast<char*>("EDC_fuel_consumption"), 0);
    return R;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::edcIPM<ReactionThermo>::Qdot() const
{
    ipm_control(1, const_cast<char*>("EDC_heat_release_rate"), 0);
    const Foam::tmp<Foam::volScalarField>& tQdot = EDC<ReactionThermo>::Qdot();
    ipm_control(-1, const_cast<char*>("EDC_heat_release_rate"), 0);
    return tQdot;
}

// ************************************************************************* //
