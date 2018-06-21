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

#include "StandardChemistryModelIPM.H"
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// ipm_control forward declare
// there's no way we can get to any of these function w/o IPM being initialized
extern "C" {
    int ipm_control(const int ctl, char *cmd, void *data);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::StandardChemistryModelIPM
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::
~StandardChemistryModelIPM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::tc() const
{
    ipm_control(1, const_cast<char*>("Chemistry_time_scale"), 0);
    tmp<volScalarField> ttc = \
        StandardChemistryModel<ReactionThermo, ThermoType>::tc();
    ipm_control(-1, const_cast<char*>("Chemistry_time_scale"), 0);
    return ttc;
}


template<class ReactionThermo, class ThermoType>
void Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::calculate()
{
    ipm_control(1, const_cast<char*>("Chemistry_reaction_rates"), 0);
    StandardChemistryModel<ReactionThermo, ThermoType>::calculate();
    ipm_control(-1, const_cast<char*>("Chemistry_reaction_rates"), 0);
}



template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    ipm_control(1, const_cast<char*>("ODE_solve"), 0);
    Foam::scalar deltaTMin = \
        Foam::StandardChemistryModel<ReactionThermo, ThermoType>::solve(deltaT);
    ipm_control(-1, const_cast<char*>("ODE_solve"), 0);
    return deltaTMin;
}


// ************************************************************************* //
