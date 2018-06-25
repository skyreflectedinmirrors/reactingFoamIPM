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
// include MPI header for PControl
#include "mpi.h"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::StandardChemistryModelIPM
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo)
{
    InfoInFunction << "Using IPM-profiled version of Chemistry." << nl;
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
    MPI_Pcontrol(1, "Chemistry_time_scale");
    tmp<volScalarField> ttc = \
        StandardChemistryModel<ReactionThermo, ThermoType>::tc();
    MPI_Pcontrol(-1, "Chemistry_time_scale");
    return ttc;
}


template<class ReactionThermo, class ThermoType>
void Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::calculate()
{
    MPI_Pcontrol(1, "Chemistry_reaction_rates");
    StandardChemistryModel<ReactionThermo, ThermoType>::calculate();
    MPI_Pcontrol(-1, "Chemistry_reaction_rates");
}



template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::StandardChemistryModelIPM<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    MPI_Pcontrol(1, "ODE_solve");
    Foam::scalar deltaTMin = \
        Foam::StandardChemistryModel<ReactionThermo, ThermoType>::solve(deltaT);
    MPI_Pcontrol(-1, "ODE_solve");
    return deltaTMin;
}


// ************************************************************************* //
