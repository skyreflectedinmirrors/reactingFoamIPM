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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "StandardChemistryModelIPM.H"
#include "BatchedChemistryModelIPM.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // base types initialized previously

     // Chemistry models based on sensibleEnthalpy
    // IPM
    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constHThermoPhysics
    );

    // Chemistry models based on sensibleInternalEnergy
    // IPM
    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        psiReactionThermo,
        constEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModelIPM,
        rhoReactionThermo,
        constEThermoPhysics
    );

    // Chemistry models based on sensibleEnthalpy
    // IPM

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );

    // Chemistry models based on sensibleInternalEnergy
    // IPM
    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );

    // Chemistry models based on sensibleEnthalpy
    // IPM

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constHThermoPhysics
    );

    // Chemistry models based on sensibleInternalEnergy
    // IPM
    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        psiReactionThermo,
        constEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        BatchedChemistryModelIPM,
        rhoReactionThermo,
        constEThermoPhysics
    );

}
// ************************************************************************* //
