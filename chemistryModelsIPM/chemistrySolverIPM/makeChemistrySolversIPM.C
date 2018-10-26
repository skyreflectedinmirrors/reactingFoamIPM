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

#include "makeChemistrySolverTypesIPM.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistrySolverTypesIPM(psiReactionThermo, constGasHThermoPhysics);
    makeChemistrySolverTypesIPM(psiReactionThermo, gasHThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesIPM
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesIPM(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistrySolverTypesIPM(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistrySolverTypesIPM(psiReactionThermo, constHThermoPhysics);

    makeChemistrySolverTypesIPM(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistrySolverTypesIPM(rhoReactionThermo, gasHThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesIPM
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypesIPM(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistrySolverTypesIPM(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistrySolverTypesIPM(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistrySolverTypesIPM(psiReactionThermo, constGasEThermoPhysics);
    makeChemistrySolverTypesIPM(psiReactionThermo, gasEThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesIPM
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesIPM(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistrySolverTypesIPM(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistrySolverTypesIPM(psiReactionThermo, constEThermoPhysics);

    makeChemistrySolverTypesIPM(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistrySolverTypesIPM(rhoReactionThermo, gasEThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesIPM
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypesIPM(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistrySolverTypesIPM(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistrySolverTypesIPM
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistrySolverTypesIPM(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
