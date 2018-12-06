/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Application
    chemFoam

Description
    Solver for chemistry problems, designed for use on single cell cases to
    provide comparison against other chemistry solvers, that uses a single cell
    mesh, and fields created from the initial conditions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "BasicChemistryModel.H"
#include "StandardChemistryModel.H"
#include "BatchedChemistryModel.H"
#include "reactingMixture.H"
#include "chemistrySolver.H"
#include "OFstream.H"
#include "thermoPhysicsTypes.H"
#include "basicSpecieMixture.H"
#include "cellModeller.H"
#include "thermoTypeFunctions.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    // create the case as usual in chemFoam
    argList::noParallel();

    #define CREATE_MESH createSingleCellMesh.H
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createSingleCellMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "readInitialConditions.H"
    #include "createControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // however, then we ditch the weird time-stepping scheme and simply solve
    // to the end time to compare to cantera

    scalarField c_(chemistry.nSpecie());
    scalar Ti = T0;

    const scalar rhoi = rho[0];
    scalar pi = p[0];

    for (label i=0; i<chemistry.nSpecie(); i++)
    {
        c_[i] = rhoi*Y[i][0]/W[i];
    }

    // Initialise time progress
    scalar dt = runTime.endTime().value();

    // Calculate the chemical source terms
    if (isA<StandardChemistryModel<rhoReactionThermo, gasHThermoPhysics>>(chemistry))
    {
        refCast<const StandardChemistryModel<rhoReactionThermo, gasHThermoPhysics>>
        (
            chemistry
        ).solve(c_, Ti, pi, dt, runTime.deltaT().value());
    }
    else
    {
        // copy into batched form
        scalarField phi(chemistry.nSpecie() + 1);
        phi[0] = Ti;
        scalar Vi = mesh.V()[0];
        phi[1] = Vi;
        scalarField pfield(1);
        scalarField dtChem_(1);
        dtChem_[0] = runTime.deltaT().value();

        pfield[0] = pi;
        for (label i=0; i<chemistry.nSpecie() - 1; i++)
        {
            phi[i + 2] = c_[i] * Vi;
        }
        refCast<const BatchedChemistryModel<rhoReactionThermo, gasHThermoPhysics>>
        (
            chemistry
        ).integrate(1, dt, phi, pfield, dtChem_);

        // and copy out
        Ti = phi[0];
        for (label i=0; i<chemistry.nSpecie() - 1; i++)
        {
            c_[i] = phi[i + 2] / Vi;
        }
    }


    // and write to the file
    post<< scientific << setprecision(16) << Ti;

    // and output all species
    forAll(c_, specieI)
    {
        post<< token::TAB;
        post<< scientific << setprecision(16) << c_[specieI];
    }
    post << nl;

    return 0;
}


// ************************************************************************* //
