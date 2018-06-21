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

Application
    reactingFoamIPM

Description
    Solver for combustion with chemical reactions.
    Instrumented with the IPM MPI-profiling library

\*---------------------------------------------------------------------------*/

// ipm_control forward declare
// here we must make sure not to call control before IPM (i.e., MPI) has been
// initialized
extern "C" {
    int ipm_control(const int ctl, char *cmd, void *data);
}



#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        ipm_control(1, const_cast<char*>("solve_density"), 0);
        #include "rhoEqn.H"
        ipm_control(-1, const_cast<char*>("solve_density"), 0);

        while (pimple.loop())
        {
            ipm_control(1, const_cast<char*>("solve_velocity"), 0);
            #include "UEqn.H"
            ipm_control(-1, const_cast<char*>("solve_velocity"), 0);
            ipm_control(1, const_cast<char*>("solve_species"), 0);
            #include "YEqn.H"
            ipm_control(-1, const_cast<char*>("solve_species"), 0);
            ipm_control(1, const_cast<char*>("solve_energy"), 0);
            #include "EEqn.H"
            ipm_control(-1, const_cast<char*>("solve_energy"), 0);

            // --- Pressure corrector loop
            ipm_control(1, const_cast<char*>("solve_pressure"), 0);
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    #include "pEqn.H"
                }
            }
            ipm_control(-1, const_cast<char*>("solve_pressure"), 0);

            ipm_control(1, const_cast<char*>("correct_turbulence"), 0);
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
            ipm_control(-1, const_cast<char*>("correct_turbulence"), 0);
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
