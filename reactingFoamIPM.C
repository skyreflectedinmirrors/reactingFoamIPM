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
#include "sigUsr1.H"
// include MPI header for PControl
#include <mpi.h>

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

    sigUsr1 _sigUsr1;
    _sigUsr1.set(true);

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

        MPI_Pcontrol(1, "solve_density");
        #include "rhoEqn.H"
        MPI_Pcontrol(-1, "solve_density");

        while (pimple.loop())
        {
            MPI_Pcontrol(1, "solve_velocity");
            #include "UEqn.H"
            MPI_Pcontrol(-1, "solve_velocity");
            MPI_Pcontrol(1, "solve_species");
            MPI_Pcontrol(1, "species_convection");
            tmp<fv::convectionScheme<scalar>> mvConvection
            (
                fv::convectionScheme<scalar>::New
                (
                    mesh,
                    fields,
                    phi,
                    mesh.divScheme("div(phi,Yi_h)")
                )
            );
            MPI_Pcontrol(-1, "species_convection");
            {
                reaction->correct();
                MPI_Pcontrol(1, "ODE_loadbalance_barrier");
                // insert barrier to avoid biased statistics due to chemistry load
                // imbalance... doesn't cause much of a slowdown
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Pcontrol(-1, "ODE_loadbalance_barrier");
                Qdot = reaction->Qdot();
                volScalarField Yt(0.0*Y[0]);

                forAll(Y, i)
                {
                    if (i != inertIndex && composition.active(i))
                    {
                        volScalarField& Yi = Y[i];

                        MPI_Pcontrol(1, "species_build_equation");
                        MPI_Pcontrol(1, "species_build_equation_ddt");
                        auto ddt = fvm::ddt(rho, Yi);
                        MPI_Pcontrol(-1, "species_build_equation_ddt");
                        MPI_Pcontrol(1, "species_build_equation_conv");
                        auto conv = mvConvection->fvmDiv(phi, Yi);
                        MPI_Pcontrol(-1, "species_build_equation_conv");
                        MPI_Pcontrol(1, "species_build_equation_laplac");
                        auto laplacian = fvm::laplacian(turbulence->muEff(), Yi);
                        MPI_Pcontrol(-1, "species_build_equation_laplac");
                        MPI_Pcontrol(1, "species_build_equation_rxn");
                        auto rxn = reaction->R(Yi);
                        MPI_Pcontrol(-1, "species_build_equation_rxn");
                        MPI_Pcontrol(1, "species_build_equation_opts");
                        auto opts = fvOptions(rho, Yi);
                        MPI_Pcontrol(-1, "species_build_equation_opts");
                        fvScalarMatrix YiEqn
                        (
                            ddt
                          + conv
                          - laplacian
                         ==
                            rxn
                          + opts
                        );
                        MPI_Pcontrol(-1, "species_build_equation");

                        MPI_Pcontrol(1, "species_solution");
                        YiEqn.relax();

                        fvOptions.constrain(YiEqn);

                        YiEqn.solve(mesh.solver("Yi"));

                        fvOptions.correct(Yi);
                        MPI_Pcontrol(-1, "species_solution");

                        Yi.max(0.0);
                        Yt += Yi;
                    }
                }

                Y[inertIndex] = scalar(1) - Yt;
                Y[inertIndex].max(0.0);
            }

            MPI_Pcontrol(-1, "solve_species");
            MPI_Pcontrol(1, "solve_energy");
            #include "EEqn.H"
            MPI_Pcontrol(-1, "solve_energy");

            // --- Pressure corrector loop
            MPI_Pcontrol(1, "solve_pressure");
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
            MPI_Pcontrol(-1, "solve_pressure");

            MPI_Pcontrol(1, "correct_turbulence");
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
            MPI_Pcontrol(-1, "correct_turbulence");
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
