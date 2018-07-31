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
// include MPI header for PControl
#include "mpi.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
{
    Info << "Using IPM-profiled version of EDC." << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::edcIPM<ReactionThermo>::~edcIPM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::edcIPM<ReactionThermo>::correct()
{
    MPI_Pcontrol(1, "EDC_correction");
    if (this->active())
    {
        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const volScalarField& epsilon = tepsilon();

        tmp<volScalarField> tmu(this->turbulence().mu());
        const volScalarField& mu = tmu();

        tmp<volScalarField> tk(this->turbulence().k());
        const volScalarField& k = tk();

        tmp<volScalarField> trho(this->rho());
        const volScalarField& rho = trho();

        scalarField tauStar(epsilon.size(), 0);

        if (version_ == EDCversions::v2016)
        {
            tmp<volScalarField> ttc(this->chemistryPtr_->tc());
            const volScalarField& tc = ttc();

            forAll(tauStar, i)
            {
                const scalar nu = mu[i]/(rho[i] + small);

                const scalar Da =
                    max(min(sqrt(nu/(epsilon[i] + small))/tc[i], 10), 1e-10);

                const scalar ReT = sqr(k[i])/(nu*epsilon[i] + small);
                const scalar CtauI = min(C1_/(Da*sqrt(ReT + 1)), 2.1377);

                const scalar CgammaI =
                    max(min(C2_*sqrt(Da*(ReT + 1)), 5), 0.4082);

                const scalar gammaL =
                    CgammaI*pow025(nu*epsilon[i]/(sqr(k[i]) + small));

                tauStar[i] = CtauI*sqrt(nu/(epsilon[i] + small));

                if (gammaL >= 1)
                {
                    kappa_[i] = 1;
                }
                else
                {
                    kappa_[i] =
                        max
                        (
                            min
                            (
                                pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                                1
                            ),
                            0
                        );
                }
            }
        }
        else
        {
            forAll(tauStar, i)
            {
                const scalar nu = mu[i]/(rho[i] + small);
                const scalar gammaL =
                    Cgamma_*pow025(nu*epsilon[i]/(sqr(k[i]) + small));

                tauStar[i] = Ctau_*sqrt(nu/(epsilon[i] + small));
                if (gammaL >= 1)
                {
                    kappa_[i] = 1;
                }
                else
                {
                    kappa_[i] =
                        max
                        (
                            min
                            (
                                pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                                1
                            ),
                            0
                        );
                }
            }
        }

        //this->chemistryPtr_->solve(tauStar);
    }
    MPI_Pcontrol(-1, "EDC_correction");
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::edcIPM<ReactionThermo>::R(volScalarField& Y) const
{
    MPI_Pcontrol(1, "EDC_fuel_consumption");
    const Foam::tmp<Foam::fvScalarMatrix>& R = EDC<ReactionThermo>::R(Y);
    MPI_Pcontrol(-1, "EDC_fuel_consumption");
    return R;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::edcIPM<ReactionThermo>::Qdot() const
{
    MPI_Pcontrol(1, "EDC_heat_release_rate");
    const Foam::tmp<Foam::volScalarField>& tQdot = EDC<ReactionThermo>::Qdot();
    MPI_Pcontrol(-1, "EDC_heat_release_rate");
    return tQdot;
}

// ************************************************************************* //
