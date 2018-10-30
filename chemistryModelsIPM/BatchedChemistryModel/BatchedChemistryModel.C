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

#include "BatchedChemistryModel.H"
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::BatchedChemistryModel
(
    ReactionThermo& thermo
)
:
    BasicChemistryModel<ReactionThermo>(thermo),
    BatchedODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),
    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_
    (
        BasicChemistryModel<ReactionThermo>::template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    RR_(nSpecie_)
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.p().mesh(),
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "BatchedChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::
~BatchedChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::tc() const
{

    NotImplemented;
    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    return ttc;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermo_[i].Hc();
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    NotImplemented
    tmp<volScalarField::Internal> tRR
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
        )
    );

    return tRR;
}


template<class ReactionThermo, class ThermoType>
void Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    NotImplemented;
}

template<class ReactionThermo, class ThermoType>
void Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::setTime
(
    const scalarField& dt,
    scalarField& dt_out,
    label celli,
    label count
)
{
    dt_out[count] = dt[celli];
}

template<class ReactionThermo, class ThermoType>
void Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::setTime
(
    const UniformField<scalar>& dt,
    UniformField<scalar>& dt_out,
    label celli,
    label count
)
{

}

template<class ReactionThermo, class ThermoType>
Foam::scalarField Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::initTime
(
    const label size,
    const scalarField& times
)
{
    return scalarField(size);
}

template<class ReactionThermo, class ThermoType>
Foam::UniformField<Foam::scalar>
Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::initTime
(
    const label size,
    const UniformField<scalar>& times
)
{
    return UniformField<scalar>(times[0]);
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c0(nSpecie_ * rho.size());
    scalarField phi((nSpecie_ + 1) * rho.size());
    DeltaTType dt = initTime(rho.size(), deltaT);
    labelField integrationMask(rho.size());
    label count = 0;

    #define concIndexACC(count, i) ((nSpecie_ + 1) * (count) + (i) + 2)
    #define concIndex(count, i) ((nSpecie_) * count + i)

    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho[celli];
            // store temperature
            phi[(nSpecie_ + 1) * count] = Ti;
            // store 'volume'
            phi[(nSpecie_ + 1) * count + 1] = 1.0;

            // convert mass fractions to concentrations
            for (label i=0; i<nSpecie_ - 1; i++)
            {
                phi[concIndexACC(count, i)] = rhoi*Y_[i][celli]/specieThermo_[i].W();
                c0[concIndex(count, i)] = phi[concIndexACC(count, i)];
            }
            c0[concIndex(count, nSpecie_ - 1)] =
                rhoi*Y_[nSpecie_ - 1][celli]/specieThermo_[nSpecie_ - 1].W();

            // Initialise time progress
            setTime(deltaT, dt, celli, count);

            // and store mask
            integrationMask[celli] = count++;
        }
        else
        {
            integrationMask[celli] = -1;
        }
    }

    this->integrate(count, dt, phi, p);


    forAll(rho, celli)
    {
        if (integrationMask[celli] >= 0)
        {
            const label mask = integrationMask[celli];
            const scalar Vinv = 1.0 / phi[(nSpecie_ + 1) * mask + 1];
            // determine concentration of last species from the others
            // C = sum(C[i], i=1...Ns) = P / RT
            //      -> C[Ns] = (P / RT) - sum(C[i], i=1...Ns-1)
            scalar cNs = p[mask] / (8314.4621 * phi[(nSpecie_ + 1) * mask]);
            for (label i=0; i<nSpecie_ - 1; i++)
            {
                scalar ci = max(Vinv * phi[concIndexACC(mask, i)], 0.0);
                RR_[i][mask] =
                    (ci - c0[concIndex(mask, i)])*specieThermo_[i].W()/deltaT[mask];
                cNs -= ci;
            }
            // and set last species reaction rate
            RR_[nSpecie_ - 1][mask] =
                (cNs - c0[concIndex(mask, nSpecie_ - 1)])*
                specieThermo_[nSpecie_ - 1].W() / deltaT[mask];
            // finally copy back the new temperature
            T[mask] = phi[(nSpecie_ + 1) * mask];
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }
    }

    #undef concIndex
    #undef concIndexACC

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::BatchedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
