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
    // get rho to get size
    tmp<volScalarField> trho(this->thermo().rho());
    const label size = trho().size();
    deltaTChem = scalarField(size);
    phi0 = scalarField(nEqns() * size);
    phi = scalarField(nEqns() * size);
    integrationMask = labelField(size);
    dt = scalarField(size);
    _p = scalarField(size);

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

#include "thermodynamicConstants.H"

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

    label count = 0;

    #define TIndex(count) (nEqns() * (count))
    #define VIndex(count) (nEqns() * (count) + 1)
    #define concIndex(count, i) (nEqns() * (count) + (i) + 2)

    scalarField nNs0(rho.size());

    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho[celli];
            const scalar Vi = this->mesh_.V()[celli];
            // store temperature
            phi[TIndex(count)] = Ti;
            // store 'volume'
            phi[VIndex(count)] = Vi;
            phi0[VIndex(count)] = Vi;
            // store pressure
            _p[count] = p[celli];

            // determine moles of last specie as:
            //      [Cns] = PV / RT - sum([Ci], i=1...Ns-1)
            nNs0[count] = _p[count] * Vi / (Foam::constant::thermodynamic::RR * Ti);

            // convert mass fractions to moles
            for (label i=0; i<nSpecie_ - 1; i++)
            {
                scalar ni = Vi*rhoi*Y_[i][celli]/specieThermo_[i].W();
                phi[concIndex(count, i)] = ni;
                phi0[concIndex(count, i)] = ni;
                nNs0[count] -= ni;
            }

            // Initialise time progress
            dt[count] = deltaT[celli];

            // and store mask
            integrationMask[celli] = count++;
        }
        else
        {
            integrationMask[celli] = -1;
        }
    }

    this->integrate(count, dt, phi, _p, deltaTChem);


    forAll(rho, celli)
    {
        if (integrationMask[celli] >= 0)
        {
            const label mask = integrationMask[celli];
            const scalar dtinv = 1.0 / deltaT[mask];
            const scalar dVinv = 1.0 / phi[VIndex(mask)];
            // determine moles of last specie as:
            //      [Cns] = PV / RT - sum([Ci], i=1...Ns-1)
            scalar nNs = _p[mask] * phi[VIndex(mask)] / (
                Foam::constant::thermodynamic::RR * phi[TIndex(mask)]);
            for (label i=0; i<nSpecie_ - 1; i++)
            {
                scalar ni = max(phi[concIndex(mask, i)], 0.0);
                // note: both the integrated (mean) moles and fine-scale (original)
                // moles are divided by the new (mean) volume, as determined by
                // accelerint -- rather than each by it's own volume
                //      see paper for details
                RR_[i][celli] = specieThermo_[i].W() * dtinv * dVinv * (
                    ni - phi0[concIndex(mask, i)]);
                nNs -= ni;
            }
            // and set last species reaction rate
            RR_[nSpecie_ - 1][celli] = specieThermo_[nSpecie_ - 1].W() * dtinv *
                dVinv * (nNs - nNs0[mask]);
            // and read back deltaTChem
            this->deltaTChem_[celli] = deltaTChem[mask];
            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);
            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

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
    #undef concIndex0
    #undef TIndex
    #undef VIndex
    #undef VIndex0

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
