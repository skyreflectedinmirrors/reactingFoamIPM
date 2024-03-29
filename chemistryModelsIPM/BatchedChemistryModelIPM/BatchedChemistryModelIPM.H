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

Class
    Foam::BatchedChemistryModelIPM

Description
    Extends the base BatchedChemistryModel to include IPM-based profiling

SourceFiles
    BatchedChemistryModelIPMI.H
    BatchedChemistryModelIPM.C

\*---------------------------------------------------------------------------*/

#ifndef BatchedChemistryModelIPM_H
#define BatchedChemistryModelIPM_H

#include "BatchedChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class BatchedChemistryModelIPM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class BatchedChemistryModelIPM
:
    public BatchedChemistryModel<ReactionThermo, ThermoType>
{
    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  [not instrumented] - we can hook onto the outer calls
        //template<class DeltaTType>
        //scalar solve(const DeltaTType& deltaT);

        //- Disallow copy constructor
        BatchedChemistryModelIPM(const BatchedChemistryModelIPM&);

        //- Disallow default bitwise assignment
        void operator=(const BatchedChemistryModelIPM&);


protected:

    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        //  [not instrumented] (inline function)
        // inline PtrList<volScalarField::Internal>& RR();


public:

    //- Runtime type information
    TypeName("batchedIPM");


    // Constructors

        //- Construct from thermo
        BatchedChemistryModelIPM(ReactionThermo& thermo);


    //- Destructor
    virtual ~BatchedChemistryModelIPM();


    // Member Functions

        //- The reactions
        //  [not instrumented] (inline function)
        // inline const PtrList<Reaction<ThermoType>>& reactions() const;

        //- Thermodynamic data of the species
        //  [not instrumented] (inline function)
        // inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        //  [not instrumented] (inline function)
        // virtual inline label nSpecie() const;

        //- The number of reactions
        //  [not instrumented] (inline function)
        // virtual inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        //  [not instrumented] (inline function)
        // inline scalar Treact() const;

        //- Temperature below which the reaction rates are assumed 0
        //  [not instrumented] (inline function)
        // inline scalar& Treact();

        //- dc/dt = omega, rate of change in concentration, for each species
        //  [not instrumented] (overhead)
        // virtual void omega
        // (
        //     const scalarField& c,
        //     const scalar T,
        //     const scalar p,
        //     scalarField& dcdt
        // ) const;


        //- Return the reaction rate for iReaction and the reference
        //  species and charateristic times
        //  [not instrumented] (overhead)
        // virtual scalar omegaI
        // (
        //     label iReaction,
        //     const scalarField& c,
        //     const scalar T,
        //     const scalar p,
        //     scalar& pf,
        //     scalar& cf,
        //     label& lRef,
        //     scalar& pr,
        //     scalar& cr,
        //     label& rRef
        // ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            //  [not instrumented] (inline function)
            // inline const volScalarField::Internal& RR
            // (
            //     const label i
            // ) const;

            //- Return non const access to chemical source terms [kg/m3/s]
            //  [not instrumented] (inline function)
            // virtual volScalarField::Internal& RR
            // (
            //     const label i
            // );

            //- Return reaction rate of the speciei in reactionI
            //  [not instrumented] (unused in normal execution)
            //virtual tmp<volScalarField::Internal> calculateRR
            //(
            //    const label reactionI,
            //    const label speciei
            //) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s3]
            //  [not instrumented] (covered by EDC Qdot)
            //virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            //  [not instrumented]
            //inline virtual label nEqns() const;

            //- [not instrumented] (too much overhead)
            // virtual void derivatives
            // (
            //     const scalar t,
            //     const scalarField& c,
            //     scalarField& dcdt
            // ) const;

            //- [not instrumented] (too much overhead)
            // virtual void jacobian
            // (
            //     const scalar t,
            //     const scalarField& c,
            //     scalarField& dcdt,
            //     scalarSquareMatrix& J
            // ) const;

            //- [not instrumented] (too much overhead)
            // virtual void solve
            // (
            //     scalarField &c,
            //     scalar& T,
            //     scalar& p,
            //     scalar& deltaT,
            //     scalar& subDeltaT
            // ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BatchedChemistryModelIPM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
