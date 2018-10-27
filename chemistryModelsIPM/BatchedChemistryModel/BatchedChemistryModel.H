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
    Foam::BatchedChemistryModel

Description
    Extends base chemistry model by adding an interface to an external batched-ODE
    solver.

SourceFiles
    BatchedChemistryModelI.H
    BatchedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef BatchedChemistryModel_H
#define BatchedChemistryModel_H

#include "BasicChemistryModel.H"
#include "Reaction.H"
#include "BatchedODESystem.H"
#include "volFields.H"
#include "UniformField.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class StandardChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class BatchedChemistryModel
:
    public BasicChemistryModel<ReactionThermo>,
    public BatchedODESystem
{
    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Disallow copy constructor
        BatchedChemistryModel(const BatchedChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const BatchedChemistryModel&);


protected:

    typedef ThermoType thermoType;


    // Protected data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction<ThermoType>>& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        //- List of reaction rate per specie [kg/m3/s]
        PtrList<volScalarField::Internal> RR_;

    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();

        // time setters / inits
        void initTime
        (
            const label size,
            const scalarField& times,
            scalarField& out_times
        );

        void initTime
        (
            const label size,
            const UniformField<scalar>& times,
            UniformField<scalar>& out_times
        );

        void setTime
        (
            const scalarField& dt,
            scalarField& dt_out,
            label celli,
            label count
        );

        void setTime
        (
            const scalarField& dt,
            UniformField<scalar>& dt_out,
            label celli,
            label count
        );


public:

    //- Runtime type information
    TypeName("batched");


    // Constructors

        //- Construct from thermo
        BatchedChemistryModel(ReactionThermo& thermo);


    //- Destructor
    virtual ~BatchedChemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<ThermoType>>& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        virtual inline label nSpecie() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar Treact() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar& Treact();

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Interface to the Batched ODE integrator
            virtual void integrate
            (
                label num,
                const scalarField& deltaT,
                scalarField& phi,
                const scalarField& p
            ) = 0;

            //- Interface to the Batched ODE integrator
            virtual void integrate
            (
                label num,
                const UniformField<scalar>& deltaT,
                scalarField& phi,
                const scalarField& p
            ) = 0;

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s3]
            virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            inline virtual label nEqns() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BatchedChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BatchedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
