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

Class
    Foam::combustionModels::EDC

Description
    Eddy Dissipation Concept (EDC) turbulent combustion model.
    This functions as a thin-wrapper around the base OpenFOAM EDC model, and
    only adds IPM profiling sections to the EDC model.

    This model considers that the reaction occurs in the regions of the flow
    where the dissipation of turbulence kinetic energy takes place (fine
    structures). The mass fraction of the fine structures and the mean residence
    time are provided by an energy cascade model.

    There are many versions and developments of the EDC model, 4 of which are
    currently supported in this implementation: v1981, v1996, v2005 and
    v2016.  The model variant is selected using the optional \c version entry in
    the \c EDCCoeffs dictionary, \eg

    \verbatim
        EDCCoeffs
        {
            version v2016;
        }
    \endverbatim

    The default version is \c v2005 if the \c version entry is not specified.

    Model versions and references:
    \verbatim
        Version v2005:

            Cgamma = 2.1377
            Ctau = 0.4083
            kappa = gammaL^exp1 / (1 - gammaL^exp2),

            where exp1 = 2, and exp2 = 2.

            Magnussen, B. F. (2005, June).
            The Eddy Dissipation Concept -
            A Bridge Between Science and Technology.
            In ECCOMAS thematic conference on computational combustion
            (pp. 21-24).

        Version v1981:

            Changes coefficients exp1 = 3 and exp2 = 3

            Magnussen, B. (1981, January).
            On the structure of turbulence and a generalized
            eddy dissipation concept for chemical reaction in turbulent flow.
            In 19th Aerospace Sciences Meeting (p. 42).

        Version v1996:

            Changes coefficients exp1 = 2 and exp2 = 3

            Gran, I. R., & Magnussen, B. F. (1996).
            A numerical study of a bluff-body stabilized diffusion flame.
            Part 2. Influence of combustion modeling and finite-rate chemistry.
            Combustion Science and Technology, 119(1-6), 191-217.

        Version v2016:

            Use local constants computed from the turbulent Da and Re numbers.

            Parente, A., Malik, M. R., Contino, F., Cuoci, A., & Dally, B. B.
            (2016).
            Extension of the Eddy Dissipation Concept for
            turbulence/chemistry interactions to MILD combustion.
            Fuel, 163, 98-111.
    \endverbatim

SourceFiles
    edcIPM.C

\*---------------------------------------------------------------------------*/

#ifndef PASR_IPM_H
#define PASR_IPM_H

#include "PaSR.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class EDC IPM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class pasrIPM
:
    public PaSR<ReactionThermo>
{
    // Private Member Functions

        //- Disallow copy construct
        pasrIPM(const pasrIPM&);

        //- Disallow default bitwise assignment
        void operator=(const pasrIPM&);


public:

    //- Runtime type information
    TypeName("pasrIPM");


    // Constructors

        //- Construct from components
        pasrIPM
        (
            const word& modelType,
            ReactionThermo& type,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~pasrIPM();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pasrIPM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
