/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::CanteraSolver

Description
    Inteface to the accelerInt ODE-integration library

SourceFiles
    AccelerInt.C

\*---------------------------------------------------------------------------*/

#ifndef Cantera_H
#define Cantera_H

#include "CanteraODESolver.H"
// cantera includes
#include "cantera/IdealGasMix.h"
#include "cantera/zerodim.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AccelerInt Declaration
\*---------------------------------------------------------------------------*/

class CanteraSolver
:
    public CanteraODESolver
{
protected:
    //- The gas object
    Cantera::IdealGasMix gas_;

    //- The reactor
    Cantera::IdealGasConstPressureReactor reac_;

    //- The reactor network
    Cantera::ReactorNet net_;


public:

    //- Runtime type information
    TypeName("cantera");


    // Constructors

        //- Construct from ODESystem
        CanteraSolver(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~CanteraSolver()
    {}


    // Member Functions

        // Override integrate from CanteraODESolver
        //- Solve num problems up to deltaT
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT
        );

        //- Resize the ODE solver
        virtual bool resize();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
