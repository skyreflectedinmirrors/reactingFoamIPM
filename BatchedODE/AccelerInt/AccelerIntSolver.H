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
    Foam::AccelerIntSolver

Description
    Inteface to the accelerInt ODE-integration library

SourceFiles
    AccelerInt.C

\*---------------------------------------------------------------------------*/

#ifndef AccelerInt_H
#define AccelerInt_H

#include <vector>
#include <string>

#include "Dictionary.H"
#include "BatchedODESolver.H"
#include "solver_interface.hpp"
// include paths library from accelerInt for covenience
#include "path.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AccelerInt Declaration
\*---------------------------------------------------------------------------*/

class AccelerIntSolver
:
    public BatchedODESolver
{
protected:
    const word platform_;
    const word accSolver_;
    const word deviceType_;
    const label vectorSize_;
    const label blockSize_;
    const word order_;
    const word pyjac_path_;
    const word our_path_;

    std::unique_ptr<opencl_solvers::SolverOptions> _options;
    std::unique_ptr<opencl_solvers::IntegratorBase> _integrator;
    std::unique_ptr<opencl_solvers::IVP> _ivp;

    inline label checkVectorSize(const label vectorSize) const;
    inline label checkBlockSize(const label blockSize) const;
    inline opencl_solvers::IntegratorType integratorType(const word accSolver_) const;
    inline opencl_solvers::DeviceType deviceType(const word deviceType_) const;


public:

    //- Runtime type information
    TypeName("AccelerInt");


    // Constructors

        //- Construct from ODESystem
        AccelerIntSolver(const BatchedODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~AccelerIntSolver()
    {}


    // Member Functions

        // Override integrate from BatchedODESolver
        //- Solve num problems up to deltaT
        virtual void integrate
        (
            label num,
            const scalarField& deltaT,
            scalarField& phi,
            const scalarField& p,
            scalarField& dtChem
        ) const;

        //- Solve num problems up to deltaT
        virtual void integrate
        (
            label num,
            const UniformField<scalar>& deltaT,
            scalarField& phi,
            const scalarField& p,
            scalarField& dtChem
        ) const;

        //- Resize the ODE solver
        virtual bool resize();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "AccelerIntSolverI.H"

#endif

// ************************************************************************* //
