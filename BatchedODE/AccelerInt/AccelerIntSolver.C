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

\*---------------------------------------------------------------------------*/

#include "AccelerIntSolver.H"
#include "addToRunTimeSelectionTable.H"
// include MPI header
#include <mpi.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AccelerIntSolver, 0);
    addToRunTimeSelectionTable(BatchedODESolver, AccelerIntSolver, dictionary);
}


class JacobianKernel
{
public:
    size_t numSpecies();
    size_t requiredMemorySize();
    size_t resize(size_t problem_size, size_t work_size, bool do_not_compile=false);
};


#define xstringify(s) (stringify(s))
#define stringify(s) (#s)

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AccelerIntSolver::AccelerIntSolver(const BatchedODESystem& ode, const dictionary& dict)
:
    BatchedODESolver(ode, dict),
    platform_(dict.lookupOrDefault<word>("platform", "")),
    accSolver_(dict.lookupOrDefault<word>("accSolver", "ROS4")),
    deviceType_(dict.lookupOrDefault<word>("deviceType", "DEFAULT")),
    vectorSize_(checkVectorSize(dict.lookupOrDefault<label>("vectorSize", 0))),
    blockSize_(checkBlockSize(dict.lookupOrDefault<label>("blockSize", 0))),
    order_("C"),
    pyjac_path_(dict.lookupOrDefault<word>("pyjacPath", "pyjac/")),
    our_path_(xstringify(WRAPPER_PATH))
{
    if (vectorSize_ && blockSize_)
    {
        FatalErrorInFunction
            << "blockSize and vectorSize may not be specified at the same time."
            << abort(FatalError);
    }
    // setup integrator
    opencl_solvers::DeviceType device_type = deviceType(deviceType_);
    opencl_solvers::IntegratorType int_type = integratorType(accSolver_);

    Info << "Selecting AccelerInt OpenCL solver " << accSolver_ << nl;
    Info << "Running on platform: " << platform_  << " with device type:"
         << deviceType_ << nl;
    Info << "Using vectorSize: " << vectorSize_ << nl;
    Info << "Using blockSize: " << blockSize_ << nl;

    _options.reset(new opencl_solvers::SolverOptions(vectorSize_, blockSize_,
                                                     absTol_[0], relTol_[0],
                                                     false, true, order_,
                                                     platform_, device_type,
                                                     1, maxSteps_,
                                                     opencl_solvers::StepperType::ADAPTIVE,
                                                     std::numeric_limits<double>::quiet_NaN(),
                                                     true, false));
    // build paths to files
    filesystem::path _p(our_path_);
    filesystem::path _pj(pyjac_path_);
    Info << "Loading OpenCL kernels from: " << pyjac_path_ << nl;
    std::vector<std::string> files = {
        (_p/filesystem::path("jac.cl")).str(),
        (_p/filesystem::path("dydt.cl")).str(),
        (_pj/filesystem::path("jacobian.ocl")).str(),
        (_pj/filesystem::path("species_rates.ocl")).str(),
        (_pj/filesystem::path("chem_utils.ocl")).str()};

    // create pyjac kernel
    int init = 0;
    MPI_Initialized(&init);
    JacobianKernel jk;

    // avoid data-races with kernel binary
    if (init)
    {
        // get rank
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0)
        {
            jk.resize(1, 1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank != 0)
        {
            // don't recompile kernel
            jk.resize(1, 1, true);
        }
    }
    else
    {
        jk.resize(1, 1);
    }
    // pyJac reports the total number of species, in order to match we should have:
    // neq == nsp - 1 + 2 for (T, V) == nsp + 1
    if (jk.numSpecies() + 1 != static_cast<size_t>(n_))
    {
        FatalErrorInFunction
            << "Number of equations reported by pyJac (" << jk.numSpecies() + 1
            << ") does not match expected number of equations from OpenFOAM ("
            << n_ << ")" << abort(FatalError);
    }
    size_t work_size = jk.requiredMemorySize();
    if (vectorSize_ || blockSize_)
    {
        //pyJac returns the complete (vectorized) memory req's, but we need
        // to pass in the unvectorized size
        if (vectorSize_)
        {
            work_size /= vectorSize_;
        }
        else
        {
            work_size /= blockSize_;
        }
    }
    // finally, create the IVP
    std::vector<std::string> include_path = {pyjac_path_};
    _ivp.reset(new opencl_solvers::IVP(files, work_size, 0, include_path));
    // and build integrator
    _integrator = std::move(opencl_solvers::init(
        int_type, n_, 1, *_ivp.get(), *_options.get()));

}

//- Solve num problems up to deltaT
void Foam::AccelerIntSolver::integrate
(
    label num,
    const scalarField& deltaT,
    scalarField& phi,
    const scalarField& p,
    scalarField& dtChem
) const
{
    if (num > 0)
        opencl_solvers::integrate_varying(
            *_integrator.get(), num, 0,
            &deltaT[0], -1, &phi[0], &p[0],
            &dtChem[0]);
}

//- Solve num problems up to deltaT
void Foam::AccelerIntSolver::integrate
(
    label num,
    const UniformField<scalar>& deltaT,
    scalarField& phi,
    const scalarField& p,
    scalarField& dtChem
) const
{
    if (num > 0)
        opencl_solvers::integrate(
            *_integrator.get(), num, 0,
            deltaT[0], -1, &phi[0], &p[0],
            &dtChem[0]);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::AccelerIntSolver::resize()
{
    return false;
}


// ************************************************************************* //
