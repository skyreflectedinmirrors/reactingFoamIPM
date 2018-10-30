/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline Foam::label Foam::AccelerIntSolver::checkVectorSize(const label vectorSize) const
{
    if (!(vectorSize >= 0))
    {
        FatalErrorInFunction
            << "Vector size " << vectorSize
            << " less than zero." << abort(FatalError);
    }
    if (vectorSize > 0 && !((((vectorSize & -vectorSize) == vectorSize) ||
                              vectorSize == 3) && vectorSize <= 16))
    {
        FatalErrorInFunction
            << "Incorrect vector size " << vectorSize
            << ". Possible values are [1, 2, 3, 4, 8, 16]." << abort(FatalError);
    }
    return vectorSize;
}

inline Foam::label Foam::AccelerIntSolver::checkBlockSize(const label blockSize) const
{
    if (!(blockSize >= 0))
    {
        FatalErrorInFunction
            << "Block size " << blockSize
            << " less than zero." << abort(FatalError);
    }
    if (blockSize > 0 && !(((blockSize & -blockSize) == blockSize) &&
                            blockSize <= 16))
    {
        FatalErrorInFunction
            << "Incorrect block size " << blockSize
            << ". Possible values are powers of two." << abort(FatalError);
    }
    return blockSize;
}

inline opencl_solvers::IntegratorType Foam::AccelerIntSolver::integratorType
(
    const Foam::word accSolver_
) const
{
    if (accSolver_.match("RKF45"))
        return opencl_solvers::IntegratorType::RKF45;
    else if (accSolver_.match("ROS3"))
        return opencl_solvers::IntegratorType::ROS3;
    else if (accSolver_.match("ROS4"))
        return opencl_solvers::IntegratorType::ROS4;
    else if (accSolver_.match("RODAS3"))
        return opencl_solvers::IntegratorType::RODAS3;
    else if (accSolver_.match("RODAS4"))
        return opencl_solvers::IntegratorType::RODAS4;
    FatalErrorInFunction
        << "Integrator type " << accSolver_ << " not recognized. "
        << "Possible devices are (RKF45, ROS3, ROS4, RODAS3, RODAS4)."
        << abort(FatalError);
    return opencl_solvers::IntegratorType::RKF45;
}

inline opencl_solvers::DeviceType Foam::AccelerIntSolver::deviceType
(
    const Foam::word deviceType_
) const
{
    if (deviceType_.match("CPU"))
        return opencl_solvers::DeviceType::CPU;
    else if (deviceType_.match("GPU"))
        return opencl_solvers::DeviceType::GPU;
    else if (deviceType_.match("ACCELERATOR"))
        return opencl_solvers::DeviceType::ACCELERATOR;
    else if (deviceType_.match("DEFAULT"))
        return opencl_solvers::DeviceType::DEFAULT;
    FatalErrorInFunction
        << "DeviceType " << deviceType_ << " not recognized. "
        << "Possible devices are (CPU, GPU, ACCELERATOR, DEFAULT)."
        << abort(FatalError);
    return opencl_solvers::DeviceType::DEFAULT;
}


// ************************************************************************* //
