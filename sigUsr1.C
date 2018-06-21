/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "sigUsr1.H"
#include "error.H"
#include "jobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// dummy constructor for  call
extern "C"
{
    void ipm_sig_handler(int sig);
}

struct sigaction Foam::sigUsr1::oldAction_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigUsr1::sigHandler(int)
{
    // Reset old handling
    if (sigaction(SIGUSR1, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGUSR1 trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo_.signalEnd();

    Info << "Terminating w/ IPM" << nl;
    ipm_sig_handler(SIGTERM);

    // Throw signal (to old handler)
    raise(SIGTERM);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigUsr1::sigUsr1()
{
    oldAction_.sa_handler = nullptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigUsr1::~sigUsr1()
{
    // Reset old handling
    if (sigaction(SIGUSR1, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGUSR1 trapping"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigUsr1::set(const bool)
{
    if (oldAction_.sa_handler)
    {
        FatalErrorInFunction
            << "Cannot call sigUsr1::set() more than once"
            << abort(FatalError);
    }

    struct sigaction newAction;
    newAction.sa_handler = sigHandler;
    newAction.sa_flags = SA_NODEFER;
    sigemptyset(&newAction.sa_mask);
    if (sigaction(SIGUSR1, &newAction, &oldAction_) < 0)
    {
        FatalErrorInFunction
            << "Cannot set SIGUSR1 trapping"
            << abort(FatalError);
    }
}


// ************************************************************************* //
