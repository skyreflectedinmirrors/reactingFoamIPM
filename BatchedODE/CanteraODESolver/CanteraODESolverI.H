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


inline Foam::label Foam::CanteraODESolver::nEqns() const
{
    return n_;
}


inline Foam::scalarField& Foam::CanteraODESolver::absTol()
{
    return absTol_;
}


inline Foam::scalarField& Foam::CanteraODESolver::relTol()
{
    return relTol_;
}

template<class Type>
inline void Foam::CanteraODESolver::resizeField(UList<Type>& f, const label n)
{
    f.shallowCopy(UList<Type>(f.begin(), n));
}


template<class Type>
inline void Foam::CanteraODESolver::resizeField(UList<Type>& f) const
{
    resizeField(f, n_);
}



// ************************************************************************* //
