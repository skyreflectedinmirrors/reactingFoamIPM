/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#ifndef makeChemistrySolverTypesIPM_H
#define makeChemistrySolverTypesIPM_H

#include "chemistrySolver.H"

#include "StandardChemistryModelIPM.H"
#include "BatchedChemistryModelIPM.H"

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"
#include "batched.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistrySolverTypeIPM(SS, Comp, Thermo)                           \
                                                                               \
    typedef SS<StandardChemistryModelIPM<Comp, Thermo>> SS##Comp##Thermo;      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(StandardChemistryModelIPM<Comp, Thermo>::typeName_()) + \
            "<" + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>"    \
            ).c_str(),                                                         \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel##Comp\
##Table_;                                                                      \

#define makeBatchedChemistrySolverTypeIPM(SS, Comp, Thermo)                    \
                                                                               \
    typedef SS<BatchedChemistryModelIPM<Comp, Thermo>> Batched##SS##Comp       \
##Thermo;                                                                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Batched##SS##Comp##Thermo,                                             \
        (#SS"<" + word(BatchedChemistryModelIPM<Comp, Thermo>::typeName_()) +  \
            "<" + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>"    \
            ).c_str(),                                                         \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<Batched##SS##Comp##Thermo>             \
        add##Batched##SS##Comp##Thermo##thermo##ConstructorTo##                \
        BasicChemistryModel##Comp##Table_;


#define makeChemistrySolverTypesIPM(Comp, Thermo)                              \
                                                                               \
    makeChemistrySolverTypeIPM                                                 \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistrySolverTypeIPM                                                 \
    (                                                                          \
        EulerImplicit,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistrySolverTypeIPM                                                 \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeBatchedChemistrySolverTypeIPM                                          \
    (                                                                          \
        batched,                                                               \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
