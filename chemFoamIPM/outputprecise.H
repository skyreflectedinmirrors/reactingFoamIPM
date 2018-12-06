    runTime.write();

    Info<< "Qdot = " << Qdot
        << ", T = " << thermo.T()[0]
        << ", p = " << thermo.p()[0]
        << ", " << Y[0].name() << " = " << Y[0][0]
        << endl;

    // and we output the state
    post<< scientific << setprecision(16) << runTime.value() << token::TAB << thermo.T()[0];

    // and output all species
    forAll(Y, specieI)
    {
        post<< token::TAB;
        post<< scientific << setprecision(16) << Y[specieI][0];
    }
    post << nl;
