import cantera as ct

gas = ct.Solution('gri30.cti')

for T in [1000, 1500]:
    for p in [1, 10, 25]:
        for phi in [0.5, 1, 2]:
            gas.TPX = T, p * ct.one_atm, 'N2:1'
            gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')
            reac = ct.IdealGasConstPressureReactor(gas)
            net = ct.ReactorNet([reac])
            while reac.T < T + 500:
                net.step()
            t_ign = net.time * 1.25
            print('T={}, P={:.0f}, phi={:.1f}, time={}'.format(
                T, p*ct.one_atm, phi, t_ign))
            # ensure ignited
            net.advance(t_ign)
            assert reac.T > 2000
