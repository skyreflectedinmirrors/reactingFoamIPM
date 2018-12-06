import cantera as ct

gas = ct.Solution('gri30.cti')

for T in [900, 1200, 1500]:
    for p in [1, 10, 25]:
        for phi in [0.6, 1, 1.4]:
            if T == 900 and p == 1:
                continue
            gas.TPX = T, p * ct.one_atm, 'N2:1'
            gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')
            reac = ct.IdealGasConstPressureReactor(gas)
            reac.volume = 1
            net = ct.ReactorNet([reac])
            net.atol = 1e-20
            net.rtol = 1e-15
            net.max_steps = int(1e9)
            while reac.T < T + 500:
                net.step()
            t_ign = net.time * 1.25
            # ensure ignited
            net.advance(t_ign)
            if reac.T <= 2000 or net.time > 1:
                continue
            print('T={}, P={:.0f}, phi={:.1f}, time={}'.format(
                T, p*ct.one_atm, phi, t_ign))
