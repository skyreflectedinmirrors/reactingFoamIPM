import cantera as ct

gas = ct.Solution('gri30.cti')

for T in [850, 1100, 1500]:
    for p in [1, 10, 25]:
        for phi in [0.5, 1, 1.5]:
            gas.TPX = T, p * ct.one_atm, 'N2:1'
            gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')
            reac = ct.IdealGasConstPressureReactor(gas)
            reac.volume = 1
            net = ct.ReactorNet([reac])
            net.atol = 1e-15
            net.rtol = 1e-10
            net.max_steps = int(1e9)
            while reac.T < T + 500:
                net.step()
            t_ign = net.time * 1.25
            # ensure ignited
            net.advance(t_ign)
            if reac.T <= 2000:
                # give it one more try
                t_ign = t_ign * 1.5
                net.advance(t_ign)

            if reac.T <= 2000 or t_ign > 5:
                print('Skipping: T={}, P={:.0f}, phi={:.1f}, time={}'.format(
                    T, p*ct.one_atm, phi, t_ign))
                continue
            print('T={}, P={:.0f}, phi={:.1f}, time={}'.format(
                T, p*ct.one_atm, phi, t_ign))
