import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.cti')

for p in [1, 10, 25]:
    for phi in [0.5, 1, 1.5]:
        temperatures = np.arange(850, 1500, 25, dtype=np.float64)
        ign_delays = np.zeros_like(temperatures)
        for i, T in enumerate(temperatures):
            gas.TPX = T, p * ct.one_atm, 'N2:1'
            gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')
            reac = ct.IdealGasConstPressureReactor(gas)
            reac.volume = 1
            net = ct.ReactorNet([reac])
            net.atol = 1e-15
            net.rtol = 1e-10
            net.max_steps = int(1e9)
            while reac.T < T + 400 and net.time < 100:
                net.step()
            if reac.T < T + 400:
                ign_delays[i] = -1
                continue
            ign_delays[i] = net.time

        masked = np.where(ign_delays > 0)
        plt.semilogy((1000 / temperatures)[masked], ign_delays[masked],
                     label='phi={}'.format(phi), linestyle='', marker='.')

    plt.legend()
    plt.title('P = {} atm'.format(p))
    plt.ylabel('Ign. delay (s)')
    plt.xlabel('1000 / T [1/K]')
    plt.show()
    plt.close()
