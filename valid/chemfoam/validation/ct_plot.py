# constant pressure ignition problem for Cantera
from argparse import ArgumentParser
import os
import cantera as ct
import numpy as np
from scipy.integrate import trapz
import matplotlib as mpl
# setup latex
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rc('text.latex',
       preamble=r'\usepackage{amsmath}'
                r'\usepackage[version=4]{mhchem}')
mpl.rc('font', family='serif')
import matplotlib.pyplot as plt # noqa

script_dir = os.path.dirname(os.path.abspath(__file__))


def ignition(pressure, temperature, phi, endtime):
    gas = ct.Solution(os.path.join(script_dir, os.pardir, 'chemkin', 'chem.cti'))
    gas.TPY = temperature, pressure, 'N2:1'
    gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')

    reac = ct.IdealGasConstPressureReactor(gas)
    net = ct.ReactorNet([reac])
    net.atol = 1e-20
    net.rtol = 1e-15
    t = [0]
    Tv = [temperature]
    while net.time <= endtime:
        net.step()
        t.append(net.time)
        Tv.append(reac.thermo.T)
    return np.array(t), np.array(Tv)


def interp_temp(time, ct_t, ct_temp):
    # find closest time
    index = np.where(ct_t >= time)[0][0]
    T1 = ct_temp[index - 1]
    T2 = ct_temp[index]
    t1 = ct_t[index - 1]
    t2 = ct_t[index]
    return T1 + (T2 - T1) * (time - t1) / (t2 - t1)


def err_norm(array, times, temps):
    # calculate err
    errs = np.zeros(array[:, 0].size)
    for i, t in enumerate(array[:, 0]):
        T = interp_temp(t, times, temps)
        errs[i] = np.power((T - array[i, 1]) / T, 2.)
    return trapz(errs, array[:, 0])


def wheel(index):
    marker_wheel = ['.', 'o', 'v', 's']
    size_wheel = [6, 10, 14, 16]
    cmap = 'inferno'
    ncolors = 4
    color_wheel = plt.get_cmap(cmap, ncolors)

    style = {'color': color_wheel(index % ncolors)}
    if index:
        style.update({
            'linestyle': '',
            'marker': marker_wheel[index % len(marker_wheel)],
            'markersize': size_wheel[index % len(size_wheel)],
            'markerfacecolor': 'none'})

    return style


def err(pressure, temperature, phi, endtime):
    acc = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.accel'.format(
                     int(pressure), int(temperature), float(phi))),
                     delimiter='\t', skiprows=1).reshape((-1, 3))
    OF = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.foam'.format(
                    int(pressure), int(temperature), float(phi))),
                    delimiter='\t', skiprows=1).reshape((-1, 3))
    times, temperatures = ignition(pressure, temperature, phi, endtime)
    # get temperature comparison
    of_err = err_norm(OF, times, temperatures)
    ai_err = err_norm(acc, times, temperatures)
    with open('{}_{}_{}.log'.format(int(pressure / ct.one_atm), int(temperature),
                                    float(phi)), 'w') as file:
        file.write('OF error-norm:{}\n'.format(of_err))
        file.write('AI error-norm:{}\n'.format(ai_err))

    thin = 25

    plt.plot(times, temperatures, label='Cantera', **wheel(0))
    plt.plot(acc[:, 0][::thin], acc[:, 1][::thin],
             label=r'\texttt{accelerInt}', **wheel(1))
    plt.plot(acc[:, 0][::thin], acc[:, 1][::thin],
             label=r'\texttt{OpenFOAM}', **wheel(2))
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(**{'loc': 0,
                  'fontsize': 16,
                  'numpoints': 1,
                  'shadow': True,
                  'fancybox': True})
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=16)
    for item in (plt.gca().title, plt.gca().xaxis.label,
                 plt.gca().yaxis.label):
        item.set_fontsize(24)
    plt.tight_layout()
    plt.savefig('{}_{}_{}.pdf'.format(int(pressure / ct.one_atm), int(temperature),
                                      phi))


if __name__ == '__main__':
    parser = ArgumentParser('ct_plot.py -- comparison constant-pressure ignition '
                            'solution of OpenFOAM and accelerInt '
                            'solvers to Cantera.')
    parser.add_argument('-pres', '--pressure',
                        type=float,
                        help='The pressure used.')
    parser.add_argument('-T', '--temperature',
                        type=float,
                        help='The temperature used.')
    parser.add_argument('-phi', '--phi',
                        type=float,
                        help='The equivalence ratio used.')
    parser.add_argument('-e', '--endtime',
                        type=float,
                        help='The simulation endtime.')
    args = parser.parse_args()
    err(args.pressure, args.temperature, args.phi, args.endtime)
