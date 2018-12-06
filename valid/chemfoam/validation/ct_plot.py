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


def get_gas():
    return ct.Solution(os.path.join(script_dir, os.pardir, 'chemkin', 'chem.cti'))


def ignition(pressure, temperature, phi, endtime):
    gas = get_gas()
    gas.TPY = temperature, pressure, 'N2:1'
    gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')

    reac = ct.IdealGasConstPressureReactor(gas)
    net = ct.ReactorNet([reac])
    lowtol = True
    if lowtol:
        print('WARNING: LOW TOLERANCES IN PLACE')
        net.atol = 1e-10
        net.rtol = 1e-06
    else:
        net.atol = 1e-20
        net.rtol = 1e-15
    phi = [np.array([0, temperature] + list(gas.Y[:]))]
    while net.time <= endtime:
        net.step()
        phi.append(np.array([net.time, reac.thermo.T] + list(reac.thermo.Y[:])))
    return np.array(phi)


def interp(time, ct_t, ct_temp):
    # find closest time
    index = np.where(ct_t >= time)[0][0]
    T1 = ct_temp[index - 1]
    T2 = ct_temp[index]
    t1 = ct_t[index - 1]
    t2 = ct_t[index]
    return T1 + (T2 - T1) * (time - t1) / (t2 - t1)


def err_norm(array, phi):
    # calculate err
    errs = np.zeros((array.shape[0], array.shape[1] - 1))
    for i, t in enumerate(array[:, 0]):
        phi_int = interp(t, phi[:, 0], phi[:, 1:])
        errs[i] = np.power((phi_int - array[i, 1:]) / (1e-300 + phi_int), 2.)
        print(errs[i, 1])
    out = np.zeros(errs.shape[1])
    for i in range(errs.shape[1]):
        out[i] = trapz(errs[:, i], array[:, 0])
    return out


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


def label(to_plot):
    if to_plot == 'T':
        return 'Temperature (K)', False
    else:
        return r'$Y_{{\ce{{{}}}}}$'.format(to_plot), True


def err(pressure, temperature, phi, endtime, to_plot):
    gas = get_gas()
    acc = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.accel'.format(
                     int(pressure), int(temperature), float(phi))),
                     delimiter='\t', skiprows=1).reshape((-1, 2 + gas.n_species))
    OF = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.foam'.format(
                    int(pressure), int(temperature), float(phi))),
                    delimiter='\t', skiprows=1).reshape((-1, 2 + gas.n_species))
    phi_ct = ignition(pressure, temperature, phi, endtime)
    # get temperature comparison
    of_err = err_norm(OF, phi_ct)
    ai_err = err_norm(acc, phi_ct)
    with open('{}_{}_{}.log'.format(int(pressure / ct.one_atm), int(temperature),
                                    float(phi)), 'w') as file:
        file.write('OF error-norm:{}\n'.format(np.array_str(of_err)))
        file.write('AI error-norm:{}\n'.format(np.array_str(ai_err)))

    thin = 30
    for val in to_plot:
        ind = (['T'] + gas.species_names).index(val)
        plt.plot(phi_ct[:, 0], phi_ct[:, 1 + ind], label='Cantera', **wheel(0))
        plt.plot(acc[:, 0][::thin], acc[:, 1 + ind][::thin],
                 label=r'\texttt{accelerInt}', **wheel(1))
        plt.plot(acc[:, 0][::thin], acc[:, 1 + ind][::thin],
                 label=r'\texttt{OpenFOAM}', **wheel(2))
        plt.xlabel('Time (s)')
        ylabel, ylog = label(val)
        plt.ylabel(ylabel)
        if ylog:
            plt.yscale('log')
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
        plt.savefig('{}_{}_{}_{}.pdf'.format(
            int(pressure / ct.one_atm), int(temperature), phi,
            val))
        plt.close()


if __name__ == '__main__':
    parser = ArgumentParser('ct_plot.py -- comparison constant-pressure ignition '
                            'solution of OpenFOAM and accelerInt '
                            'solvers to Cantera.')
    parser.add_argument('-pres', '--pressure',
                        type=float,
                        help='The pressure used.',
                        required=True)
    parser.add_argument('-T', '--temperature',
                        type=float,
                        help='The temperature used.',
                        required=True)
    parser.add_argument('-phi', '--phi',
                        type=float,
                        help='The equivalence ratio used.',
                        required=True)
    parser.add_argument('-e', '--endtime',
                        type=float,
                        help='The simulation endtime.',
                        required=True)
    parser.add_argument('-t', '--to_plot',
                        nargs='+',
                        default=['T', 'CH4', 'OH', 'HO2', 'NO'])
    args = parser.parse_args()
    err(args.pressure, args.temperature, args.phi, args.endtime, args.to_plot)
