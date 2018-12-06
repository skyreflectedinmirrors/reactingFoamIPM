# constant pressure ignition problem for Cantera
from argparse import ArgumentParser
import os
import cantera as ct
import numpy as np
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

    # and the comparison state
    comp = np.array([reac.thermo.T] + list(reac.thermo.concentrations))
    return np.array(phi), comp


def interp(time, ct_t, ct_temp):
    # find closest time
    index = np.where(ct_t >= time)[0][0]
    T1 = ct_temp[index - 1]
    T2 = ct_temp[index]
    t1 = ct_t[index - 1]
    t2 = ct_t[index]
    return T1 + (T2 - T1) * (time - t1) / (t2 - t1)


def err_norm(test, ct):
    # calculate err
    return np.abs(test - ct) / np.abs(1e-300 + ct)


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


def sample(times, arr, dt=5e-03, percent_diff=75):
    """ return a reasonably thinned data array for plotting """
    indicies = [0]
    last = arr[0]

    def percent(val):
        return percent_diff <= 100 * (np.abs(
                arr[i + 1] - last) / np.abs(1e-300 + last))

    for i in range(arr.size - 1):
        if times[i + 1] - times[indicies[-1]] >= dt:
            indicies.append(i)
            last = arr[i + 1]
        elif percent(arr[i + 1]):
            percent(arr[i + 1])
            indicies.append(i)
            last = arr[i + 1]

    return times[np.array(indicies)], arr[np.array(indicies)]


def err(pressure, temperature, phi, endtime, to_plot):
    gas = get_gas()
    acc = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.accelode'.format(
                     int(pressure), int(temperature), float(phi))),
                     delimiter='\t', skiprows=1)
    OF = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.foamode'.format(
                    int(pressure), int(temperature), float(phi))),
                    delimiter='\t', skiprows=1)
    phi_ct, comp_ct = ignition(pressure, temperature, phi, endtime)
    # get temperature comparison
    of_err = err_norm(OF, comp_ct)
    ai_err = err_norm(acc, comp_ct)
    with open('{}_{}_{}.log'.format(int(pressure / ct.one_atm), int(temperature),
                                    float(phi)), 'w') as file:
        file.write('OF error-norm:{}\n'.format(np.array_str(of_err)))
        file.write('AI error-norm:{}\n'.format(np.array_str(ai_err)))

    # load plotter's
    acc = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.accel'.format(
                     int(pressure), int(temperature), float(phi))),
                     delimiter='\t', skiprows=1).reshape((-1, gas.n_species + 2))
    OF = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.foam'.format(
                    int(pressure), int(temperature), float(phi))),
                    delimiter='\t', skiprows=1).reshape((-1, gas.n_species + 2))

    for val in to_plot:
        ind = (['T'] + gas.species_names).index(val)

        sample_args = {}
        if val == 'T':
            sample_args['percent_diff'] = 5
        elif val in ['NO', 'OH']:
            sample_args['percent_diff'] = 100

        ylabel, ylog = label(val)
        plt.plot(phi_ct[:, 0], phi_ct[:, 1 + ind],
                 label='Cantera', **wheel(0))
        plt.plot(*sample(acc[:, 0], acc[:, 1 + ind], **sample_args),
                 label=r'\texttt{accelerInt}', **wheel(1))
        plt.plot(*sample(OF[:, 0], OF[:, 1 + ind], **sample_args),
                 label=r'\texttt{OpenFOAM}', **wheel(2))
        plt.xlabel('Time (s)')
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
