# constant pressure ignition problem for Cantera
from argparse import ArgumentParser
import os
import cantera as ct
import numpy as np
import matplotlib as mpl
import sys
# setup latex
mpl.use('Agg')
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rc('text.latex',
       preamble=r'\usepackage{amsmath}'
                r'\usepackage[version=4]{mhchem}')
mpl.rc('font', family='serif')
import matplotlib.pyplot as plt  # noqa

script_dir = os.path.dirname(os.path.abspath(__file__))


def get_gas():
    return ct.Solution(os.path.join(script_dir, os.pardir, 'chemkin', 'chem.cti'))


def ignition(pressure, temperature, phi, endtime=None, sample_times=None):

    def init():
        gas = get_gas()
        gas.TPY = temperature, pressure, 'N2:1'
        gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')

        reac = ct.IdealGasConstPressureReactor(gas)
        net = ct.ReactorNet([reac])
        net.max_steps = int(1e9)
        lowtol = False
        if lowtol:
            print('WARNING: LOW TOLERANCES IN PLACE')
            net.atol = 1e-10
            net.rtol = 1e-06
        else:
            net.atol = 1e-20
            net.rtol = 1e-15
        return gas, reac, net

    state = []
    if endtime is not None:
        gas, reac, net = init()
        state = [np.array([0, temperature] + list(gas.Y[:]))]
        while net.time <= endtime:
            net.step()
            state.append(np.array([net.time, reac.thermo.T] + list(
                reac.thermo.Y[:])))
    elif sample_times is not None:
        for time in sample_times:
            gas, reac, net = init()
            net.advance(time)
            assert np.isclose(time, net.time)
            state.append(np.array([net.time, reac.thermo.T] + list(
                reac.thermo.Y[:])))

    return np.array(state)


def interp(time, ct_t, ct_temp):
    # find closest time
    index = np.where(ct_t >= time)[0][0]
    T1 = ct_temp[index - 1]
    T2 = ct_temp[index]
    t1 = ct_t[index - 1]
    t2 = ct_t[index]
    return T1 + (T2 - T1) * (time - t1) / (t2 - t1)


def err_norm(test, ct, norm=2):
    # calculate err
    err = np.abs(test[:, 1:] - ct[:, 1:]) / np.abs(1e-10 + ct[:, 1:])
    err = np.linalg.norm(err, ord=norm, axis=0)
    return err


def truncate_cmap(cmap, n_min=0, n_max=256, n=4):
    """ Generate a truncated colormap

    Colormaps have a certain number of colors (usually 255) which is
    available via cmap.N

    This is a simple function which returns a new color map from a
    subset of the input colormap. For example, if you only want the jet
    colormap up to green colors you may use
    tcmap = truncate_cmap(plt.cm.jet,n_max=150)
    This function is especially useful when you want to remove the white
    bounds of a colormap

    Parameters
    cmap : plt.matplotlib.colors.Colormap
    n_min : int
    n_max : int

    Return
    truncated_cmap : plt.matplotlib.colors.Colormap

    https://gist.github.com/astrodsg/09bfac1b68748967ed8b#file-mpl_colormap_tools
    """
    color_index = (n_min + (n_max - n_min) * np.linspace(0, 1, n)).astype(int)
    colors = cmap(color_index)
    name = "truncated_{}".format(cmap.name)
    return plt.matplotlib.colors.ListedColormap(colors, name=name)


def get_color(index):
    ncolors = 4
    color_wheel = plt.cm.inferno
    crange = 210
    color_wheel = truncate_cmap(color_wheel, n_max=crange, n=ncolors)
    return color_wheel(index % ncolors)


def wheel(index):
    marker_wheel = ['.', 'o', 'v', 's']
    size_wheel = [12, 15]
    style = {'color': get_color(index)}
    if index:
        style.update({
            'linestyle': '',
            'marker': marker_wheel[(index - 1) % len(marker_wheel)],
            'markersize': size_wheel[(index - 1) % len(size_wheel)],
            'markerfacecolor': 'none'})

    return style


def label(to_plot):
    if to_plot == 'T':
        return 'Temperature (K)', False
    else:
        return r'$Y_{{\ce{{{}}}}}$'.format(to_plot), True


def arc_length(x, y, aspect=None, interval=None, ylog=False):
    """
    Based on simple arc length algorithm (see citation below).

    Modified to:
        - take an optional argument `interval` that if specified, will prompt
        the method to return the indicies most closely matching arc-lengths in evenly
        spaced intervals
        - Consider the normalized arc-length (that-is, by the max / min y & x vals)
        - Consider the figure aspect ratio, if specified (width / height), and scale
          the normalized arc-length by this if supplied
        - Consider the y-log scale, if specified

    @MISC {19385,
    TITLE = {Calculate contour line length},
    AUTHOR = {nicoguaro (https://scicomp.stackexchange.com/users/9667/nicoguaro)},
    HOWPUBLISHED = {Computational Science Stack Exchange},
    NOTE = {URL:https://scicomp.stackexchange.com/q/19385 (version: 2015-04-21)},
    EPRINT = {https://scicomp.stackexchange.com/q/19385},
    URL = {https://scicomp.stackexchange.com/q/19385}
    }
    """

    xscale = 1 / np.abs(np.max(x) - np.min(x))
    yscale = 1 / np.abs(np.max(y) - np.min(y))
    if aspect:
        xscale *= aspect

    def get_arc(point):
        ax = ((x[point] - x[point - 1]) * xscale)**2
        ay = ((y[point] - y[point - 1]) * yscale)**2
        return np.sqrt(ax + ay**2)

    indicies = []
    last_arc = 0
    npts = x.size
    arc = 0
    for k in range(0, npts):
        arc = arc + get_arc(k)
        if interval is not None and arc >= interval + last_arc:
            indicies.append(k)
            last_arc = arc

    return arc if interval is None else np.array(indicies)


def find_nearest(array, values):
    # https://stackoverflow.com/a/53111759/1667311
    array = np.asarray(array)

    # the last dim must be 1 to broadcast in (array - values) below.
    values = np.expand_dims(values, axis=-1)

    indicies = np.abs(array - values).argmin(axis=-1)

    return indicies


def slice(times, arr, skip, percent_diff=10):
    last = arr[0]
    indicies = [0]
    for i, point in enumerate(arr):
        if i - indicies[-1] >= skip or (100. * np.abs(
                last - point) / last) > percent_diff:
            indicies.append(i)
            last = point

    indicies = np.array(indicies)
    return times[indicies], arr[indicies]


def sample(times, arr, base_times, base, npoints=75, aspect_ratio=8.0 / 6.0,
           thin_mode='slice'):

    if thin_mode == 'slice':
        skip = np.maximum(int(times.size / npoints), 1)
        return slice(times, arr, skip)
    else:
        interval = arc_length(base_times, base, aspect=aspect_ratio) / npoints
        indicies = arc_length(base_times, base, interval=interval,
                              aspect=aspect_ratio)
        bvals = base[indicies]
        # find the points in array that are closest to the base vals
        indicies = find_nearest(arr, bvals)
        return times[np.array(indicies)], arr[np.array(indicies)]


def err(pressure, temperature, phi, endtime, to_plot,
        draw_samples, thin_mode, do_err=True):
    gas = get_gas()
    acc = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.accelode'.format(
                     int(pressure), int(temperature), float(phi))),
                     delimiter='\t', skiprows=1).reshape((-1, gas.n_species + 2))
    OF = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.foamode'.format(
                    int(pressure), int(temperature), float(phi))),
                    delimiter='\t', skiprows=1).reshape((-1, gas.n_species + 2))

    sample_times = acc[:, 0]
    assert np.allclose(sample_times, OF[:, 0])
    if do_err:
        # get temperature comparison
        comp_ct = ignition(pressure, temperature, phi,
                           sample_times=sample_times)
        errs = {}
        errs['n_samples'] = OF[:, 0].size
        errs['of_err_mean'] = err_norm(OF, comp_ct, norm=2)
        errs['ai_err_mean'] = err_norm(acc, comp_ct, norm=2)
        errs['of_err_inf'] = err_norm(OF, comp_ct, norm=np.inf)
        errs['ai_err_inf'] = err_norm(acc, comp_ct, norm=np.inf)
        filename = '{}_{}_{}.npz'.format(int(pressure / ct.one_atm),
                                         int(temperature), float(phi))
        np.savez(filename, **errs)

    # load plotter's
    phi_ct = ignition(pressure, temperature, phi, endtime=endtime)
    acc = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.accel'.format(
                     int(pressure), int(temperature), float(phi))),
                     delimiter='\t', skiprows=1).reshape((-1, gas.n_species + 2))
    OF = np.loadtxt(os.path.join(script_dir, '{}_{}_{}.foam'.format(
                    int(pressure), int(temperature), float(phi))),
                    delimiter='\t', skiprows=1).reshape((-1, gas.n_species + 2))

    for val in to_plot:
        if draw_samples:
            for x in sample_times:
                plt.axvline(x=x, linestyle='--', color=get_color(3))
        ind = (['T'] + gas.species_names).index(val)
        ylabel, ylog = label(val)
        plt.plot(phi_ct[:, 0], phi_ct[:, 1 + ind],
                 label='Cantera', **wheel(0))
        plt.plot(*sample(acc[:, 0], acc[:, 1 + ind],
                         phi_ct[:, 0], phi_ct[:, 1 + ind],
                         thin_mode=thin_mode),
                 label=r'\texttt{accelerInt}', **wheel(1),
                 zorder=2)
        plt.plot(*sample(OF[:, 0], OF[:, 1 + ind],
                         phi_ct[:, 0], phi_ct[:, 1 + ind],
                         thin_mode=thin_mode),
                 label=r'\texttt{OpenFOAM}', **wheel(2),
                 zorder=1)
        plt.xlabel('Time (s)')
        plt.ylabel(ylabel)
        if ylog:
            plt.yscale('log')
        plt.legend(**{'loc': 0,
                      'fontsize': 16,
                      'numpoints': 1,
                      'shadow': True,
                      'fancybox': True})
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.tick_params(axis='both', which='minor', labelsize=14)
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
    parser.add_argument('-tp', '--to_plot',
                        nargs='+',
                        default=['T', 'CH4', 'OH', 'HO2', 'NO'])
    parser.add_argument('-s', '--draw_sample_times',
                        action='store_true',
                        required=False,
                        default=False,
                        help='Draw the samples tested for the error-norm.')
    parser.add_argument('-th', '--thin_mode',
                        choices=['arclen', 'slice'],
                        required=False,
                        default='slice',
                        help='The method used to thin the data such that the plot '
                             'is more readable.')
    parser.add_argument('-n', '--no_error',
                        action='store_true',
                        required=False,
                        default=False,
                        help='Plot only.')
    args = parser.parse_args()
    err(args.pressure, args.temperature, args.phi, args.endtime, args.to_plot,
        args.draw_sample_times, args.thin_mode, not args.no_error)

    sys.exit(0)
