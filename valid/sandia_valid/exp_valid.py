import os
import subprocess
import argparse
import numpy as np
import matplotlib as mpl
from string import Template
# setup latex
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rc('text.latex',
       preamble=r'\usepackage{amsmath},\usepackage{siunitx},'
                r'\usepackage[version=4]{mhchem}')
mpl.rc('font', family='serif')
import matplotlib.pyplot as plt  # noqa


axial0 = "0 0 0.0"
axial1 = "0 0 0.5"

vertical0 = "0    0  {z}"
vertical1 = "0.15 0  {z}"

skeleton = Template(r"""
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Web:      www.OpenFOAM.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Writes graph data for specified fields along a line, specified by start
    and end points.

\*---------------------------------------------------------------------------*/

// Extract state of solution along a slice of the flame
start   (${start});
end     (${end});
fields  ${fields};
setFormat   csv;


// Sampling and I/O settings
#includeEtc "caseDicts/postProcessing/graphs/sampleDict.cfg"

setConfig
{
    nPoints         100;
}

// Must be last entry
#includeEtc "caseDicts/postProcessing/graphs/graph.cfg"

// ************************************************************************* //
""".strip())


def valid(home=os.getcwd(), caselist=[]):
    if not caselist:
        caselist = os.listdir(home)
    for case in caselist:
            if not os.path.isdir(os.path.join(os.getcwd(), case)):
                continue
            yield case


def times(case, xslice='extractAxial'):
    path = os.path.join(case, 'postProcessing', xslice)
    for time in os.listdir(path):
        if not os.path.isdir(os.path.join(path, time)):
            continue
        yield os.path.join(path, time)


# required fields to use the app for post-processing:
req_fields = set(['U', 'alphat', 'nut', 'k', 'epsilon', 'G'])


def _make_full_fields(fields, for_extract=True):
    subtract = set()
    return sorted((req_fields | set(fields)) - subtract)


def _make_fields(fields, for_extract=True):
    f = _make_full_fields(fields, for_extract=for_extract)
    if for_extract:
        return '({})'.format(' '.join(f))
    else:
        return '_'.join(f)


def _field_iter(fields, for_extract=True):
    """
    An iterable to parse the fields into file names that aren't too long
    for extraction
    """
    field_list = fields[:]
    while field_list:
        fl_len = min(10, len(field_list))
        yield field_list[:fl_len]
        field_list = field_list[fl_len:]


def _field_index(field, fields, for_extract=True):
    return _make_full_fields(fields, for_extract=for_extract).index(field)


def _num_fields(fields, for_extract=True):
    return len(_make_full_fields(fields, for_extract=for_extract))


def _get_slices(zlist=[7.5, 30, 45]):
    return ['extractAxialVelocity'] + ['extractVerticalVelocity_{:2g}'.format(
        z) for z in zlist]


def extract(fields, r, xlist, timelist=[], caselist=[], force=False):
    home = os.getcwd()

    def _make_times():
        return "{}".format(','.join([str(x) for x in timelist]))

    slices = _get_slices()
    for case in valid(home, caselist=caselist):
        # try to extract
        os.chdir(case)
        try:
            # reconstruct our desired times
            call = ['reconstructPar']
            if not force:
                call += ['-newTimes']
            if timelist:
                call += ['-time', _make_times()]
            subprocess.check_call(call)

            for f in _field_iter(fields):
                # write the radial slices
                with open(os.path.join(
                        'system', 'extractAxialVelocity'), 'w') as file:
                    assert os.path.basename(file.name) in _get_slices()
                    file.write(skeleton.substitute(start=axial0,
                                                   end=axial1,
                                                   fields=_make_fields(f)))

                # write the vertical slices
                for z in zlist:
                    z_real = z * r
                    with open(os.path.join('system',
                                           'extractVerticalVelocity_{:2g}'.format(
                                            z)),
                              'w') as file:
                        assert os.path.basename(file.name) in _get_slices()
                        file.write(skeleton.substitute(
                            start=vertical0.format(z=z_real),
                            end=vertical1.format(z=z_real),
                            fields=_make_fields(f)))

                includeFunc = ''
                for file in sorted(slices):
                    if includeFunc:
                        includeFunc += ' '
                    includeFunc += '#includeFunc {}'.format(file)

                # set the extractor
                subprocess.check_call(['foamDictionary', '-entry', 'functions',
                                       '-set', '{{{}}}'.format(includeFunc),
                                       'system/controlDict'])

                app = subprocess.check_output([
                    'foamDictionary', '-entry', 'application',
                    '-value', 'system/controlDict']).strip()

                # chdir
                call = [app, '-postProcess']
                if timelist:
                    call += ['-time', _make_times()]
                subprocess.check_call(call)
        except FileNotFoundError:
            pass
        except subprocess.CalledProcessError:
            pass
        finally:
            os.chdir(home)


name_map = {'SandiaD_LTS': r'OF (\texttt{ROS4})',
            'SandiaD_LTS_seulex': r'OF (\texttt{Seulex})',
            'SandiaD_LTS_accelerint': r'AI (\texttt{ROS4})'}


def fieldnames(field):
    defaults = {'U': r'$\frac{U}{U_\infty}$'}
    return defaults[field]


def yscale():
    U0 = 49.6  # m / s
    return 1. / U0


def xscale(inv=True):
    r = 7.2e-3  # mm
    if inv:
        return 1. / r
    return r


def limits(xslice):
    xmax = 70
    if '7.5' in xslice:
        xmax = 3
    elif '30' in xslice:
        xmax = 7
    elif '45' in xslice:
        xmax = 8
    return (0, xmax)


def islog(field):
    defaults = {'U': False}
    if field in defaults:
        return defaults[field]
    return True


def _timecheck(t1, t2):
    return np.isclose(t1, t2, atol=1e-10, rtol=1e-10)


def load(fields, timelist, cases):
    results = {}
    home = os.getcwd()
    timev = set(timelist)
    for case in valid(home, caselist=cases):
        nicecase = os.path.basename(case)
        results[nicecase] = {}
        try:
            for xslice in _get_slices():
                results[nicecase][xslice] = {}
                # load all data
                for time in times(case, xslice):
                    t = float(os.path.basename(time))
                    if timelist and not any(_timecheck(x, t) for x in timelist):
                        continue

                    composite_fields = None
                    vals = np.fromfile(os.path.join(time, 'line_U.xy'), sep='\n')
                    # have three velocity components
                    vals = vals.reshape((-1, 4))
                    # want Uz
                    indicies = np.array([3])
                    svals = vals[:, indicies]
                    if composite_fields is None:
                        composite_fields = np.hstack((
                            np.expand_dims(vals[:, 0], 1), svals))
                    else:
                        composite_fields = np.hstack((composite_fields, svals))

                    nice_t = next(x for x in timelist if _timecheck(x, t))
                    results[nicecase][xslice][nice_t] = composite_fields
                    timev.add(float(nice_t))
                if not results[nicecase][xslice]:
                    del results[nicecase][xslice]
        except FileNotFoundError:
            pass
        if not results[nicecase]:
            del results[nicecase]

    return timev, results


def plot(timev, results, grey=False, base='SandiaD_LTS_accelerint'):
    marker_wheel = ['.', 'o', 'v', 's']
    size_wheel = [6, 10, 14, 16]
    cmap = 'Greys' if grey else 'inferno'
    color_wheel = plt.get_cmap(cmap, len(size_wheel) + 1)
    try:
        os.mkdir('figs')
    except OSError:
        pass
    for index, field in enumerate(sorted(fields)):
        for time in sorted(timev):
            fig = plt.figure(figsize=(16, 4))
            lines = []
            last = []
            axes = []
            for i, xslice in enumerate(_get_slices()):
                ax0 = fig.add_subplot(1, len(_get_slices()), i + 1)
                axes.append(ax0)

                for j, case in enumerate(results):
                    if xslice not in results[case]:
                        continue
                    if time not in results[case][xslice]:
                        continue
                    vals = results[case][xslice][time]
                    if not vals.size:
                        continue
                    x = vals[:, 0] * xscale()
                    y = vals[:, 1 + index] * yscale()
                    label = case
                    if case in name_map:
                        label = name_map[case]
                    plotter = ax0.semilogy if islog(field) else ax0.plot
                    line = plotter(x, y, label=label,
                                   linestyle='',
                                   marker=marker_wheel[j % len(marker_wheel)],
                                   markersize=size_wheel[j % len(size_wheel)],
                                   markerfacecolor='none',
                                   color=color_wheel(j % len(size_wheel)))
                    if not i:
                        lines.extend(line)

                plt.ylim((0, 1.3))
                plt.xlim(*limits(xslice))
                labels = [l.get_label() for l in lines + last]
                plt.legend(lines + last, labels,
                           **{'loc': 0,
                              'fontsize': 16,
                              'numpoints': 1,
                              'shadow': True,
                              'fancybox': True})
            for ax in axes:
                ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
                ax.tick_params(axis='both', which='major', labelsize=20)
                ax.tick_params(axis='both', which='minor', labelsize=16)
                for item in (ax.title, ax.xaxis.label,
                             ax.yaxis.label):
                    item.set_fontsize(24)
                if 'Vertical' in xslice:
                    ax.set_xlabel(r'$\frac{r}{d}$')
                else:
                    ax.set_xlabel(r'$\frac{x}{d}$')
                ax.set_ylabel(fieldnames(field))
            plt.tight_layout()
            plt.savefig(os.path.join('figs', 'exp_velocity_comp_{time}.pdf'.format(
                time=time)), bbox_inches='tight')
            plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        'exp_valid.py - extract / plot velocity data for comparison to '
        'experimental results for Sandia Flame-D for valdation.')

    parser.add_argument('-p', '--plot',
                        action='store_true',
                        default=False,
                        help='Plot the extracted fields for the various solutions.')

    parser.add_argument('-e', '--extract',
                        action='store_true',
                        default=False,
                        help='Plot the extracted fields for the various solutions.')

    parser.add_argument('-t', '--times',
                        nargs='+',
                        default=[1500, 5000, 5000.01],
                        help='The times to plot / extract.')

    parser.add_argument('-c', '--cases',
                        nargs='+',
                        default=['SandiaD_LTS', 'SandiaD_LTS_accelerint',
                                 'SandiaD_LTS_seulex'],
                        help='The cases to process.')

    parser.add_argument('-r', '--force_reextraction',
                        action='store_true',
                        default=False,
                        help='If specified, re-extract the post-processing data '
                             'regardless of whether it already exists or not.')

    parser.add_argument('-v', '--validate',
                        action='store_true',
                        default=False,
                        help='If specified, emit a norm-based evaluation of the '
                             'differences between the various solutions.')

    args = parser.parse_args()

    fields = ['U']
    zlist = [7.5, 30, 45]
    r = xscale(inv=False)
    if args.extract:
        extract(fields, r, zlist, args.times, args.cases, args.force_reextraction)
    t = None
    results = None
    if args.plot or args.validate:
        t, results = load(fields, args.times, args.cases)
    if args.plot:
        plot(t, results)
