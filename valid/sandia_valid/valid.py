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
import matplotlib.pyplot as plt


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

// axial velocity varying along y-direction, useful for monitoring solution
start   (0.0091 0 0.0);
end     (0.0091 0 0.5);
fields  ${fields};
setFormat   csv;


// Sampling and I/O settings
#includeEtc "caseDicts/postProcessing/graphs/sampleDict.cfg"

setConfig
{
    nPoints         50;
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


def times(case):
    path = os.path.join(case, 'postProcessing', 'extractAxial')
    for time in os.listdir(path):
        if not os.path.isdir(os.path.join(path, time)):
            continue
        yield os.path.join(path, time)


# required fields to use the app for post-processing:
req_fields = set(['U', 'alphat', 'nut', 'k', 'epsilon', 'G'])


def _make_full_fields(fields, for_extract=True):
    subtract = set()
    if not for_extract:
        subtract = set(['U'])
    return sorted((req_fields | set(fields)) - subtract)


def _make_fields(fields, for_extract=True):
    f = _make_full_fields(fields, for_extract=for_extract)
    if for_extract:
        return '({})'.format(' '.join(f))
    else:
        return '_'.join(f)


def _field_index(field, fields, for_extract=True):
    return _make_full_fields(fields, for_extract=for_extract).index(field)


def _num_fields(fields, for_extract=True):
    return len(_make_full_fields(fields, for_extract=for_extract))


def extract(fields, timelist=[], caselist=[], force=False):
    home = os.getcwd()

    def _make_times():
        return "{}".format(','.join([str(x) for x in timelist]))

    for case in valid(home, caselist=caselist):
        # try to extract
        os.chdir(case)
        try:
            with open(os.path.join('system', 'extractAxial'), 'w') as file:
                file.write(skeleton.substitute(fields=_make_fields(fields)))

            # set the extractor
            subprocess.check_call(['foamDictionary', '-entry', 'functions',
                                   '-set', '{#includeFunc extractAxial}',
                                   'system/controlDict'])

            # reconstruct our desired times
            call = ['reconstructPar']
            if not force:
                call += ['-newTimes']
            if timelist:
                call += ['-time', _make_times()]
            subprocess.check_call(call)

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
        finally:
            os.chdir(home)


name_map = {'SandiaD_LTS': 'OF (ROS4)',
            'SandiaD_LTS_seulex': 'OF (seulex)',
            'SandiaD_LTS_accelerint': 'AI (ROS4)'}


def fieldnames(field):
    defaults = {'p': 'Pressure (Pa)',
                'T': 'Temperature (K)'}
    if field in defaults:
        return defaults[field]

    return Template(r'$$\text{Y}_{\ce{${field}}}$$').substitute(field=field)


def limits(field):
    defaults = {}
    if field in defaults:
        return defaults[field]

    return (None, None)


def islog(field):
    defaults = {'T': False,
                'p': False,
                'N2': False}
    if field in defaults:
        return defaults[field]
    return True


def plot(fields, timelist, show, grey=False):

    def _timecheck(t1, t2):
        return np.isclose(t1, t2, atol=1e-10, rtol=1e-10)
    marker_wheel = ['.', 'o', 'v', 's']
    size_wheel = [6, 10, 14, 16]
    cmap = 'Greys' if grey else 'inferno'
    color_wheel = plt.get_cmap(cmap, len(size_wheel) + 1)
    home = os.getcwd()
    results = {}
    timev = set(timelist)
    for case in valid(home):
        nicecase = os.path.basename(case)
        results[nicecase] = {}
        try:
            # load all data
            for time in times(case):
                t = float(os.path.basename(time))
                if timelist and not any(_timecheck(x, t) for x in timelist):
                    continue
                vals = np.fromfile(os.path.join(time, 'line_{}.xy'.format(
                    _make_fields(fields, for_extract=False))), sep='\n')
                vals = vals.reshape((-1, _num_fields(fields, for_extract=False) + 1))
                nice_t = next(x for x in timelist if _timecheck(x, t))
                results[nicecase][nice_t] = vals
                timev.add(float(nice_t))
        except FileNotFoundError:
            pass
        if not results[nicecase]:
            del results[nicecase]

    try:
        os.mkdir('figs')
    except OSError:
        pass
    for field in fields:
        for time in sorted(timev):
            ymin = np.finfo(np.float64).max
            ymax = -np.finfo(np.float64).max
            for j, case in enumerate(results):
                if time not in results[case]:
                    continue
                vals = results[case][time]
                if not vals.size:
                    continue
                index = _field_index(field, fields, for_extract=False)
                if case in name_map:
                    case = name_map[case]

                ymin = np.minimum(ymin, np.min(vals[:, 1 + index]))
                ymax = np.maximum(ymax, np.max(vals[:, 1 + index]))
                plotter = plt.semilogy if islog(field) else plt.plot
                plotter(vals[:, 0], vals[:, 1 + index], label=case,
                        linestyle='',
                        marker=marker_wheel[j % len(marker_wheel)],
                        markersize=size_wheel[j % len(size_wheel)],
                        markerfacecolor='none',
                        color=color_wheel(j % len(size_wheel)))

            plt.ylim(*limits(field))
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
            plt.xlabel('Axial-position (m)')
            plt.ylabel(fieldnames(field))
            plt.tight_layout()
            plt.savefig(os.path.join('figs', '{field}_{time}.pdf'.format(
                time=time, field=field)))
            if show:
                plt.show()
            plt.close()
    return timev, results


def validate(times, results, fields, base='SandiaD_LTS_seulex'):
    indicies = []
    for field in fields:
        index = _field_index(field, fields, for_extract=False)
        indicies.append(index)
    indicies = np.array(indicies, dtype=np.int32)
    for time in times:
        comp = results[base][time]
        for case in results:
            if case == base:
                continue
            if time not in results[case]:
                continue
            test = results[case][time]
            diff = np.abs(comp[:, 1 + indicies] - test[:, 1 + indicies]) / (
                1e-30 + np.abs(comp[:, 1 + indicies]))
            diff = np.linalg.norm(diff, ord=2, axis=1)
            rel_err = np.linalg.norm(diff, ord=np.inf, axis=0) * 100.

            diff = np.abs(comp[:, 1 + indicies] - test[:, 1 + indicies]) / (
                1e-10 + 1e-06 * np.abs(comp[:, 1 + indicies]))
            diff = np.linalg.norm(diff, ord=2, axis=1)
            weighted_err = np.linalg.norm(diff, ord=np.inf, axis=0)
            print(case, time, "rel:", rel_err, "%", "weighted:", weighted_err)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('valid.py - extract / plot data from Sandia '
                                     'Flame-D for valdation.')
    parser.add_argument('-f', '--fields',
                        nargs='+',
                        default=['T', 'p', 'CH4', 'CO', 'CO2', 'N2', 'O',
                                 'H', 'OH', 'HO2', 'H2O', 'NO', 'CH2O', 'HCO',
                                 'N2O'],
                        type=str,
                        required=False,
                        help='The fields to extract / plot')

    parser.add_argument('-p', '--plot',
                        action='store_true',
                        default=False,
                        help='Plot the extracted fields for the various solutions.')

    parser.add_argument('-e', '--extract',
                        action='store_true',
                        default=False,
                        help='Plot the extracted fields for the various solutions.')

    parser.add_argument('-s', '--show',
                        action='store_true',
                        default=False,
                        help='If true, show the plots before saving to file '
                             '[warning: this will generate _many_ plots]')

    parser.add_argument('-t', '--times',
                        nargs='+',
                        default=[5000, 5000.01],
                        type=int,
                        help='The times to plot / extract.')

    parser.add_argument('-c', '--cases',
                        nargs='+',
                        default=None,
                        type=str,
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

    if args.extract:
        extract(args.fields, args.times, args.cases, args.force_reextraction)
    t = None
    results = None
    if args.plot:
        t, results = plot(args.fields, args.times, args.show)
    if args.validate:
        assert args.plot, "must plot to validate"
        validate(t, results, args.fields)
