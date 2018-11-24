import os
import subprocess
import argparse
import numpy as np
import matplotlib.pyplot as plt


skeleton = (r"""
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
start   (0 0 -0.1);
end     (0 0 0.5);
fields  ({fields});
setFormat   csv;


// Sampling and I/O settings
#includeEtc "caseDicts/postProcessing/graphs/sampleDict.cfg"

setConfig
{{
    nPoints         50;
}}

// Must be last entry
#includeEtc "caseDicts/postProcessing/graphs/graph.cfg"

// ************************************************************************* //
""").strip()


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


def extract(fields, timelist=[], caselist=[]):
    home = os.getcwd()

    def _make_times():
        return "{}".format(','.join([str(x) for x in timelist]))
    for case in valid(home, caselist=caselist):
        # try to extract
        os.chdir(case)
        try:
            with open(os.path.join('system', 'extractAxial'), 'w') as file:
                file.write(skeleton.format(fields=' '.join(fields)))

            # set the extractor
            subprocess.check_call(['foamDictionary', '-entry', 'functions',
                                   '-set', '{#includeFunc extractAxial}',
                                   'system/controlDict'])

            # reconstruct our desired times
            call = ['reconstructPar', '-fields', '({})'.format(' '.join(fields)),
                    '-newTimes']
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
                    '_'.join(fields))), sep='\n')
                vals = vals.reshape((-1, len(fields) + 1))
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
    for i, field in enumerate(fields):
        for time in sorted(timev):
            for j, case in enumerate(results):
                if time not in results[case]:
                    continue
                vals = results[case][time]
                if not vals.size:
                    continue
                plt.semilogy(vals[:, 0], vals[:, 1 + i], label=case,
                             linestyle='',
                             marker=marker_wheel[j % len(marker_wheel)],
                             markersize=size_wheel[j % len(size_wheel)],
                             markerfacecolor='none',
                             color=color_wheel(j % len(size_wheel)))
            plt.title(field + ' comparison at time {}s'.format(time))
            plt.legend(loc=0)
            plt.xlabel('Z-position (m)')
            plt.ylabel(field)
            plt.savefig(os.path.join('figs', '{field}_{time}.pdf'.format(
                time=time, field=field)))
            if show:
                plt.show()
            plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser('valid.py - extract / plot data from Sandia '
                                     'Flame-D for valdation.')
    parser.add_argument('-f', '--fields',
                        nargs='+',
                        default=['T', 'p', 'CH4', 'CO', 'CO2', 'N2', 'O',
                                 'H', 'OH', 'HO2', 'H2O'],
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

    args = parser.parse_args()

    if args.extract:
        extract(args.fields, args.times, args.cases)
    if args.plot:
        plot(args.fields, args.times, args.show)
