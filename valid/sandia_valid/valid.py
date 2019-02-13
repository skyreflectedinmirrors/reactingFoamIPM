import os
import subprocess
import argparse
import cantera as ct
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


axial0 = "0.0091 0 0.0"
axial1 = "0.0091 0 0.5"

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
    if not for_extract:
        subtract = set(['U'])
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


def _get_slices():
    return ['extractAxial'] + ['extractVertical_{:.2}'.format(
        z) for z in np.linspace(0.1, 0.4, 4, endpoint=True)]


def extract(fields, timelist=[], caselist=[], force=False):
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
                # write the axial slice
                with open(os.path.join('system', 'extractAxial'), 'w') as file:
                    assert os.path.basename(file.name) in _get_slices()
                    file.write(skeleton.substitute(start=axial0,
                                                   end=axial1,
                                                   fields=_make_fields(f)))

                # write the vertical slices
                for z in np.linspace(0.1, 0.4, 4, endpoint=True):
                    with open(os.path.join('system',
                                           'extractVertical_{:.2}'.format(z)),
                              'w') as file:
                        assert os.path.basename(file.name) in _get_slices()
                        file.write(skeleton.substitute(start=vertical0.format(z=z),
                                                       end=vertical1.format(z=z),
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
                    for f in _field_iter(fields, for_extract=False):
                        vals = np.fromfile(os.path.join(time, 'line_{}.xy'.format(
                            _make_fields(f, for_extract=False))), sep='\n')
                        vals = vals.reshape((
                            -1, _num_fields(f, for_extract=False) + 1))
                        indicies = np.array([1 + _field_index(
                                             x, f, for_extract=False) for x in f])
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
        for xslice in _get_slices():
            for time in sorted(timev):

                fig, ax0 = plt.subplots()
                axes = [ax0]
                lines = []
                last = []

                for j, case in enumerate(results):
                    if xslice not in results[case]:
                        continue
                    if time not in results[case][xslice]:
                        continue
                    vals = results[case][xslice][time]
                    if not vals.size:
                        continue
                    label = case
                    if case in name_map:
                        label = name_map[case]
                    plotter = ax0.semilogy if islog(field) else ax0.plot
                    line = plotter(vals[:, 0], vals[:, 1 + index], label=label,
                                   linestyle='',
                                   marker=marker_wheel[j % len(marker_wheel)],
                                   markersize=size_wheel[j % len(size_wheel)],
                                   markerfacecolor='none',
                                   color=color_wheel(j % len(size_wheel)))
                    lines.extend(line)

                    if field != 'T' and case == base:
                        # plot the temperature on the other axis
                        ax1 = ax0.twinx()
                        axes.append(ax1)
                        T_ind = fields.index('T')
                        assert base in results and xslice in results[base] and \
                            time in results[base][xslice], 'Missing base-temperature'
                        base_vals = results[base][xslice][time]
                        T = ax1.plot(base_vals[:, 0], base_vals[:, 1 + T_ind],
                                     label='Temp.',
                                     linestyle='-',
                                     color='k')
                        last.extend(T)
                        ax1.set_ylabel(fieldnames('T'))
                        ax1.set_ylim(*limits('T'))

                plt.ylim(*limits(field))
                labels = [l.get_label() for l in lines + last]
                plt.legend(lines + last, labels,
                           **{'loc': 'lower left',
                              'fontsize': 16,
                              'numpoints': 1,
                              'shadow': True,
                              'fancybox': True,
                              'bbox_to_anchor': (0, 1.02, 1, 0.2),
                              'mode': 'expand',
                              'ncol': 2})
                for ax in axes:
                    ax.tick_params(axis='both', which='major', labelsize=20)
                    ax.tick_params(axis='both', which='minor', labelsize=16)
                    for item in (ax.title, ax.xaxis.label,
                                 ax.yaxis.label):
                        item.set_fontsize(24)
                if 'Axial' in xslice:
                    ax0.set_xlabel('Axial-position (m)')
                else:
                    assert 'Vertical' in xslice
                    ax0.set_xlabel('Distance from centerline (m)')
                ax0.set_ylabel(fieldnames(field))
                plt.tight_layout()
                plt.savefig(os.path.join('figs', '{field}_{slice}_{time}.pdf'.format(
                    time=time, field=field, slice=xslice)),
                    bbox_inches='tight')
                plt.close()


def validate(times, results, fields, base='SandiaD_LTS_accelerint',
             reacting_cutoff=500, axial_pad=0):

    nfields = len(fields)
    for time in times:
        for case in results:
            npoints = 0
            if case == base:
                continue
            # result store for non per-slice
            validation = {}

            def __update(key, val, inds, xslice):
                if key == 'inf':
                    # first, take the infinum over the sample axis, filtered by
                    # the supplied indicies
                    val = np.linalg.norm(val[inds, :], ord=np.inf, axis=0)

                    # compute infinum norm
                    temp = val
                    if 'inf' in validation:
                        temp = np.maximum(validation[key], val)
                    # store maximum error locations / slices
                    locs = np.where(temp == val)
                    if 'locs' not in validation:
                        validation['locs'] = np.array([xslice] * temp.size,
                                                      dtype=object)
                    else:
                        validation['locs'][locs] = xslice
                    validation['inf'] = temp
                else:
                    assert key == 'mean'
                    # first, take the squared sum over the sample axis,
                    # filtered by the supplied indicies
                    mean = np.sum(np.power(val[inds, :], 2), axis=0)
                    # and add to stored
                    if 'mean' not in validation:
                        validation['mean'] = mean
                    else:
                        validation['mean'] += mean

            def __finalize():
                # divide by total sampled points to find mean
                validation['mean'] /= npoints
                # and take sqrt to complete L2 norm
                validation['mean_scalar'] = np.sqrt(np.sum(
                    validation['mean'], axis=0) / nfields)
                validation['mean'] = np.sqrt(validation['mean'])
                # similarly, do the infimum
                validation['inf_scalar'] = np.linalg.norm(
                    validation['inf'], ord=np.inf, axis=0)
                # and store location of max err
                max_loc = np.where(
                    validation['inf_scalar'] == validation['inf'])[0][0]
                validation['locs'] = (max_loc, validation['locs'][max_loc])
                # and convert to percent
                validation['mean'] *= 100
                validation['inf'] *= 100
                validation['mean_scalar'] *= 100
                validation['inf_scalar'] *= 100

            for xslice in _get_slices():
                comp = results[base][xslice][time]
                if xslice not in results[case]:
                    continue
                if time not in results[case][xslice]:
                    continue
                test = results[case][xslice][time]
                # first, find where the temperature of the base solution >= 500
                inds = np.where(comp[:, 1 + fields.index('T')] >= reacting_cutoff)[0]
                if axial_pad and xslice == 'extractAxial':
                    inds = inds[np.where(np.logical_and(
                        inds >= axial_pad,
                        inds < inds.size - axial_pad))]
                diff = np.abs(comp[:, 1:] - test[:, 1:]) / (1e-10 + np.abs(comp[
                    :, 1:]))

                # update total sampled points
                npoints += inds.size

                # and update the norms
                __update('inf', diff, inds, xslice)
                __update('mean', diff, inds, xslice)

            # finally, finalize the mean and output
            __finalize()
            if 'T' in fields:
                ind = fields.index('T')
                print('Maximum temperature error:', validation['inf'][ind])
            if 'p' in fields:
                ind = fields.index('p')
                print('Maximum pressure error:', validation['inf'][ind])

            max_err_loc, max_err_slice = validation['locs']
            print(case, time, "rel(max, mean)", validation['inf_scalar'],
                  validation['mean_scalar'], "%",
                  'Max error for: ', fields[max_err_loc], 'in slice',
                  max_err_slice)
            print()
            print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser('valid.py - extract / plot data from Sandia '
                                     'Flame-D for valdation.')

    species = ct.Solution('gri30.cti').species_names[:]
    species[species.index('CH2(S)')] = 'CH2S'
    fields = ['T', 'p', 'Qdot'] + species
    parser.add_argument('-f', '--fields',
                        nargs='+',
                        default=fields,
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

    parser.add_argument('-t', '--times',
                        nargs='+',
                        default=[5000, 5000.01],
                        type=float,
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

    parser.add_argument('-a', '--axial_pad',
                        type=int,
                        help='The number of points to discard near the walls of the '
                             'combustion chamber for the axial validation slice',
                        default=0)

    parser.add_argument('-b', '--validation_base',
                        type=str,
                        help='The case to use as the reference answer for '
                             'validation',
                        default='SandiaD_LTS_accelerint')

    parser.add_argument('-T', '--reacting_cutoff',
                        type=float,
                        help='The temperature below which reations are considered '
                             'inactive',
                        default=500)

    args = parser.parse_args()

    if args.validation_base not in args.cases:
        parser.error('--validation_base must be one of the supplied --cases.')

    fields = sorted(args.fields)
    if args.extract:
        extract(fields, args.times, args.cases, args.force_reextraction)
    t = None
    results = None
    if args.plot or args.validate:
        t, results = load(fields, args.times, args.cases)
    if args.plot:
        plot(t, results)
    if args.validate:
        validate(t, results, fields, axial_pad=args.axial_pad,
                 base=args.validation_base, reacting_cutoff=args.reacting_cutoff)
