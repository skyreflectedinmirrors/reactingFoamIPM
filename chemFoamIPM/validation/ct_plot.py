# constant pressure ignition problem for Cantera
import cantera as ct
import numpy as np
from scipy.integrate import trapz
import os
import re


def read_ics():
    constant = None
    re_con = re.compile(r'constantProperty\s+(\w+);')
    basis = None
    re_basis = re.compile(r'fractionBasis\s+(\w+);')
    pressure = None
    re_p = re.compile(r'p\s+([\d\.e+-]+);')
    temperature = None
    re_T = re.compile(r'T\s+([\d\.e+-]+);')
    re_fracs = re.compile(r'fractions\s*$')
    # species name identifier, not complete but works for us
    re_frac = re.compile(r'([A-Za-z0-9]+)\s+([\d\.e+-]+);')
    in_fracs = False
    fracs = {}

    with open(os.path.join(os.pardir, 'constant', 'initialConditions'), 'r') as file:
        ics = file.readlines()
    for line in ics:
        if re_con.search(line):
            constant = re_con.search(line).groups(1)[0]
        elif re_basis.search(line):
            basis = re_basis.search(line).groups(1)[0]
        elif re_p.search(line):
            pressure = float(re_p.search(line).groups(1)[0])
        elif re_T.search(line):
            temperature = float(re_T.search(line).groups(1)[0])
        if re_fracs.search(line):
            # search next lines for mass / mole fractions
            in_fracs = True
        if in_fracs and re_frac.search(line):
            spec, val = re_frac.search(line).groups()
            fracs[spec] = float(val)
        if re.search('}', line) and in_fracs:
            in_fracs = False

    return temperature, pressure, fracs, basis, constant


def ignition():
    gas = ct.Solution('chem.cti')
    # get initial conditions
    T, P, v, basis, con = read_ics()
    if basis == 'mole':
        gas.TPX = T, P, v
    elif basis == 'mass':
        gas.TPY = T, P, v
    else:
        raise Exception('Unknown basis {}'.format(basis))

    assert con == 'pressure' or con == 'volume', (
        "Unknown constant parameter {}".format(con))

    reac = ct.IdealGasReactor if con == 'volume' else ct.IdealGasConstPressureReactor
    reac = reac(gas)
    net = ct.ReactorNet([reac])
    net.atol = 1e-20
    net.rtol = 1e-15
    t = []
    Tv = []

    if not os.path.exists('ct.log'):
        with open('ct.log', 'w') as file:
            file.write("# Time [s]    Temperature [K]\n")

            def __write():
                file.write('{:.10e}\t{:.10e}\n'.format(net.time, reac.thermo.T))

            __write()
            t.append(0)
            Tv.append(T)
            while net.time <= 0.07:
                net.step()
                t.append(net.time)
                Tv.append(reac.thermo.T)
                __write()
    else:
        arry = np.loadtxt('ct.log', delimiter='\t', skiprows=1).reshape((-1, 2))
        t = arry[:, 0]
        Tv = arry[:, 1]
    return np.array(t), np.array(Tv)


def interp_temp(time, ct_t, ct_temp):
    # find closest time
    index = np.where(ct_t >= time)[0][0]
    T1 = ct_temp[index - 1]
    T2 = ct_temp[index]
    t1 = ct_t[index - 1]
    t2 = ct_t[index]
    return T1 + (T2 - T1) * (time - t1) / (t2 - t1)


def err_norm(array, times, temps, label):
    # calculate err
    errs = np.zeros(array[:, 0].size)
    for i, t in enumerate(array[:, 0]):
        T = interp_temp(t, times, temps)
        errs[i] = np.power((T - array[i, 1]) / T, 2.)
    print('{}_err:'.format(label), trapz(errs, array[:, 0]))


def err():
    times, temperatures = ignition()
    acc = np.loadtxt('../chemFoam.out', delimiter='\t',
                     skiprows=1).reshape((-1, 3))
    OF = np.loadtxt('../chemFoam_base.out', delimiter='\t',
                    skiprows=1).reshape((-1, 3))
    # get temperature comparison
    err_norm(OF, times, temperatures, 'OF')
    err_norm(acc, times, temperatures, 'ACC')


if __name__ == '__main__':
    err()
