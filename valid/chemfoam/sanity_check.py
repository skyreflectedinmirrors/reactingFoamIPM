from argparse import ArgumentParser
import sys
import re
import cantera as ct
import numpy as np


if __name__ == '__main__':
    parser = ArgumentParser('sanity-check.py')
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
    parser.add_argument('-log', '--log_file',
                        type=str,
                        help='The log file to check the state of.',
                        required=True)

    args = parser.parse_args()
    gas = ct.Solution('chemkin/chem.cti')
    gas.TP = args.temperature, args.pressure
    gas.set_equivalence_ratio(args.phi, 'CH4', 'O2:1, N2:3.76')

    # open log file
    with open(args.log_file, 'r') as file:
        file = file.readlines()

    pres_re = re.compile(r'^Pressure=([^\n]+)$')
    temp_re = re.compile(r'^Temperature=([^\n]+)$')
    dens_re = re.compile(r'^rho=([^\n]+)$')
    conc_re = re.compile(r'^concentration\[(\d+)\]=([^\n]+)$')

    checked = 0
    for line in file:
        if pres_re.search(line):
            checked += 1
            pressure = float(pres_re.search(line).group(1))
            assert np.isclose(pressure, gas.P)
        if temp_re.search(line):
            checked += 1
            temperature = float(temp_re.search(line).group(1))
            assert np.isclose(temperature, gas.T)
        if dens_re.search(line):
            checked += 1
            density = float(dens_re.search(line).group(1))
            assert np.isclose(density, gas.density)
        if conc_re.search(line):
            checked += 1
            match = conc_re.search(line)
            index, conc = int(match.group(1)), float(match.group(2))
            assert np.isclose(conc, gas.concentrations[index])

    assert checked == 3 + gas.n_species
    sys.exit(0)
