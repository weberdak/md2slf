#! /usr/bin/env python

import nmrglue as ng
import numpy as np
import argparse

# Written by: Daniel K. Weber, University of Minnesota, MN, USA
# Script has been adapted from NMRGlue documentation.

def parse_args():
    parser = argparse.ArgumentParser(description='Simulate 15N-1H SLF spectrum from SPARKY peak list file.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i','--infile', type=str,
        help='Peak list file in SPARKY format.'
    )
    parser.add_argument(
        '-o','--outfile', type=str,
        help='Filename of output spectrum.'
    )
    parser.add_argument(
        '--sw_n', type=float,
        help='Spectral width of 15N chemical shift dimension (Hz). Default: 20000.0.',
        default=20000.0
    )
    parser.add_argument(
        '--sw_nh', type=float,
        help='Spectral width of 15N-1H dipolar coupling dimension (Hz). Default: 31250.0',
        default=31250.0
    )
    parser.add_argument(
        '--size_n', type=int,
        help='Number of points in 15N chemical shift dimension. Default: 2048.',
        default=2048
    )
    parser.add_argument(
        '--size_nh', type=int,
        help='Number of points in 15N-1H dipolar coupling dimension. Default: 1024.',
        default=1024
    )
    parser.add_argument(
        '--lw_n', type=float,
        help='Peak linewidth in 15N chemical shift dimension (Hz). Default: 100.0.',
        default=100.0
    )
    parser.add_argument(
        '--lw_nh', type=float,
        help='Peak linewidth in 15N-1H dipolar coupling dimension (Hz). Default: 100.0.',
        default=100.0
    )
    parser.add_argument(
        '--freq_n', type=float,
        help='Observe frequency of 15N (MHz). Default: 60.7639142.',
        default=60.7639142
    )
    parser.add_argument(
        '--carr_n', type=float,
        help='Carrier offset frequency of 15N (Hz). Default: 5600.0.',
        default=5600.0
    )
    args = parser.parse_args()
    return args


def main():
    
    args = parse_args()
    peaks = args.infile
    outfile = args.outfile
    
    sw_N = args.sw_n
    sw_H = args.sw_nh
    size_N = args.size_n
    size_H = args.size_nh
    lw_N_Hz = args.lw_n
    lw_H_Hz = args.lw_nh
    freq_n = args.freq_n
    carr_n = args.carr_n

    # Make NMRGlue universal dictionary
    udic = {
        'ndim': 2,
        0: {'car': 0.0,
            'complex': False,
            'encoding': 'states',
            'freq': True,
            'label': '1H',
            'obs': 1000.0,
            'size': size_H,
            'sw': sw_H,
            'time': False},
        1: {'car': carr_n,
            'complex': False,
            'encoding': 'direct',
            'freq': True,
            'label': '15N',
            'obs': freq_n,
            'size': size_N,
            'sw': sw_N,
            'time': False}
    }
    dic = ng.sparky.create_dic(udic)
    data = np.empty((size_N, size_H), dtype='float32')

    # read in the peak list
    peak_list = np.recfromtxt(peaks, names=True)
    npeaks = len(peak_list)
    #print(peak_list)

    # convert the peak list from PPM to points
    uc_H = ng.sparky.make_uc(dic, None, 0)
    uc_N = ng.sparky.make_uc(dic, None, 1)

    lw_N = (lw_N_Hz/sw_N)*size_N    # 15N dimension linewidth in points
    lw_H = (lw_H_Hz/sw_H)*size_H    # 1H dimension linewidth in points

    params = []
    for name, ppm_H, ppm_N in peak_list:
        pts_H = uc_H.f(ppm_H, 'ppm')
        pts_N = uc_N.f(ppm_N, 'ppm')
        params.append([(pts_H, lw_H), (pts_N, lw_N)])

    # simulate the spectrum
    shape = (size_H, size_N)      # size should match the dictionary size
    lineshapes = ('g', 'g')  # gaussian in both dimensions
    amps = [100.0] * npeaks
    data = ng.linesh.sim_NDregion(shape, lineshapes, params, amps)

    # save the spectrum
    ng.sparky.write(outfile, dic, data.astype('float32'), overwrite=True)


if __name__ == '__main__':
    main()
