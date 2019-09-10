import numpy as np
import argparse


# Written by: Daniel K. Weber, University of Minnesota, MN, USA

def parse_args():
    parser = argparse.ArgumentParser(description='Convert outputs of slf.tcl to SPARKY peak list.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i', '--infile', type=str,
        help='Input file name. An output from the slf.tcl VMD script.'
    )
    parser.add_argument(
        '-o', '--outfile', type=str,
        help='Filename to output SPARKY peak list file.'
    )
    parser.add_argument(
        '-n', '--negative',
        help='Convert to negative dipolar couplings.',
        action='store_true'
    )
    parser.add_argument(
        '-r', '--renumber', type=int,
        help='Add this number to shift residue numbers.',
    )
    args = parser.parse_args()
    return args
    
    
def main():
    
    args = parse_args()
    f = args.infile
    fo = open(args.outfile, 'w')
    fo.write('      Assignment         w1         w2  \n\n')
    shift = args.renumber
    
    resnames = np.genfromtxt(f,usecols=(0),dtype=str)
    resids = np.genfromtxt(f,usecols=(1),dtype=int)
    shifts = np.genfromtxt(f,usecols=(2),dtype=float)
    couplings = np.genfromtxt(f,usecols=(3),dtype=float)

    conversion = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "HSE": "H",
        "HID": "H",
        "HSD": "H",
        "HSP": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "ASX": "B",
        "GLX": "Z"
    }

    for resname,resid,cs,dc in zip(resnames,resids,shifts,couplings):

        if args.negative:
            dc = abs(dc)*-1
        
        assignment = conversion[resname]+str(resid+shift)+'H-N'
        assignment = assignment.rjust(17)
        dc = str(dc)
        dc = dc.rjust(11,' ')
        cs = str(cs)
        cs = cs.rjust(11,' ')
        fo.write('{}{}{} \n'.format(assignment,dc,cs))
        
    fo.close()

if __name__ == '__main__':
    main()
