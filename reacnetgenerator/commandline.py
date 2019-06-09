import argparse
from . import __version__

def _commandline():
    parser = argparse.ArgumentParser(
        description=f'ReacNetGenerator {__version__}')
    parser.add_argument(
        '-i', '--inputfilename', nargs='*',
        help='Input trajectory file, e.g. bonds.reaxc', required=True)
    parser.add_argument('-a', '--atomname',
                        help='Atomic names in the trajectory, e.g. C H O',
                        nargs='*', required=True)
    parser.add_argument(
        '--nohmm', help='Process trajectory without Hidden Markov Model (HMM)',
        action="store_true")
    parser.add_argument(
        '--dump', help='Process the LAMMPS dump file', action="store_true")
    parser.add_argument(
        '-n', '-np', '--nproc', help='Number of processes', type=int)
    parser.add_argument(
        '-s', '--selectatoms',
        help='Select atoms in the reaction network, e.g. C', nargs='*')
    parser.add_argument(
        '--stepinterval', help='Step interval', type=int, default=1)
    parser.add_argument(
        '--split', help='Split number for the time axis', type=int, default=1)
    parser.add_argument(
        '--maxspecies', help='Max number of nodes (species) in the network', type=int, default=20)
    args = parser.parse_args()
    from reacnetgen import ReacNetGenerator
    ReacNetGenerator(
        inputfilename=args.inputfilename, atomname=args.atomname,
        runHMM=not args.nohmm,
        inputfiletype=('lammpsdumpfile' if args.dump else 'lammpsbondfile'),
        nproc=args.nproc, selectatoms=args.selectatoms,
        stepinterval=args.stepinterval,
        split=args.split,
        maxspecies=args.maxspecies,
    ).runanddraw()
