# Generates reads using art_illumina

import argparse
import logging
import subprocess
import time

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


def execute_art(cur_idx, art, coverage, out):
    # generating paired reads (-p) of length (-l) 100, with average insert length (-m) 350 and standard deviation
    # (-s) 50
    return subprocess.Popen([art, '-p', '-l', '100', '-m', '350', '-s', '50', '-i', args.fasta, '-f',
                             str(coverage / 2.), '-rs', str(cur_idx), '-o', f'{out}{cur_idx}'])


if __name__ == '__main__':
    args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    args.add_argument('--start', help='First cell id', default=0, type=int)
    args.add_argument('--stop', help='Last cell id, exclusive', default=100, type=int)
    args.add_argument('-p', help='Number of cells to generate in parallel', default=100, type=int)
    args.add_argument('--coverage', help='Coverage for this cell (2 pairs with half coverage will be generated)',
                      default=0.05, type=float)
    args.add_argument('--out', help='The prefix of the output file', default='./art_')
    args.add_argument('--art', help='Location of the art_illumina binary', default='art_illumina')
    args.add_argument('--fasta', help='Location of the template fasta file from which reads are generated',
                      default=None)

    args = args.parse_args()

    if args.stop - args.start < 0:
        logger.exception('Difference between start and stop must be > 0')
        exit(1)

    cur_idx = args.start
    process_list = []
    while cur_idx - args.start < args.p and cur_idx < args.stop:
        process_list.append((execute_art(cur_idx, args.art, args.coverage, args.out), cur_idx))
        cur_idx += 1

    while process_list:
        for p, idx in list(process_list):
            return_code = p.poll()
            if return_code is not None:
                process_list.remove((p, idx))
                if return_code == 0:
                    logger.info(f'File {idx} generated successfully')
                else:
                    logger.info(f'Generation for {idx} failed. Re-trying')
                if cur_idx < args.stop:
                    process_list.append((execute_art(cur_idx, args.art, args.coverage, args.out), cur_idx))
                    cur_idx += 1
        logger.info('Back to sleep')
        time.sleep(1)
