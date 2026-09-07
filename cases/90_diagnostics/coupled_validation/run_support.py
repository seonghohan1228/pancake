"""Shared launch controls and machine provenance for validation campaigns."""
import os
import platform
import socket
import subprocess


def add_launch_arguments(parser):
    parser.add_argument('--ranks', type=int, default=1, help='MPI ranks per case')
    parser.add_argument('--mpi-launcher', default='mpiexec')
    parser.add_argument('--timeout', type=float, default=1800, help='Seconds per case')


def machine_provenance(args):
    if args.ranks < 1 or args.timeout <= 0:
        raise ValueError('ranks and timeout must be positive')
    return dict(hostname=socket.gethostname(), platform=platform.platform(),
                cpu_count=os.cpu_count(), ranks=args.ranks,
                mpi_launcher=args.mpi_launcher if args.ranks > 1 else None,
                omp_threads=1, openblas_threads=1)


def run_case(args, config, log_path):
    command = [str(args.binary.resolve()), '-c', str(config)]
    if args.ranks > 1:
        command = [args.mpi_launcher, '-n', str(args.ranks)] + command
    with log_path.open('w') as log:
        try:
            return subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'),
                timeout=args.timeout).returncode
        except subprocess.TimeoutExpired:
            log.write(f'\nCampaign timeout after {args.timeout:g} seconds.\n')
            return 124
