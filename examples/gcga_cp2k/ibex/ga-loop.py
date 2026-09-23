# Ibex GCGA worker loop: ONE sbatch job = ONE worker (16 cores on a shared Turin node, never a whole node).
# Submit N copies of ibex/run-worker.sbatch for N concurrent workers (e.g. popSize per mu search, <= 1300 cores/account).
# Kid numbering is race-free across jobs: a directory is claimed by os.mkdir, which fails if another job got it first.
#   python -u ibex/ga-loop.py <totConf> <minConf>        (from ibex/run-worker.sbatch; touch STOP to end)
# Run ibex/ga-prep.py ONCE before the first worker (initializes gmid / natural selection).
import datetime
import os
import subprocess
import sys
import input
from gocia.ga.popGrandCanon import PopulationGrandCanonical

totConf = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
minConf = int(sys.argv[2]) if len(sys.argv) > 2 else 200
pop = PopulationGrandCanonical(
    gadb='gcga.db', substrate='substrate.vasp', popSize=input.popSize,
    convergeCrit=input.popSize * 10, subsPot=input.subsPot,
    chemPotDict=input.chemPotDict, zLim=input.zLim)


def claim():
    k = max([int(d) for d in os.listdir('.') if d.isdigit()] + [0]) + 1
    while True:
        try:
            os.mkdir('%06d' % k)
            return k
        except FileExistsError:
            k += 1


while 'STOP' not in os.listdir('.'):
    nkid = sum(d.isdigit() for d in os.listdir('.'))
    if nkid >= totConf or (pop.is_converged2() and nkid >= minConf):
        break
    k = claim()
    with open('%06d/worker.log' % k, 'w') as log:
        subprocess.run([sys.executable, '-u', '../ga-worker.py'], cwd='%06d' % k, stdout=log, stderr=subprocess.STDOUT)
    print('Job %i done @%s' % (k, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')), flush=True)
print('CONVERGED!' if pop.is_converged() else 'TERMINATED!', flush=True)
