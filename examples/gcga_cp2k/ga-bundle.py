# Bundled GCGA master for Shaheen (paper Fig. 3b scheme): ONE sbatch job owns N nodes,
# this script keeps up to `nworker` ga-worker.py processes running, each launching CP2K
# through `srun --exact` inside the allocation. Replaces ga-slurm.py (one job per kid), which
# the shared-budget rule (<=10 concurrent 01__ jobs) forbids.
#
#   nohup python -u ga-bundle.py <nworker> <totConf> <minConf> &      (from run-bundle.sbatch)
#   touch STOP to end gracefully.
import os
import sys
import time
import subprocess
import datetime
import input
from gocia.ga.popGrandCanon import PopulationGrandCanonical

nworker = int(sys.argv[1]) if len(sys.argv) > 1 else 4
totConf = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
minConf = int(sys.argv[3]) if len(sys.argv) > 3 else 200

pop = PopulationGrandCanonical(
    gadb='gcga.db', substrate='substrate.vasp', popSize=input.popSize,
    convergeCrit=input.popSize * 10, subsPot=input.subsPot,
    chemPotDict=input.chemPotDict, zLim=input.zLim)

if not os.path.isfile('gmid'):          # fresh start; comment out to restart
    pop.initializeDB()
    pop.natural_selection()

# resume numbering after existing kid directories
kidnum = max([int(d) for d in os.listdir('.') if d.isdigit()] + [0])
running = {}

def spawn(k):
    d = '%06d' % k
    os.makedirs(d, exist_ok=True)
    log = open(os.path.join(d, 'worker.log'), 'w')
    p = subprocess.Popen([sys.executable, '-u', '../ga-worker.py'], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    running[k] = (p, log)
    print('Job %i\t@%s' % (k, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')), flush=True)

while True:
    for k in list(running):
        p, log = running[k]
        if p.poll() is not None:
            log.close()
            del running[k]
    stop = 'STOP' in os.listdir('.')
    done_enough = kidnum >= totConf or (pop.is_converged2() and kidnum >= minConf)
    if (stop or done_enough) and not running:
        break
    while not stop and not done_enough and len(running) < nworker:
        kidnum += 1
        spawn(kidnum)
        time.sleep(5)
    time.sleep(60)

print('CONVERGED!' if pop.is_converged() else 'TERMINATED!', flush=True)
