# Run once (login node is fine, database only) before submitting Ibex workers: same start-up as ga-bundle.py.
import os
import input
from gocia.ga.popGrandCanon import PopulationGrandCanonical

pop = PopulationGrandCanonical(
    gadb='gcga.db', substrate='substrate.vasp', popSize=input.popSize,
    convergeCrit=input.popSize * 10, subsPot=input.subsPot,
    chemPotDict=input.chemPotDict, zLim=input.zLim)
if not os.path.isfile('gmid'):
    pop.initializeDB()
    pop.natural_selection()
    print('initialized: gmid written')
else:
    print('gmid exists: nothing to do (restart)')
