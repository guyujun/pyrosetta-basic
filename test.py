#!/usr/bin/env python

from pyrosetta import pose_from_pdb,init,create_score_function
init()
pdb_pose = pose_from_pdb('./data/pose_demo.pdb')

# 发送pose到pymol中
from pyrosetta.rosetta.protocols.moves import PyMOLMover
pymover = PyMOLMover()
pymover.update_energy(True)
pymover.apply(pdb_pose)

scoring = create_score_function('ref2015')
scoring(pose)
