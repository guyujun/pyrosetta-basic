#!/usr/bin/env python
import pyrosetta
import random
import math
import seaborn as sns
from pyrosetta import init, pose_from_sequence, create_score_function
from pyrosetta.rosetta.core.pose import Pose
import pandas as pd
import matplotlib as plt
plt.style.use('ggplot')
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300


def sample(pose):
    delta_phi = random.uniform(-5, 5)
    delta_psi = random.uniform(-5, 5)
    cur_psi = pose.psi(1)
    cur_phi = pose.phi(2)

    new_phi = cur_phi+delta_phi
    new_psi = cur_psi+delta_psi

    if new_phi > 180:
        new_phi = -180 + (new_phi-180)
    elif new_phi < -180:
        new_phi = 180 - (new_phi+180)

    if new_psi > 180:
        new_psi = -180 + (new_psi-180)
    elif new_psi < -180:
        new_psi = 180 - (new_psi+180)

    pose.set_phi(2, new_phi)
    pose.set_psi(1, new_psi)


def run_mcmc(pose):
    #
    old_score = score(pose)
    old_pose = Pose()
    old_pose.assign(pose)

    # sample;
    sample(pose)
    new_score = score(pose)

    # accept?
    P = math.exp(-1*(new_score-old_score)/1.5)
    print(P)

    if P > 1:
        return pose
    else:
        p = random.uniform(0, 1.0)
        if p < P:
            return pose
        else:
            return old_pose


if __name__ == '__main__':
    init()
    pose = pose_from_sequence('AA')
    score = create_score_function('ref2015')
    pose.set_phi(2, random.uniform(-180, 180))
    pose.set_psi(1, random.uniform(-180, 180))
    print(pose.psi(1), pose.phi(2))
    pose.dump_pdb('mc.pdb')
    df = pd.DataFrame()
    energy = score(pose)
    df.loc[0, 'psi'] = pose.psi(1)
    df.loc[0, 'phi'] = pose.phi(2)
    # df.loc[0, 'energy'] = energy

    for i in range(1, 10000):
        pose = run_mcmc(pose)
        cur_psi = pose.psi(1)
        cur_phi = pose.phi(2)
        energy = score(pose)
        df.loc[i, 'psi'] = cur_psi
        df.loc[i, 'phi'] = cur_phi
        # df.loc[i, 'energy'] = energy

    # df = df.sort_values('energy').head(250)
    df
    ax = sns.jointplot('psi', 'phi', data=df, xlim=(-180,180), ylim=(-180, 180), kind="hex", color='b')
    ax.savefig('mcmc.png')
