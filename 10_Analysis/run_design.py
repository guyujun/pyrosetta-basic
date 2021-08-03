# 举例使用FastDesign快速设计一些序列和结构:
from pyrosetta import pose_from_pdb, init, create_score_function
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.io import poses_to_silent

# init
init('')

# load pose
starting_pose = pose_from_pdb('./data/EHEE_rd4_0976.pdb')
ref2015 = create_score_function('ref2015')
design_tf = TaskFactory()

# setup FastDesign
fastdesign = FastDesign(ref2015, 1)
fastdesign.set_default_movemap() #使用默认的Movemap()
fastdesign.set_task_factory(design_tf)


# design for 10 times:
for i in range(10):
    design_pose = Pose()
    design_pose.assign(starting_pose)  # assign pose
    fastdesign.apply(design_pose)  ## apply design
    # output to silent file;
    poses_to_silent(design_pose, './data/design_result.silent')
