# Pyrosetta Basic

Pyrosetta Basic中文教程，本教程由浅入深，讲解Rosetta的基本原理以及在PyRosetta中的应用实例。

@文档贡献者：

1. 吴炜坤 @晶泰人工智能研发中心
2. 黄健 @晶泰人工智能研发中心
3. 张博文 @晶泰人工智能研发中心
4. 槐喆 @晶泰人工智能研发中心 实习生

@校对：

1. 王天元 @晶泰人工智能研发中心
2. 郭宁 @晶泰人工智能研发中心
3. 张晨虹 @晶泰人工智能研发中心



@外援支持:

1. 刘源 博士后 北京大学王初课题组



## 大纲内容:

零、安装与入门介绍

0.0 [Installation](https://github.com/guyujun/chinese-pyrosetta/blob/master/0_0_Installation.ipynb)

0.1 [Python_Basic](https://github.com/guyujun/chinese-pyrosetta/blob/master/0_1_Python_Basic.ipynb)

0.2 Utils



一、Pose与Structure IO: 负责介绍PyRosetta对结构文件的处理，以及Pose对象的重要作用

负责人:@吴炜坤  进度: 100% 

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/02.00-Introduction-to-PyRosetta.ipynb

- 1.0 [Pose Object Abstract](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_0_Pose_Abstract.ipynb)

- 1.1 [Pose IO](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_1_Pose_IO.ipynb)

- 1.2 [PymolMover](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_2_PyMover_PyRosetta.ipynb)

- 1.3 [Pose & PDBinfo](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_3_Pose_PDBinfo.ipynb)

- 1.4 [Atom & Residue](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_4_Atom_Residue.ipynb)

- 1.5 [Conformation & Protein Geometry](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_5_Conformation_Geometry.ipynb)

- 1.6 [Pose Operation](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/1_6_Pose_Operating.ipynb)



二、Energy Function与Constraint: 介绍Rosetta的能量函数与物理约束

负责人: @黄健 进度: 60% 

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/03.00-Rosetta-Energy-Score-Functions.ipynb

Constraint的API总结: https://zhuanlan.zhihu.com/p/58897635

- 2.0 [Atom Model](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/2_Energy/2_1_Atom_Model.ipynb)

- 2.1 [Energy Terms and Score Function](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/2_Energy/2_2_Energy_Function.ipynb)

- 2.2 [Constraints](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/2_Energy/2_3_Constraint.ipynb)



三、Kinematics与MoveMap: 介绍Rosetta的自由度控制

负责人:@张博文 进度: 10% 

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/04.00-Introduction-to-Folding.ipynb

Foldtree的概念: https://zhuanlan.zhihu.com/p/59863638

- 3.0 FoldTree
- 3.1 Docking Tree & Jumps
- 3.2 MoveMap



四、Monte Carlo与Folding: 介绍Rosetta中的Foldtree与蒙特卡洛算法

负责人:@吴炜坤  进度: 100% 

- 4.0 [Metropolis & Simulated annealing](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/4_MCMC/Metropolis_Monte_Carlo.ipynb)

- 4.1 [Movers & MC object ](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/4_MCMC/Movers_MC_object.ipynb)



五、Residue Selector: 介绍残基选择器

负责人:@槐喆  进度: 50% 。校对:@吴炜坤 

中文总结：https://zhuanlan.zhihu.com/p/58348980

官网资料: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/ResidueSelectors/ResidueSelectors

进度: 10%, 需要差缺补漏

[residue selector preview](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/ResidueSelectors.ipynb) 

- 5.0 Residue Selector的逻辑

- 5.1 Residue Selector的API查询

- 5.2 Residue Selector在pymol中的显示



六、Packer与TaskOperation: 介绍Packer与氨基酸侧链自由度控制

负责人:@吴炜坤 进度: 50% 

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.00-Introduction-to-Packing-and-Design.ipynb

TaskOperation: [Pack和Design用法.pdf](https://xtalpi.feishu.cn/file/boxcnb4h8Gl8QNLmRgJikidqN9c) 

- 6.1 Rotamers & RotamerLib
- 6.2 Packer & Design
- 6.3 TaskOperation与TaskFactory（Rotamer自由度控制）
- 6.4 NCAA(调色板)



七、SimpleMetric: 新一代特征计算和记录工具

负责人:@槐喆 @黄健 进度: 30% 

SimpleMetric的API总结 https://zhuanlan.zhihu.com/p/58383955

- 7.0 SimpleMetric



八、Filters: 过滤器

负责人: @黄健 @张博文 进度: 0% 

官网资料: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/Filters-RosettaScripts

- 8.0 Filters



九、xmlObject & RosettaScript: xmlObject如何解决Rosetta历史遗留问题

负责人:@黄健 进度: 0% 

xmlObject的API总结: https://zhuanlan.zhihu.com/p/58381573

官网资料: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts

- 9.0 RosettaScript
- 9.1 XmlObject
- 9.2 自定义Mover
