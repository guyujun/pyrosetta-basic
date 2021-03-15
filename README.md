# Pyrosetta Core

PyRosetta Core中文教程，本教程由浅入深，讲解Rosetta的基本原理以及在PyRosetta中的应用实例。

贡献者：

1. 吴炜坤 @晶泰人工智能研发中心
2. 翟珂 @晶泰人工智能研发中心
3. 胡志运 @晶泰人工智能研发中心
4. 张博文 @晶泰人工智能研发中心



## 3 中文章节分配

成员可自行领取@all

#### 大纲内容:

零、安装与入门介绍

0.0 [Installation](https://github.com/guyujun/chinese-pyrosetta/blob/master/0_0_Installation.ipynb)

0.1 [Python_Basic](https://github.com/guyujun/chinese-pyrosetta/blob/master/0_1_Python_Basic.ipynb)

0.2 Utils



一、Pose与Structure IO: 负责介绍PyRosetta对结构文件的处理，以及Pose对象的重要作用@吴炜坤 

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/02.00-Introduction-to-PyRosetta.ipynb

进度: 100% @吴炜坤 

1 [Pose Object Abstract](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_0_Pose_Abstract.ipynb)

1.1 [Pose IO](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_1_Pose_IO.ipynb)

1.2 [PymolMover](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_2_PyMover_PyRosetta.ipynb)

1.3 [Pose & PDBinfo](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_3_Pose_PDBinfo.ipynb)

1.4 [Atom & Residue](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_4_Atom_Residue.ipynb)

1.5 [Conformation & Protein Geometry](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_5_Conformation_Geometry.ipynb)

1.6 [Pose Operation](https://github.com/guyujun/chinese-pyrosetta/blob/master/1_6_Pose_Operating.ipynb)



二、Energy Function与Constraint: 介绍Rosetta的能量函数与物理约束

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/03.00-Rosetta-Energy-Score-Functions.ipynb

Constraint的API总结: https://zhuanlan.zhihu.com/p/58897635

2.0 Atom Model

2.1 Energy Terms and Score Function

2.2 Constraints



三、Kinematics与MoveMap: 介绍Rosetta的自由度控制

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/04.00-Introduction-to-Folding.ipynb

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/05.00-Structure-Refinement.ipynb

3.0 FoldTree

3.1 MoveMap



四、Monte Carlo与Folding: 介绍Rosetta中的Foldtree与蒙特卡洛算法@吴炜坤 

Foldtree的概念: https://zhuanlan.zhihu.com/p/59863638

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/04.00-Introduction-to-Folding.ipynb

4.0 Monte Carlo object

4.1 Fragment & Folding

4.2 Refinement



五、Residue Selector: 介绍残基选择器

中文总结：https://zhuanlan.zhihu.com/p/58348980

官网资料: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/ResidueSelectors/ResidueSelectors

进度: 10%, 需要差缺补漏

5.0 [residue selector preview](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/ResidueSelectors.ipynb)



六、Packer与TaskOperation: 介绍Packer与氨基酸侧链自由度控制

相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.00-Introduction-to-Packing-and-Design.ipynb

6.0 Rotamers & Rotamer Lib

6.1 Packer & Repack

6.2 TaskOperation & TaskFactory

6.3 PackerPalettes & D-amino acids



七、SimpleMetric: 新一代的打分系统

SimpleMetric的API总结 https://zhuanlan.zhihu.com/p/58383955

7.0 SimpleMetric



八、Filters: 过滤器，大过滤器！

官网资料: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/Filters-RosettaScripts

8.0 Filters



九、xmlObject & RosettaScript: xmlObject如何解决Rosetta历史遗留问题

xmlObject的API总结: https://zhuanlan.zhihu.com/p/58381573

官网资料: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts

9.0 RosettaScript

9.1 xmlObject



专题:

1. Protein Docking 
2. Loop Modelling
3. Antibody Modelling
4. Denovo Design
5. trRosetta
6. 对称性处理