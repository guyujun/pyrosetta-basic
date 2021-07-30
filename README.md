# PyRosetta Basic

PyRosetta Basic中文教程，讲解Rosetta的基本原理以及在PyRosetta中的应用实例。

@文档贡献者：

1. 吴炜坤 @晶泰人工智能研发中心
2. 黄健 @晶泰人工智能研发中心
3. 张博文 @晶泰人工智能研发中心
4. 槐喆 @晶泰人工智能研发中心

@校对：

1. 王天元 @晶泰人工智能研发中心
2. 郭宁 @晶泰人工智能研发中心
3. 张晨虹 @晶泰人工智能研发中心



@外援支持:

1. 刘源 博士后 北京大学王初课题组



## 大纲内容:

### 零、安装与入门介绍

0.0 [Installation](https://github.com/guyujun/chinese-pyrosetta/blob/master/0_0_Installation.ipynb)

0.1 [Python_Basic](https://github.com/guyujun/chinese-pyrosetta/blob/master/0_1_Python_Basic.ipynb)

0.2 Utils



### 一、Pose与Structure IO

> 介绍PyRosetta对结构文件的处理，以及Pose对象的重要作用 

负责人:@吴炜坤  进度: 100% 

- 1.0 [Pose Object Abstract](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_0_Pose_Abstract.ipynb)

- 1.1 [Pose IO](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_1_Pose_IO.ipynb)

- 1.2 [PymolMover](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_2_PyMover_PyRosetta.ipynb)

- 1.3 [Pose & PDBinfo](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_3_Pose_PDBinfo.ipynb)

- 1.4 [Atom & Residue](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_4_Atom_Residue.ipynb)

- 1.5 [Conformation & Protein Geometry](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_5_Conformation_Geometry.ipynb)

- 1.6 [Pose Operation](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_6_Pose_Operating.ipynb)



### 二、Energy Function与Constraint: 介绍Rosetta的能量函数与物理约束

负责人: @黄健 校对: @吴炜坤 进度: 100% 

> 相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/03.00-Rosetta-Energy-Score-Functions.ipynb

> Constraint的API总结: https://zhuanlan.zhihu.com/p/58897635

- 2.0 [Atom Model](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_1_Atom_Model.ipynb)
- 2.1 [Energy Terms and Score Function](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_2_Energy_Function.ipynb)
- 2.2 [Constraints](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_3_Constraint.ipynb)
- 2.3 [Constraint_API](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_4_Contsraint_API.ipynb)



### 三、Kinematics与Trees: 介绍Rosetta的骨架自由度控制

负责人:@张博文 校对: @吴炜坤 进度: 100% 

> 相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/04.00-Introduction-to-Folding.ipynb

> Foldtree的概念: https://zhuanlan.zhihu.com/p/59863638

- 3.0 [FoldTree与顺序性](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/3_Kinematics/3_0_FoldTree.ipynb)

- 3.1 [Jump_Cutpoint](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/3_Kinematics/3_1_Jump_Cutpoint.ipynb)




### 四、Monte Carlo: 介绍Rosetta中的蒙特卡洛算法【核心】

负责人:@吴炜坤  进度: 100% 

> 相关的官方章节:https://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/4.02-Low-Res-Scoring-and-Fragments.ipynb

> 相关的官方章节: http://nbviewer.jupyter.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/05.00-Structure-Refinement.ipynb

- 4.0 [Metropolis & Simulated annealing](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/4_0_Metropolis_Monte_Carlo.ipynb)
- 4.1 [Movers & MC object ](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/4_1_Movers_MC_object.ipynb)
- 4.2 [Fragment_Folding](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/4_2_Fragment_Folding.ipynb)
- 4.3 [Extend_Reading](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/Extended_Reading_Metropolis_Monte_Carlo.ipynb)



### 五、Residue Selector: 介绍残基选择器

> 介绍残基选择器，自定义选择范围。

负责人: @槐喆  校对: @吴炜坤 进度: 100% 

- 5.0 [Residue Selector的逻辑](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/5_Residue_Selector/5_0_ResidueSelectors_Logic.ipynb)
- 5.1 [Residue Selector的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/5_Residue_Selector/5_1_ResidueSelector_ApiSearch.ipynb)



### 六、Packer与TaskOperation: 介绍Packer与氨基酸侧链自由度控制

> 介绍Packer与氨基酸侧链自由度控制，如何使用PyRosetta进行设计

负责人:@吴炜坤 进度: 100% 

- 6.0 [Rotamers & Packer](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_0_Rotamer_Packer.ipynb)
- 6.1 [TaskOperation & PackTask](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_1_PackTask_TaskOP.ipynb)
- 6.2 [Resfile System](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_2_Resfile System.ipynb)
- 6.3 [TaskOperation API](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_3_TaskOperation_API.ipynb)



### 七、SimpleMetric: 新一代特征计算和记录工具

> 新一代特征计算和记录工具

负责人:@槐喆  校对: @吴炜坤 

- 7.0 [SimpleMetric逻辑](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/7_Simple_Metrics/7_0_Simple_Metrics_Logic.ipynb)
- 7.1 [SimpleMetric的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/7_Simple_Metrics/7_1_Simple_Metrics_ApiSearch.ipynb)



### 八、Filters: 过滤器

> 过滤器，大过滤器！

负责人: @黄健 进度: 50% 校对: @吴炜坤 

- 8.0 [Filters的逻辑]([ 8_1_Filter_Introduction.ipynb](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/8_Filter/8_1_Filter_Introduction.ipynb))

- 8.1 [Filters的API](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/8_Filter/8_2_Filter_api.ipynb)



### 九、xmlObject & RosettaScript: xmlObject如何解决Rosetta历史遗留问题

负责人:@黄健 进度: 90% 校对: @吴炜坤 

> xmlObject如何解决Rosetta C++ to Python历史接口

xmlObject的API总结: https://zhuanlan.zhihu.com/p/58381573

- 9.0 [RosettaScript基础](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/9_xmlObject_RosettaScript/9_1_RS_basis.ipynb)
- 9.1  [RosettaScript进阶](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/9_xmlObject_RosettaScript/9_2_RS_advanced.ipynb)



## 参考资料:

中文开源计划的地址: https://github.com/guyujun/chinese-pyrosetta

PyRosetta Notebook开源地址: https://github.com/RosettaCommons/PyRosetta.notebooks

PyRosetta API查询: https://graylab.jhpytu.edu/PyRosetta.documentation/search.html?q=cdr

Rosetta中文知乎: https://www.zhihu.com/column/rosettastudy