# PyRosetta基础中文教程

PyRosetta的官方Notebook在很久之前就正式发布了，这次由晶泰科技团队带来的PyRosetta基础中文教程具有更好的易读性、更全面的用法介绍以及更多的实例展示。在这个教程中，读者通过动手实操一步步了解Rosetta的底层逻辑与建模思想，可以从零了解如何使用API组件去一步步搭建设计蛋白质、多肽的计算流程，同时教程中从零翻译了大量的PyRosetta API接口以备不时之需。

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



## 大纲内容

### 零、PyRosetta与Pymol服务器的安装配置:

负责人:@吴炜坤  

-  0.0 [安装与配置](https://nbviewer.jupyter.org/github/guyujun/chinese-pyrosetta/blob/master/1_PoseIO/0_0_Installation.ipynb)



### 一、Pose与Structure IO: 

负责介绍PyRosetta对结构文件的处理，以及Pose对象的重要作用 

负责人:@吴炜坤  

- 1.0 [Pose IO介绍](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_1_Pose_IO.ipynb)
- 1.1 [PymolMover与可视化](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_2_PyMover_PyRosetta.ipynb)
- 1.2 [Pose与PDBinfo的交互](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_3_Pose_PDBinfo.ipynb)
- 1.3 [Atom与Residue的层级](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_4_Atom_Residue.ipynb)
- 1.4 [Conformation层级](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_5_Conformation_Geometry.ipynb)
- 1.5 [Pose的操作变换一览](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/1_PoseIO/1_6_Pose_Operating.ipynb)



### 二、Energy Function与Constraint

介绍Rosetta的能量函数与物理约束

负责人: @黄健 

- 2.0 [Atom与Structure的表示方式](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_1_Atom_Model.ipynb)
- 2.1 [Energy Terms and Score Function介绍](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_2_Energy_Function.ipynb)
- 2.2 [Constraints的介绍与用法](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_3_Constraint.ipynb)
- 2.3 [Contsraint的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/2_Energy/2_4_Contsraint_API.ipynb)



### 三、Kinematics与Trees

介绍Rosetta的骨架自由度控制

负责人:@张博文 

- 3.0 [FoldTree的连续性](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/3_Kinematics/3_0_FoldTree.ipynb)
- 3.1 [FoldTree的不连续性](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/3_Kinematics/3_1_Jump_Cutpoint.ipynb)



### 四、Monte Carlo

介绍Rosetta中的蒙特卡洛算法【核心】

负责人:@吴炜坤  

- 4.0 [Metropolis Monte Carlo采样](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/4_0_Metropolis_Monte_Carlo.ipynb)
- 4.1 [Movers & MC object运作机理](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/4_1_Movers_MC_object.ipynb)
- 4.2 [Fragment如何加速Folding采样](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/4_2_Fragment_Folding.ipynb)
- 4.3 [Extend_Reading(选读)](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/4_MCMC/Extended_Reading_Metropolis_Monte_Carlo.ipynb)



### 五、Residue Selector

介绍残基选择器，自定义选择范围。

负责人:@槐喆

- 5.0 [Residue Selector简介](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/5_Residue_Selector/5_0_ResidueSelectors_Logic.ipynb)
- 5.1 [Residue Selector的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/5_Residue_Selector/5_1_ResidueSelector_ApiSearch.ipynb)



### 六、Packer与TaskOperation

介绍Packer与氨基酸侧链自由度控制，如何使用PyRosetta进行设计

负责人:@吴炜坤 

- 6.0 [Rotamers & Packer简介](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_0_Rotamer_Packer.ipynb)
- 6.1 [TaskOperation & PackTask任务机制](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_1_PackTask_TaskOP.ipynb)
- 6.2 [Resfile System与Rotamer控制](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_2_Resfile_System.ipynb)
- 6.3 [TaskOperation的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/6_Packer_TaskOperation/6_3_TaskOperation_API.ipynb)



### 七、SimpleMetric

新一代特征计算和记录工具

负责人:@槐喆 

- 7.0 [SimpleMetric逻辑](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/7_Simple_Metrics/7_0_Simple_Metrics_Logic.ipynb)
- 7.1 [SimpleMetric的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/7_Simple_Metrics/7_1_Simple_Metrics_ApiSearch.ipynb)



### 八、Filters

过滤器也是计算器。

负责人: @吴炜坤@黄健 

- 8.0 [Filters的逻辑](https://xtalpi.feishu.cn/docs/[ 8_1_Filter_Introduction.ipynb](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/8_Filter/8_1_Filter_Introduction.ipynb))
- 8.1 [Filters的API查询](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/8_Filter/8_2_RosettaFilter_API.ipynb)



### 九、RosettaScript & XmlObject

xmlObject如何解决Rosetta历史遗留问题

负责人:@黄健@吴炜坤 

- 9.0 [RosettaScript基础](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/9_xmlObject_RosettaScript/9_1_RS_basis.ipynb)
- 9.1  [RosettaScript进阶](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/9_xmlObject_RosettaScript/9_2_RS_advanced.ipynb)
- 9.2 [XmlObject的用法](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/9_XmlObject/9_3_XmlObject.ipynb)



### 十、进阶分析

Silent与rstoolbox的完美结合。

负责人:@吴炜坤 

- 10.0 [rstoolbox对Silent文件解析](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/10_Analysis/10_0_rstoolbox.ipynb)
- 10.1 [rstoolbox更多有用的API](https://nbviewer.jupyter.org/github/guyujun/pyrosetta-basic/blob/master/10_Analysis/10_1_more_api.ipynb)



### 参考资料:

PyRosetta中文基础计划的地址: https://github.com/guyujun/chinese-pyrosetta

官方PyRosetta Notebook开源地址: https://github.com/RosettaCommons/PyRosetta.notebooks

官方PyRosetta API查询: https://graylab.jhpytu.edu/PyRosetta.documentation/search.html?q=cdr

Rosetta中文社区知乎: https://www.zhihu.com/column/rosettastudy

