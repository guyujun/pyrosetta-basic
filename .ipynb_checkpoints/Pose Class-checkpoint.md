# 第一讲: Pose与Structure IO

蛋白质是大型生物分子，它由一个或多个由α-氨基酸残基组成的长链条组成。α-氨基酸分子呈线性排列，相邻α-氨基酸残基的羧基和氨基通过肽键连接在一起。
蛋白质的分子结构可划分为四级，以描述其不同的方面：</br>
蛋白质一级结构：组成蛋白质多肽链的线性氨基酸序列。</br>
蛋白质二级结构：依靠不同氨基酸之间的C=O和N-H基团间的氢键形成的稳定结构，主要为α螺旋和β折叠。</br>
蛋白质三级结构：通过多个二级结构元素在三维空间的排列所形成的一个蛋白质分子的三维结构。</br>
蛋白质四级结构：用于描述由不同多肽链（亚基）间相互作用形成具有功能的蛋白质复合物分子。</br>
![title](./img/pose.png)

### 1.1 Pose的组织构架
因此如果要在计算机中建立一个蛋白质的结构模型，就清楚地描述每一个原子的信息。在Rosetta中，Pose是管理蛋白质信息的中心，可以描述蛋白质一到四级结构所有的信息。而且这些信息是分层管理的比如:

Conformation: 负责管理原子类型(AtomType)、氨基酸类型(ResidueType)、氨基酸的原子坐标(xyz)、氨基酸连接方式的定义(FoldTree/AtomTree)</br>
Energy: 负责管理氨基酸直接的能量计算所需的信息(EnergyGraph/energies)</br>
ConstraintSet: 负责管理原子间的约束信息(constraints)</br>
DataCache: 负责管理用户自定义的信息</br>
分层式管理使得Pose的信息修改和更新变得容易。

除此以外，还有些外部对象如PDBinfo也负责转换和储存Pose与PDB之间的信息变换。

以下是一个Pose中的架构的示意图：
![title](./img/PoseObject.png)

### 1.2 Pose的生成与输出
Rosetta兼容最常规的几种两种记录结构格式：PDB和Silent文件：
- PDB文件可以从https://www.rcsb.org/ 数据库中获取；
- Silent文件为Rosetta开发的pose压缩文件（其功能也是储存结构等信息，但其体积比PDB小了10倍之多，非常适合在超算中心进行的数据文件的传输）

Rosetta为PDB的结构信息读入提供了非常丰富的接口，此处我们介绍主要4种结构信息读入相关的函数:
- pose_from_pdb：从pdb文件读入并生成pose
- pose_from_sequence：从一级序列信息生成pose
- poses_from_silent：从silent文件读入并生成pose
- pose_from_rcsb：从rcsb数据库远程获取PDB code数据并读入和生成pose

一般而言，Rosetta输入的结构信息大多来源于PDB结构文件，经过处理后生成对应的Pose。以下将逐一介绍如何使用这些外部文件来生成Pose对象。


```python
# 初始化（必须执行的步骤）:
from pyrosetta import init
init()
```

    PyRosetta-4 2020 [Rosetta PyRosetta4.conda.mac.python37.Release 2020.02+release.22ef835b4a2647af94fcd6421a85720f07eddf12 2020-01-05T17:31:56] retrieved from: http://www.pyrosetta.org
    (C) Copyright Rosetta Commons Member Institutions. Created in JHU by Sergey Lyskov and PyRosetta Team.
    [0mcore.init: {0} [0mChecking for fconfig files in pwd and ./rosetta/flags
    [0mcore.init: {0} [0mRosetta version: PyRosetta4.conda.mac.python37.Release r242 2020.02+release.22ef835b4a2 22ef835b4a2647af94fcd6421a85720f07eddf12 http://www.pyrosetta.org 2020-01-05T17:31:56
    [0mcore.init: {0} [0mcommand: PyRosetta -ex1 -ex2aro -database /opt/miniconda3/lib/python3.7/site-packages/pyrosetta/database
    [0mbasic.random.init_random_generator: {0} [0m'RNG device' seed mode, using '/dev/urandom', seed=-1854535783 seed_offset=0 real_seed=-1854535783 thread_index=0
    [0mbasic.random.init_random_generator: {0} [0mRandomGenerator:init: Normal mode, seed=-1854535783 RG_type=mt19937


#### 1.2.1 从PDB文件生成Pose


```python
# 此处以denovo的小蛋白为例作为读入
from pyrosetta import pose_from_pdb
pose = pose_from_pdb('./data/pose_demo.pdb')
print(pose)
```

    [0mcore.chemical.GlobalResidueTypeSet: {0} [0mFinished initializing fa_standard residue type set.  Created 980 residue types
    [0mcore.chemical.GlobalResidueTypeSet: {0} [0mTotal time to initialize 0.729505 seconds.
    [0mcore.import_pose.import_pose: {0} [0mFile './data/pose_demo.pdb' automatically determined to be of type PDB
    [0mcore.conformation.Conformation: {0} [0mFound disulfide between residues 2 11
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 2 CYS
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 11 CYS
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 2 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 11 CYD
    [0mcore.conformation.Conformation: {0} [0mFound disulfide between residues 5 25
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 5 CYS
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 25 CYS
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 5 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 25 CYD
    PDB file name: ./data/pose_demo.pdb
    Total residues: 26
    Sequence: HCFHCRNIRFCSEDEEELRRAREECK
    Fold tree:
    FOLD_TREE  EDGE 1 26 -1 


##### 结果解读:
表明，PDB文件自动读入成功，并发现2对二硫键，位于2、11、5、25位。总氨基酸数目是26、氨基酸的序列等信息都可以直接打印出来。

#### 1.2.2 从序列文件生成Pose


```python
# 从一级序列生成full-atom Pose
from pyrosetta import pose_from_sequence
seq_pose = pose_from_sequence('AAAAA', "fa_standard")
seq_pose
```




    <pyrosetta.rosetta.core.pose.Pose at 0x112954f70>



##### 结果解读:
<pyrosetta.rosetta.core.pose.Pose at ***********> 表明新的pose已经生成并储存在内存当中，该序列是线性的，可以加上后缀指定生成的是centriod("centroid")还是full-atom("fa_standard")模型的pose

#### 1.2.3 从PDB代号生成Pose


```python
# 我们还可以从rscb的PDB数据库中直接下载并读入生成Pose
from pyrosetta.toolbox.rcsb import pose_from_rcsb
rcsb_pose = pose_from_rcsb('4R80')
rcsb_pose
```

    [0mcore.import_pose.import_pose: {0} [0mFile '4R80.clean.pdb' automatically determined to be of type PDB
    [0mcore.conformation.Conformation: {0} [0m[1m[ WARNING ][0m missing heavyatom:  OXT on residue SER:CtermProteinFull 76
    [0mcore.conformation.Conformation: {0} [0m[1m[ WARNING ][0m missing heavyatom:  OXT on residue SER:CtermProteinFull 152





    <pyrosetta.rosetta.core.pose.Pose at 0x1128a03f0>



##### 结果解读:
File '4R80.clean.pdb' automatically determined to be of type PDB
运行结束时间根据网络情况而定，下载的代号为4R80的PDB结构被自动清洗和读入。
<pyrosetta.rosetta.core.pose.Pose at ***********> 表明新的pose已经生成并储存在内存当中

#### 1.2.3 从Silent文件读取生成Pose


```python
# 读入Silent文件
from pyrosetta.io import poses_from_silent
poses = poses_from_silent('./data/demo.silent')
print(poses)
for pose in poses:
    print(pose)
```

    <generator object poses_from_silent at 0x14cf29950>
    [0mcore.io.silent.SilentFileData: {0} [0mReading all structures from ./data/demo.silent
    [0mcore.io.silent.SilentFileData: {0} [0mFinished reading 2 structures from ./data/demo.silent
    [0mcore.conformation.Conformation: {0} [0mFound disulfide between residues 2 11
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 2 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 11 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 2 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 11 CYD
    [0mcore.conformation.Conformation: {0} [0mFound disulfide between residues 5 25
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 5 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 25 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 5 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 25 CYD
    PDB file name: /Users/kunkun/Desktop/PyRosetta培训资料/data/pose_demo.pdb
    Total residues: 26
    Sequence: HCFHCRNIRFCSEDEEELRRAREECK
    Fold tree:
    FOLD_TREE  EDGE 1 26 -1 
    [0mcore.conformation.Conformation: {0} [0mFound disulfide between residues 2 11
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 2 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 11 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 2 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 11 CYD
    [0mcore.conformation.Conformation: {0} [0mFound disulfide between residues 5 25
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 5 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 25 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 5 CYD
    [0mcore.conformation.Conformation: {0} [0mcurrent variant for 25 CYD
    PDB file name: /Users/kunkun/Desktop/PyRosetta培训资料/data/pose_demo1.pdb
    Total residues: 26
    Sequence: HCFHCRNIRFCSEDEEELRRAREECK
    Fold tree:
    FOLD_TREE  EDGE 1 26 -1 


##### 结果解读:
poses_from_silent返回的是一个迭代生成器对象, 其中包含了所有的pose。
<generator object poses_from_silent at ***> 表明新的poses object已经生成并储存在内存当中.
只需要通过for循环即可遍历获取pose并操作。此处可见，返回了silent文件中的两个Pose

#### 1.2.4 从Pose生成输出Silent/PDB文件
从已经生成的pose输出PDB结构或Silent文件，仅需要调用Pose类中dump_pdb的方法函数即可或poses_to_silent函数即可。


```python
# 输出pdb
seq_pose.dump_pdb('./data/AAAAA.pdb')
```




    True




```python
# 输出5个重复的pose到silent文件中
from pyrosetta.io import poses_to_silent
for i in range(5):
    poses_to_silent(seq_pose, './data/multi_AAAAA.silent')
```

### 1.3 Pose中的PDB信息

Pose是通常是从PDB文件中衍生出来的，通常除了原子的坐标信息以外，PDB文件中含包含了许多额外的信息，而这些信息是储存在PDBinfo中。比如温度因子数据(bfactor)、晶体解析数据(crystinfo)、原子的占用率(occupancy)、Pose编号与PDB编号的转换以及Pose的序列信息等。**如果Pose中的氨基酸发生了插入和删除，那么这部分的信息就需要重新从Pose中进行更新**。以下列举PDBinfo的一些重要功能:

##### **1.3.1 Rosetta编号与PDB编号**
在PDB编号中，氨基酸的编号是依赖于其所在的链，如1A，2A...120A，1B，2B...130B等。
而在Pose类中，氨基酸的编号是忽略链的分隔，按照链的顺序，**均从1开始递增**，因此在Pose中的氨基酸和PDB编号的对应是棘手的，但是我们可以通过PDB_info这个类中的pdb2pose和pose2pdb来解决转换问题


```python
# 获取PDB号为24A的氨基酸残基所在的Pose残基编号
pose = pose_from_pdb('./data/4R80.pdb')
pose_pdbinfo = pose.pdb_info()
pose_number = pose_pdbinfo.pdb2pose('A', 24)
print(pose_number)
```

    [0mcore.import_pose.import_pose: {0} [0mFile './data/4R80.pdb' automatically determined to be of type PDB
    [0mcore.conformation.Conformation: {0} [0m[1m[ WARNING ][0m missing heavyatom:  OXT on residue SER:CtermProteinFull 76
    [0mcore.conformation.Conformation: {0} [0m[1m[ WARNING ][0m missing heavyatom:  OXT on residue SER:CtermProteinFull 152
    24



```python
# 获取Pose残基编号为24的氨基酸残基所在的PDB号
pdb_number = pose_pdbinfo.pose2pdb(24)
print(pdb_number)
```

    24 A 



```python
# 获取pose number为24的氨基酸残基所在的PDB链号
pose_pdbinfo.chain(24)
```




    'A'



##### **1.3.2 PDB中晶体解析信息的提取**
除了基本的编号信息以外，一些晶体相关的信息也可以轻松进行提取:


```python
# 提取bfactor信息
pose.pdb_info().bfactor(1,1)  # 返回温度因子信息
```




    49.13




```python
# 获取PDB的晶体信息
crystinfo = pose.pdb_info().crystinfo()
print(crystinfo)
```

    <CrystInfo>{0,0,0,90,90,90 : P 1}



```python
# 获取原子的occupancy
pose.pdb_info().occupancy(1, 1)
```




    1.0



### 1.4 Pose中的结构信息
PDB的信息可以通过PDBinfo获取，除此以外，Pose中还有大量的结构信息，我们可以轻松通过各类函数来获取更多的信息，以下将逐一地列举。

#### 1.4.1 一级与二级结构信息的提取


```python
# Pose中的基本信息，如残基数量，序列，FoldTree(后续讲述)
print(pose)
```

    PDB file name: ./data/4R80.pdb
    Total residues: 152
    Sequence: PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    Fold tree:
    FOLD_TREE  EDGE 1 76 -1  EDGE 1 77 1  EDGE 77 152 -1 


**Tips: 除此以外，还有一些具体的方法可以获取更多的细节:**


```python
# 返回pose中链的数目
pose.num_chains()
```




    2




```python
# 获取Pose的序列信息的方法
print(f'所有的氨基酸:{pose.sequence()}\n')
print(f'前5个氨基酸:{pose.sequence(1,5)}\n')
print(f'1号链的所有氨基酸:{pose.chain_sequence(1)}\n')
print(f'2号链的所有氨基酸:{pose.chain_sequence(2)}\n')
```

    所有的氨基酸:PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    
    前5个氨基酸:PSEEE
    
    1号链的所有氨基酸:PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    
    2号链的所有氨基酸:PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    



```python
# 获取Pose的氨基酸总数量方法
pose.total_residue()
```




    152




```python
# 通过DSSP获取二级结构信息
from pyrosetta.rosetta.protocols.membrane import get_secstruct
get_secstruct(pose)
```

    [0mprotocols.DsspMover: {0} [0mLLHHHHHHHHHHHHHHHHHHHLLLLLEEEEEEEEELLEEEEEEEEEELLEEEEEEEEEEEELLEEEEEEEEEEEELLLHHHHHHHHHHHHHHHHHHHLLLLLEEEEEEEEELLEEEEEEEEEELLEEEEEEEEEEEELLEEEEEEEEEEEEL





    vector1_char[L, L, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, L, L, L, L, L, E, E, E, E, E, E, E, E, E, L, L, E, E, E, E, E, E, E, E, E, E, L, L, E, E, E, E, E, E, E, E, E, E, E, E, L, L, E, E, E, E, E, E, E, E, E, E, E, E, L, L, L, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, H, L, L, L, L, L, E, E, E, E, E, E, E, E, E, L, L, E, E, E, E, E, E, E, E, E, E, L, L, E, E, E, E, E, E, E, E, E, E, E, E, L, L, E, E, E, E, E, E, E, E, E, E, E, E, L]



#### 1.4.2 氨基酸信息
Residue Object是Pose的重要组成部分，整个Pose的Conformation是由一个个具体的氨基酸的具体构象组成，每个氨基酸有一个单独的object来描述，其中记录了所有的氨基酸信息。
同样，通过Pose类，我们可以轻松地访问每个氨基酸对象，并提取其中的信息。


```python
# 获取的第24个氨基酸残基对象(Residue objects)
residue24 = pose.residue(24)
residue24
```




    <pyrosetta.rosetta.core.conformation.Residue at 0x1127b6eb0>




```python
# 获取残基的Rosetta残基名、单字母缩写、三字母缩写、标注名
print(pose.residue(1).name())
print(pose.residue(2).name())
```

    PRO:NtermProteinFull
    SER


##### 结果解读:
Rosetta中N段和C段以及形成了二硫键的半胱氨酸等都是以"被修饰"的状态存在，因为这些末端或成环氨基酸与那些处于肽链中的氨基酸的原子数目不同，因此描述他们拓扑结构的params文件不能直接复用，因此在Rosetta中提出了Patch系统，给这些特殊的氨基酸打上补丁来描述他们的拓扑结构。</br>
因此为了区分，他们的命名也带上了补丁字样，比如上述的例子中，1号氨基酸名称为PRO:NtermProteinFull，其中NtermProteinFull就是他的"补丁名"。而第二号位于肽链中的氨基酸就是正常的三字母缩写的丝氨酸。


```python
# 获取氨基酸其他的缩写信息
print(pose.residue(1).name1())
print(pose.residue(1).name3())
print(pose.residue(1).annotated_name())
```

    P
    PRO
    P[PRO:NtermProteinFull]


##### 结果解读:
氨基酸的名称不止一种，还可以通过residue中的一些属性获得单字母缩写、三字母缩写以及标注名。</br>
以1号氨基酸为例，它的单字母缩写为P，三字母缩写为PRO，而标注名为单字母缩写外加Rosetta残基名。


```python
# 还可以直接打印这氨基酸的所有信息:
print(residue24)
```

    Residue 24: SER (SER, S):
    Base: SER
     Properties: POLYMER PROTEIN CANONICAL_AA SC_ORBITALS POLAR METALBINDING ALPHA_AA L_AA
     Variant types:
     Main-chain atoms:  N    CA   C  
     Backbone atoms:    N    CA   C    O    H    HA 
     Side-chain atoms:  CB   OG  1HB  2HB   HG 
    Atom Coordinates:
       N  : 64.583, -3.559, 23.925
       CA : 65.325, -2.363, 24.296
       C  : 65.522, -1.459, 23.088
       O  : 66.646, -1.048, 22.794
       CB : 66.679, -2.743, 24.897
       OG : 67.576, -3.19, 23.894
       H  : 65.0884, -4.36939, 23.5965
       HA : 64.7499, -1.81747, 25.0456
      1HB : 67.1044, -1.88054, 25.4094
      2HB : 66.54, -3.52835, 25.6388
       HG : 67.093, -3.14001, 23.0656
    Mirrored relative to coordinates in ResidueType: FALSE
    


##### 结果解读:
使用print打印信息后，可以获取这个对象中所有的信息:</br>
如: 氨基酸的类型为丝氨酸(SER), 骨架原子和侧链的原子组成、残基性质、以及每个单独原子的三维坐标信息


```python
# 通过Residue Object还可以判断氨基酸化学性质
print(pose.residue(5).is_polar())
print(pose.residue(5).is_aromatic())
print(pose.residue(5).is_charged())
```

    True
    False
    True


##### 结果解读:
可见半胱氨酸残基并不是极性残基、不含有疏水芳香环、也不带电。

#### 1.4.3 原子信息
每个Residue Object都是由Atom Object组成，其中记录了所有的原子基本属性信息:</br>
如原子名、元素类型、原子的一些化学和物理属性等。（**但不包括坐标信息**）


```python
# 获取第24号残基的原子的信息, 可见24号氨基酸共有11个原子
residue24.natoms()
```




    11




```python
# 获取每个原子的信息:
for atom_id in range(1, 11):
    atom_name = residue24.atom_name(atom_id)  # 获取原子名称
    print(f'atom_id is:{atom_id}, atom_name is: {atom_name}')
```

    atom_id is:1, atom_name is:  N  
    atom_id is:2, atom_name is:  CA 
    atom_id is:3, atom_name is:  C  
    atom_id is:4, atom_name is:  O  
    atom_id is:5, atom_name is:  CB 
    atom_id is:6, atom_name is:  OG 
    atom_id is:7, atom_name is:  H  
    atom_id is:8, atom_name is:  HA 
    atom_id is:9, atom_name is: 1HB 
    atom_id is:10, atom_name is: 2HB 


##### 结果解读:
此处我们获取了氨基酸内部的原子编号以及原子的PDB原子名


```python
# 获取前3个原子的属性信息
for atom_id in range(1, 4):
    atom_type = residue24.atom_type(atom_id)  # 获取原子名称
    print(atom_type)
```

    Atom Type: Nbb
    	element: N
    	Lennard Jones: radius=1.80245 wdepth=0.161725
    	Lazaridis Karplus: lambda=3.5 volume=15.992 dgfree=-9.96949
    	properties: DONOR 
    Extra Parameters: 1.8725 1.55 0.79 1.55 1.44 1.65 1.55 -5.95 -4.293 -1.145 -5 -0.47 2 1 0.65 1.85 8.52379 0.025 0.01 0.005 -289.292 -0.697267 -1933.88 -1.56243 -93.2613 93.2593 0.00202205 715.165 74.6559 -74.6539 0.00268963 -1282.36 0 0 0 0 0 0
    
    Atom Type: CAbb
    	element: C
    	Lennard Jones: radius=2.01176 wdepth=0.062642
    	Lazaridis Karplus: lambda=3.5 volume=12.137 dgfree=2.53379
    	properties: 
    Extra Parameters: 2.14 1.7 0.72 1.7 2.1285 1.87 1.75 -0.187 -0.487 -0.645 1 0.07 0 0 0 2.275 9.52518 0.025 0.01 0.005 746.028 -1.30263 -639.876 -1.99623 -73.5409 73.1524 0.0018258 1212.99 178.3 -178.298 0.00302069 -1035.77 0 0 0 0 0 0
    
    Atom Type: CObb
    	element: C
    	Lennard Jones: radius=1.91666 wdepth=0.141799
    	Lazaridis Karplus: lambda=3.5 volume=13.221 dgfree=3.10425
    	properties: 
    Extra Parameters: 2.14 1.7 0.72 1.7 1.89 1.76 1.65 0 0 0 1 0.51 0 0 0 2 8.81363 0.025 0.01 0.005 147.227 -0.811304 -8117.41 -2.17625 -85.8924 85.8904 0.00196363 900.14 168.481 -168.287 0.00113765 -6725.43 0 0 0 0 0 0
    


##### 结果解读:
通过atom_type方法，我们可以获取每个原子的细节的信息，如原子的Rosetta类型、原子的元素名、范德华半径、Lazaridis Karplus溶剂化参数等


```python
# 也可以通过原子名反向查找原子的残基内原子的ID编号:
ca_id = residue24.atom_index('CA')
N_id = residue24.atom_index('N')
print(ca_id, N_id)
```

    2 1


##### 结果解读:
24氨基酸中, N原子为1号原子, Cα为2号原子。


```python
# 获取原子的坐标的方式:
x, y, z = residue24.xyz("CA")
print(f'x: {x}, y:{y}, z:{z}')
```

    x: 65.325, y:-2.363, z:24.296



```python
# 通过原子序号原子的坐标的方式:
x, y, z = residue24.xyz(2)
print(f'x: {x}, y:{y}, z:{z}')
```

    x: 65.325, y:-2.363, z:24.296


#### 1.4.4 蛋白质几何信息
Pose中描述具体的蛋白质构象，除了氨基酸类型以外，更是由原子间的键长、键角，二面角等一系列的具体参数构成。Pose中的Conformation对象负责记录了这些连接信息。
![title](./img/phipsiomega.png)

为了定位原子的信息，首先需要构建atom identifier对象，相当于创建一个ID卡，让Rosetta知道我们指定的原子是位于哪个氨基酸中的。通过AtomID，提供残基号，原子号，就可以创建atom identifier对象


```python
# 获取原子间的键长、键角信息前需要构建atom identifier objects
from pyrosetta.rosetta.core.id import AtomID
atom1 = AtomID(1, 24)  # 24号残基的第一个原子
atom2 = AtomID(2, 24)  # 24号残基的第二个原子
atom3 = AtomID(3, 24)  # 24号残基的第三个原子
atom4 = AtomID(4, 24)  # 24号残基的第四个原子
print(atom1, atom2, atom3, atom4)
```

     atomno= 1 rsd= 24   atomno= 2 rsd= 24   atomno= 3 rsd= 24   atomno= 4 rsd= 24 


知道原子的ID后，就可以轻松的通过conformation对象来获取键长、键角等数据了。一般在Rosetta中键长和键角都设定为**理想值**。


```python
# 通过conformation层获取键长数据
bond_length = pose.conformation().bond_length(atom1, atom2)

# 通过conformation层获取键角数据(弧度)
bond_angle = pose.conformation().bond_angle(atom1, atom2, atom3)

print(f'原始键长:{bond_length}, 原始键角:{bond_angle}')
```

    原始键长:1.45554835027903, 原始键角:1.930305491631243



```python
# 通过pose获取氨基酸的骨架二面角数据
phi = pose.phi(24)
psi = pose.psi(24)
omega = pose.omega(24)
print(f'原始phi角:{phi}, 原始psi角:{psi}, 原始omega角:{omega}')
```

    原始phi角:-91.24043946940483, 原始psi角:51.92030782373845, 原始omega角:-174.71026990318242


### 1.5 Pose的操作
Pose中储存了非常多的信息，同时还提供了接口可以让用户方便的对其中的信息进行修改（采样）。

#### 1.5.1 Pose的创建和复制
前几节中提及过，pose是一个容器，理所当然我们可以创建一个空的容器，里面什么都不放。</br>
很多时候，我们需要创建空的Pose对象，便于保存当前pose的实例化状态，可作为可回溯点或初始状态，方便反复调用。


```python
# 通过创建一个新的Pose
from pyrosetta import Pose
starting_pose = Pose()
starting_pose2 = None
```

此处通过两种方法，将已有的Pose信息放入新的容器里，一种是通过assign函数复制，一种是通过python赋值符号来赋值。


```python
# 方法1：通过assign复制新的构象
starting_pose.assign(pose)

# 方法2：通过python的直接赋值符号
starting_pose2 = pose

print(starting_pose)
print('\n')
print(starting_pose2)
```

    PDB file name: ./data/4R80.pdb
    Total residues: 152
    Sequence: PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    Fold tree:
    FOLD_TREE  EDGE 1 76 -1  EDGE 1 77 1  EDGE 77 152 -1 
    
    
    PDB file name: ./data/4R80.pdb
    Total residues: 152
    Sequence: PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    Fold tree:
    FOLD_TREE  EDGE 1 76 -1  EDGE 1 77 1  EDGE 77 152 -1 


**结果解读:</br>
可见，两种方式“看”起来里面都有了新的pose信息。但真的如此么？**


```python
# 尝试调整mypose中的24号氨基酸phi值:
pose.set_phi(24, 170.0)

# 查看对starting_pose以及starting_pose2的影响:
print(f'starting_pose 24 residue phi:{starting_pose.phi(24)}')
print(f'starting_pose2 24 residue phi:{starting_pose2.phi(24)}')
```

    starting_pose 24 residue phi:-91.24043946940483
    starting_pose2 24 residue phi:170.0


##### 结果解读:
结果可见，starting_pose2是用过python直接赋值的Pose并没有复制，而只是pose的一个"超链接"符，并没有进行"复制"的操作。
而通过assign的复制，原始的pose的任何调整都没有对starting_pose造成任何的影响，可见其独立性。

#### 1.5.2 链与氨基酸的增减操作
除了对整体信息的迁移，用户还可以对Pose中的链以及氨基酸进行操作。

##### **1.5.2.1 链的切割处理**

尽管pose的氨基酸编号是忽略多肽链的分隔的，但是pose中的链依然是根据多肽链的物理结构进行编号的，同理也是从1开始编号。</br>
如一个蛋白中有2条链A和B，那么链编号结果即为1和2。其中A对应1号链，B对应2号链，与PDB的链顺序有关（当然A链的顺序如果在后面，那么B链就是1号链）。</br>
Pose类中许多的方法可以很方便对链的增减进行操作, 以下2个举例进行说明:


```python
# 将Pose按照链的数量进行切割
pose_list = pose.split_by_chain()
pose_list
```




    vector1_std_shared_ptr_core_pose_Pose_t[0x7fc4a287d1a0, 0x7fc4a2681d80]



此处的pose_list中存放了2个数据，说明链已经被切割成2个独立的pose对象了。</br>
通过python的索引，可以获得具体的pose:


```python
# 获取只含有第一个链的pose
chain1_pose = pose.split_by_chain()[1]  # 直接切片索引链号。
chain2_pose = pose.split_by_chain()[2]  # 直接切片索引链号。

# check
print(f'chain1_pose中的氨基酸总数:{chain1_pose.total_residue()}')
print(f'chain2_pose中的氨基酸总数:{chain2_pose.total_residue()}')
print(f'原始pose中的氨基酸总数:{pose.total_residue()}')
```

    chain1_pose中的氨基酸总数:76
    chain2_pose中的氨基酸总数:76
    原始pose中的氨基酸总数:152


##### **1.5.2.2 链的合并处理**

除了分隔操作，用户还可以通过一些简单的方式把链合并到一个pose中，此处使用append_pose_to_pose函数就可以达到目的。但需要注意，pose中的氨基酸、链的数量变化后，都需要对PDBinfo进行更新。否则PDBinfo的信息与Pose信息不对称。


```python
# 两条链的合并;
# add binder to pose;
from pyrosetta.rosetta.core.pose import append_pose_to_pose
print(f'原始chain1_pose中的氨基酸总数:{chain1_pose.total_residue()}')
append_pose_to_pose(chain1_pose, chain2_pose)
chain1_pose.update_residue_neighbors()
chain1_pose.update_pose_chains_from_pdb_chains()
chain1_pose.conformation().detect_disulfides()

# update pdbinfo; 别忘了更新pdbinfo;
# 更新pdb_info; [别忘了]
from pyrosetta.rosetta.core.pose import renumber_pdbinfo_based_on_conf_chains
renumber_pdbinfo_based_on_conf_chains(pose)

print(f'append之后的chain1_pose中的氨基酸总数:{chain1_pose.total_residue()}')
chain1_pose.sequence()

# 检查PDBinfo是否正确: Returns true if PDBInfo is obsolete and needs updating
print(f'PDBinfo是否需要被更新:{pose.pdb_info().obsolete()}')
```

    原始chain1_pose中的氨基酸总数:76
    append之后的chain1_pose中的氨基酸总数:152
    PDBinfo是否需要被更新:False


##### **1.5.3.3 氨基酸的删减操作**
除了对链的合并之外，我们还可以对链中的氨基酸进行添加、删除的操作！具体的过程是用户需要创建一个独立的氨基酸(residue object)，并将这个氨基酸加载到现有的构像中。</br>
加载的方式可以是前置后后置，根据使用的函数不同而定。


```python
# 在链的前端添加新的氨基酸或删除氨基酸
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.chemical import ChemicalManager

print(f'原始氨基酸总数:{pose.total_residue()}')
print(f'原始氨基酸序列:{pose.sequence()}\n')

# 向前添加氨基酸
chm = ChemicalManager.get_instance()
rts = chm.residue_type_set("fa_standard")
new_rsd = ResidueFactory.create_residue(rts.name_map('ALA')) # 创建一个residue object
pose.prepend_polymer_residue_before_seqpos(new_rsd, 1, True)  # 在第一个氨基酸前添加一个ALA

print(f'向前添加之后氨基酸总数:{pose.total_residue()}')
print(f'向前添加之后氨基酸序列:{pose.sequence()}\n')
```

    原始氨基酸总数:152
    原始氨基酸序列:PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    
    [0mcore.conformation.Residue: {0} [0m[1m[ WARNING ][0m Residue connection id changed when creating a new residue at seqpos 1
    [0mcore.conformation.Residue: {0} [0m[1m[ WARNING ][0m ResConnID info stored on the connected residue (residue 2) is now out of date!
    [0mcore.conformation.Residue: {0} [0m[1m[ WARNING ][0m Connection atom name (in src):  C
    向前添加之后氨基酸总数:153
    向前添加之后氨基酸序列:APSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGS
    



```python
# 向后添加氨基酸
last_residue = pose.total_residue()
pose.append_polymer_residue_after_seqpos(new_rsd, last_residue, True)  # 在第一个氨基酸前添加一个ALA

print(f'向后添加之后氨基酸总数:{pose.total_residue()}')
print(f'向后添加之后氨基酸序列:{pose.sequence()}\n')
```

    向后添加之后氨基酸总数:154
    向后添加之后氨基酸序列:APSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSA
    



```python
# 删除氨基酸
pose.delete_polymer_residue(1)  # 删除第一个氨基酸

print(f'删除之后氨基酸总数:{pose.total_residue()}')
print(f'删除之后氨基酸序列:{pose.sequence()}\n')
```

    删除之后氨基酸总数:153
    删除之后氨基酸序列:PSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSA
    



```python
# 还可以范围性的删除氨基酸
pose.delete_residue_range_slow(1,5) # 删除第一个至第五个氨基酸

print(f'删除之后氨基酸总数:{pose.total_residue()}')
print(f'删除之后氨基酸序列:{pose.sequence()}\n')
```

    删除之后氨基酸总数:148
    删除之后氨基酸序列:EKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSPSEEEEKRRAKQVAKEKILEQNPSSKVQVRRVQKQGNTIRVELEITENGKKTNITVEVEKQGNTFTVKRITETVGSA
    


##### **1.5.3.4 PBDinfo更新**


```python
# 更新pdb_info; [别忘了!]
from pyrosetta.rosetta.core.pose import renumber_pdbinfo_based_on_conf_chains

renumber_pdbinfo_based_on_conf_chains(pose)  # 更新PDBinfo.

# 检查PDBinfo是否正确: Returns true if PDBInfo is obsolete and needs updating
print(f'PDBinfo是否需要被更新:{pose.pdb_info().obsolete()}')
```

    PDBinfo是否需要被更新:False


#### 1.5.3 构象的调整

除了对多肽链的氨基酸数量的调整，我们还可以通过Pose中的一些函数来调整蛋白质的具体构象，如主链的phi/psi角、化学键中的键长与键角数据等。

##### **1.5.3.1 化学键的数据调整**


```python
# 修改键长键角必须通过conformation层进行处理:
print(f'原始键长:{bond_angle}, 原始键角:{bond_length}')

pose.conformation().set_bond_angle(atom1, atom2, atom3, 0.66666 * 3.14)
new_bond_angle = pose.conformation().bond_angle(atom1, atom2, atom3)

pose.conformation().set_bond_length(atom1, atom2, 1.44)
new_bond_length = pose.conformation().bond_length(atom1, atom2)

print(f'新的键长:{new_bond_length}, 新的键角:{new_bond_angle}')
```

    原始键长:1.930305491631243, 原始键角:1.45554835027903
    新的键长:1.44, 新的键角:2.0933124000000003



```python
# 修改phi、psi、chi、omega角可以直接通过pose的函数:
# 通过pose获取氨基酸的骨架二面角数据
print(f'原始phi角:{pose.phi(24)}, 原始psi角:{pose.psi(24)}, 原始omega角:{pose.omega(24)}')
pose.set_phi(24, 66.0)
pose.set_psi(24, 55.0)
pose.set_omega(24, 180.0)

print(f'调整后phi角:{pose.phi(24)}, 调整后psi角:{pose.psi(24)}, 调整后omega角:{pose.omega(24)}')
```

    原始phi角:-97.41405901201394, 原始psi角:125.84230362614217, 原始omega角:-174.63774486370116
    调整后phi角:66.0, 调整后psi角:55.0, 调整后omega角:180.0


##### **1.5.3.2 氨基酸类型的调整(突变)**
除了具体的化学键数据的调整，在PyRosetta中进行氨基酸的类型调整也是很方便的


```python
# 调整氨基酸的类型
from pyrosetta.toolbox import mutate_residue
print(f'原始氨基酸类型:{pose.residue(1).name()}')
print('突变氨基酸中...')
mutate_residue(pose, 1, 'A', 9.0)  # 1 代表氨基酸突变的pose编号，9.0代表对氨基酸附近9埃范围内的氨基酸进行侧链优化，适应新的突变。
print(f'突变后氨基酸类型:{pose.residue(1).name()}')
```

    原始氨基酸类型:PRO:NtermProteinFull
    突变氨基酸中...
    [0mcore.scoring.ScoreFunctionFactory: {0} [0mSCOREFUNCTION: [32mref2015[0m
    [0mcore.scoring.etable: {0} [0mStarting energy table calculation
    [0mcore.scoring.etable: {0} [0msmooth_etable: changing atr/rep split to bottom of energy well
    [0mcore.scoring.etable: {0} [0msmooth_etable: spline smoothing lj etables (maxdis = 6)
    [0mcore.scoring.etable: {0} [0msmooth_etable: spline smoothing solvation etables (max_dis = 6)
    [0mcore.scoring.etable: {0} [0mFinished calculating energy tables.
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/hbonds/ref2015_params/HBPoly1D.csv
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/hbonds/ref2015_params/HBFadeIntervals.csv
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/hbonds/ref2015_params/HBEval.csv
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/hbonds/ref2015_params/DonStrength.csv
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/hbonds/ref2015_params/AccStrength.csv
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/rama/fd/all.ramaProb
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/rama/fd/prepro.ramaProb
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/omega/omega_ppdep.all.txt
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/omega/omega_ppdep.gly.txt
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/omega/omega_ppdep.pro.txt
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/omega/omega_ppdep.valile.txt
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/P_AA_pp/P_AA
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/P_AA_pp/P_AA_n
    [0mcore.scoring.P_AA: {0} [0mshapovalov_lib::shap_p_aa_pp_smooth_level of 1( aka low_smooth ) got activated.
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/P_AA_pp/shapovalov/10deg/kappa131/a20.prop
    [0mcore.pack.task: {0} [0mPacker task: initialize from command line()
    [0mbasic.io.database: {0} [0mDatabase file opened: scoring/score_functions/elec_cp_reps.dat
    [0mcore.scoring.elec.util: {0} [0mRead 40 countpair representative atoms
    [0mcore.pack.dunbrack.RotamerLibrary: {0} [0mshapovalov_lib_fixes_enable option is true.
    [0mcore.pack.dunbrack.RotamerLibrary: {0} [0mshapovalov_lib::shap_dun10_smooth_level of 1( aka lowest_smooth ) got activated.
    [0mcore.pack.dunbrack.RotamerLibrary: {0} [0mBinary rotamer library selected: /opt/miniconda3/lib/python3.7/site-packages/pyrosetta/database/rotamer/shapovalov/StpDwn_0-0-0/Dunbrack10.lib.bin
    [0mcore.pack.dunbrack.RotamerLibrary: {0} [0mUsing Dunbrack library binary file '/opt/miniconda3/lib/python3.7/site-packages/pyrosetta/database/rotamer/shapovalov/StpDwn_0-0-0/Dunbrack10.lib.bin'.
    [0mcore.pack.dunbrack.RotamerLibrary: {0} [0mDunbrack 2010 library took 0.180545 seconds to load from binary
    [0mcore.pack.pack_rotamers: {0} [0mbuilt 156 rotamers at 8 positions.
    [0mcore.pack.pack_rotamers: {0} [0mRequesting all available threads for interaction graph computation.
    [0mcore.pack.interaction_graph.interaction_graph_factory: {0} [0mInstantiating PDInteractionGraph
    [0mbasic.thread_manager.RosettaThreadManager: {?} [0mCreating a thread pool of 16 threads.
    [0mbasic.thread_manager.RosettaThread: {?} [0mLaunching thread 3.
    [0mbasic.thread_manager.RosettaThread: {?} [0mLaunching thread 2.
    [0mbasic.thread_manager.RosettaThread: {?} [0mLaunching thread 9.
    [0mbasic.thread_manager.RosettaThread: {?} [0mLaunching thread 4.
    [0mbasic.thread_manager.RosettaThread: {?} [0mLaunching thread 7.
    [0mbasic.thread_manager.RosettaThreadPool: {?} [0mLaunched 15 new threads.
    [0mbasic.thread_manager.RosettaThread: {10} [0mLaunching thread 10.
    [0mbasic.thread_manager.RosettaThread: {6} [0mLaunching thread 6.
    [0mbasic.thread_manager.RosettaThread: {5} [0mLaunching thread 5.
    [0mbasic.random.init_random_generator: {4} [0m'RNG device' seed mode, using '/dev/urandom', seed=529492867 seed_offset=0 real_seed=529492871 thread_index=4
    [0mbasic.random.init_random_generator: {4} [0mRandomGenerator:init: Normal mode, seed=529492871 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThread: {13} [0mLaunching thread 13.
    [0mbasic.thread_manager.RosettaThread: {1} [0mLaunching thread 1.
    [0mbasic.random.init_random_generator: {3} [0m'RNG device' seed mode, using '/dev/urandom', seed=-137995808 seed_offset=0 real_seed=-137995805 thread_index=3
    [0mbasic.random.init_random_generator: {3} [0mRandomGenerator:init: Normal mode, seed=-137995805 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThread: {11} [0mLaunching thread 11.
    [0mbasic.thread_manager.RosettaThread: {8} [0mLaunching thread 8.
    [0mbasic.random.init_random_generator: {2} [0m'RNG device' seed mode, using '/dev/urandom', seed=1311789193 seed_offset=0 real_seed=1311789195 thread_index=2
    [0mbasic.random.init_random_generator: {2} [0mRandomGenerator:init: Normal mode, seed=1311789195 RG_type=mt19937
    [0mbasic.random.init_random_generator: {7} [0m'RNG device' seed mode, using '/dev/urandom', seed=1028457733 seed_offset=0 real_seed=1028457740 thread_index=7
    [0mbasic.random.init_random_generator: {7} [0mRandomGenerator:init: Normal mode, seed=1028457740 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThread: {14} [0mLaunching thread 14.
    [0mbasic.random.init_random_generator: {10} [0m'RNG device' seed mode, using '/dev/urandom', seed=-1085986616 seed_offset=0 real_seed=-1085986606 thread_index=10
    [0mbasic.random.init_random_generator: {10} [0mRandomGenerator:init: Normal mode, seed=-1085986606 RG_type=mt19937
    [0mbasic.random.init_random_generator: {1} [0m'RNG device' seed mode, using '/dev/urandom', seed=-1541843351 seed_offset=0 real_seed=-1541843350 thread_index=1
    [0mbasic.random.init_random_generator: {1} [0mRandomGenerator:init: Normal mode, seed=-1541843350 RG_type=mt19937
    [0mbasic.random.init_random_generator: {13} [0m'RNG device' seed mode, using '/dev/urandom', seed=164600221 seed_offset=0 real_seed=164600234 thread_index=13
    [0mbasic.random.init_random_generator: {6} [0m'RNG device' seed mode, using '/dev/urandom', seed=-928086000 seed_offset=0 real_seed=-928085994 thread_index=6
    [0mbasic.random.init_random_generator: {6} [0mRandomGenerator:init: Normal mode, seed=-928085994 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThread: {12} [0mLaunching thread 12.
    [0mbasic.random.init_random_generator: {14} [0m'RNG device' seed mode, using '/dev/urandom', seed=542029705 seed_offset=0 real_seed=542029719 thread_index=14
    [0mbasic.random.init_random_generator: {13} [0mRandomGenerator:init: Normal mode, seed=164600234 RG_type=mt19937
    [0mbasic.random.init_random_generator: {5} [0m'RNG device' seed mode, using '/dev/urandom', seed=1503161305 seed_offset=0 real_seed=1503161310 thread_index=5
    [0mbasic.random.init_random_generator: {5} [0mRandomGenerator:init: Normal mode, seed=1503161310 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThread: {15} [0mLaunching thread 15.
    [0mbasic.random.init_random_generator: {8} [0m'RNG device' seed mode, using '/dev/urandom', seed=-1441949229 seed_offset=0 real_seed=-1441949221 thread_index=8
    [0mbasic.random.init_random_generator: {8} [0mRandomGenerator:init: Normal mode, seed=-1441949221 RG_type=mt19937
    [0mbasic.random.init_random_generator: {14} [0mRandomGenerator:init: Normal mode, seed=542029719 RG_type=mt19937
    [0mbasic.random.init_random_generator: {11} [0m'RNG device' seed mode, using '/dev/urandom', seed=-1420986851 seed_offset=0 real_seed=-1420986840 thread_index=11
    [0mbasic.random.init_random_generator: {11} [0mRandomGenerator:init: Normal mode, seed=-1420986840 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThreadManager: {4} [0mThread 4 completed 22 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {0} [0mThread 0 completed 7 of 39 work units.
    [0mbasic.random.init_random_generator: {15} [0m'RNG device' seed mode, using '/dev/urandom', seed=-400430038 seed_offset=0 real_seed=-400430023 thread_index=15
    [0mbasic.random.init_random_generator: {15} [0mRandomGenerator:init: Normal mode, seed=-400430023 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThreadManager: {1} [0mThread 1 completed 1 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {5} [0mThread 5 completed 1 of 39 work units.
    [0mbasic.random.init_random_generator: {9} [0m'RNG device' seed mode, using '/dev/urandom', seed=1363095172 seed_offset=0 real_seed=1363095181 thread_index=9
    [0mbasic.thread_manager.RosettaThreadManager: {3} [0mThread 3 completed 1 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {10} [0mThread 10 completed 2 of 39 work units.
    [0mbasic.random.init_random_generator: {9} [0mRandomGenerator:init: Normal mode, seed=1363095181 RG_type=mt19937
    [0mbasic.random.init_random_generator: {12} [0m'RNG device' seed mode, using '/dev/urandom', seed=1010233423 seed_offset=0 real_seed=1010233435 thread_index=12
    [0mbasic.random.init_random_generator: {12} [0mRandomGenerator:init: Normal mode, seed=1010233435 RG_type=mt19937
    [0mbasic.thread_manager.RosettaThreadManager: {7} [0mThread 7 completed 1 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {11} [0mThread 11 completed 0 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {15} [0mThread 15 completed 0 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {9} [0mThread 9 completed 0 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {14} [0mThread 14 completed 0 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {13} [0mThread 13 completed 1 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {6} [0mThread 6 completed 1 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {8} [0mThread 8 completed 1 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {12} [0mThread 12 completed 0 of 39 work units.
    [0mbasic.thread_manager.RosettaThreadManager: {2} [0mThread 2 completed 1 of 39 work units.
    [0mcore.pack.rotamer_set.RotamerSets: {0} [0mCompleted interaction graph pre-calculation in 16 available threads (16 had been requested).
    突变后氨基酸类型:ALA:NtermProteinFull


##### **1.5.3.3 原子坐标的修改**
原子坐标的修改需要获取residue对象，并获取原子ID(atom identifier objects)。通过pose.set_xyz函数设定新的xyz坐标, 但用户一般不需要”显式“地修改原子坐标, 除非你明白这样操作的意义。</br>
此处以创建一个镜像原子进行说明:


```python
# 原子坐标的修改（一般不需要这样操作）
from pyrosetta.rosetta.numeric import xyzVector_double_t

# 对第24个氨基酸的所有原子的x坐标乘上一个负号:
residue24 = pose.residue(24)  # 获取residue对象
for atom_id, atom in enumerate(residue24.atoms()):
    x, y, z = atom.xyz()
    print(f'坐标进行修改前信息: 原子号:{atom_id+1}, x:{x}, y:{y}, z:{z}')
    
    mirror_xyz = xyzVector_double_t(-x, y, z)  # 乘上负号.
    atom_index = AtomID(atom_id+1, 24)   # 24号氨基酸的第x个原子的id
    pose.set_xyz(atom_index, mirror_xyz) # 设置xyz坐标

print('\n')
    
for atom_id, atom in enumerate(residue24.atoms()):
    x, y, z = atom.xyz()
    print(f'坐标进行修改后信息:  原子号:{atom_id+1}, x:{x}, y:{y}, z:{z}')
```

    坐标进行修改前信息: 原子号:1, x:76.92159052534156, y:-8.174496228716757, z:15.98103320685106
    坐标进行修改前信息: 原子号:2, x:77.36438744744312, y:-9.473088221114633, z:15.543780736237686
    坐标进行修改前信息: 原子号:3, x:76.41027860252882, y:-10.656849105413853, z:15.59304005880104
    坐标进行修改前信息: 原子号:4, x:76.72016199264006, y:-11.68632780177099, z:16.190650214571495
    坐标进行修改前信息: 原子号:5, x:78.5769299962636, y:-9.959911626299197, z:16.357609673708236
    坐标进行修改前信息: 原子号:6, x:79.01274945042509, y:-11.33892387835475, z:15.887363342457212
    坐标进行修改前信息: 原子号:7, x:79.72344508837577, y:-8.968056999916671, z:16.245312996630993
    坐标进行修改前信息: 原子号:8, x:76.84663984966622, y:-7.427148755840927, z:15.305790931069307
    坐标进行修改前信息: 原子号:9, x:77.66402417738198, y:-9.402369247935077, z:14.49748949847452
    坐标进行修改前信息: 原子号:10, x:78.28577338752454, y:-10.056135282169324, z:17.403697589934474
    坐标进行修改前信息: 原子号:11, x:79.87110134338924, y:-11.66882632897991, z:16.473256188412318
    坐标进行修改前信息: 原子号:12, x:78.19218533052184, y:-12.044084133499098, z:16.01828274769628
    坐标进行修改前信息: 原子号:13, x:79.28925340858027, y:-11.293585810152619, z:14.8341891608537
    坐标进行修改前信息: 原子号:14, x:80.57253875808294, y:-9.326278855525196, z:16.826662284346718
    坐标进行修改前信息: 原子号:15, x:80.0157555951842, y:-8.86655858792907, z:15.199972692734988
    坐标进行修改前信息: 原子号:16, x:79.40476467527382, y:-7.998544441100375, z:16.628662238398416
    
    
    坐标进行修改后信息:  原子号:1, x:-76.92159052534156, y:-8.174496228716757, z:15.98103320685106
    坐标进行修改后信息:  原子号:2, x:-77.36438744744312, y:-9.473088221114633, z:15.543780736237686
    坐标进行修改后信息:  原子号:3, x:-76.41027860252882, y:-10.656849105413853, z:15.59304005880104
    坐标进行修改后信息:  原子号:4, x:-76.72016199264006, y:-11.68632780177099, z:16.190650214571495
    坐标进行修改后信息:  原子号:5, x:-78.5769299962636, y:-9.959911626299197, z:16.357609673708236
    坐标进行修改后信息:  原子号:6, x:-79.01274945042509, y:-11.33892387835475, z:15.887363342457212
    坐标进行修改后信息:  原子号:7, x:-79.72344508837577, y:-8.968056999916671, z:16.245312996630993
    坐标进行修改后信息:  原子号:8, x:-76.84663984966622, y:-7.427148755840927, z:15.305790931069307
    坐标进行修改后信息:  原子号:9, x:-77.66402417738198, y:-9.402369247935077, z:14.49748949847452
    坐标进行修改后信息:  原子号:10, x:-78.28577338752454, y:-10.056135282169324, z:17.403697589934474
    坐标进行修改后信息:  原子号:11, x:-79.87110134338924, y:-11.66882632897991, z:16.473256188412318
    坐标进行修改后信息:  原子号:12, x:-78.19218533052184, y:-12.044084133499098, z:16.01828274769628
    坐标进行修改后信息:  原子号:13, x:-79.28925340858027, y:-11.293585810152619, z:14.8341891608537
    坐标进行修改后信息:  原子号:14, x:-80.57253875808294, y:-9.326278855525196, z:16.826662284346718
    坐标进行修改后信息:  原子号:15, x:-80.0157555951842, y:-8.86655858792907, z:15.199972692734988
    坐标进行修改后信息:  原子号:16, x:-79.40476467527382, y:-7.998544441100375, z:16.628662238398416


### 1.6 Pose的能量
如果现在我们已有一个Pose，我们想评估这个蛋白质的能量评分，可以直接通过创建一个打分函数并对Pose的能量进行计算，</br>
再可以通过pose中的energies对象将氨基酸残基的One-body, Two-body的能量信息列出，达到残基级别能量"分解"的目的。

#### 1.6.1 对结构进行能量计算
create_score_function函数可以用于快速创建一个打分函数。


```python
## 创建标准打分函数
from pyrosetta import create_score_function
scorefxn = create_score_function('ref2015')

# 对当前Pose中的构象进行能量计算
weighted_total_score = scorefxn(pose)
print(weighted_total_score)
```

    -53.972765295810895


#### 1.6.2 能量信息
Rosetta的能量是加权后的能量，每个能量项有自己的权重，通过energies获取的是能量项的原始结果（unweighted）。


```python
# 获取能量对象
scores = pose.energies()

# 获取1号残基的所有能量项的信息:
print(scores.show(1))
```

    [0mcore.scoring.Energies: {0} [0mE               fa_atr        fa_rep        fa_sol  fa_intra_repfa_intra_sol_x   lk_ball_wtd       fa_elec     pro_close   hbond_sr_bb   hbond_lr_bb   hbond_bb_sc      hbond_sc     dslf_fa13         omega        fa_dun       p_aa_pp yhh_planarity           ref   rama_prepro
    [0mcore.scoring.Energies: {0} [0mE(i)   1         -1.79          0.08          1.44          0.28          0.00         -0.23         -0.35          0.00          0.00          0.00          0.00          0.00          0.00          0.02          0.00          0.00          0.00          1.32          0.00
    None


Rosetta中的Score项有许多，如fa_atr代表范德华吸引势力, fa_rep代表范德华排斥项, fa_elec代表静电项等。</br>
比如通过ScoreType类下的属性，即可搜索到对应的能量项。</br>
获取总能中的某一个项的值可以直接使用python的索引功能，十分方便:


```python
# 获取总fa_atr项能量项的得分结果:
from pyrosetta.rosetta.core.scoring import ScoreType
pose.energies().total_energies()[ScoreType.fa_atr]
```




    -793.2904073710969




```python
# 单独获得第5号氨基酸的fa_atr项能量得分:
pose.energies().residue_total_energies(5)[ScoreType.fa_atr]
```




    -5.85223523327472



### 1.7 自定义信息
Pose中含有让用户自定义写入任何信息的功能，比如在程序设计过程中，中间生成的临时数值或字符都可以写入到PoseExtraScore中，这些信息会随着Pose一并输出到PDB或则Silent文件中，在后续的分析和处理的过程中非常方便。


```python
# 给pose加入额外的信息: 比如filter计算的值就可以储存.
from pyrosetta.rosetta.core.pose import setPoseExtraScore, getPoseExtraScore

setPoseExtraScore(pose, "distance", 1.0)
setPoseExtraScore(pose, "angle", '120.5')

# 提取信息
print(getPoseExtraScore(pose, 'distance'))
print(getPoseExtraScore(pose, 'angle'))  # 目前有bug，但是信息已经储存在pose中了
```

    1.0



    ---------------------------------------------------------------------------

    RuntimeError                              Traceback (most recent call last)

    <ipython-input-46-dacdb1758aa1> in <module>
          7 # 提取信息
          8 print(getPoseExtraScore(pose, 'distance'))
    ----> 9 print(getPoseExtraScore(pose, 'angle'))  # 目前有bug，但是信息已经储存在pose中了
    

    RuntimeError: 
    
    File: /Volumes/MacintoshHD3/benchmark/W.fujii.release/rosetta.Fujii.release/_commits_/main/source/src/core/pose/extra_pose_info_util.cc:297
    [ ERROR ] UtilityExitException
    ERROR: Assertion `getPoseExtraScore( pose, name, value )` failed.
    



## 课后作业

1. 通过Pose对象，生成骨架多肽链的ContactMap（提示:xyz坐标）

2. 尝试将一个四聚体进行2等份的切割并输出PDB结构。


```python

```
