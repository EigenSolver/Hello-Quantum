# 量子渐进优化算法（QAOA）阅读笔记


## 摘要
量子近似优化算法（Quantum Approximate Optimization Algorithm）是一种经典量子混合算法（Classical-Quantum Hybrid），用于解决经典计算机难以处理的组合优化问题（Combinatorial Optimization Problem）。

## 定义问题
组合优化问题（Combinatorial Optimization Problem）是一类NP-Complete问题的统称，因而无法被经典计算机在多项式时间内解决。然而，很多组合优化问题在应用中具有重要价值，为人熟知的旅行商（TSP，Travelling Salesman Problem）问题，最小生成树问题（Minimum Spanning Tree Problem， MST），最大切割问题（Max Cut Problem）和布尔可满足性问题（Boolean Satisﬁability Problem
）等都属于组合优化问题。参见Wiki$^{[1]}$。

大部分组合优化问题都可以抽象为图论问题，其形式化表述取决于具体问题。以最大切割问题（MaxCut）为例，最大切割问题要求我们在给定的图中找到一条曲线，使得该曲线与图中相交的边的数量取得最大值，其中每条边只能被切割一次。

![Wiki: MaxCut](https://upload.wikimedia.org/wikipedia/commons/c/cf/Max-cut.svg)

假设给定图
$G=(E,V,W)$
其中$E$图中包含的边，$V$为图的节点，$W$为边的权重。

假设图中包含$n$个节点，即$||V||=n$，那么MaxCut问题等效于给出对这$n$个节点的分割映射

$f: V\mapsto \{0,1\}$

其中如果两个节点$i,j\in V$有$f(i)\neq f(j)$则说明$i,j$间的边被切割，$f$使得目标函数
$$O=\sum_{f(i)\neq f(j)}^{ij\in E}w_{ij}$$
取最大值。

## 全局最优解-绝热量子计算（Quantum Adiabatic Computation）
实际上对于优化问题，量子渐进优化算法并不能找到全局最优解（Global Optimal），其仅仅是全局求解的绝热算法（Quantum Adiabatic Algorithm, QAA）的近似形式$^{[2]}$。旨在通过离散化的角度演化（Trotterization）近似QAA中的时间演化使得绝热退火（Adiabatic Annealing）过程$^{[3]}$可以在电路模型（Circuit Model）中执行，同时取得更高的效率。

### 绝热理论（Adiabatic Theorem）
绝热理论（Adiabatic Theorem）描述量子系统在绝热过程下的演化过程（Adiabatic Evolution）。

绝热过程通常是外界环境变化极其缓慢，以致于哈密顿量的变化率近似为零$\dot H \simeq0$。当系统哈密顿量通过绝热过程演化时，系统所处的量子态同样从前初态哈密顿量的本征态演化到末态哈密顿量的对应本征态。

例如系统哈密顿量在绝热过程下演化

$H_i \mapsto H_j$

那么系统所处态应演化如下

$|\psi_i^n \rangle \mapsto |\psi_j^n \rangle$。

绝热过程的具体理论可见格里菲斯的量子力学导论[5]。

### 使用量子系统解决传统问题-构造哈密顿量（Hamiltonian）
以前文提到的MaxCut问题为例，我们可以针对这一问题构造一个Hamiltonian，使得系统基态给出MaxCut问题的最优解。事实上这也是应用绝热量子计算或者量子退火算法解决优化问题的根本方法$^{[3][4]}$。

例如MaxCut问题中的分割映射$f: V\mapsto \{0,1\}$，等效于对$V$中的n个节点各自标记赋值0或1。其组合可以用n个bit的字符串（bit string）表示，这一字符串可以用包含n个qubit的量子态表示。例如$|0101110101\rangle$表示一个具有10个节点的图，其中标记分别为0和1的两个节点间的边被切割。

目标函数可以写为

$$C(z)=\frac{1}{2}\sum _{\langle i,j \rangle \in E}(z_i\oplus z_j)w_{ij}$$

其中$z_i$表示字符串的$i$位，而$\oplus$为XOR算符。在不考虑权重的情况下，如果采用$f: V\mapsto \{-1,1\}$标记，则原函数可写为

$$C(z)=\frac{1}{2}\sum _{\langle i,j \rangle \in E}(1-z_i z_j)$$

可以看出当且仅当$z_i,z_j$异号（分别为1和-1时）目标函数才不会抵消，优化得到$C(z)$的最大值即为最大切割数。

由此定义目标函数的哈密顿量

$$H_z=\sum _{\langle i,j \rangle \in E} (\sigma ^z_i \sigma ^z_j -1)$$

其中
$\sigma ^z_i, \sigma ^z_j$分别表示作用于第$i,j$个qubit的Pauli Z算符。注意到$|0\rangle$和$|1\rangle$均为PauliZ的本征态且本征值分别为1和-1，因而构造的Hamiltonian与目标函数等价。

使该Hamiltonian处于能量最小值的基态即对应MaxCut问题的最优解。

### 量子绝热算法（Quantum Adiabatic Algorithm）- 量子退火算法（Quantum Annealing）

考虑如下的绝热演化过程

$$\hat H(s)=s\hat H_z+(1-s)\hat H_x$$

其中$s=s(t)$随时间变化，且$s(0)=0,s(T)=1$

显然，$t: 0\rightarrow T$ 对应 $\hat H(s): \hat H_x \rightarrow \hat H_z$

在绝热量子计算中，通常系统初态为较易制备的哈密顿量，例如Pauli-X，对于MaxCut问题中的多qubit系统 

$$\hat H_x=\sum_{j=1}^n \sigma^x_j$$

其中 $\sigma^x_j$ 作用于第 $j$ 个qubit上的Pauli-X算符 $\hat\sigma_x$。此处的 $n$ 是对应组合优化问题中的节点数量。

$\hat H_x$ 对应系统基态为 $|s\rangle=|-\rangle_1|-\rangle_2 \dots|-\rangle_n$



**如前文绝热理论所述，如果我们可以制备一个处于基态 $|s\rangle$ 的量子系统，并设法使其哈密顿量绝从易制备的初态 $\hat H_x$ 绝热演化至我们想要的目标函数对应的目标哈密顿量 $\hat H_z$，那么演化结束后系统仍处于基态——$\hat H_z$ 的基态$|s\rangle'$，就等效于完成了这一优化的计算过程，这就是量子绝热计算（也称量子退火算法，Quantum Annealing）的基本原理。$^{[3][4]}$**


## 近似解-量子渐进优化算法（Quantum Approximate Optimization Algorithm）


### 时间演化算符离散化近似（Trotter Splitting）

如前文所述，实际上QAOA是对QAA的一种近似。

在量子力学$^{[5]}$中，时变哈密顿量带的的态的变化可以用时间演化算符描述。
$$
\psi(t)=\hat U(t) \psi(0)
$$
对于QAA中的绝热演化过程，考虑时间演化算符
$$
\hat U(t)=exp\left(-\frac{i}{\hbar}\int_0^t \hat H(s(\tau)) d\tau \right)
$$

其指数项上的连续积分可以用离散求和近似（类比数值积分算法）。


$$
\int_0^t\hat H(s(\tau )) d\tau \simeq \sum_{n=1}^N \hat H(s_n)\Delta t_n
$$

![Wiki: Numerical Integration](https://upload.wikimedia.org/wikipedia/commons/f/f2/Integral_as_region_under_curve.svg)

其中，$N$ 为离散区间的总数， $\sum_{n=1}^N \Delta t_n=T$，$s_n=s(\sum_{m=1}^n \Delta t_m)$

经过离散近似的时间演化算符可以写作
$$U(t)=\prod_{n=1}^N e^{-\frac{i}{\hbar} \hat H(s_n)\Delta t_n}$$

带入绝热演化过程的时变哈密顿量 

$\hat H(s_n)=s_n\hat H_z+(1-s_n)\hat H_x$

可得
$$
\begin{aligned}
    U(t)
    &=\prod_{n=1}^N e^{-\frac{i}{\hbar} s_n\hat H_z\Delta t_n} e^{-\frac{i}{\hbar} (1-s_n)\hat H_x\Delta t_n}\\
    &=\prod_{n=1}^N e^{-i \beta_n H_x} e^{-i \gamma_n H_z}
\end{aligned}
$$

其中

$$
\begin{cases}
\beta_n=(1-s_n)\Delta t_n/\hbar\\
\gamma_n= s_n \Delta t_n/\hbar
\end{cases}
$$
如果记
$$
\begin{cases}
U(\hat H_x,\beta_n)=e^{-i \beta_n H_x}\\
U(\hat H_z,\gamma_n)=e^{-i \gamma_n H_z}
\end{cases}
$$
那么显然
$$
U(t)=U(\hat H_x,\beta_N)U(\hat H_z,\gamma_N)\dots
U(\hat H_x,\beta_1)U(\hat H_z,\gamma_1)
$$
这一近似过程也被称为 *Trotter Splitting* 或者 *Trotterization*。
$$
\hbar \sum_{n=1}^N (\beta_n+\gamma_n)=T
$$

### 渐进优化算法

QAOA中，通常将处态哈密顿量 $\hat H_x$ 记作 $B$，末态目标哈密顿量 $\hat H_z$ 记作 $C$，总离散区间数量 $N$ 记作参数 $p$。

$$B=\hat H_x=\sum_{j=1}^n \sigma^x_j$$

$$ C=\hat H_z=\frac{1}{2}\sum _{\langle i,j \rangle \in E} (1-\sigma ^z_i \sigma ^z_j)$$


系统初态 $|s\rangle=|+\rangle_1|+\rangle_2 \dots|+\rangle_n$


值得注意的是此处的 $\hat H_z$ 与前文相比符号取反并且加入 1/2 系数，其本征值即MaxCut问题中最大切割数量。且初态和末态对应的量子态均为系统最高激发态而非基态。这样做的目的是为了复现14年MIT的QAOA算法论文$^{[2]}$中描述的结构。


当然我们也可以沿用QAA中的描述，并给出损失函数的最小值。Anyway，最终的结果和核心思路是相同的。

得到时间演化算符后，带入
$|s\rangle'=U(t)|s\rangle$
可以得到
$$
|s\rangle'=|\boldsymbol{\gamma},\boldsymbol{\beta}\rangle=U(B,\beta_p)U(C,\gamma_p)\dots
U(B,\beta_1)U(C,\gamma_1)|s\rangle
$$
其中向量
$\boldsymbol{\gamma}=(\gamma_1,\gamma_2\dots\gamma_p)$，
$\boldsymbol{\beta}=(\beta_1,\gamma_2\dots\beta_p)$

理想情况下，适当选取这两组的值可以令系统哈密顿量演化至理想目标函数，此时末态量子态 $|s\rangle'$ 即为问题的解——给出最大切割的映射$f: V\mapsto \{0,1\}$。但是实际上我们并不知道参数的取值，即我们并不知道系统是否演化到了目标末态。我们需要调整参数取值以使目标函数达到极值，这是QAOA与QAA的最大不同之处。

而实际上，参数调整的过程是反馈过程，我们可以使用类似梯度递降（Gradient Descent）的经典优化算法来实现参数的调整，使目标函数的期望取极值。

在参数递降过程中，评估目标函数的值需要对当下量子态进行观测。

对于MaxCut问题，目标函数的期望是
$$
F_p(\boldsymbol{\gamma},\boldsymbol{\beta})=\langle\boldsymbol{\gamma},\boldsymbol{\beta}|C|\boldsymbol{\gamma},\boldsymbol{\beta}\rangle
$$
在优化过程中，量子计算系统的用于快速评估哈密顿量的期望值，从而得到当前参数下目标函数的值。而优化过程即参数调整过程由传统计算机完成，最终对于给定的参数p，优化过程应使目标函数接近最大值
$$
M_p=\max\limits_{\boldsymbol{\gamma},\boldsymbol{\beta}} F_{p}(\boldsymbol{\gamma},\boldsymbol{\beta})
$$

## 总结和展望
QAOA是典型的经典量子混合算法，虽然QAOA来源于量子退火算法，相比传统QAA以及经典SDP算法，其具有更高的效率。由于组合优化问题的重要性和困难度（NP-Complete），学界对于在优化方向应用量子计算寄予厚望。而且相对QAA，QAOA可以解决的问题范围更广，研究已经证明QAOA可以解决某些QAA无法解决的问题$^{[1]}$。但同时QAOA也存在某些局限性，例如有研究指出QAOA在处理某些问题时并不能在误差范围达到最优解。$^{[6]}$。相关领域研究进展很快，进一步的研究和工作仍然需要深入的阅读。

## 附录
QAOA代码-from Github

1. Tensorflow: https://github.com/rigetti/quantumflow-qaoa
2. Scipy: https://github.com/GiggleLiu/QAOA
3. ProjectQ: https://github.com/jiosue/QAOAPython
4. ProjectQ: https://github.com/gavarela/projectq_qaoa
5. PyQuil: https://github.com/sanand13/QAOA (Application-Landscape Search)
6. PyQuil: https://github.com/Tommy-Moffat/Maxcut-QAOA (Application MaxCut)
7. Qiskit: https://github.com/bernovie/QAOA-MaxClique (Application-MaxClique @ Princeton University)

## 参考文献


[1] https://en.wikipedia.org/wiki/Combinatorial_optimization

[2] https://arxiv.org/abs/1411.4028v1

[3] https://arxiv.org/abs/quant-ph/0001106v1

[4] https://arxiv.org/abs/1906.08948

[5] Introduction to Quantum Mechanics, David Griffiths, 3rd Edition

[6] https://arxiv.org/abs/1906.11259v1 
