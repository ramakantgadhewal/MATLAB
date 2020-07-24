# Algorithm

### Vorticity-Stream function Method

##### 适用范围

一般只用于二维不可压缩流动的求解，这里介绍有限差分形式

##### 基本方程及边界条件

1. 涡量对流扩散方程（从速度到涡量）：$\omega_t+u\omega_x+v\omega_y=\frac{1}{\mathrm{Re}}\nabla^2\omega$
   * 给定速度边界：用$\psi_{x_ix_j}$结合虚拟网格点计算
   * 给定速度梯度边界：用相邻内点计算

2. 涡量-流函数关联方程（从涡量到流函数）：$\nabla^2\psi=\omega$
   * 给定速度边界：直接采用虚拟网格点计算
   * 给定速度梯度边界：用相邻内点计算

3. 流函数定义（从流函数到速度）：$\psi_x=-v,\ \psi_y=u$
4. （额外）压力泊松方程：$\nabla^2p=-(u_x^2+2u_yv_x+v_y^2)$
   * 给定法向梯度边界：${\bf{n}}\cdot\nabla p=\frac{1}{\mathrm{Re}}{\bf{n}}\cdot{\Delta}{\bf{u}}$

##### 差分格式

直接采用有限差分格式即可

### SIMPLE Method

##### 适用范围

可用于二维或三维不可压缩流动的求解，这里介绍有限体积形式

##### 基本方程

$u_x+v_y=0$

$u_t+E_x+F_y=0$

$v_t+G_x+H_y=0$

其中，

$E=p+u^2-\frac{1}{\mathrm{Re}}u_x,\ F=uv-\frac{1}{\mathrm{Re}}u_y$

$G=uv-\frac{1}{\mathrm{Re}}v_x,\ H=p+v^2-\frac{1}{\mathrm{Re}}v_y$

##### 离散格式

1. u,v-动量方程（从压力估计值到速度估计值）
   * 给定速度边界：直接给定边界处的预估速度场或通过虚拟网格点保证速度边界条件
2. 压力修正方程（从速度估计值更新压力修正值）
   * 给定法向速度：利用速度修正方程结合相邻内点计算
   * 给定压力边界：通过虚拟网格点保证边界处压力修正量为零
3. 速度修正方程（从压力修正值到速度修正值并更新压力和速度）

###### u-动量方程（下式中的 $u,\ p$ 均实为估计值 $\tilde{u},\ \tilde{p}$）

$\begin{aligned} a_{i+1/2,j} u_{i+1 / 2, j}=\ &a_{i-1 / 2, j} u_{i-1 / 2, j}+a_{i+3 / 2, j} u_{i+3 / 2, j}+a_{i+1 / 2, j-1} u_{i+1 / 2, j-1}+a_{i+1 / 2, j+1} u_{i+1 / 2, j+1}\\ &-\Delta y(p_{i+1, j}-p_{i, j})+\frac{\Delta x \Delta y}{\Delta t} u_{i+1 / 2, j}^{n}+\hat{a}_{i+1 / 2, j}\end{aligned}$

其中，

$a_{i+1 / 2, j}=\frac{\Delta x \Delta y}{\Delta t}+\Delta y\left(\alpha_{i+1, j}^{+}-\alpha_{i, j}^{-}+\frac{2}{\operatorname{Re} \Delta x}\right)+\Delta x\left(\alpha_{i+1 / 2, j+1 / 2}^{+}-\alpha_{i+1 / 2, j-1 / 2}^{-}+\frac{2}{\operatorname{Re} \Delta y}\right)$
$a_{i-1 / 2, j}=\Delta y\left(\alpha_{i, j}^{+}+\frac{1}{\operatorname{Re} \Delta x}\right), a_{i+3 / 2, j}=\Delta y\left(-\alpha_{i+1, j}^{-}+\frac{1}{\operatorname{Re} \Delta x}\right)$
$a_{i+1 / 2, j-1}=\Delta x\left(\alpha_{i+1 / 2, j-1 / 2}^{+}+\frac{1}{\operatorname{Re} \Delta y}\right), a_{i+1 / 2, j+1}=\Delta x\left(-\alpha_{i+1 / 2, j+1 / 2}^{-}+\frac{1}{\operatorname{Re} \Delta y}\right)$
$\hat{a}_{i+1 / 2, j}=\Delta y\left(\gamma_{i+1, j}-\gamma_{i, j}\right)+\Delta x\left(\gamma_{i+1 / 2, j+1 / 2}-\gamma_{i+1 / 2, j-1 / 2}\right)$

这里，

$\alpha^+_{i,j}=\max(0,\bar{u}_{i,j}),\ \alpha^-_{i,j}=\min(0,\bar{u}_{i,j})$

$\alpha^+_{i+1/2,j+1/2}=\max(0,\bar{v}_{i+1/2,j+1/2}),\ \alpha^-_{i+1/2,j+1/2}=\min(0,\bar{v}_{i+1/2,j+1/2})$

$\bar{u}_{i,j}=\left(u_{i+1/2,j}^n+u_{i-1/2,j}^n\right)/2,\ \bar{v}_{i+1/2,j+1/2}=\left(v_{i,j+1/2}^n+v_{i+1,j+1/2}^n\right)/2$

$\gamma$ 表达式过于复杂，此处略去（目前只实现了 $\gamma=0$ 的一阶迎风格式）

###### v-动量方程（下式中的 $v,\ p$ 均实为估计值 $\tilde{v},\ \tilde{p}$）

$\begin{aligned} b_{i, j+1 / 2} v_{i, j+1 / 2}=\ & b_{i-1, j+1 / 2} v_{i-1, j+1 / 2}+b_{i+1, j+1 / 2} v_{i+1, j+1 / 2}+b_{i, j-1 / 2} v_{i, j-1 / 2}+b_{i, j+3 / 2} v_{i, j+3 / 2} \\ &-\Delta x\left(p_{i, j+1}-p_{i, j}\right)+\frac{\Delta x \Delta y}{\Delta t} v_{i, j+1 / 2}^n+\hat{b}_{i, j+1 / 2} \end{aligned}$

其中，

$b_{i, j+1/2}=\frac{\Delta x \Delta y}{\Delta t}+\Delta y\left(\beta_{i+1/2, j+1/2}^{+}-\beta_{i-1/2, j+1/2}^{-}+\frac{2}{\operatorname{Re} \Delta y}\right)+\Delta x\left(\beta_{i, j+1}^{+}-\beta_{i, j}^{-}+\frac{2}{\operatorname{Re} \Delta x}\right)$
$b_{i-1, j+1/2}=\Delta y\left(\beta_{i-1/2, j+1/2}^{+}+\frac{1}{\operatorname{Re} \Delta x}\right), b_{i+1, j+1/2}=\Delta y\left(-\beta_{i+1/2, j+1/2}^{-}+\frac{1}{\operatorname{Re} \Delta x}\right)$
$b_{i, j-1/2}=\Delta x\left(\beta_{i, j}^{+}+\frac{1}{\operatorname{Re} \Delta y}\right), b_{i, j+3/2}=\Delta x\left(-\beta_{i, j+1}^{-}+\frac{1}{\operatorname{Re} \Delta y}\right)$
$\hat{b}_{i, j+1/2}=\Delta y\left(\delta_{i+1/2, j+1/2}-\delta_{i+1/2, j-1/2}\right) + \Delta x\left(\delta_{i+1, j}-\delta_{i, j}\right)$

这里，

$\beta^+_{i,j}=\max(0,\bar{v}_{i,j}),\ \beta^-_{i,j}=\min(0,\bar{v}_{i,j})$

$\beta^+_{i+1/2,j+1/2}=\max(0,\bar{u}_{i+1/2,j+1/2}),\ \beta^-_{i+1/2,j+1/2}=\min(0,\bar{u}_{i+1/2,j+1/2})$

$\bar{v}_{i,j}=\left(v_{i,j+1/2}^n+v_{i,j-1/2}^n\right)/2,\ \bar{u}_{i+1/2,j+1/2}=\left(u_{i+1/2,j}^n+u_{i+1/2,j+1}^n\right)/2$

$\delta$ 表达式过于复杂，此处略去（目前只实现了 $\delta=0$ 的一阶迎风格式）

###### 压力修正方程：$p_{i,j}=\tilde{p}_{i,j}+p_{i,j}^\prime$

$c_{i,j}p_{i,j}^\prime=c_{i-1,j}p_{i-1,j}^\prime+c_{i+1,j}p_{i+1,j}^\prime+c_{i,j-1}p_{i,j-1}^\prime+c_{i,j+1}p_{i,j+1}^\prime+\hat{c}_{i,j}$

其中，

$c_{i,j}=(\Delta y)^2\left(\frac{1}{a_{i+1/2,j}}+\frac{1}{a_{i-1/2,j}}\right)+(\Delta x)^2\left(\frac{1}{b_{i,j+1/2}}+\frac{1}{b_{i,j-1/2}}\right)$

$c_{i-1,j}={(\Delta y)^2}/{a_{i-1/2,j}},\ c_{i+1,j}={(\Delta y)^2}/{a_{i+1/2,j}}$

$c_{i,j-1}={(\Delta x)^2}/{b_{i,j-1/2}},\ c_{i,j+1}={(\Delta x)^2}/{b_{i,j+1/2}}$

$\hat{c}_{i,j}=\Delta y \left(\tilde{u}_{i-1/2,j}-\tilde{u}_{i+1/2,j}\right)+\Delta x \left(\tilde{v}_{i,j-1/2}-\tilde{v}_{i,j+1/2}\right)$

###### 速度修正方程：$u_{i+1/2,j}=\tilde{u}_{i+1/2,j}+u_{i+1/2,j}^\prime,\ v_{i,j+1/2}=\tilde{v}_{i,j+1/2}+v_{i,j+1/2}^\prime$

$u_{i+1/2,j}^\prime=-\frac{\Delta y}{a_{i+1/2,j}}(p_{i+1,j}^\prime-p_{i,j}^\prime),\ v_{i,j+1/2}^\prime=-\frac{\Delta x}{b_{i,j+1/2}}(p_{i,j+1}^\prime-p_{i,j}^\prime)$

