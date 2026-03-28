

**版本**：3.0 (最终完整版)
**适用对象**：托卡马克物理初学者、科学计算开发者
**核心框架**：FyTok / TORAX / FreeGS
**关键词**：Grad-Shafranov方程，1.5D输运，有限体积法(FVM)，Newton-Raphson，Picard迭代，多物理场耦合

---

## 目录

1.  [第一章 磁流体平衡基础：Grad-Shafranov 方程](#第一章-磁流体平衡基础grad-shafranov-方程)
2.  [第二章 平衡方程的数值求解 (FreeGS)](#第二章-平衡方程的数值求解-freegs)
3.  [第三章 多组分含时演化输运方程组](#第三章-多组分含时演化输运方程组)
4.  [第四章 输运方程的数值离散与求解算法体系（详解）](#第四章-输运方程的数值离散与求解算法体系详解)
5.  [第五章 输运与平衡的闭环耦合实现](#第五章-输运与平衡的闭环耦合实现)
6.  [第六章 总结](#第六章-总结)

---

## 第一章 磁流体平衡基础：Grad-Shafranov 方程

在轴对称托卡马克装置中，等离子体宏观上处于力学平衡状态（$\nabla p = \mathbf{J} \times \mathbf{B}$）。描述这种二维平衡状态的核心方程是 **Grad-Shafranov (G-S) 方程**。

### 1.1 物理方程
在柱坐标系 $(R, \phi, Z)$ 下，G-S 方程建立了极向磁通量 $\psi$ 与等离子体电流源项之间的关系：

$$
\Delta^* \psi = -\mu_0 R^2 \frac{dp(\psi)}{d\psi} - F(\psi)\frac{dF(\psi)}{d\psi}
$$

**符号定义**：
* **$\psi(R,Z)$**：极向磁通量（Poloidal Flux），待求解的二维标量场。磁面由 $\psi = \text{const}$ 定义。
* **$\Delta^*$**：Grad-Shafranov 椭圆算子：
    $$
    \Delta^* = R \frac{\partial}{\partial R} \left( \frac{1}{R} \frac{\partial}{\partial R} \right) + \frac{\partial^2}{\partial Z^2}
    $$
* **$p(\psi)$**：等离子体热压强剖面（由输运方程计算提供）。
* **$F(\psi) = R B_\phi$**：极向电流函数，与环向磁场有关。
* **$J_\phi(R, Z)$**：环向电流密度（方程右端项），是平衡求解的“源”。

### 1.2 边界条件类型
* **固定边界 (Fixed Boundary)**：给定最后闭合磁面 (LCFS) 的几何形状，仅求解内部 $\psi$ 分布。适用于 CHEASE 代码。
* **自由边界 (Free Boundary)**：同时求解内部和外部真空区。边界条件由外部极向场线圈 (PF Coils) 电流决定，等离子体形状和位置是自洽演化的。适用于 FreeGS 代码。

---

## 第二章 平衡方程的数值求解 (FreeGS)

FreeGS 采用 **Picard 迭代法** 结合 **格林函数法 (Green's Function)** 求解自由边界问题。

### 2.1 求解算法流程
1.  **初始化**：设定网格 $(R_{grid}, Z_{grid})$，初始化 $\psi^{(0)}$（通常设为真空场）。
2.  **迭代循环 ($k=0, 1, 2...$)**：
    * **A. 更新源项 $J_\phi$**：
        根据当前 $\psi^{(k)}$ 找到磁轴和边界，归一化得到 $\psi_{norm}$。代入输运求解器给出的 $p'(\psi)$ 和 $FF'(\psi)$ 计算空间各点的 $J_\phi(R, Z)$。
    * **B. 更新边界 (格林函数)**：
        计算等离子体电流 $J_\phi$ 在边界产生的磁通 $\psi_{plasma}$，叠加外部线圈产生的 $\psi_{ext}$，作为边界条件。
    * **C. 求解 PDE**：
        在已知源项和边界条件下，求解线性椭圆方程 $\Delta^* \psi^{(k+1)} = -\mu_0 R J_\phi^{(k)}$ 得到新的 $\psi^{(k+1)}$。
    * **D. 收敛检查**：
        若 $||\psi^{(k+1)} - \psi^{(k)}|| < \epsilon$，停止迭代。
3.  **后处理**：
    计算一维输运所需的磁面平均几何量：$V(\rho)$（体积），$S(\rho)$（面积），$\langle |\nabla \rho|^2 \rangle$（几何因子）等。

---

## 第三章 多组分含时演化输运方程组

输运求解器（如 TORAX）负责演化等离子体的热力学状态。这是一组 **1.5D 抛物型偏微分方程**。

### 3.1 统一的输运方程形式
为了推导清晰，我们将所有输运方程统一写为关于物理量 $u(\rho, t)$ 的守恒律形式：

$$
\frac{\partial u}{\partial t} = \frac{1}{V'} \frac{\partial}{\partial \rho} \left[ V' \left( D \langle |\nabla \rho|^2 \rangle \frac{\partial u}{\partial \rho} - v_{conv} \langle |\nabla \rho| \rangle u \right) \right] + S
$$

**符号定义**：
* **$u$**：待求解的物理量（如电子密度 $n_e$，电子温度 $T_e$，极向磁通 $\psi$）。
* **$\rho$**：归一化环向磁通半径 $\hat{\rho} \in [0, 1]$。
* **$V' = \frac{\partial V}{\partial \rho}$**：几何雅可比行列式（体积随半径的变化率）。
* **$D$**：扩散系数（Diffusivity，如热扩散率 $\chi$）。
* **$v_{conv}$**：对流速度（Convection velocity，如粒子箍缩）。
* **$S$**：源项（Source/Sink，如加热功率、欧姆加热）。

### 3.2 物理含义
* **非线性**：系数 $D, v_{conv}, S$ 通常强依赖于 $u$ 及其梯度 $\nabla u$。例如，湍流热扩散系数 $\chi \propto T^{3/2} \cdot \nabla T$。
* **耦合**：不同方程通过源项和系数相互耦合（例如 $n_e$ 影响 $T_e$ 的源项，$T_e$ 影响 $\psi$ 的电阻率）。

---
这是一个非常扎实的需求。为了让您彻底理解每种求解器的底层运作机制，我将以**统一的数学框架**（从连续方程到代数方程组），对第四章中列出的所有求解方法进行详细的推导。

我们将统一使用以下符号：

- **物理量**：$u$（待求解向量，如 $[n_e, T_e]$）
    
- **时间步**：$n$（当前已知时刻），$n+1$（下一时刻）
    
- **迭代步**：$k$（同一时间步内的第 $k$ 次猜测）
    
- 通用输运方程（算子形式）：
    
    $$\frac{\partial u}{\partial t} = \nabla \cdot (D(u) \nabla u) + S(u)$$
    

---
### 3.3**边界条件与台基物理接口**

输运方程是二阶偏微分方程，需要两个边界条件。

- **轴心处 ($\rho=0$)**：必须满足自然边界条件，即梯度为零（$\nabla T = 0, \nabla n = 0$），以保证物理量的连续性和光滑性。
- **边缘处 ($\rho=\rho_{sep}$ 或 $\rho_{ped}$) **：这是最棘手的地方。
    - **L 模**：可以直接给定最外层闭合磁面（LCFS）处的温度和密度值（通常很低）。
    - **H 模**：输运求解器通常只求解到台基顶（$\rho \approx 0.9 - 0.95$）。台基区的温度和密度通常由 **EPED 模型** 或 **经验定标率** 给出，作为芯部输运方程的狄利克雷（Dirichlet）边界条件。
    - **耦合逻辑**：TORAX 需要留出接口，允许根据当前的 $\beta_p$ 和整形参数，动态调用 EPED 计算台基高度，从而更新边界条件。
# 第四章 输运方程求解算法体系详解
##### 基本知识补充

######  步骤 1：写出积分形式的守恒方程

**原始微分方程**：

$$
\frac{\partial n}{\partial t} = -\frac{1}{V'} \frac{\partial}{\partial \psi} (V' \Gamma) + S
$$

两边同乘 $V'$（即 $dV/d\psi$）：

$$
V' \frac{\partial n}{\partial t} = - \frac{\partial}{\partial \psi} (V' \Gamma) + S V'
$$

> 💡 注意：$V' d\psi = dV$，因此左边为 $\frac{\partial n}{\partial t} dV$。

现在，对第 $i$ 个控制体（磁面范围 $\psi \in [\psi_{i-1/2}, \psi_{i+1/2}]$）积分：

$$
\int_{\psi_{i-1/2}}^{\psi_{i+1/2}} V' \frac{\partial n}{\partial t} \, d\psi 
= - \int_{\psi_{i-1/2}}^{\psi_{i+1/2}} \frac{\partial}{\partial \psi} (V' \Gamma) \, d\psi 
+ \int_{\psi_{i-1/2}}^{\psi_{i+1/2}} S V' \, d\psi
$$

---

###### 步骤 2：简化每一项

  **左边：时间变化率 × 体积**

$$
\int_{\psi_{i-1/2}}^{\psi_{i+1/2}} V' \frac{\partial n}{\partial t} \, d\psi 
= \frac{\partial}{\partial t} \int_{\psi_{i-1/2}}^{\psi_{i+1/2}} n V' \, d\psi 
\approx \frac{\partial n_i}{\partial t} \int_{\psi_{i-1/2}}^{\psi_{i+1/2}} V' \, d\psi 
= \frac{\partial n_i}{\partial t} \Delta V_i
$$

> 📝 **说明**：假设单元内密度均匀（$n \approx n_i$），可将$\n_i$ 提出积分号。

 **中间项：通量的净流出（用微积分基本定理）**

$$
- \int_{\psi_{i-1/2}}^{\psi_{i+1/2}} \frac{\partial}{\partial \psi} (V' \Gamma) \, d\psi 
= - \left[ (V' \Gamma) \Big|_{\psi_{i+1/2}} - (V' \Gamma) \Big|_{\psi_{i-1/2}} \right] 
= - \left( A_{i+1/2} \Gamma_{i+1/2} - A_{i-1/2} \Gamma_{i-1/2} \right)
$$

其中：

- $A_{i\pm1/2} = V'(\psi_{i\pm1/2})$：界面面积（磁面“周长”因子）；
- $\Gamma_{i\pm1/2}$：界面处的粒子通量。

> 📝 **说明**：**这是 FVM 的核心——将空间导数转化为“界面通量差”。**

###### ✅ 右边：源项积分

$$
\int_{\psi_{i-1/2}}^{\psi_{i+1/2}} S V' \, d\psi 
\approx S_i \int_{\psi_{i-1/2}}^{\psi_{i+1/2}} V' \, d\psi 
= S_i \Delta V_i
$$

> 📝 **说明**：同样假设源项在单元内均匀（$S \approx S_i$）。

---

######  步骤 3：合并得到积分守恒方程
将径向 $\rho \in [0, 1]$ 划分为 $N$ 个网格。为了数值稳定，使用 **参差网格 (Staggered Grid)**：
* **网格中心 ($i$)**：存储状态变量 $u_i$ ($n, T, \psi$)。
* **网格界面 ($i \pm 1/2$)**：存储通量 $\Gamma_{i \pm 1/2}$。
将三部分代入，得到第 $i$ 个单元的**连续时间守恒方程**：

$$
\frac{\partial n_i}{\partial t} \Delta V_i = - \left( A_{i+1/2} \Gamma_{i+1/2} - A_{i-1/2} \Gamma_{i-1/2} \right) + S_i \Delta V_i
$$

✅ **这就是 FVM 的核心离散方程**（空间已离散，时间仍连续）。
### 4.1 空间离散：有限体积法 (FVM)

将径向 $\rho \in [0, 1]$ 划分为 $N$ 个网格。为了数值稳定，使用 **参差网格 (Staggered Grid)**：
* **网格中心 ($i$)**：存储状态变量 $u_i$ ($n, T, \psi$)。
* **网格界面 ($i \pm 1/2$)**：存储通量 $\Gamma_{i \pm 1/2}$。

**半离散方程推导**：
在第 $i$ 个网格单元上对通用方程积分：

$$
\frac{d u_i}{dt} \Delta V_i = \Phi_{i+1/2} - \Phi_{i-1/2} + S_i \Delta V_i
$$

其中 $\Phi$ 是通过界面的总通量（包含扩散和对流）：
$$
\Phi_{i+1/2} \approx \left[ V' \left( D \langle |\nabla \rho|^2 \rangle \frac{u_{i+1} - u_i}{\Delta \rho} - v_{conv} \langle |\nabla \rho| \rangle u_{interface} \right) \right]_{i+1/2}
$$

### 步骤 4.2：引入时间离散（Theta 方法）

现在对时间导数进行离散。使用 Theta 格式（$\theta \in [0,1]$）：

$$
\frac{n_i^{n+1} - n_i^n}{\Delta t} \Delta V_i = 
- \left[ \theta \left( A \Gamma^{n+1} \right)_{\text{net}} + (1 - \theta) \left( A \Gamma^{n} \right)_{\text{net}} \right] 
+ \left( \theta S_i^{n+1} + (1 - \theta) S_i^n \right) \Delta V_i
$$

其中净通量：
$$
\left( A \Gamma \right)_{\text{net}} = A_{i+1/2} \Gamma_{i+1/2} - A_{i-1/2} \Gamma_{i-1/2}
$$

展开后即为您看到的完整方程：

$$
\frac{n_i^{n+1} - n_i^n}{\Delta t} \Delta V_i = 
- \theta \left( \Gamma_{i+1/2}^{n+1} A_{i+1/2} - \Gamma_{i-1/2}^{n+1} A_{i-1/2} \right) 
- (1 - \theta) \left( \Gamma_{i+1/2}^{n} A_{i+1/2} - \Gamma_{i-1/2}^{n} A_{i-1/2} \right) 
+ \left( \theta S_i^{n+1} + (1 - \theta) S_i^n \right) \Delta V_i
$$
其中 $A_{i\pm1/2} = V'(\psi_{i\pm1/2})$ 是界面面积（磁面周长相关）。
直观看如下：

对于任意一个网格单元 $i$，体积为 $\Delta V_i$，物理量的守恒方程为：

$$ \frac{n_i^{new} - n_i^{old}}{\Delta t} \Delta V_i = - (\text{Flux}_{out} A_{out} - \text{Flux}_{in} A_{in}) + S_i \Delta V_i $$

我们将所有 $n^{new}$ 项移到等式左边（构建矩阵 $A$），已知项移到右边（构建向量 $b$）。

$$ \underbrace{\frac{\Delta V_i}{\Delta t} n_i^{new}}_{\text{时间项}} + \underbrace{(\text{Flux项})}_{\text{输运项}} = \underbrace{\frac{\Delta V_i}{\Delta t} n_i^{old} + S_i \Delta V_i}_{\text{历史项 + 源项}} $$
---

 **物理意义再解读**

| 项                                               | 物理含义                 |
| ----------------------------------------------- | -------------------- |
| $\frac{n_i^{n+1} - n_i^n}{\Delta t} \Delta V_i$ | 单元 $i$ 内粒子总数的变化量     |
| $- (\text{右通量} - \text{左通量})$                   | **净流出**：右界面流出减去左界面流入 |
| $S_i \Delta V_i$                                | 单元内源项产生的总粒子数         |

> ✅ **负号关键**：如果右通量 > 左通量，净流出为正，导致单元内粒子减少，符合物理直觉。
> ⚠️ **关键**：通量 $\Gamma$ 在 $n+1$ 时刻是非线性的（因 $D = D(n^{n+1})$, $V = V(n^{n+1})$），需迭代求解。
我们将严格遵循 **Theta 时间离散 + FVM 空间离散** 框架，逐步展示如何从物理通量 $\Gamma$ 最终得到三对角线性系统。整个过程适用于 **Picard 迭代**（固定系数）或 **牛顿法**（线性化后）。

---

##### # 1. 起点：半离散守恒方程（已时间离散）

对第 $i$ 个单元，Theta 方法给出：

$$
\frac{n_i^{n+1} - n_i^n}{\Delta t} \Delta V_i = 
- \theta \left( \Gamma_{i+1/2}^{n+1} A_{i+1/2} - \Gamma_{i-1/2}^{n+1} A_{i-1/2} \right) 
- (1 - \theta) \left( \Gamma_{i+1/2}^{n} A_{i+1/2} - \Gamma_{i-1/2}^{n} A_{i-1/2} \right) 
+ \left( \theta S_i^{n+1} + (1 - \theta) S_i^n \right) \Delta V_i
$$

> ✅ **目标**：将所有含 $n^{n+1}$ 的项移到左边，其余移到右边。

---
**磁轴处的数值处理**

在网格中心 $i=1$（最靠近磁轴的点），通量 $\Phi_{1/2}$（即 $\rho=0$ 处）必须强制为 0。但在计算径向电场或扩散系数时，经常遇到 $\frac{1}{\rho}\frac{\partial}{\partial \rho}$ 形式的算子。
- **物理约束**：必须利用洛必达法则（L'Hôpital's rule）。例如，在 $\rho \to 0$ 时，体积 $V \propto \rho^2$，面积 $S \propto \rho$，因此 $V' \propto \rho$。
- **代码实现**：在计算几何系数数组时，必须显式定义 `g[0]` 的极限值，而不是让计算机去算 $0/0$。例如，$\langle |\nabla \rho|^2 \rangle$ 在轴心处通常是有限值，需通过抛物线拟合邻近点来外推。
##### 2. 通量离散：将 $\Gamma^{n+1}$ 表示为 $n^{n+1}$ 的线性函数

#### **(a) 扩散项离散（中心差分）**
$$
\left( -D \frac{\partial n}{\partial \psi} \right)_{i+1/2}^{n+1} 
\approx -D_{i+1/2}^{(k-1)} \frac{n_{i+1}^{(k)} - n_i^{(k)}}{\Delta \psi_{i+1/2}}
$$
- $D_{i+1/2}^{(k-1)}$：在 Picard 第 $k$ 步，用上一次迭代的 $n^{(k-1)}$ 计算 $D$；
- $n_{i}^{(k)}$：当前迭代待求解的未知量。

#### **(b) 对流项离散（迎风格式）**
假设 $V_{i+1/2} > 0$（向外流）：
$$
\left( V n \right)_{i+1/2}^{n+1} \approx V_{i+1/2}^{(k-1)} n_i^{(k)}
$$

#### **(c) 合并通量**
$$
\Gamma_{i+1/2}^{n+1} = 
\underbrace{\left( \frac{D_{i+1/2}^{(k-1)}}{\Delta \psi_{i+1/2}} + V_{i+1/2}^{(k-1)} \right)}_{c_1^{(k-1)}} n_i^{(k)} 
+ \underbrace{\left( -\frac{D_{i+1/2}^{(k-1)}}{\Delta \psi_{i+1/2}} \right)}_{c_2^{(k-1)}} n_{i+1}^{(k)}
$$

> 🔑 **关键**：在 Picard 迭代中，$c_1, c_2$ 被视为**已知常数**（由上一轮迭代确定），因此 $\Gamma^{n+1}$ 是 $n^{(k)}$ 的**线性函数**。

---

### 3. 代入守恒方程并移项

将 $\Gamma^{n+1}$ 表达式代入守恒方程，并将所有 $n^{(k)}$ 项移到左边：

#### 左边（未知量）：
$$
\frac{\Delta V_i}{\Delta t} n_i^{(k)} 
+ \theta \left[ \left(c_1^{(k-1)} n_i^{(k)} + c_2^{(k-1)} n_{i+1}^{(k)} \right) A_{i+1/2} 
- \left(c_1^{\text{left}} n_{i-1}^{(k)} + c_2^{\text{left}} n_i^{(k)} \right) A_{i-1/2} \right]
$$

> 💡 其中 $c_1^{\text{left}}, c_2^{\text{left}}$ 是左界面 $i-1/2$ 的系数。

展开后，按 $n_{i-1}^{(k)}, n_i^{(k)}, n_{i+1}^{(k)}$ 合并：

- **$n_{i-1}^{(k)}$ 的系数**（下对角线）：
  $$
  a_{i,i-1} = - \theta \cdot c_1^{\text{left}} \cdot A_{i-1/2}
  $$

- **$n_i^{(k)}$ 的系数**（主对角线）：
  $$
  a_{i,i} = \frac{\Delta V_i}{\Delta t} + \theta \left( c_1^{(k-1)} A_{i+1/2} - c_2^{\text{left}} A_{i-1/2} \right)
  $$

- **$n_{i+1}^{(k)}$ 的系数**（上对角线）：
  $$
  a_{i,i+1} = \theta \cdot c_2^{(k-1)} \cdot A_{i+1/2}
  $$

#### 右边（已知量）：
$$
b_i = \frac{\Delta V_i}{\Delta t} n_i^n 
- (1 - \theta) \left( \Gamma_{i+1/2}^{n} A_{i+1/2} - \Gamma_{i-1/2}^{n} A_{i-1/2} \right) 
+ \left( \theta S_i^{(k-1)} + (1 - \theta) S_i^n \right) \Delta V_i
$$

> 💡 $S_i^{(k-1)}$ 同样由上一轮迭代的 $n^{(k-1)}$ 计算。

---

### 4. 矩阵组装：基于面的贡献分配（Face-based Assembly）

上述系数可通过遍历**每个界面**高效计算：

``` python
	# 伪代码：基于面的矩阵组装
	for i in range(N-1):  # 遍历界面 i+1/2
	    # 1. 计算界面系数 (c1, c2)
	    D_face = 0.5 * (D[i] + D[i+1])      # 算术平均
	    dr = psi[i+1] - psi[i]
	    c1 = D_face / dr + V_face
	    c2 = -D_face / dr
	    
	    # 2. 对左侧单元 i 的贡献
	    main_diag[i]   += theta * c1 * A_face
	    upper_diag[i]  += theta * c2 * A_face  # 影响 n_{i+1}
	    
	    # 3. 对右侧单元 i+1 的贡献
	    lower_diag[i]  += -theta * c1 * A_face  # 影响 n_i
	    main_diag[i+1] += -theta * c2 * A_face
```

### 4. 通量离散化（Flux Discretization）

以界面 $i + 1/2$ 为例，通量：

$$
\Gamma_{i+1/2} = -D_{i+1/2} \frac{n_{i+1} - n_i}{\Delta \psi_{i+1/2}} + V_{i+1/2} n_{\text{upwind}}
$$

- **扩散项**：中心差分；
- **对流项**：迎风格式（$V > 0$ 时 $n_{\text{upwind}} = n_i$）；
- **输运系数**：通常取算术平均 $D_{i+1/2} = (D_i + D_{i+1})/2$。

于是 $\Gamma_{i+1/2}$ 可表示为 $n_i$ 和 $n_{i+1}$ 的线性组合：

$$
\Gamma_{i+1/2} = c_1 n_i + c_2 n_{i+1}
$$

其中：
$$
c_1 = \frac{D_{i+1/2}}{\Delta \psi} + V_{i+1/2}, \quad c_2 = -\frac{D_{i+1/2}}{\Delta \psi} \quad \text{（若 } V > 0 \text{）}
$$

### 2. 矩阵系数 ($A \cdot x = b$) 的构成

方程最终被整理为三对角形式：
$$ a_{i,i-1} n_{i-1} + a_{i,i} n_i + a_{i,i+1} n_{i+1} = b_i $$

#### **A. 主对角线 (Main Diagonal, $a_{i,i}$)**
这是 $n_i$ 的系数，它代表了“**本单元的惯性 + 流出损失**”。
*   **时间惯性 (Time Inertia)**: $$self.dV[i] / dt$$
    *   代码对应: $main_diag[i] += self.dV[i] / dt$
    *   物理意义: 这一项越大（$\Delta t$ 越小），$n_i$ 越倾向于保持原来的值，矩阵越稳定（对角占优）。
*   **扩散流出 (Diffusion Loss)**: $\frac{D}{\Delta r} A$
    *   代码对应: `main_diag[i] += (coeff_i) * A_face` (其中 `coeff_i` 为正)
    *   物理意义: 浓度越高，向外扩散越快，这是一种“损失”机制，所以系数为正（移项到左边）。
*   **对流流出 (Convection Loss)**: $V A$ (如果 $V > 0$)
    *   代码对应: `main_diag[i] += (conv_i) * A_face`
    *   物理意义: 如果速度向外，粒子被带走。

#### **B. 非对角线 (Off-Diagonals, $a_{i,i-1}$ 和 $a_{i,i+1}$)**
这是邻居单元 $n_{i-1}$ 和 $n_{i+1}$ 的系数，代表“**邻居对本单元的影响（流入/耦合）**”。
*   **扩散流入 (Diffusion Gain)**: $-\frac{D}{\Delta r} A$
    *   代码对应: `upper_diag[i]` 或 `lower_diag[i]`
    *   物理意义: 邻居浓度越高，扩散进入本单元的粒子越多。因为在等式左边，所以通常为**负值**。
*   **对流流入 (Convection Gain)**: $-V A$ (如果 $V$ 指向本单元)
    *   代码对应: Upwind 格式处理。
    *   物理意义: 上游单元的粒子随流体进入本单元。

#### **C. 右端项 (RHS, $b_i$)**
这是方程的驱动源。
*   **历史记忆 (History)**: $n_i^{old} \frac{\Delta V_i}{\Delta t}$
    *   代码对应: `b[i] += self.n[i] * self.dV[i] / dt`
    *   物理意义: 如果没有输运和源，下一时刻密度应该等于上一时刻。
*   **外部源 (Source)**: $S_i \Delta V_i$
    *   代码对应: `b[i] += self.S[i] * self.dV[i]`
    *   物理意义: 外部注入的粒子（如中性束注入、电离等）。
#### D. 对角物理图像如下
![[00-输运方程中的多物理场耦合关系.png]]

---
**总结对应关系：**
*   **$\Delta t$ 变小** $\rightarrow$ `dV/dt` 变大 $\rightarrow$ 主对角线 $a_{i,i}$ 变大 $\rightarrow$ 矩阵更稳定，变化更慢。
*   **$D$ 变大** $\rightarrow$ `D/dr` 变大 $\rightarrow$ 主对角线和非对角线绝对值都变大 $\rightarrow$ 空间耦合更强，平滑作用更强。
*   **$S$ 变大** $\rightarrow$ RHS $b_i$ 变大 $\rightarrow$ 密度 $n$ 整体上升。
###### *关键的非对角耦合项**

雅可比矩阵的复杂性主要来自以下物理机制的强耦合，必须确保 JAX 能够正确捕捉这些导数：
1. **自举电流 (Bootstrap Current)**：$J_{BS} \propto \frac{\partial n}{\partial \rho} + \frac{\partial T}{\partial \rho}$。这意味着磁扩散方程（求解 $\psi$）对温度梯度和密度梯度极其敏感。这是最强的刚性来源之一。
2. **能量交换项 ($Q_{ei}$)**：$Q_{ei} \propto n^2 (T_e - T_i) / T_e^{1.5}$。当 $T_e$ 和 $T_i$ 接近时，该项对温度的导数非常大，直接耦合了电子和离子通道。
3. **磁几何反馈**：$D \propto q^2$（安全因子）。由于 $q \sim \partial \psi / \partial \rho$，这意味着扩散系数本身依赖于磁通的梯度。
## 4.3 基础求解策略：线性化与解耦

### 4.3.1 算子分裂法 (Operator Splitting)

**核心思想：**
将复杂的多物理场耦合方程 $\frac{\partial u}{\partial t} = \mathcal{L}_{total}(u)$ 拆解为多个简单的子过程，按顺序依次求解。
数学推导：
假设系统包含粒子输运 $\mathcal{L}_n$ 和能量输运 $\mathcal{L}_T$。方程写作：
$$\frac{\partial \mathbf{u}}{\partial t} = \mathcal{L}_n(\mathbf{u}) + \mathcal{L}_T(\mathbf{u})$$
我们将从 $t^n$ 到 $t^{n+1}$ 的求解过程拆分为两步（Godunov Splitting）：
**Step 1: 求解密度方程（冻结温度 $T$）**
构造中间状态 $u^*$。此时只演化密度算子：
$$\frac{n^* - n^n}{\Delta t} = \nabla \cdot (D_n(T^n) \nabla n^*) + S_n$$
这是一个关于 $n^*$ 的线性方程，求解得到 $n^*$（$T$ 保持为 $T^n$）。
**Step 2: 求解温度方程（使用最新的密度 $n^*$）**
从中间状态 $u^*$ 演化到 $u^{n+1}$。此时只演化能量算子：
$$\frac{T^{n+1} - T^n}{\Delta t} = \frac{1}{n^*} \nabla \cdot (n^* \chi(T^{n+1}) \nabla T^{n+1}) + \dots$$
通常这里会对 $\chi$ 进行线性化处理（取 $\chi(T^n)$）。
**算法流程**：
1. 构建密度扩散矩阵 $\mathbf{A}_n$。
2. 解线性方程 $\mathbf{A}_n \mathbf{n}^* = \mathbf{b}_n$，得中间密度。
3. 构建温度扩散矩阵 $\mathbf{A}_T$（利用 $n^*$）。
4. 解线性方程 $\mathbf{A}_T \mathbf{T}^{n+1} = \mathbf{b}_T$，得最终温度。

---
### 4.3.2 线性 Theta 方法 (Linear Theta Method)

核心思想：
不对物理过程进行拆分，而是联立求解。但为了避免迭代，将所有非线性系数（$D, \chi, S$）**“冻结”**在上一时刻 $t^n$。
数学推导：
离散化方程（全隐式 $\theta=1$）：
$$\frac{u^{n+1} - u^n}{\Delta t} = \nabla \cdot (D(u^{?}) \nabla u^{n+1}) + S(u^{?})$$
线性化假设：
令 $D(u^{?}) \approx D(u^n)$， $S(u^{?}) \approx S(u^n) + S'(u^n)(u^{n+1}-u^n)$。
离散化与矩阵构建：
对于空间点 $i$，扩散项展开为三对角形式：
$$\nabla \cdot (D(u^n) \nabla u^{n+1}) \approx \alpha_i(u^n) u_{i-1}^{n+1} - \beta_i(u^n) u_{i}^{n+1} + \gamma_i(u^n) u_{i+1}^{n+1}$$
代入并整理，得到标准线性方程组 $\mathbf{A} \mathbf{u}^{n+1} = \mathbf{b}$：
$$\left[ \frac{1}{\Delta t}\mathbf{I} - \mathbf{L}(u^n) \right] \mathbf{u}^{n+1} = \frac{1}{\Delta t}\mathbf{u}^n + S(u^n)$$
其中 $\mathbf{L}(u^n)$ 是基于上一时刻系数构建的拉普拉斯矩阵。
**算法流程**：
1. 利用 $u^n$ 计算全场 $D, \chi, v_{conv}$。
2. 组装稀疏矩阵 $\mathbf{A}$ 和右端向量 $\mathbf{b}$。
3. 调用线性求解器（如 `scipy.sparse.linalg.spsolve` 或 Thomas 算法）一次性解出 $u^{n+1}。
---
## 4.4 进阶求解策略：非线性迭代

### 4.4.1 预测-校正法 (Predictor-Corrector)

核心思想：
通过两步计算来更新非线性系数，从而减少线性化带来的滞后误差。
**数学推导**：
Step 1: 预测 (Predictor)
使用线性 Theta 方法求得一个估计值 $u^*$：
$$\frac{u^* - u^n}{\Delta t} = \nabla \cdot (D(u^n) \nabla u^*) + S(u^n)$$

解得 $u^*$。

Step 2: 更新系数

利用预测值 $u^*$ 计算 $t^{n+1/2}$ 或 $t^{n+1}$ 时刻的物理系数：

$$D^* = D(u^*) \quad (\text{或者 } D(0.5(u^n + u^*)))$$

Step 3: 校正 (Corrector)
使用新系数 $D^*$ 重新构建方程求解最终值 $u^{n+1}$：
$$\frac{u^{n+1} - u^n}{\Delta t} = \nabla \cdot (D^* \nabla u^{n+1}) + S(u^*)$$

**算法流程**：
1. **Predict**: 用 $u^n$ 的系数解线性方程 $\to u^*$。
2. **Evaluate**: 计算 $u^*$ 下的新系数 $D^*, S^*$。
3. **Correct**: 用 $D^*, S^*$ 再次构建并解线性方程 $\to u^{n+1}$
---
### 4.4.2 皮卡迭代法 (Picard Iteration)
核心思想：
构造不动点迭代序列 $u^{(k+1)} = G(u^{(k)})$。本质上是反复执行“更新系数 -> 解线性方程”的过程，直到系数和解不再变化。
数学推导：
将原方程写为：
$$\mathcal{A}(u) \cdot u = b(u)$$

其中算子 $\mathcal{A}(u)$包含了依赖于 $u$ 的扩散系数，$\mathcal{A}(u) \sim \frac{1}{\Delta t} - \nabla \cdot D(u) \nabla$。
迭代格式：
在第 $k$ 次迭代中，冻结系数于 $u^{(k)}$，求解 $u^{(k+1)}$：
$$\mathcal{A}(u^{(k)}) \cdot u^{(k+1)} = b(u^{(k)})$$
具体离散形式（第 $i$ 网格）：
$$\frac{u_i^{(k+1)} - u_i^n}{\Delta t} = \alpha_i(u^{(k)}) u_{i-1}^{(k+1)} - \beta_i(u^{(k)}) u_{i}^{(k+1)} + \gamma_i(u^{(k)}) u_{i+1}^{(k+1)} + S_i(u^{(k)})$$

**算法流程**：

1. 令初值 $u^{(0)} = u^n$。
2. **Loop** $k = 0, 1, 2 \dots$:
    - 计算系数 $D^{(k)} = D(u^{(k)}), S^{(k)} = S(u^{(k)})$。
    - 构建矩阵 $\mathbf{M}^{(k)}$ 和向量 $\mathbf{RHS}^{(k)}$。
    - 求解线性系统 $\mathbf{M}^{(k)} u^* = \mathbf{RHS}^{(k)}$。
    - **欠松弛**：$u^{(k+1)} = \omega u^* + (1-\omega) u^{(k)}$ （防止震荡）。
    - **检查收敛**：若 $||u^{(k+1)} - u^{(k)}|| < \epsilon$，退出。
3. 令 $u^{n+1} = u^{(k+1)}$。
    

---
### 4.4.3 全隐式牛顿-拉夫逊法 (Newton-Raphson) —— TORAX 核心

核心思想：
寻找非线性残差函数 $F(u^{n+1}) = 0$ 的根。利用泰勒展开保留一阶导数项（雅可比矩阵），实现二次收敛。
数学推导：
定义全隐式残差向量 $F(u)$（省略时间上标 $n+1$）：
$$F(u) = \frac{u - u^n}{\Delta t} - \nabla \cdot (D(u) \nabla u) - S(u) = 0$$
在第 $k$ 次迭代值 $u^{(k)}$ 处进行泰勒展开：
$$F(u^{(k)} + \delta u) \approx F(u^{(k)}) + \frac{\partial F}{\partial u}\bigg|_{u^{(k)}} \cdot \delta u = 0$$
由此得到关于修正量 $\delta u$ 的线性方程：
$$\mathbf{J}^{(k)} \cdot \delta u = -F(u^{(k)})$$
其中 $\mathbf{J}$ 是雅可比矩阵。

雅可比矩阵的项 (The Hard Part vs. JAX)：
$J_{ij} = \frac{\partial F_i}{\partial u_j}$。对于扩散项 $\nabla \cdot (D(u)\nabla u)$，其导数包含两部分：
$$\frac{\partial}{\partial u} [\nabla \cdot (D(u)\nabla u)] = \underbrace{\nabla \cdot (D(u) \nabla (\delta u))}_{\text{线性部分}} + \underbrace{\nabla \cdot (\frac{\partial D}{\partial u} \delta u \cdot \nabla u)}_{\text{高度非线性部分}}$$

- **传统方法**：如果 $D$ 是 $T$ 的复杂函数（如湍流模型），第二项极难手算推导。
- **TORAX (JAX)**：`J = jax.jacobian(residual_fn)(u_k)`。计算机自动应用链式法则，精确计算包含上述两项的 $\mathbf{J}$。

**算法流程**：
1. 初值 $u^{(0)} = u^n$。
2. **Loop** $k = 0, 1 \dots$:
    - 计算残差向量 $F(u^{(k)})$。
    - 计算雅可比矩阵 $\mathbf{J}(u^{(k)})$ (JAX 自动完成)。
    - 求解线性系统 $\mathbf{J} \delta u = -F$。
    - 更新 $u^{(k+1)} = u^{(k)} + \delta u$。
    - 检查 $||F||$ 或 $||\delta u||$ 是否收敛。
---
### 4.4.4 无雅可比牛顿-克雷洛夫法 (JFNK)

核心思想：
在牛顿法中，求解 $\mathbf{J} \delta u = -F$ 时，使用 Krylov 子空间迭代（GMRES）。GMRES 不需要显式的 $\mathbf{J}$，只需要计算矩阵-向量乘积 $\mathbf{J} \cdot v$。
数学推导：
我们需要计算 $\mathbf{J}(u) \cdot v$。根据导数定义：
$$\mathbf{J}(u) \cdot v = \lim_{\epsilon \to 0} \frac{F(u + \epsilon v) - F(u)}{\epsilon}$$
在数值上，取一个极小的 $\epsilon$ 进行有限差分近似：
$$\mathbf{J} v \approx \frac{F(u + \epsilon v) - F(u)}{\epsilon}$$

**算法流程**：
1. **外层 (Newton)**：
    - 计算当前残差 $F_k = F(u^{(k)})$。
    - 调用 GMRES 求解器寻找 $\delta u$。 
2. **内层 (GMRES)**：
    - GMRES 算法给出一个搜索方向向量 $v$。  
    - **Matrix-Free 操作**：不查矩阵，而是运行两次物理残差函数：计算 $F(u^{(k)} + \epsilon v)$。  
    - 返回 $(F_{pert} - F_k)/\epsilon$ 给 GMRES。  
    - GMRES 迭代直到满足线性收敛标准。  
3. 更新 $u^{(k+1)} = u^{(k)} + \delta u$。
---

### 4.4.5 优化器求解法 (Optimizer Method)

核心思想：
将求根问题 $F(u)=0$ 转化为最小化问题（最小二乘法）。
数学推导：
定义标量损失函数（Loss Function）：
$$L(u) = \frac{1}{2} || F(u) ||^2 = \frac{1}{2} \sum_i F_i(u)^2$$

求 $L(u)$ 的极小值等价于求 $F(u)=0$。
为了使用梯度下降，我们需要 $L$ 对 $u$ 的梯度：
$$\nabla L = \left( \frac{\partial L}{\partial u} \right)^T = \left( \frac{\partial F}{\partial u} \right)^T F = \mathbf{J}^T \cdot F$$

**算法流程**：

1. 初值 $u^{(0)}$。
2. **Loop** $k = 0, 1 \dots$:
    - 计算残差 $F(u^{(k)})$。
    - 计算梯度 $g^{(k)} = \text{grad}(L)(u^{(k)})$ (利用 JAX `jax.grad`)。
    - 更新：利用优化器规则（如 Adam）：
$$m_k = \beta_1 m_{k-1} + (1-\beta_1) g^{(k)} \quad (\text{动量})$$
$$u^{(k+1)} = u^{(k)} - \eta \frac{m_k}{\sqrt{v_k} + \epsilon}$$
3. 当 $L(u)$ 足够小时停止。


---
### 4.4 算法选择指南
### 方法选择指南（基于推导的总结）

1. **Linear / Operator Splitting**: 省去了构建 $\mathbf{J}$ 和迭代的过程，直接解 $\mathbf{A}x=b$。$\to$ **快，但不准**。
2. **Picard**: 迭代更新 $\mathbf{A}(u)$，忽略了 $\frac{\partial D}{\partial u}$ 项。$\to$ **简单，但对强非线性收敛慢**。
3. **Newton**: 构建包含 $\frac{\partial D}{\partial u}$ 的完整 $\mathbf{J}$。$\to$ **最难实现（JAX除外），但收敛最快最稳**。
4. **Optimizer**: 避开矩阵求逆，只走下坡路。$\to$ **在牛顿法走投无路（矩阵奇异）时的救命稻草**。

| 求解方法 | 非线性处理 | 耦合强度 | 收敛速度 | 稳定性 | 适用场景 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **算子分裂** | 解耦求解 | 弱 | N/A | 差 | 教学、弱耦合问题 |
| **线性 Theta** | 冻结系数 | 弱 | N/A | 中 | 快速扫描、初值猜测 |
| **Picard** | 迭代更新系数 | 中/强 | 线性 (慢) | 高 | 模块化代码、中等非线性 |
| **Newton-Raphson** | 全耦合导数 | **极强** | **二次 (极快)** | **极高** | **TORAX 推荐**，刚性问题 |
| **Optimizer** | 梯度下降 | 强 | 线性 | **极高 (鲁棒)** | 初值恶劣、牛顿法发散时 |

---

## 第五章 输运与平衡的闭环耦合实现

本章阐述如何在为什么需要平衡和输运的自洽迭代，如何实现几何与剖面的自洽演化。
关键字：**解决粒子/能量输运与磁通扩散（电流演化）的强耦合求解问题**
![[00-三个方程耦合.png]]
### 5.0 为什么要耦合一起计算


 一、从方程出发：三大核心方程的强耦合

我们需要同时求解以下三个方程，它们通过多个物理量**深度耦合**，形成一个闭环系统：

 1. 磁扩散方程（决定电流如何演化）
描述极向磁通 $\psi$ 随时间演化的方程：
$$
\sigma_\parallel \frac{\partial \psi}{\partial t} = \frac{1}{\mu_0 R^2} \nabla \cdot \left( \frac{R^2}{\sigma_\parallel} \nabla \psi \right) - j_{\text{ni}}
$$

- **耦合点 A（电阻率 $\sigma_\parallel$）**：  
  $\sigma_\parallel \propto T_e^{3/2}$  
  **物理含义**：电子温度 $T_e$ 决定了等离子体电阻。温度越高，电阻越小，电流扩散越慢（趋肤效应）。  
  **若分开求解**：在计算电流扩散时使用“上一时刻”的 $T_e$，在加热突变（如 NBI 开启）或电流爬升阶段将导致**巨大误差**——因为 $T_e$ 的瞬时变化会立即改变电流扩散速率。

---

 2. 输运方程（决定温度/密度如何演化）
以电子能量守恒为例：
$$
\frac{3}{2} \frac{\partial (n_e T_e)}{\partial t} = \nabla \cdot (n_e \chi_e \nabla T_e) + P_{\text{OH}} + P_{\text{fus}} + \cdots
$$

- **耦合点 B（欧姆加热 $P_{\text{OH}}$）**：  
  $P_{\text{OH}} = \eta j^2 \propto T_e^{-3/2} (\nabla^2 \psi)^2$  
  **物理含义**：欧姆加热同时依赖于温度（通过电阻率 $\eta$）和电流密度（通过 $\psi$ 的二阶导）。  
  **关键问题**：$T_e$ 和 $\psi$ **互为因果**——电流分布变化 → 加热分布变化 → 温度变化 → 电阻变化 → 电流再变化。  
  **必须同步更新**：只有通过统一雅可比矩阵（Jacobian）在同一时间步内联立求解，才能保证自洽性。

- **耦合点 C（输运系数 $\chi_e$）**：  
  $\chi_e = \chi_e(s, q, \dots)$，其中磁剪切 $s$ 和安全因子 $q$ 由 $\psi$ 决定。  
  **物理含义**：磁面结构控制微观湍流强度。  
  **若冻结 $\chi_e$**：无法捕捉如**反剪切（Reverse Shear）形成内部输运垒（ITB）** 时的输运系数突变，导致模拟失真。

---

 3. Grad-Shafranov (GS) 方程（决定磁面形状）
$$
\Delta^* \psi = -\mu_0 R^2 p'(\psi) - F F'(\psi)
$$

- **耦合点 D（几何反馈）**：  
  **物理含义**：输运方程在磁面坐标 $\rho(\psi)$ 上求解。压强 $p$（由 $n, T$ 决定）和电流的变化会改变磁面形状（如 Shafranov 位移 $\Delta_{\text{sh}}$），进而改变几何因子 $\langle |\nabla \rho|^2 \rangle$，**直接修正输运方程中的有效扩散项**。  
  **自洽性痛点**：在电流爬升段（Ramp-up），磁面快速变形。若采用“先输运 → 后平衡”的算子分裂法，将**忽略磁面挤压对输运的瞬时调制效应**，极易导致数值发散或物理失真。

---
必须联立求解，是因为在托卡马克全放电周期中，‘热-电-磁’三个物理过程的耦合具有极强的**刚性（Stiffness）**，分开求解会切断关键的物理反馈链条。具体体现在三个方面：  
>   
> 1. **电阻率的瞬时反馈（热-电耦合）**：  
>    在电流爬升阶段，电子温度 $T_e$ 的升高会瞬间降低电阻率（$\eta \propto T_e^{-3/2}$），从而‘冻结’电流分布。如果我们分开求解，电流扩散方程就会使用‘过时’的电阻率，导致无法正确模拟出**反剪切（Reverse Shear）结构**的形成——而这正是高性能放电的关键。  
>    
> 2. **源项的高度非线性（源-场耦合）**：  
>    欧姆加热 $P_{\text{OH}} \propto T_e^{-3/2} j^2$ 同时依赖于温度和电流。在我们的统一框架中，雅可比矩阵（Jacobian）包含了交叉导数项 $\partial P_{\text{OH}} / \partial \psi$ 和 $\partial P_{\text{OH}} / \partial T_e$。这意味着我们在数学上保证了能量方程和磁扩散方程在同一个时间步内是**严格自洽的**，避免了数值震荡。  
>    
> 3. **几何度规的动态响应（输运-平衡耦合）**：  
>    特别是在高 $\beta_p$ 运行模式下，压强梯度的变化会瞬间引起 Shafranov 位移，改变磁通面的几何因子 $\langle |\nabla \rho|^2 \rangle$。我们的方案通过全隐式迭代，确保输运过程能感知到磁面形状的**实时变化**，这对于模拟 VDE（垂直位移事件）或 L-H 转换等瞬态过程是至关重要的。  
>    
> 所以，联立求解不是为了炫技，而是为了在数学上守住这些物理方程的**本征耦合特性**。”

---
三、总结图示（逻辑闭环）

建议在答辩时心中构建如下关联图，逻辑将极为清晰：

- 温度 $T_e$ $\xrightarrow{\text{决定}}$ 电阻 $\eta$ $\xrightarrow{\text{控制}}$ 电流 $j$  
- 电流 $j$ $\xrightarrow{\text{决定}}$ 欧姆热 $P_{\text{OH}}$ $\xrightarrow{\text{加热}}$ 温度 $T_e$  
- 压强 $p$ $\xrightarrow{\text{推挤}}$ 磁面形状 $\psi$ $\xrightarrow{\text{改变}}$ 几何因子 $\langle |\nabla \rho|^2 \rangle$ $\xrightarrow{\text{修正}}$ 输运率  

> **结论**：这是一个**物理闭环**。  

    
### 5.1 耦合逻辑
耦合的核心在于**时间尺度的分离**：输运是快过程，平衡变化是慢过程。
![[00-平衡-输运迭代时间步.png]]
1.  **输运 $\to$ 平衡 (Transport to Equilibrium)**：
    * TORAX 计算出 $n_e, T_e, \psi$。
    * 导出压强 $p(\rho) = \sum n T$ 和环向电流 $J_{tot}(\rho)$。
    * 转换为 G-S 源项：$p'(\psi)$ 和 $FF'(\psi)$。
2.  **平衡 $\to$ 输运 (Equilibrium to Transport)**：
    * FreeGS 求解 G-S 方程，得到二维磁面 $\psi(R,Z)$。
    * 积分计算一维几何系数：$V(\rho), S(\rho), \langle |\nabla \rho|^2 \rangle$ 等。
    * 这些系数直接修正输运方程中的通量项。

### 5.2 **高精度网格重映射**

- **问题**：FreeGS 计算出的 $\psi(R,Z)$ 需要映射到 TORAX 的 $\rho_{tor}$ 网格上。
    
- **方案**：
    
    1. **磁面追踪**：在 FreeGS 中沿 $\psi$ 等值线进行高精度追踪，计算环向磁通 $\Phi(\psi)$，从而建立 $\rho_{tor} \leftrightarrow \psi_{pol}$ 的映射关系。
        
    2. **样条插值 (Spline)**：使用三次样条插值将几何系数（如 $V', \langle \nabla \rho \rangle$）插值到 TORAX 的节点上。
        
    3. **守恒性检查**：在插值过程中，必须保证总等离子体体积和总感应磁通守恒，否则长时间运行会导致“虚假”的粒子或能量源。

### **平滑处理 (Smoothing)**

> 
> 在耦合迭代中，如果 FreeGS 更新后的几何形状变化过大（例如 X 点位置跳变），会导致输运求解器在下一个时间步崩溃。
> - **松弛技术**：在更新几何时，不要直接使用 $G_{new}$，而是使用 $G_{input} = \alpha G_{new} + (1-\alpha) G_{old}$（例如 $\alpha=0.5$）进行平滑过渡。
>
---
### 第五部分：验证与基准测试 (V & V)**
> 
> 1. **静态验证**：固定所有源项和几何，将时间步设为无限大（或跑至稳态），对比 TORAX 的结果与解析解（如圆柱近似下的 Bessel 函数解）。   
> 2. **代码对标 (Code-to-Code Benchmark)**：
>     - **输运对标**：使用标准的 ITER 算例，将 TORAX 结果与 **ASTRA** 或 **RAPTOR** 进行对比。  
>     - **平衡对标**：将 FreeGS 的结果与 **EFIT** 或 **CHEASE** 对比。    
> 3. **守恒性验证**：在模拟过程中实时监控 $\int n dV$（总粒子数）和 $\int \frac{3}{2}nT dV$（总能量）。除了物理源项和边界损失外，数值误差导致的漂移应低于 $10^{-5}$。
## 第六章 总结

本报告建立了一套完整的托卡马克自洽模拟技术框架。

1.  从 **Grad-Shafranov 平衡** 和 **1.5D 输运方程** 出发，确保了物理模型的完备性。
2.  详细推导了 **FVM 空间离散** 和 **全隐式时间离散**，奠定了数值基础。
3.  系统阐述了 **Picard、Newton-Raphson、Optimizer** 等求解器的数学原理及适用场景，特别是利用 **JAX 自动微分** 解决牛顿法实现难点。
4.  设计了 **输运-平衡闭环耦合** 流程，实现了等离子体形状与内部剖面的自洽演化。

这套框架既具有极高的数值稳定性（适合长脉冲模拟），又具备高度的灵活性（可更换物理模型），是现代聚变模拟研究的有力工具。