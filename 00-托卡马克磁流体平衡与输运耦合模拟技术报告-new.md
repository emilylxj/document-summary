# 托卡马克磁流体平衡与输运耦合模拟技术报告

## 1. 引言

...[existing content]...

## 2. 背景

...[existing content]...

## 3. 方法

### 3.1 Picard迭代

在Picard迭代中，我们通过以下方式更新解：

A(u^(k-1))·u^(k) = b(u^(k-1))，而没有使用雅可比矩阵。

### 3.2 牛顿法

...[existing content]...

## 4. 离散化流

### 4.1 离散化过程

...[existing content]...

### 4.2 形状函数

...[existing content]...

### 4.3 离散化流总结

| Continuous PDE | [FVM space] | Semi-discrete ODE | [Theta time] | Fully discrete | [Flux linearization] | Linear system Ax=b |
|----------------|--------------|-------------------|---------------|----------------|----------------------|---------------------|

...(discretization flow summary table)...

### 4.4 其他部分

...[existing content]...

### 4.5 结论

...[existing content]...

### 4.8.2 Picard迭代

在此部分，我们将会去除错误的泰勒展开及雅可比矩阵描述，替换为正确的Picard定点迭代。

...[additional content updated]...

## 5. 结论

...[existing content]...

## 版本号

当前版本为 3.2 (Plan A Fix)