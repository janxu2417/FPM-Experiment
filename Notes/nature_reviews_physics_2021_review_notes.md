# Nature Reviews Physics 2021: Fourier ptychography Review 笔记

原文 PDF：
[Nature Reviews Physics 2021 Concept, implementations and applications of Fourier ptychography.pdf](/D:/1_徐前/物理/傅里叶叠层显微术/原始论文/Nature%20Reviews%20Physics%202021%20Concept,%20implementations%20and%20applications%20of%20Fourier%20ptychography.pdf)

提取文本：
[clean.txt](/D:/1_徐前/物理/傅里叶叠层显微术/pdf_output/nature_reviews_physics_2021/Nature%20Reviews%20Physics%202021%20Concept,%20implementations%20and%20applications%20of%20Fourier%20ptychography_clean.txt)

## 1. 文章定位

- 这是一篇 `Nature Reviews Physics` 的技术综述，不是单一实验论文。
- 主线不是“报告一个新系统参数”，而是系统梳理：
  - Fourier ptychography (FP/FPM) 的基本概念
  - 与 ptychography、synthetic aperture、phase retrieval 的关系
  - 不同硬件实现
  - 典型应用
  - 系统误差校正与未来方向

对实验复现最有价值的部分有三块：

1. FP 的 forward model 和 alternating projection 逻辑
2. `Box 2` 给出的实验搭建与 LED 对准方法
3. system correction 一节中关于 LED 强度、LED 入射角、pupil aberration 的处理

## 2. 核心问题意识

文章开头强调的核心矛盾是：

- 传统显微镜通常要在 `分辨率` 和 `视场` 之间折中
- 物理光学系统的通量可用 `space-bandwidth product (SBP)` 描述
- 大多数现成物镜的 SBP 只有 `~10 megapixels`

FP 的意义在于：

- 不依赖大幅改造硬件
- 利用计算，把“提升 SBP”的问题从“复杂光学设计”转移到“计算重建”

## 3. FPM 核心原理

### 3.1 两个基础支柱

FP = `synthetic aperture imaging + phase retrieval`

- synthetic aperture：通过多次采样合成更大的有效孔径
- phase retrieval：只测强度，不直接测相位，通过迭代恢复复振幅

### 3.2 成像模型

对于第 `i` 个 LED：

- LED 提供一个斜入射平面波
- 该照明把样品频谱平移到新的位置
- 物镜 pupil 只截取其中一个圆形低通区域
- 相机记录对应的低分辨率强度图

因此：

- 每张图像对应傅里叶域中的一个子孔径
- 多个入射角对应多个 shifted pupils
- 这些子孔径在傅里叶域有重叠时，就可以做稳定相位恢复

### 3.3 迭代恢复逻辑

FP 的标准流程可概括为：

1. 初始化高分辨率复数物体
2. 取某个 LED 对应的频谱子区域
3. 逆变换得到当前低分辨率复场估计
4. 用实测强度平方根替换振幅，保留相位
5. 正变换回傅里叶域，把更新后的子区域写回全局频谱
6. 对所有 LED 重复
7. 全局再迭代至收敛

这就是文章里反复强调的 alternating projection 思想：

- 空间域：modulus constraint
- 傅里叶域：pupil support constraint

### 3.4 为什么 overlap 很关键

文章明确指出：

- 相邻采样的子孔径必须有重叠
- 没有重叠时，每个子问题会彼此独立，回到 phase problem 的固有歧义

因此，对实验来说：

- `频谱重叠率` 是成败关键量
- 它由 `NA_obj`、`LED 几何`、`样品到光源距离` 共同决定

## 4. FP 与 real-space ptychography 的关系

文章把两者关系讲得很清楚：

- real-space ptychography：
  - 空间域用 probe support 约束
  - 傅里叶域用衍射测量模值约束
- FP：
  - 傅里叶域用 pupil support 约束
  - 空间域用图像强度模值约束

可以理解为：

- FP 用了物镜
- 因此在一定意义上“交换了” real-space ptychography 里的空间域/傅里叶域角色

这个类比很重要，因为后面很多 correction 方法，就是从 ptychography 那边借过来的。

## 5. 实现层面的关键点

### 5.1 最原始的 LED array microscope

文章回顾的原始 FP 实现：

- 用低 NA 物镜保持大视场
- 用 LED 阵列提供不同入射角
- 每个 LED 对应一张低分辨率图像
- 重建后获得高分辨率复数图像，同时保留大 FOV

### 5.2 非均匀采样

文章指出：

- 原始实现中，很多采样点分布是准周期网格
- 可改用非均匀采样减少采集数量
- 例如从 `137` 幅减少到约 `68` 幅

动机有两个：

1. 中心频谱能量高，可以更密采
2. 边缘 darkfield 区域可以降低 overlap，减少采样数

实验启发：

- 不是所有 LED 都必须等权采
- 采样可以围绕“信息密度”和“重叠需求”来设计

### 5.3 多种硬件实现

review 讲了很多变体：

- 传统 LED 阵列
- annular illumination
- multiplexed coded illumination
- laser scanning / mirror array / DMD
- aperture-scanning FP
- camera-scanning FP
- reflection-mode FP
- single-shot FP
- translated speckle illumination

对当前显微实验复现来说，最相关的仍然是：

- LED 阵列方案
- annular illumination
- multiplexed illumination

## 6. Box 2：实验搭建与 LED 几何标定

这是最值得反复读的部分。

### 6.1 背景与杂散光校正

文章建议：

- 在空白载玻片上采集 darkfield reference
- 用这些 reference 去减掉 ambient light 和 stray light

实验意义：

- 这一步很便宜，但很重要
- 尤其对高角度暗场 LED，背景和杂散光会显著污染约束

### 6.2 用 brightfield-to-darkfield transition 做 LED 阵列对准

文章的实操方法是：

1. 选一个 LED 作为中心参考点
2. 选 4 个彼此中心对称、且入射角接近物镜接受角边界的 LED
3. 观察这 4 张图里的 brightfield-to-darkfield transition
4. 用这些 transition 的位置来估计：
   - LED 阵列相对样品/光轴的位置
   - 阵列方向角

这是一个非常实用的方法，因为：

- brightfield/darkfield 分界本质上对应 `NA_illumination ≈ NA_objective`
- 对几何偏移非常敏感
- 不需要复杂样品

### 6.3 分块重建的意义

文章提到：

- 原始图像通常会被分成小块，例如 `256 × 256`

原因：

- 每个 tile 内，pupil aberration 可以近似为空间不变
- LED illumination 可以近似为平面波

这对你实验的直接启发是：

- 先在小块上调好参数，再推广到整幅图
- 不要一开始全视场同时估 LED 几何和 pupil

### 6.4 mask 掉不可信区域

文章还明确说：

- brightfield-to-darkfield transition 区域
- stray-light 亮斑区域

这些区域可以用 binary mask 标出来，在更新 object 时不参与 modulus constraint。

这非常重要，因为：

- 不干净的像素如果强行进入迭代，会把模型拉偏

## 7. LED 几何标定：对实验真正可执行的流程

结合文章内容，可以整理成下面这套流程。

### 7.1 第一层：全局几何粗标定

要标的全局参数：

- LED 阵列中心相对光轴偏移
- 阵列旋转角
- LED 到样品的有效高度 `H`

操作：

1. 放空白片或均匀弱散射样品
2. 采 darkfield 背景
3. 采中心 LED
4. 采 4 个中心对称、接近 brightfield/darkfield 分界的 LED
5. 拟合 transition 边界位置
6. 解出中心、旋转、H

### 7.2 第二层：wavevector refinement

如果阵列 pitch 不够准，或不是规则 LED 阵列，只做全局几何不够。

文章建议：

- 在迭代重建中进一步 refine illumination wavevector
- 这与 ptychography 的 positional correction 是同类问题

实验上建议：

1. 先固定全局参数
2. 只允许每个 LED 的 `(k_x^i, k_y^i)` 做小范围微调
3. 用重建误差、高频清晰度、频谱重叠一致性作为优化目标

### 7.3 实验上最稳的参数优化顺序

建议顺序：

1. 背景校正
2. LED 全局几何标定
3. 小块中心视场重建
4. 局部 wavevector refinement
5. pupil / wavefront 校正
6. 全视场重建

而不是：

- 一开始就全 LED 自由优化
- 同时自由恢复几何和高维 pupil

后者很容易参数耦合。

## 8. 波前校正 / pupil correction

这是 review 另一块非常重要的内容。

### 8.1 三类典型系统误差

文章明确写了三个最常见的系统误差：

1. `LED intensity variations`
2. `incident angles of the LED elements`
3. `pupil aberrations`

这其实就给了实验排错优先级。

### 8.2 LED 强度不一致

文章指出：

- 可用 image-quality metric 在迭代中更新 LED intensity

启发：

- 每个 LED 并不一定同亮
- 不先做强度归一化，重建会把亮度误差错当成样品/波前结构

实验上可做：

- 先做独立平场标定
- 或在重建中把每个 LED 的强度系数作为待估参数

### 8.3 pupil aberration 的两条路线

文章给了两条主路线：

#### 路线 A：用标定靶直接测

- 用 calibration target 测 `spatially varying aberrations`
- 对不同视场位置恢复 pupil
- 拟合 field-dependent aberration model

优点：

- 可解释性强
- 易于工程化

#### 路线 B：联合恢复 object + pupil

- 在迭代中同时恢复样品和 pupil
- 文中对应 `embedded pupil function recovery`

优点：

- 不一定需要单独标定样品

缺点：

- 更容易与 LED 几何误差耦合
- 更依赖初值和正则化

### 8.4 field-dependent aberration

文章特别强调：

- 不应简单把每个 tile 的 pupil 完全独立恢复
- 更好的方法是建模 `pupil` 随视场位置变化

这样可以：

- 降低自由度
- 提高全视场鲁棒性
- 做 full-field aberration metrology

这对大视场 FPM 尤其重要。

### 8.5 defocus 是最先该处理的像差

文章明确说：

- defocus 本质上也是一种 aberration
- 因此 FP 能做 computational refocus

实验建议：

- 在高阶像差前，先单独扫 defocus
- 这通常最稳、收益最大

## 9. 对当前复现实验的直接指导

结合你的参数：

- `NA_obj=0.15`
- `wlength=628nm`
- `H≈4cm`
- `LEDp=4mm`

### 9.1 先关心什么

优先级建议：

1. LED 阵列中心与旋转是否正确
2. 有效高度 `H` 是否准确
3. LED 强度是否一致
4. 中心视场的 defocus 是否明显
5. 再看高阶 pupil 像差

### 9.2 为什么你的系统对 LED 几何更敏感

因为：

- `LEDp=4 mm`
- `H≈4 cm`

意味着相邻 LED 的角度步进不算很小。  
一旦 `H`、中心偏移、阵列旋转有误，映射到 `(k_x, k_y)` 的误差会直接影响频谱拼接位置。

### 9.3 对你最建议的实际操作顺序

建议你按下面做：

1. 用空白片先采 darkfield reference
2. 用 4 个近 brightfield/darkfield 分界的中心对称 LED 做几何标定
3. 只在中心小块重建
4. 先扫 defocus
5. 再做小范围 LED wavevector refinement
6. 最后再引入 pupil 恢复

### 9.4 不建议一开始做的事

- 不建议一开始全视场同时恢复几何和 pupil
- 不建议在几何还没稳时就放开高阶 Zernike
- 不建议忽略 darkfield 背景和 stray-light mask

## 10. 文章里值得记住的几句结论

- FP 的本质是用 computation 扩展系统 SBP。
- phase retrieval 的稳定性依赖数据冗余，也就是频谱 overlap。
- LED 几何误差、LED 强度误差、pupil 像差都应被看作 forward model 的一部分。
- aberration correction 在 FP 中是“后校正问题”，不只是“光学设计问题”。
- field-dependent aberration 建模是大视场高质量重建的关键方向。

## 11. 后续精读建议

建议你后续重点细读这几部分：

1. 页 1–3：基本原理与 Box 1
2. 页 5：`Box 2` 的实验搭建和 LED 对准
3. system correction 那一段：LED intensity / incident angle / pupil aberration
4. aberration metrology 那一段：field-dependent pupil model

如果继续做实验复现，这篇 review 最有价值的不是“给某组固定参数”，而是给出了：

- 一套完整的系统误差分类框架
- 一套合理的参数校正顺序
- 多种可迁移到你系统里的 correction 思路
