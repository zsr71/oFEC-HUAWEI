# Chase-Pyndiah 解码器工作原理

实现了 Pyndiah 1998 提出的 SISO Chase 组件。以下以 plain 版本为例，说明各阶段的计算公式。

## 1. 输入与符号

对每个 Decoder 传入 256 维 LLR：  
$$
\mathbf{y} = (y_0, y_1, \dotsc, y_{255}),\quad y_j \in \mathbb{R}
$$
其中前 255 位对应 BCH(255,239) 码芯，最后一位是扩展偶校验。为方便描述，定义硬判决：
$$
\hat{c}_j = \begin{cases}
0, & y_j \ge 0 \\
1, & y_j < 0
\end{cases}
$$

## 2. 选取最不可靠的位置

按照 \(|y_j|\) 从小到大排序，取得 `Params::CHASE_L = L = 6` 个索引集 \(\mathcal{P} = \{p_0,\dots,p_{L-1}\}\)，其对应的幅度 \(a_\ell = |y_{p_\ell}|\)。这些位置将被翻转以生成测试向量。

## 3. 生成测试模式

`gen_test_patterns()` 通过枚举翻转模式生成 \(N_\text{test}=64\) 个布尔向量 \(\mathbf{b}^{(c)} \in \{0,1\}^L\)。据此构造候选码字  
\[
\mathbf{u}^{(c)} = \hat{\mathbf{c}} \oplus \sum_{\ell=0}^{L-1} b^{(c)}_\ell \mathbf{e}_{p_\ell},
\]
其中 \(\mathbf{e}_{p_\ell}\) 为在索引 \(p_\ell\) 处为 1、其余位置为 0 的标准基向量。随后，将 \(\mathbf{u}^{(c)}\) 输入 BCH(255,239) 硬判决译码器；若译码成功，则为其补上 overall parity 比特，得到合法码字 \(\mathbf{C}^{(c)}\)。


## 4. 候选度量与 ML 选择

对每个有效候选求内积度量：
$$
S^{(c)} = -\sum_{j=0}^{255} |y_j| \cdot \mathbb{1}\{ C^{(c)}_j \ne \hat{c}_j \}
$$
得分越大表示越接近接收向量。取得分最大的候选作为 ML 码字 \(\mathbf{C}^{\text{ML}}\)。若所有候选均译码失败，则退化为 \(\mathbf{C}^{\text{ML}} = \hat{\mathbf{c}}\)。

## 5. 外信息计算（Pyndiah 式）

对每个比特 \(j\)，在所有译码成功的候选码字中，分别寻找第 \(j\) 位为 0/1 时的最佳得分：
\[
S_j^{(+)} = \max_{c : C^{(c)}_j = 0} S^{(c)}, \quad
S_j^{(-)} = \max_{c : C^{(c)}_j = 1} S^{(c)}.
\]
若二者都存在，则根据 Pyndiah (14)–(17) 将“竞争码字”在其余分量上的差异投影为：
\[
\omega_j = \sum_{l \ne j} r_l \, c_l^{(+)} \, \mathbb{1}\{ c_l^{(+)} \ne c_l^{(-)} \},
\]
其中 \(r_l = y_l\)、\(c_l^{(+)} = (-1)^{C^{(+)}_l}\)、\(c_l^{(-)} = (-1)^{C^{(-)}_l}\)。若 \(S_j^{(+)}\) 或 \(S_j^{(-)}\) 之一不存在，则按照式 (20) 采用常量回退：
\[
\omega_j = \beta \cdot (-1)^{C^{\text{ML}}_j},
\]
由 `Params::beta` 控制该常量外信息的幅值。


## 6. Extrinsic LLR 输出

每个行的软输出仅包含外信息：
$$
\Lambda_j^\text{ext} = \omega_j
$$

