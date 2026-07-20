# 以 Tile 为中心的 C++ 类划分设计

## 文档范围

这篇文档只讨论后续 C++ 重构里的类职责划分，以及类之间的大体协作关系。

这篇文档暂时不讨论：

- 具体文件怎么拆
- 头文件和源文件怎么落地
- 命名空间怎么组织
- 迁移顺序怎么安排
- 模板和继承体系怎么设计

这一阶段只先回答一个问题：

> 如果以后以 `Tile` 作为最核心的对象，那么围绕它应该有哪些类，每个类分别负责什么？

## 总体方向

这个方向我认为是对的，而且很适合当前仓库的实际情况。

原因是现在的 tile 主流程本身就已经隐含了几个比较清楚的阶段：

1. 准备 tile 输入
2. 做 early-stop / hybrid 分类
3. 构建调度计划
4. 对允许执行的行做解码
5. 把结果写回 tile 输出

也就是说，当前代码不是没有结构，而是结构更多还停留在“函数阶段”和“隐式约定”上。  
这次重构最合适的做法，不是重新发明一套新的算法流程，而是把已经存在的流程显式地落到类职责上。

## 核心原则

`Tile` 应该是总控类，但不应该是“万能大类”。

也就是说：

- `Tile` 负责组织一次 tile 的完整生命周期
- `Tile` 负责协调各个子模块
- `Tile` 负责定义处理顺序
- 但 `Tile` 不应该把存储、分类、调度、解码、统计全部自己吞进去

如果把所有逻辑都塞进 `Tile` 里，那么最后只是把“大量自由函数”换成了“一个超大类”，结构问题并没有真正解决。

所以更合理的方向是：

- `Tile` 做编排
- 各个子类分别承担清晰、单一的职责
- 阶段之间通过明确的数据对象交接

## 推荐的类划分

## 1. `Tile`

`Tile` 是 tile 级处理流程的核心对象。

它的职责应该是：

- 接收一次 tile 处理所需的输入和配置快照
- 创建或持有本次 tile 执行需要的协作对象
- 按固定顺序驱动 tile 主流程
- 最终返回一个 `TileResult`

`Tile` 主要回答的问题是：

- 一个 tile 要经历哪些处理阶段？
- 这些阶段按什么顺序执行？
- 每个阶段之间传递什么中间结果？

`Tile` 不应该直接负责：

- early-stop 的具体判定细节
- hybrid 分类器的具体实现
- soft / hard decode 的底层实现
- 存储结构内部怎么组织

换句话说，`Tile` 更像“总调度者”，而不是“亲自做所有事情的人”。

## 2. `TileStorage`

`TileStorage` 负责 tile 局部数据的持有、组织和访问。

它的职责应该是：

- 保存 tile 的输入矩阵和相关参考数据
- 提供分类器和解码器需要的行视图、线性化视图或中间缓冲区
- 保存本轮 tile 处理中可变的工作数据
- 在解码结束后承接结果并形成最终 tile 输出

`TileStorage` 的重点是“管理数据”，而不是“做决策”。

它应该知道：

- 原始 tile 输入是什么
- channel tile 输入是什么
- 行级线性化后的数据长什么样
- 中间工作缓冲区放在哪里
- 最终 tile 输出如何组织

它不应该知道：

- 某一行属于 `HardFinish` 还是 `SoftDecode`
- 某一行为什么进入 early-stop
- 哪些行要参与 MUX 竞争

所以 `TileStorage` 的角色是“数据所有者”，不是“算法决策者”。

## 3. `TileClassifier`

`TileClassifier` 负责 tile 内部行级语义的判定。

它的职责应该是：

- 执行 early-stop 条件判断
- 执行 hybrid prepass 分类
- 给每一行打上算法意义上的分类标签
- 生成结构化的分类结果，供后续调度阶段使用

这类逻辑包括：

- 这一行是否命中 early-stop
- 这一行属于 `Clean`、`ParityOnly`、`OneMain`、`TwoMain`、`HardFail` 中的哪一类
- 这一行是否具备 hard-finish 候选资格

这里最重要的边界是：

`TileClassifier` 只负责说明“这一行是什么”，不负责决定“这一行这一轮最终怎么处理”。

也就是说：

- 分类器产出的是分类结果
- 不是最终调度动作

这个边界一定要守住，否则分类和调度又会重新混在一起。

## 4. `TileDispatchPlanner`

虽然你最初的想法主要是 `Tile + 存储类 + 分类器类 + 解码器类`，但我非常建议把“调度规划”单独拆成一类。

原因很简单：

如果没有这个类，那么分类器和解码器之间的那一层 MUX / SISO / 行动作决策，最后还是会重新混到别的类里。

`TileDispatchPlanner` 的职责应该是：

- 把分类结果转换成这一轮真正的执行计划
- 应用 MUX / SISO budget 规则
- 决定每一行最终属于哪种执行动作
- 生成稳定的 `TileDispatchPlan`

它最终要决定的内容通常包括：

- 哪些行走 `EarlyStopAction`
- 哪些行走 `HardFinish`
- 哪些行走 `SoftDecode`
- 哪些行在本轮变成 `Unscheduled`

这个类存在的意义非常大，因为：

- 分类说的是“行的算法属性”
- 调度说的是“系统这一轮准备对它做什么”

这两个问题本质上不是一回事。

## 5. `TileDecoder`

`TileDecoder` 负责执行真正的解码工作。

它的职责应该是：

- 接收准备好的行数据和调度计划
- 对需要执行的行运行 hard decode 或 soft decode
- 产出解码结果
- 尽量不掺入上层的实验控制逻辑

`TileDecoder` 应该只关心一件事：

> 已经决定要解的那些行，怎么解？

它不应该负责：

- 为什么某行拿到了 SISO 资源
- 为什么某行被 early-stop 了
- sweep / probe / 双流实验如何组织
- 上层总统计怎么汇总

也就是说，`TileDecoder` 是“执行器”，不是“策略制定者”。

## 6. `TileResult`

`TileResult` 是一次 tile 处理的统一输出对象。

它的职责应该是：

- 承接最终 tile 输出
- 承接 tile 级统计信息
- 必要时承接行级决策结果，用于调试和观测
- 为 window / frame 上层提供稳定接口

`TileResult` 很重要，因为它能避免上层直接依赖 tile 内部那些临时工作状态和底层缓冲区。

上层只应该依赖一个清晰、稳定的结果对象，而不是去碰 tile 内部细节。

## 建议配套的中间数据对象

除了上面的几个主类，我建议把几个阶段边界上的数据也显式做成对象。

这样做的目的不是“为了面向对象而面向对象”，而是为了让每个阶段的输入输出清楚可见。

## 1. `TileInput`

`TileInput` 表示一次 tile 处理的原始输入。

它可以包含：

- tile 输入矩阵
- channel tile 矩阵
- tile 顶部全局行号
- 可选的 tx 参考
- 本次 tile 相关的配置快照

它强调的是“输入是输入”，不要一开始就和工作状态混在一起。

## 2. `TileClassificationResult`

`TileClassificationResult` 表示分类阶段的输出。

它可以包含：

- early-stop 标志
- row detail
- hybrid 分类结果
- hard-finish 候选信息
- 分类阶段的摘要统计

这个对象的意义是让分类结果在进入调度前就定型，不要靠多个零散数组隐式配合。

## 3. `TileDispatchPlan`

`TileDispatchPlan` 表示这一轮 tile 的最终执行计划。

它可以包含：

- 每一行的最终动作标签
- scheduled / unscheduled 状态
- hard-finish 标记
- soft candidate 列表
- soft scheduled 列表
- 调度阶段统计

这个对象最好成为 tile 主流程里最关键的中间结果之一。

## 4. `TileDecodeResult`

`TileDecodeResult` 表示解码执行阶段产出的原始结果。

它可以包含：

- hard-finish 输出
- soft decode 输出
- 每行的执行结果
- 解码阶段附带的局部统计

有了这个对象以后，`TileDecoder` 和 `TileStorage` 的边界会更清楚：

- `TileDecoder` 负责算
- `TileStorage` 负责接收并组织这些结果

## 推荐的协作流程

如果按这个设计走，一次 tile 执行的协作关系大致应该是：

1. `Tile` 接收 `TileInput`
2. `Tile` 创建或持有 `TileStorage`
3. `TileStorage` 准备后续阶段所需的数据视图
4. `TileClassifier` 基于存储内容做分类，返回 `TileClassificationResult`
5. `TileDispatchPlanner` 基于分类结果生成 `TileDispatchPlan`
6. `TileDecoder` 按计划执行解码，返回 `TileDecodeResult`
7. `TileStorage` 接收解码结果并整理最终 tile 输出
8. `Tile` 组装并返回 `TileResult`

这个流程有一个很大的好处：

每个阶段只回答一个主要问题。

- 存储：我有哪些数据？
- 分类：这些行分别是什么？
- 调度：这一轮该怎么处理这些行？
- 解码：真正执行的行会产出什么结果？
- Tile：这些阶段怎么串成一次完整 tile 处理？

## 为什么这样比“一个大 Tile 类”更好

如果所有逻辑都直接塞进 `Tile`，那么这个类很快就会同时拥有：

- tile 输入输出矩阵
- 行级中间缓冲
- early-stop 判定逻辑
- hybrid 分类逻辑
- MUX 逻辑
- soft / hard decode 逻辑
- writeback 逻辑
- 统计逻辑

这样一来，看起来用了类，但本质上还是“一个中心化的大过程对象”。

而这里建议的方案更好，是因为它把复杂度拆成了四种不同性质的问题：

- 数据组织复杂度
- 算法分类复杂度
- 调度规划复杂度
- 解码执行复杂度

这四种复杂度分开之后，后面继续加策略、加实验、加观测，都不会那么容易把结构再次拖乱。

## 必须坚持的几个边界

## 1. 分类不等于调度

这是最关键的边界。

例如：

- `HardFail` 是分类结果
- `Unscheduled` 是调度结果
- `ParityOnly` 是分类标签
- `SoftDecode` 是最终执行动作

这些概念如果再混回去，设计很快就会退回到现在这种“看得懂，但不够一眼清楚”的状态。

## 2. 存储不等于业务逻辑

`TileStorage` 不应该变成一个隐藏的算法引擎。

它可以提供访问接口、准备视图、管理缓冲区，但不应该替代分类器、调度器和解码器做业务决策。

## 3. 解码器不等于策略层

`TileDecoder` 应该执行策略，而不是偷偷决定策略。

也就是说，尽量不要让解码器内部暗含：

- MUX 决策
- hybrid 重分类
- 写回归属不清的状态修改

否则类虽然拆了，策略仍然会藏在难追踪的地方。

## 对上层 window / frame 的意义

如果 tile 这一层抽象做得好，那么上层的 window 和 frame 逻辑会自然变简单。

未来 window / frame 层 ideally 只需要做这些事：

- 生成 tile 输入
- 调用 `Tile`
- 收集 `TileResult`
- 向上继续聚合统计

这说明一个好的 tile 抽象，不只是让 tile 自己更清楚，也能反过来降低上层复杂度。

## 建议优先显式类型化的概念

为了让类设计在后面保持稳定，我建议下面这些概念尽量做成明确类型，而不是散落成若干布尔值和数组约定：

- `HybridClassification`
- `RowDispatchTag`
- `SoftSchedulingState`
- `TileClassificationResult`
- `TileDispatchPlan`
- `TileResult`

一旦这些概念都有了明确类型，代码里很多边界就不容易再被无意中混掉。

## 最终建议

你的方向是对的，但我建议做一个小升级：

不是只停在：

- `Tile`
- `TileStorage`
- `TileClassifier`
- `TileDecoder`

而是进一步明确为：

- `Tile`
- `TileStorage`
- `TileClassifier`
- `TileDispatchPlanner`
- `TileDecoder`
- `TileResult`

其中：

- `Tile` 是总控
- `TileStorage` 管数据
- `TileClassifier` 管“这行是什么”
- `TileDispatchPlanner` 管“这一轮怎么处理这行”
- `TileDecoder` 管“真正怎么解”
- `TileResult` 管最终输出

这个划分我认为是比较平衡的：

- 足够面向对象，职责边界明确
- 又不会一上来过度设计太多抽象层
- 和当前仓库已有的算法阶段天然对齐

## 这一阶段先不要急着定的事情

在类职责还没有完全达成共识前，下面这些事情先不要过早定死：

- 每个类具体放在哪个文件
- 是否所有类都要做接口抽象
- 是否一开始就全模板化
- 旧的 `.ipp` 要保留多少
- 哪些实现要先搬、哪些实现后搬

这些都属于下一阶段的问题。

当前最重要的是先统一：

> 每一类复杂度到底由谁来负责。

只要这个问题先定清楚，后面的落地就会顺很多。
