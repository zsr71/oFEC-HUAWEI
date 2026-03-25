# VS Code 调试时查看 `std::vector` 完整内容说明

## 1. 问题现象

在 VS Code 里调试 C++ 程序时，像 `std::vector` 这样的容器，常见会遇到两个问题：

- 变量窗口里只能看到前一部分元素，例如前 `999` 个
- 展开变量时，看不到完整大小，或者不方便直接定位某个下标

这通常不是程序本身有问题，而是调试器默认的显示限制。

---

## 2. 先看 `vector` 的大小

如果想看一个 `std::vector` 的实际大小，可以在 VS Code 的 `Debug Console` 或 `Watch` 里直接输入：

```cpp
inInfoBits.size()
```

如果当前表达式求值不稳定，也可以直接用 gdb 命令：

```gdb
-exec p inInfoBits.size()
```

如果还想看容量，也可以输入：

```gdb
-exec p inInfoBits.capacity()
```

---

## 3. 看某一个下标的值

如果只想确认某个特定位置的内容，可以直接输入：

```cpp
inInfoBits[1200]
```

或者：

```gdb
-exec p inInfoBits[1200]
```

这种方式适合快速检查：

- 某个位置是不是 0 / 1
- 某个异常比特附近是不是已经错了

---

## 4. 看一段连续元素

如果想看从某个下标开始的一小段，可以在 `Debug Console` 里输入：

```gdb
-exec p inInfoBits[1000]@20
```

含义是：

- 从 `inInfoBits[1000]` 开始
- 连续打印 `20` 个元素

例如，它可以用来查看：

- 第 `1000 ~ 1019` 个元素

这种方式比直接展开整个 vector 更适合定位局部问题。

---

## 5. 取消“只显示前若干个元素”的限制

gdb 默认会限制容器打印的元素个数，所以经常只显示前 `999` 个元素。

可以在 `Debug Console` 里输入：

```gdb
-exec set print elements 0
```

这里的 `0` 表示：

- 不限制打印元素数量

执行之后，再展开变量，通常就能看到完整容器，而不再只停在前 `999` 个元素。

---

## 6. 如果变量值仍然被截断

有时候除了元素个数限制，还会遇到“值太大被截断”的问题。  
这时可以继续输入：

```gdb
-exec set max-value-size unlimited
```

这样可以放宽 gdb 对大对象显示大小的限制。

---

## 7. 建议写进 `launch.json`

如果每次调试都要手工输入这些命令，会比较麻烦。  
建议把它们加到 `.vscode/launch.json` 的 `setupCommands` 里。

示例：

```json
{
  "description": "Enable pretty-printing",
  "text": "-enable-pretty-printing",
  "ignoreFailures": true
},
{
  "description": "Show all container elements",
  "text": "-gdb-set print elements 0",
  "ignoreFailures": true
},
{
  "description": "Unlimited value size",
  "text": "-gdb-set max-value-size unlimited",
  "ignoreFailures": true
}
```

这样以后每次启动调试时，gdb 都会自动应用这些设置。

---

## 8. 实际调试时最常用的几条命令

如果只是想快速查看一个 `vector`，通常最常用的是下面这几条：

```gdb
-exec set print elements 0
-exec p inInfoBits.size()
-exec p inInfoBits[1200]
-exec p inInfoBits[1000]@20
```

可以把它们理解成：

- 看完整元素数量
- 看某个特定位置
- 看某一小段连续内容

---

## 9. 总结

当 VS Code 调试时看不到 `std::vector` 的完整内容，不是程序数据少了，而是调试器默认限制了显示。

最实用的处理方式是：

1. 用 `inInfoBits.size()` 看大小
2. 用 `inInfoBits[idx]` 看单个元素
3. 用 `-exec p inInfoBits[start]@count` 看一段
4. 用 `-exec set print elements 0` 取消元素显示上限
5. 最好把这些设置写进 `launch.json`
