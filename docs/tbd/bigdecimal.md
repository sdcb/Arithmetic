# BigDecimal（待定）

> 状态：讨论中，尚未决定是否实现。

## 是否有必要

结论：**有必要，优先级高于 Float128**。它几乎完全复用 `BigInt`，成本低，填补的是 .NET 生态里真实存在的空缺。

| | 精确十进制 | 精度 |
|---|---|---|
| `decimal`（.NET 10 及以前） | ✓ | 固定 28–29 位 |
| `Decimal128`（.NET 11，IEEE 754 decimal128） | ✓ | 固定 34 位 |
| `Real` | ✗（二进制，0.1 不精确） | 任意 |
| `BigRational` | ✓ | 任意，但没有"小数位数"概念，1/3 这种结果会让分母不断膨胀 |
| **`BigDecimal`** | **✓** | **任意** |

真实需求：

- 金融/记账：超过 28 位，或需要按指定小数位做舍入（HalfUp、HalfEven 等）。
- 数据库 `DECIMAL(38..65)`（Oracle NUMBER、PostgreSQL numeric）。
- 区块链金额（uint256 + 18 位小数，约 78 位数字）。
- JSON/协议中超出 `decimal` 范围的数字需要原样往返。

## 实现成本

- 表示：`BigInt` 未缩放值 + `int` 小数位数（scale）。加、减、乘、比较直接调用 `BigInt`。
- 除法/舍入：`BigInt.DivRem` + 十进制舍入规则（HalfUp、HalfEven、HalfDown、Up、Down、Ceiling、Floor，以及"不精确就抛异常"）。
- 乘以 10^k 复用已缓存的幂表（`PowerTable`）；Parse/ToString 复用现有的分治实现。
- 与 `decimal`、`double`、`BigRational`、`Real`、`Decimal128` 的互转都是现成能力的组合。

## 建议的 MVP 范围

1. 不可变的 `BigDecimal` 结构体。
2. 精确的加、减、乘。
3. 除法和 sqrt 必须显式指定上下文：有效位数或小数位数 + 舍入模式（类似 Java 的 `MathContext`）。
4. `Round`、`SetScale`、`StripTrailingZeros`；普通写法与科学计数法的 `ToString` / `Parse`。
5. 接入泛型数学（`INumber<BigDecimal>`），与现有三种类型互转。
6. **相等语义按值比较**，即 `1.0 == 1.00`，避开 Java `equals` 同时比较 scale 的老坑；scale 单独暴露为属性。
7. 超越函数（exp/log/sin 等）暂不做（Java 到 9 才加 sqrt）。以后可以用 `Real` 计算再按十进制正确舍入，但要处理二进制转十进制恰好落在舍入边界的情况，需要 Ziv 式的方法，不简单。

## 对照基准

### 正确性

| 对照对象 | 用途 | 说明 |
|---|---|---|
| **General Decimal Arithmetic 规范的 `.decTest` 测试集**（Mike Cowlishaw / IBM） | **主要语义测试** | 十进制算术事实上的标准测试集，覆盖加、减、乘、除、quantize、各种舍入、精度上下文、特殊值，每条用例都有期望结果。CPython `decimal`、IBM decNumber、libmpdec 都用它。拷进仓库前需确认许可证。 |
| **libmpdec / mpdecimal**（C，BSD，CPython `decimal` 的内核） | 逐位比对 | 类似现有的 `gmp-exact`：P/Invoke 跑随机输入，比较系数、指数和舍入结果。任意精度、正确舍入。 |
| **`BigRational`**（库内现成） | 进程内精确对照 | BigDecimal 的值都是有理数：先用 `BigRational` 精确计算，再按十进制规则舍入，与结果比较。成本低，适合放进 xunit 和 CI。 |
| **Java `BigDecimal`** | 对齐 API 语义 | `setScale`、`stripTrailingZeros`、`MathContext` 等行为。可以在 CI 里用 setup-java 生成测试向量文件，比较输出字符串，不需要跨进程调用。 |
| **`decimal`**（.NET 10） | 小数值范围冒烟测试 | 28 位以内比对。注意它的舍入和缩放规则有自己的特点，不是规范定义的语义。 |
| **`Decimal128`**（.NET 11） | 34 位精度比对 | IEEE 754 decimal128 等价于规范中"精度 34、emax 6144"的上下文，可以直接比。.NET 11 发布后再加，或用条件编译先接 preview。 |

### 性能

| 对照对象 | 意义 |
|---|---|
| **libmpdec（原生）** | **主要标杆**，作用相当于 GMP。复用 `perf` 表格格式：我们的耗时 / 原生耗时。 |
| **Java `BigDecimal`（JMH）** | 业界最常见的参照。跨运行时比较有噪声，但用户最关心这个数字。 |
| **.NET 社区库**（例如 ExtendedNumerics.BigDecimal、Singulink.Numerics.BigDecimal，包名待核实） | 大多基于 `BigInteger`，放进 BenchmarkDotNet 做同运行时对比，类似现有的 BigInteger 那一列。 |
| **`decimal` / `Decimal128`** | 只比小数值、固定精度的场景，说明"小数值时慢多少"。 |

## 待定的架构问题：系数表示

| 方案 | 优点 | 缺点 |
|---|---|---|
| **二进制系数**（Java 的做法，复用 `BigInt`） | 加、减、乘很快，直接用上 AVX2 内核和 NTT | "按精度 p 舍入"要先知道十进制位数，再除以 10^k，是一次真正的大数除法；BigDecimal 每次运算都要舍入，会成为热点 |
| **以 10^19 为基的 limb**（libmpdec 的做法） | 舍入、取十进制位数、移小数点都是 O(1) 或按 limb 截断 | 乘法要另写一套（libmpdec 用 NTT），加法进位规则不同，复用不了 `BigInt` |

倾向：**先用二进制系数**，实现成本最低。缓解办法：

- 十进制位数用"bit 长度 × log10(2) 估算，再查一次 10^k 表修正"。
- 除以 10^k 复用 `PowerTable`，必要时加倒数缓存。

再用 libmpdec 的 benchmark 判断舍入开销是否可以接受，不行再考虑换表示。

## 如果开搞，第一步

先只搭对照框架，再写 BigDecimal 本体：

1. libmpdec 的 P/Invoke 接入（需要原生运行时包或自行编译）。
2. `.decTest` 解析器。
3. 跑通 `BigRational` 精确对照。
