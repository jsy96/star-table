# 卫星在轨标准星表与视方向计算（星敏感器姿态确定支持）

从 Gaia DR3 生成的**标准星表**（定长二进制，C/C++ 直接 `fread`）+ 一个
**C/C++ 函数**：输入观测时刻（UTC）与卫星状态矢量（GCRS 或 **ECEF/ITRS**，后者经
IERS EOP C04 极移+DUT1 转换），输出该时刻该位置的整表恒星**视方向**（自行传播、
周年视差、太阳/地球引力光线偏折、光行差全部改正）。算法链路与 IAU SOFA `atciq`
完全一致（逐微角秒校验通过），但不施加 bias-岁差-章动旋转——该共模旋转由姿态
解算吸收，逐星物理改正才需要逐星施加。EOP 数据由 `update_eop.py` 流水线自动
下载、更新、重编译、全量测试并生成 Markdown 报告。

## 更新记录（最新在前）

- **2026-11-15 v1.1.261115-5**：快速产品实装落地 + **一键脚本**。
  - `run_eop_pipeline.bat`：双击即可执行全流程（下载→快速合并→编译→测试→报告；
    自动探测 `eop/finals.daily`，纯英文脚本避免 GBK 乱码）；
  - finals.daily 解析改**固定列宽切片**（月日在日<10 时 token 数漂移、标志字母
    与数值粘连 `I-0.0581`/`P-0.159`，按列切才可靠），预测行 LOD 缺省置 0；
  - **实测合并效果**：181 行快速数据并入 C04 之后至 2027-02-13，有效滞后 0 天
    （预测段覆盖未来 89 天）；重叠段 30 天实测 |快速−最终|：极移 ≤0.13 mas、
    DUT1 ≤0.059 ms——快速产品真实精度实测确认（此前 30 ms 假差异系 token
    错位）；DUT1 拼接偏置校正自动拟合（本次 a=−0.015 ms ≈ 0，残差 0.022 ms）。
- **2026-11-15 v1.1.261115-4**：EOP 策略升级为**双模自动**——快速产品
  （USNO `finals.daily`）可达时：C04 历史基底 + 快速产品补充末端（合并为同名
  同格式文件，**C 库零改动**），重叠段保留 C04 最终值，滞后从 ~30 天降至 ~1 天；
  不可达时自动退化为纯 C04 + 外推（v3 策略）。依据（实测对比）：C04 外推 31 天
  残差 DUT1 ≤21 ms / 极移 ≤44 mas / 旋转 ≤0.32″，而快速产品典型 ~0.02–0.05 ms /
  ~0.2–0.5 mas，优 1~2 个数量级。手动启用：下载
  `https://maia.usno.navy.mil/ser7/finals.daily` 放 `eop/`，加参数
  `--input-rapid eop/finals.daily`。报告新增合并详情与重叠段 |快速−最终| 实测节。
- **2026-11-15 v1.1.261115-3**：EOP 策略定版——**固定使用长期产品 C04，不使用 IERS
  快速产品**。数据末端越界时 `ac_eop_get` 自动外推（DUT1 按 dut1−LOD·Δ 解析关系、
  极移按末端 30 天线性趋势，dX/dY 端点保持），仍返回 `AC_W_EOP_RANGE` 警告；
  `update_eop.py` 新增**留出回测**（C04 自身历史 6 截断点×逐日评估），实测：
  31 天跨度 DUT1 外推残差 ≤21 ms（端点保持 30 ms）、极移 ≤0.05″，折合地固旋转
  ≤0.32″；C/Python 外推逻辑同步。修复报告生成器变量复用与历史条目误捞两处 bug。
- **2026-11-15 v1.1.261115-2**：接入 **IERS EOP C04（IAU2000）**。
  - 新增 `update_eop.py` 流水线：**下载（hpiers 主源+存档回退）→ 更新（原子替换+备份）→
    编译（vswhere 定位 MSVC，优先 v145）→ 测试 → 生成 `EOP_report.md`**；
  - C 库新增 `c/eop_c04.c`：EOP C04 解析/插值、`ac_itrs_to_gcrs()`/`ac_gcrs_to_itrs()`
    （EOP 极移+DUT1+IAU2006/2000A CIO 岁差章动，erfa c2t06a 链路）、
    `ac_apparent_catalog_itrs()`（**ECEF 状态矢量直接输入**）；
  - 测试：T1–T7 全过（解析/插值 1e-12/越界警告/往返 1e-9/极移效应量级/整表）；
    ITRS→GCRS 旋转 C vs Python 交叉验证 max|Δr|=2.0e-5 km、max|Δv|=2.1e-8 km/s PASS；
  - 修复：测试程序 T5 浅拷贝清零极移污染后续用例（C 库本身无恙）。
- **2026-11-15 v1.0.261115-1**：初始交付。
  - `build_std_catalog.py`：Gaia DR3（J2016.0）→ 标准星表 `.bin`(C) + `.npz`(Python)，
    含 `pmsafe` 严格历元传播 J2016.0→J2000.0；G≤9 共 177,426 星；
    Barnard(10.39″/yr)、Kapteyn、Groombridge 1830、比邻星(π=768 mas) 等极端星确认在表。
  - C/C++ 库 `c/apparent_catalog.{h,c}`（零第三方依赖，内置 ERFA 子集）。
  - 验证：Python 链路 vs SOFA `apci13+atciq` 最大差 **0.000 μas**；
    C vs Python 逐星交叉验证（2000 星）**全部 <0.0001″ PASS**；
    LEO 示例统计与理论预测吻合（自行+视差中位数 0.1776″ vs 理论 0.1783″，
    最大 112.90″ = Barnard 星 10.87 年自行）。
  - 修复记录：ERFA 视差单位约定（arcsec 非 mas）、ERFA 赤经自行约定
    （dRA/dt，需除 cosδ）——两处约定错误曾在统计中表现为 12%~10³ 倍偏差。

## 文件结构

```
louyifan/
├─ gaia_dr3_13.csv            # 源数据: Gaia DR3 摘要 7,369,632 星 (G≤13, J2016.0)
├─ build_std_catalog.py       # 标准星表生成器 (uv 运行)
├─ apparent_catalog.py        # Python 参考实现 + 自检 + ITRS/EOP 支持
├─ validate_c_vs_python.py    # C ↔ Python 逐星交叉验证
├─ update_eop.py              # IERS EOP C04 流水线: 下载→更新→编译→测试→报告
├─ EOP_report.md              # 流水线生成的最新报告
├─ std_catalog_G9p0.bin       # 生成的标准星表 (C/C++, 80B 头 + 80B/星)
├─ std_catalog_G9p0.npz       # 生成的标准星表 (Python)
├─ eop/
│  └─ eopc04_IAU2000.62-now   # IERS EOP C04 数据 (update_eop.py 维护, 旧版 .prev)
└─ c/
   ├─ apparent_catalog.h      # C/C++ 公开 API (extern "C"), 含 EOP/ITRS 接口
   ├─ apparent_catalog.c      # 视方向核心实现
   ├─ eop_c04.c               # EOP C04 加载/插值 + ITRS<->GCRS + ECEF 直入接口
   ├─ test_apparent_catalog.c # 星表链路测试 (统计与 CSV 导出)
   ├─ test_eop.c              # EOP 测试 T1-T7 (含旋转样本导出)
   ├─ CMakeLists.txt          # MSVC /utf-8, 工具集 v145
   └─ third_party/erfa/       # ERFA 子集 (LGPL-3.0, 原样收录)
```

## 快速开始

### 1. 生成标准星表（只需一次，或换星等截断时重跑）

```bash
cd louyifan
uv run --no-project --with numpy --with pandas --with pyerfa \
    python -I build_std_catalog.py --input gaia_dr3_13.csv --gmax 9.0
```

### 2. 编译 C 库（MSVC v145，经 vswhere 自动定位）

```bash
cd c
cmake -S . -B build -T v145
cmake --build build --config Release
# 测试 + 示例输出（LEO 550 km 示例状态）
./build/Release/test_apparent_catalog.exe ../std_catalog_G9p0.bin \
    2026 11 15 12.0  6928.137 0 0  0 -1.0065 7.5219  2000
```

### 3. C/C++ 调用

```c
#include "apparent_catalog.h"

ac_catalog_t cat;
if (ac_load_catalog("std_catalog_G9p0.bin", &cat) != AC_OK) { /* ... */ }

double jd = ac_utc_to_jd(2026, 11, 15, 12.0);        /* UTC → JD */

/* 方式 A: GCRS 惯性系状态矢量 */
double r[3] = {6928.137, 0, 0};                      /* GCRS, km  */
double v[3] = {0, -1.0065, 7.5219};                  /* GCRS, km/s */
static double ra[N], dec[N];
static ac_breakdown_t brk[N];
ac_apparent_catalog(&cat, jd, r, v, AC_DEFL_SUN_EARTH, ra, dec, NULL, brk);

/* 方式 B: ECEF/ITRS 状态矢量 + IERS EOP（推荐, OD 产品常为地固系） */
ac_eop_t eop;
ac_eop_load("eop/eopc04_IAU2000.62-now", &eop);
double ri[3] = { /* GPS/ECEF 位置 km */ }, vi[3] = { /* ECEF 速度 km/s */ };
ac_apparent_catalog_itrs(&cat, &eop, jd, ri, vi,
                         AC_DEFL_SUN_EARTH, ra, dec, NULL, brk);
/* 或仅做坐标转换: ac_itrs_to_gcrs() / ac_gcrs_to_itrs() */

ac_free_catalog(&cat);
ac_eop_free(&eop);
```

性能（Release, 单线程）：整表 177,426 星约 67–132 ms，≈0.4–0.7 μs/星。

### 3b. IERS EOP 流水线（下载/更新/编译/测试/报告）

**一键执行（推荐，双击或命令行）**：`run_eop_pipeline.bat`
（自动探测 `eop\finals.daily`，存在即合并快速产品；参数透传给 update_eop.py，
如 `run_eop_pipeline.bat --skip-download`；脚本纯英文防 GBK 乱码）。

手动方式：

```bash
uv run --no-project --with numpy --with pyerfa python -I update_eop.py
#   可选: --skip-download / --skip-build / --skip-test
#         --input 本地C04文件 / --input-rapid eop/finals.daily
#   产物: eop/eopc04_IAU2000.62-now (旧版自动备份 .prev) + EOP_report.md
```

数据源：IERS EOP 产品中心（巴黎天文台 hpiers，主源+存档回退）。**固定使用长期
产品 C04**：约月更、末端滞后约 30 天属正常，滞期间自动外推（策略见"注意事项"
第 2 条），外推精度留出回测随报告输出。每次运行自动重编译、全量测试
（含 C↔Python 旋转交叉验证与外推回测）并刷新报告。

### 4. Python 参考实现（含自检与 ITRS 支持）

```bash
uv run --no-project --with numpy --with pyerfa \
    python -I apparent_catalog.py --catalog std_catalog_G9p0.npz
```

```python
import apparent_catalog as ac
cat = ac.load_catalog("std_catalog_G9p0.npz")
out = ac.apparent_catalog(cat, "2026-11-15 12:00:00",
                          r_sat=[6928.137, 0, 0], v_sat=[0, -1.0065, 7.5219],
                          include_breakdown=True)   # 各项改正量
```

## 标准星表格式

`.bin`（小端）＝ 80 字节文件头 + N×80 字节记录：

| 头字段 | 类型 | 说明 |
|---|---|---|
| magic | char[8] | `"STDCAT01"` |
| format_version | int32 | 1 |
| n_stars | int32 | 星数 |
| g_max | double | 生成时 G 截断 |
| jd_epoch_source | double | 2457389.0（Gaia J2016.0） |
| jd_epoch_propagated | double | 2451545.0（J2000.0） |
| reserved | int32[10] | 0 |

记录（与 `ac_star_record_t` 一一对应，无填充）：
`source_id`(i8), `ra_j2000`/`dec_j2000`(rad), `pmra_j2000`(**dRA/dt**, rad/yr,
不含 cosδ——ERFA 约定), `pmdec_j2000`(rad/yr), `parallax_as`(**arcsec**——ERFA 约定),
`rv_kms`(km/s, 缺测=0), `g_mag`, `ra_j2016_deg`/`dec_j2016_deg`(原始 Gaia 值,
误差分析基准)。

`.npz` 字段同名（`ra_j2000`…`dec_j2016_deg`），供 `apparent_catalog.py` 使用。

## 物理模型与精度

链路（每星）：`eraPmpx`（自行+视差，含 Roemer 光行时近似）→ `eraLd`
（太阳+地球引力偏折，dlim 与 SOFA `ldsun` 一致）→ `eraAb`（光行差，
观测者质心系总速度＝地球公转+卫星轨道）。
星历/时间：`eraEpv00`（1900–2100，速度精度 42 mm/s → 光行差误差 ~29 μas）、
`eraDat` 闰秒表（截至 2017-01-01 的 37 s，此后按 37 s 外推）、UTC→TT。
GCRS≈BCRS（差 <0.1 μas）；TT≈TDB（差 <0.1 μas）。

改正量量级（G≤9 星表，2026-11-15 LEO 550 km 示例，全表统计）：

| 项 | 均值 | 中位数 | P99 | 最大 |
|---|---|---|---|---|
| 自行+视差 | 0.3616″ | 0.1776″ | 2.9532″ | 112.90″(Barnard) |
| 引力偏折 | 6.48 mas | 4.43 mas | 36.3 mas | 2168 mas |
| 光行差 | 18.16″ | 19.79″ | 22.17″ | 22.17″ |
| 总改正 | 18.21″ | 19.81″ | 22.77″ | 126.11″ |

## 验证记录（最新在前）

- 2026-11-15：C ↔ Python 交叉验证 2000 星：max|Δra·cosδ|=max|Δdec| <0.0001″，PASS；
  Python 链路 vs SOFA `apci13+atciq`（地心、仅太阳偏折，bpn 逆旋转后比对）：
  max|Δ|=0.000 μas；C 测试程序全表统计与 Gaia 星表直接统计的理论预测一致
  （自行+视差 median 0.1776″ vs 0.1783″，P99 2.953″ vs 2.953″，max 112.90″ vs 112.98″）。

## 注意事项与已知限制

1. **输入坐标系**：`ac_apparent_catalog` 接受 GCRS（≈J2000 惯性系）状态矢量；
   ECEF/ITRS 状态矢量请用 `ac_apparent_catalog_itrs`（需 EOP，极移+DUT1 自动插值）。
2. **EOP 双模策略**：优先**快速产品**（USNO `finals.daily`，滞后 ~1 天，精度
   ~0.02–0.05 ms / ~0.2–0.5 mas），合并进 C04 历史基底（同名同格式，C 库零改动）；
   快速源不可达（如当前网络）自动退化为纯 C04 + 外推并返回 `AC_W_EOP_RANGE` 警告
   ——DUT1 按 `dut1(t0+Δ)=dut1(t0)−LOD·Δ`（解析关系）、极移按末端 30 天线性趋势。
   退化模式留出回测（`EOP_report.md`）实测 31 天跨度残差：DUT1 ≤21 ms（端点保持
   30 ms）、极移 ≤0.05″、折合地固旋转 ≤0.32″——比快速产品差 1~2 个数量级，
   网络恢复或手动提供 finals.daily 后自动升级。
3. **时间范围**：1900–2100（`eraEpv00` 适用域）；新增闰秒后需更新 ERFA `dat.c` 表。
4. **不做岁差章动旋转**：恒星视方向输出为 GCRS；`ac_itrs_to_gcrs` 内部使用
   IAU2006/2000A 完整岁差章动（这是坐标转换所必需的，与视方向"不旋转"不矛盾）。
5. **负视差星**：`pmsafe` 对视差 ≤0.1 mas 的 ~2000 星做了距离保护（ERFA 标准行为，
   μas 级影响）。
6. **视向速度缺测**置 0（丢失的透视加速度 μas 级）。
7. **许可**：`c/third_party/erfa/` 为 LGPL-3.0（ERFA，源自 IAU SOFA），原样收录并保留
   许可文本；对外发布产品时需遵守 LGPL 条款。
8. 星表二进制为小端格式（x86/ARM-LE 通用；大端平台需转换）。
