# 模残差格基复原

计量室从一台旧扫描器拿到八至十组已知顺序的二维系数和两到三个质数模下的线性残差。原始整数格基没有随票据保存，只剩这些同余观测。审计程序要从观测恢复格，再输出唯一方向规范基。

命令入口固定为 `/app/lattice_audit.py`，只接收输入和结果两个路径：

```text
python3 /app/lattice_audit.py --input <输入路径> --output <结果路径>
```

## 输入票据

票据最外层只保留以下五个键：

```text
schema_version, case_id, coefficients, moduli, limits
```

`schema_version` 必须是整数 `1`。`case_id` 是非空字符串。

`coefficients` 恰好有 6 个有序的二维整数系数，每个成员恰好有两个元素。两个分量都必须位于闭区间 `[-2,2]`，且不能同时为零。第 `i` 个系数与每个模数对象中的第 `i` 条残差相对应。

`moduli` 恰好有 2 个模数视角。每个对象的键恰好为：

```text
modulus, residuals
```

`modulus` 必须是互不相同的质数。`residuals` 的长度必须等于系数数；第 `i` 项恰好是两个整数 `[left, right]`，都在 `[0, modulus-1]` 内。

设待恢复基的列向量为 `b1`、`b2`，第 `i` 个系数为 `(z1,z2)`。每条残差满足：

```text
left  = z1*b1[0] + z2*b1[1]  (mod modulus)
right = z1*b2[0] + z2*b2[1]  (mod modulus)
```

`limits` 只保留八个规模与界值字段：

```text
min_observations, max_observations, max_dimension, min_moduli,
max_moduli, max_coefficient, max_basis_entry, max_determinant
```

这些值依次固定为 `6 / 6 / 2 / 2 / 2 / 2 / 8 / 32`。边界值本身允许。

## 规范基

恢复出的两个基向量必须使格为 primitive：四个基元素与行列式绝对值的最大公因数为 `1`。四个基元素都位于闭区间 `[-max_basis_entry,max_basis_entry]`。行列式必须大于零且不超过 `max_determinant`。

所有基向量分量都是整数。报告采用二维 Gauss 约化后的唯一方向规范基：

1. `det(b1,b2)` 必须为正；
2. `b1[0]` 必须为正；
3. Gauss 约化后 `b2[1]` 必须为正；
4. 若有多个候选，先取最短基向量的平方范数最小值；
5. 平方范数相同时，再按 `b1`、`b2` 的字典序决胜。

## 输出与失败

结果最外层只保留四个键：

```text
schema_version, case_id, basis, check
```

`schema_version` 和 `case_id` 原样沿用输入。`basis` 必须是 `[[b1x,b1y],[b2x,b2y]]`，四个元素都是 JSON 整数。

`check` 只保留三个字段：

```text
observation_count, modulus_count, satisfied_equations
```

`observation_count` 等于系数数，`modulus_count` 等于模数对象数。`satisfied_equations` 等于 `2*observation_count*modulus_count`，表示每个系数在两个基向量上的全部同余式均被验证。

结果使用无空格 JSON，键名升序，文件只以一个换行收尾。成功完成时 `stdout` 与 `stderr` 都必须为空字符串。

坏票据、不可解同余系统、越界基和方向规范不唯一统一返回状态码 `3`，控制台保持沉默。失败不得改写既有结果。结果经同目录临时文件安装；若最终换名失败，既有文件按旧内容保留。
