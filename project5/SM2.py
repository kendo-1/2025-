import time
import random
from typing import Tuple, List, Dict, Optional
import sys

# SM2曲线参数
P = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
A = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC
B = 0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93
N = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123
Gx = 0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7
Gy = 0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0


class Point:
    """椭圆曲线点类 (支持仿射和雅可比坐标)"""
    __slots__ = ('x', 'y', 'z')  # 优化内存使用

    def __init__(self, x: int, y: int, z: int = 1):
        self.x = x % P if x != 0 else 0
        self.y = y % P if y != 0 else 0
        self.z = z % P if z != 0 else 0

    def to_affine(self) -> Tuple[int, int]:
        """雅可比坐标转仿射坐标"""
        if self.is_infinity():
            return (0, 0)
        z_inv = pow(self.z, P - 2, P)  # 费马小定理求逆
        z_inv_sq = (z_inv * z_inv) % P
        x_aff = (self.x * z_inv_sq) % P
        y_aff = (self.y * z_inv_sq * z_inv) % P
        return (x_aff, y_aff)

    def __eq__(self, other) -> bool:
        if not isinstance(other, Point):
            return False
        if self.is_infinity() and other.is_infinity():
            return True
        if self.is_infinity() or other.is_infinity():
            return False

        # 使用雅可比坐标直接比较（避免转换开销）
        if self.z == 0 or other.z == 0:
            return self.z == 0 and other.z == 0

        # 转换为仿射坐标比较
        x1, y1 = self.to_affine()
        x2, y2 = other.to_affine()
        return x1 == x2 and y1 == y2

    def __str__(self):
        if self.is_infinity():
            return "Point(Infinity)"
        x, y = self.to_affine()
        return f"Point(0x{x:064x}, 0x{y:064x})"

    def is_infinity(self) -> bool:
        return self.z == 0

    def copy(self) -> 'Point':
        return Point(self.x, self.y, self.z)

    def negate(self) -> 'Point':
        """返回该点的负点"""
        if self.is_infinity():
            return self.copy()
        return Point(self.x, P - self.y, self.z)


# 基础大数运算
def mod_add(a: int, b: int, mod: int = P) -> int:
    return (a + b) % mod


def mod_sub(a: int, b: int, mod: int = P) -> int:
    return (a - b) % mod


def mod_mul(a: int, b: int, mod: int = P) -> int:
    """优化模乘 - 使用Python内置大整数运算"""
    return (a * b) % mod


def mod_inv(a: int, mod: int = P) -> int:
    """模逆 (费马小定理)"""
    return pow(a, mod - 2, mod)


# 点运算优化 (雅可比坐标)
def point_double_jacobian(p: Point) -> Point:
    """雅可比坐标点倍乘 (优化实现)"""
    if p.is_infinity() or p.y == 0:
        return Point(0, 0, 0)

    # 雅可比坐标公式 (4M+4S)
    XX = mod_mul(p.x, p.x)
    YY = mod_mul(p.y, p.y)
    YYYY = mod_mul(YY, YY)
    S = mod_mul(4, mod_mul(p.x, YY))
    M = mod_mul(3, XX)

    # 处理曲线参数A
    if A != 0:
        ZZ = mod_mul(p.z, p.z)
        ZZ_sq = mod_mul(ZZ, ZZ)
        M = mod_add(M, mod_mul(A, ZZ_sq))

    x3 = mod_sub(mod_mul(M, M), mod_mul(2, S))
    y3 = mod_sub(mod_mul(M, mod_sub(S, x3)), mod_mul(8, YYYY))
    z3 = mod_mul(2, mod_mul(p.y, p.z))

    return Point(x3, y3, z3)


def point_add_jacobian(p: Point, q: Point) -> Point:
    """雅可比坐标点加 (优化实现)"""
    if p.is_infinity():
        return q.copy()
    if q.is_infinity():
        return p.copy()

    # 自动处理z=1的特殊情况
    if q.z == 1:
        return point_add_jacobian_z1(p, q)

    # 雅可比坐标公式 (12M+4S)
    Z1Z1 = mod_mul(p.z, p.z)
    Z2Z2 = mod_mul(q.z, q.z)
    U1 = mod_mul(p.x, Z2Z2)
    U2 = mod_mul(q.x, Z1Z1)
    S1 = mod_mul(p.y, mod_mul(q.z, Z2Z2))
    S2 = mod_mul(q.y, mod_mul(p.z, Z1Z1))

    H = mod_sub(U2, U1)
    R = mod_sub(S2, S1)

    if H == 0:
        if R == 0:
            return point_double_jacobian(p)
        return Point(0, 0, 0)

    HH = mod_mul(H, H)
    HHH = mod_mul(H, HH)
    V = mod_mul(U1, HH)

    x3 = mod_sub(mod_sub(mod_mul(R, R), HHH), mod_mul(2, V))
    y3 = mod_sub(mod_mul(R, mod_sub(V, x3)), mod_mul(S1, HHH))
    z3 = mod_mul(H, mod_mul(p.z, q.z))

    return Point(x3, y3, z3)


def point_add_jacobian_z1(p: Point, q: Point) -> Point:
    """优化点加 (当q.z=1时)"""
    if p.is_infinity():
        return q.copy()
    if q.is_infinity():
        return p.copy()

    # 优化点加公式 (8M+3S)
    Z1Z1 = mod_mul(p.z, p.z)
    U2 = mod_mul(q.x, Z1Z1)
    S2 = mod_mul(q.y, mod_mul(p.z, Z1Z1))

    H = mod_sub(U2, p.x)
    R = mod_sub(S2, p.y)

    if H == 0:
        if R == 0:
            return point_double_jacobian(p)
        return Point(0, 0, 0)

    HH = mod_mul(H, H)
    HHH = mod_mul(H, HH)
    V = mod_mul(p.x, HH)

    x3 = mod_sub(mod_sub(mod_mul(R, R), HHH), mod_mul(2, V))
    y3 = mod_sub(mod_mul(R, mod_sub(V, x3)), mod_mul(p.y, HHH))
    z3 = mod_mul(H, p.z)

    return Point(x3, y3, z3)


# 点乘算法优化
def point_mul_binary(k: int, p: Point) -> Point:
    """二进制点乘算法 (优化实现)"""
    # 处理特殊情况
    if k == 0 or p.is_infinity():
        return Point(0, 0, 0)
    if k == 1:
        return p.copy()

    result = Point(0, 0, 0)  # 无穷远点
    temp = p.copy()

    # 二进制展开 (从最低位开始)
    while k > 0:
        if k & 1:
            result = point_add_jacobian(result, temp)
        temp = point_double_jacobian(temp)
        k >>= 1

    return result


def precompute_fixed_point(p: Point, window_size: int = 8) -> List[Point]:
    """预计算固定点表 (优化实现)"""
    table_size = 1 << window_size
    table = [Point(0, 0, 0)] * table_size

    # 基础点
    table[0] = Point(0, 0, 0)  # 无穷远点
    table[1] = p.copy()  # 1倍点

    # 计算2^i倍点
    current = p.copy()
    for i in range(2, table_size):
        # 直接累加基点
        table[i] = point_add_jacobian(table[i - 1], table[1])

    return table


def point_mul_fixed_precompute(k: int, precomputed: List[Point], window_size: int = 8) -> Point:
    """固定点预计算点乘 (优化实现)"""
    # 处理特殊情况
    if k == 0:
        return Point(0, 0, 0)

    result = Point(0, 0, 0)
    num_bits = 256  # SM2曲线是256位

    # 从最高位开始处理
    for i in range(num_bits - 1, -1, -1):
        result = point_double_jacobian(result)

        # 检查当前位
        if (k >> i) & 1:
            # 找到对应的预计算点
            result = point_add_jacobian(result, precomputed[1])

    return result


def naf_encode(k: int, w: int = 5) -> List[int]:
    """NAF编码 (非相邻形式) - 优化实现"""
    naf = []
    half_w = 1 << (w - 1)
    mask = (1 << w) - 1

    while k > 0:
        if k & 1:
            # 取低w位
            digit = k & mask
            if digit > half_w:
                digit -= 1 << w
            k -= digit
            naf.append(digit)
        else:
            naf.append(0)
        k //= 2

    return naf[::-1]


def point_mul_naf(k: int, p: Point, w: int = 5) -> Point:
    """NAF点乘算法 (优化实现)"""
    # 生成NAF编码
    naf_rep = naf_encode(k, w)

    # 主循环
    result = Point(0, 0, 0)
    for digit in naf_rep:
        result = point_double_jacobian(result)
        if digit > 0:
            # 直接计算正倍数
            temp = p.copy()
            for _ in range(digit - 1):
                temp = point_add_jacobian(temp, p)
            result = point_add_jacobian(result, temp)
        elif digit < 0:
            # 计算负倍数
            abs_digit = -digit
            temp = p.copy()
            for _ in range(abs_digit - 1):
                temp = point_add_jacobian(temp, p)
            result = point_add_jacobian(result, temp.negate())

    return result


# 性能测试框架
def benchmark(func, *args, runs: int = 100, warmup: int = 10, name: str = "") -> float:
    """优化的基准测试函数"""
    # 预热
    for _ in range(warmup):
        func(*args)

    # 精确计时
    min_time = float('inf')
    for _ in range(3):  # 多次测试取最小值
        start = time.perf_counter()
        for _ in range(runs):
            result = func(*args)
            # 验证结果不是无穷远点
            if isinstance(result, Point) and result.is_infinity():
                raise ValueError(f"{name} produced infinity point")
        elapsed = (time.perf_counter() - start) / runs
        min_time = min(min_time, elapsed)

    return min_time


def compare_optimizations():
    """比较不同优化技术的性能 (改进版)"""
    # 准备测试点
    G = Point(Gx, Gy)
    k = random.randint(1, N - 1)
    P_rand = point_mul_binary(random.randint(1, N - 1), G)

    # 固定点预计算表 (w=8)
    window_size_fixed = 8
    precomputed_fixed = precompute_fixed_point(G, window_size_fixed)

    # 测试不同点乘算法
    results = {}

    # 固定点乘
    results["固定点-基础二进制"] = benchmark(point_mul_binary, k, G, runs=50, name="point_mul_binary (fixed)")
    results["固定点-预计算(w=8)"] = benchmark(
        point_mul_fixed_precompute, k, precomputed_fixed, window_size_fixed,
        runs=100, name="point_mul_fixed_precompute"
    )

    # 非固定点乘
    results["非固定点-基础二进制"] = benchmark(point_mul_binary, k, P_rand, runs=50, name="point_mul_binary (random)")
    results["非固定点-NAF(w=5)"] = benchmark(point_mul_naf, k, P_rand, 5, runs=100, name="point_mul_naf")

    # 点运算优化
    G_z1 = Point(Gx, Gy, 1)  # z=1的点
    results["点加-雅可比坐标"] = benchmark(
        point_add_jacobian, G, G, runs=1000, name="point_add_jacobian"
    )
    results["点加-雅可比Z1优化"] = benchmark(
        point_add_jacobian_z1, G, G_z1, runs=1000, name="point_add_jacobian_z1"
    )

    # 模运算优化 (批量测试)
    a = random.randint(1, P - 1)
    b = random.randint(1, P - 1)
    results["模乘-基础(1000次)"] = benchmark(
        lambda: mod_mul(a, b), runs=1000, name="mod_mul"
    )

    # 打印结果
    print("\nSM2算法优化性能对比 (单次操作时间, 秒):")
    print("-" * 75)
    print(f"{'优化技术':<25} | {'操作':<30} | {'时间':<15} | {'加速比'}")
    print("-" * 75)

    # 计算并显示加速比
    baseline_point_mul_fixed = results.get("固定点-基础二进制", 1e-9)
    baseline_point_mul_var = results.get("非固定点-基础二进制", 1e-9)
    baseline_point_add = results.get("点加-雅可比坐标", 1e-9)
    baseline_mod_mul = results.get("模乘-基础(1000次)", 1e-9)

    # 设置参考基准
    ref_fixed = baseline_point_mul_fixed
    ref_var = baseline_point_mul_var
    ref_point_add = baseline_point_add
    ref_mod_mul = baseline_mod_mul

    for name, time_val in results.items():
        if time_val <= 0:
            time_val = 1e-9

        if "固定点" in name:
            speedup = ref_fixed / time_val
            op = "点乘 (固定点)"
        elif "非固定点" in name:
            speedup = ref_var / time_val
            op = "点乘 (随机点)"
        elif "点加" in name:
            speedup = ref_point_add / time_val
            op = "点加"
        elif "模乘" in name:
            speedup = ref_mod_mul / time_val
            op = "模乘 (1000次平均)"
        else:
            op = "其他"
            speedup = 1.0

        print(f"{name:<25} | {op:<30} | {time_val:.8f} | {speedup:.2f}x")

    # 总结关键优化效果
    print("\n关键优化技术效果总结:")
    if "固定点-预计算(w=8)" in results:
        ratio = ref_fixed / results["固定点-预计算(w=8)"]
        print(f"固定点预计算 (w=8) 加速比: {ratio:.2f}x")

    if "非固定点-NAF(w=5)" in results:
        ratio = ref_var / results["非固定点-NAF(w=5)"]
        print(f"NAF编码 (w=5) 加速比: {ratio:.2f}x")

    if "点加-雅可比Z1优化" in results:
        ratio = ref_point_add / results["点加-雅可比Z1优化"]
        print(f"雅可比Z1优化加速比: {ratio:.2f}x")


def verify_correctness():
    """验证算法正确性"""
    G = Point(Gx, Gy)

    # 验证点运算
    G2_double = point_double_jacobian(G)
    G2_add = point_add_jacobian(G, G)

    print("正确性验证:")
    print(f"2G (点倍): {G2_double.to_affine()}")
    print(f"2G (点加): {G2_add.to_affine()}")
    print(f"点倍 == 点加: {G2_double == G2_add}")

    # 验证点乘一致性
    k = 0x123456789ABCDEF
    k2 = 0x9876543210FEDCBA

    # 二进制点乘
    P1_bin = point_mul_binary(k, G)
    P2_bin = point_mul_binary(k2, G)

    # 预计算点乘
    precomputed = precompute_fixed_point(G, 8)
    P1_pre = point_mul_fixed_precompute(k, precomputed, 8)
    P2_pre = point_mul_fixed_precompute(k2, precomputed, 8)

    # NAF点乘
    P1_naf = point_mul_naf(k, G, 5)
    P2_naf = point_mul_naf(k2, G, 5)

    # 验证结果一致性
    print(f"二进制点乘 == 预计算点乘 (k1): {P1_bin == P1_pre}")
    print(f"二进制点乘 == 预计算点乘 (k2): {P2_bin == P2_pre}")
    print(f"二进制点乘 == NAF点乘 (k1): {P1_bin == P1_naf}")
    print(f"二进制点乘 == NAF点乘 (k2): {P2_bin == P2_naf}")

    # 验证点加Z1优化
    G_z1 = Point(Gx, Gy, 1)
    P_add_normal = point_add_jacobian(G, G)
    P_add_z1 = point_add_jacobian_z1(G, G_z1)
    print(f"点加 (正常) == 点加 (Z1优化): {P_add_normal == P_add_z1}")


if __name__ == "__main__":
    print(f"Python版本: {sys.version}")
    print(f"SM2优化测试 | 曲线阶数: {N.bit_length()}位")

    # 验证算法正确性
    verify_correctness()

    # 运行性能对比
    compare_optimizations()