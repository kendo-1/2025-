# poseidon2_hash.py
import sys

# 定义 BN128 曲线的素数
p = 21888242871839275222246405745257275088548364400416034343698204186575808495617


def mod_p(x):
    """在素数域 p 上进行模运算"""
    return x % p


def pow5(x):
    """计算 x^5 mod p，使用高效方法"""
    x2 = mod_p(x * x)
    x4 = mod_p(x2 * x2)
    return mod_p(x4 * x)


def poseidon2_hash(in0, in1):
    """计算 Poseidon2 哈希值 (t=2, d=5)"""
    # 轮常数（与 Circom 电路完全一致）
    round_constants = [
        0x0e7d0d2c7e5c5d57, 0x0c850c0c0c0c0c0c,
        0x0d2d2d2d2d2d2d2d, 0x0a5a5a5a5a5a5a5a,
        0x0b6b6b6b6b6b6b6b, 0x0c7c7c7c7c7c7c7c,
        0x0d8d8d8d8d8d8d8d, 0x0e9e9e9e9e9e9e9e,
        0x0fafafafafafafaf, 0x0101010101010101
    ]

    # MDS 矩阵（与 Circom 电路完全一致）
    mds_matrix = [
        [2, 1],
        [1, 3]
    ]

    # 初始化状态
    state = [mod_p(in0), mod_p(in1)]

    # 输出中间状态（调试用）
    print(f"初始状态: {state}")

    # 轮函数处理 (5轮)
    for r in range(5):
        print(f"\n--- 轮 {r} ---")

        # 1. 添加轮常数
        state[0] = mod_p(state[0] + mod_p(round_constants[r * 2]))
        state[1] = mod_p(state[1] + mod_p(round_constants[r * 2 + 1]))
        print(f"添加轮常数后: {state}")

        # 2. S-box应用
        # 完整轮：所有元素应用S-box (x^5)
        if r < 2 or r >= 3:  # 第0,1,3,4轮
            state[0] = pow5(state[0])
            state[1] = pow5(state[1])
            print(f"完整轮 S-box 后: {state}")
        # 部分轮：仅第一个元素应用S-box
        else:  # 第2轮
            state[0] = pow5(state[0])
            print(f"部分轮 S-box 后: {state}")

        # 3. MDS矩阵混合
        new_state0 = mod_p(mds_matrix[0][0] * state[0] + mds_matrix[0][1] * state[1])
        new_state1 = mod_p(mds_matrix[1][0] * state[0] + mds_matrix[1][1] * state[1])
        state = [new_state0, new_state1]
        print(f"MDS混合后: {state}")

    # 返回第一个状态元素作为哈希值
    return state[0]


if __name__ == "__main__":
    # 解析命令行参数
    if len(sys.argv) != 3:
        print("用法: python poseidon2_hash.py <in0> <in1>")
        print("示例: python poseidon2_hash.py 10 20")
        sys.exit(1)

    try:
        in0 = int(sys.argv[1])
        in1 = int(sys.argv[2])
    except ValueError:
        print("错误: 输入必须是整数")
        sys.exit(1)

    # 计算哈希值
    hash_value = poseidon2_hash(in0, in1)

    # 输出结果
    print("\n" + "=" * 50)
    print(f"输入: [{in0}, {in1}]")
    print(f"Poseidon2 哈希值 (t=2, d=5): {hash_value}")
    print("=" * 50)