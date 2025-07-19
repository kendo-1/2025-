import hashlib
import random
from cryptography.hazmat.primitives.asymmetric import ec
from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat, load_der_public_key
from cryptography.hazmat.backends import default_backend
from phe import paillier
import binascii
from ecdsa import NIST256p, ellipticcurve
from ecdsa.util import number_to_string, string_to_number

# 预定义常见椭圆曲线的阶
CURVE_ORDERS = {
    "secp256r1": 115792089210356248762697446949407573529996955224135760342422259061068512044369
}


class PrivateIntersectionSum:
    def __init__(self, curve_name="secp256r1"):
        # 创建曲线对象
        self.curve = ec.SECP256R1()
        self.order = CURVE_ORDERS[curve_name]

        # 为ecdsa库创建曲线对象
        self.ecdsa_curve = NIST256p.curve

        # 生成随机私钥标量 (k1 或 k2)
        self.private_scalar = random.randint(1, self.order - 1)

        # 用于存储中间结果
        self.identifier_map = {}
        self.point_bytes_map = {}

    def hash_to_point(self, identifier):
        """将标识符哈希到椭圆曲线点"""
        hash_func = hashlib.sha256()
        hash_func.update(identifier.encode())
        hash_bytes = hash_func.digest()

        # 直接使用哈希值作为私钥值
        private_value = int.from_bytes(hash_bytes, 'big') % self.order
        if private_value == 0:
            private_value = 1  # 避免零值

        private_key = ec.derive_private_key(
            private_value,
            self.curve,
            default_backend()
        )
        return private_key.public_key()

    def point_to_bytes(self, point):
        """将点转换为DER编码的字节串"""
        return point.public_bytes(
            encoding=Encoding.DER,
            format=PublicFormat.SubjectPublicKeyInfo
        )

    def bytes_to_point(self, point_bytes):
        """将字节串转换回椭圆曲线点"""
        return load_der_public_key(point_bytes, default_backend())

    def extract_point_coordinates(self, public_key):
        """从公钥中提取坐标 (x, y)"""
        public_numbers = public_key.public_numbers()
        return public_numbers.x, public_numbers.y

    def create_ecdsa_point(self, x, y):
        """从坐标创建ecdsa点对象"""
        return ellipticcurve.Point(self.ecdsa_curve, x, y)

    def scalar_multiply_point(self, point, scalar):
        """使用ecdsa库实现正确的标量乘法：scalar * point"""
        # 提取点的坐标
        x, y = self.extract_point_coordinates(point)

        # 创建ecdsa点对象
        ecdsa_point = self.create_ecdsa_point(x, y)

        # 执行标量乘法
        result_point = scalar * ecdsa_point

        # 将结果点转换回cryptography公钥
        x_int = result_point.x()
        y_int = result_point.y()

        # 创建公钥数字
        public_numbers = ec.EllipticCurvePublicNumbers(
            x_int, y_int, self.curve
        )

        return public_numbers.public_key(default_backend())

    def protocol_p1_step1(self, identifiers):
        """P1 - 步骤1: 哈希标识符并存储映射"""
        self.original_points = []
        self.identifier_map = {}
        self.raw_points = {}  # 存储原始点用于调试

        for id in identifiers:
            # 哈希标识符到点
            raw_point = self.hash_to_point(id)
            self.raw_points[id] = raw_point

            # 乘以k1: [k1] * H(id)
            encrypted_point = self.scalar_multiply_point(raw_point, self.private_scalar)
            point_bytes = self.point_to_bytes(encrypted_point)
            self.original_points.append(point_bytes)
            self.identifier_map[point_bytes] = id
            self.point_bytes_map[point_bytes] = encrypted_point

        # 随机打乱顺序增强隐私
        random.shuffle(self.original_points)
        return self.original_points

    def protocol_p2_step2(self, p1_points, id_value_pairs):
        """P2 - 步骤2: 处理标识符和关联值"""
        # 生成Paillier密钥
        self.paillier_public_key, self.paillier_private_key = paillier.generate_paillier_keypair()

        # 双加密P1的点: [k2] * ([k1] * H(id))
        self.double_encrypted = []
        for point_bytes in p1_points:
            point = self.bytes_to_point(point_bytes)
            # 乘以k2: [k2] * point
            encrypted_point = self.scalar_multiply_point(point, self.private_scalar)
            encrypted_bytes = self.point_to_bytes(encrypted_point)
            self.double_encrypted.append(encrypted_bytes)

        # 处理P2自己的数据: [k2] * H(w_j)
        self.paired_data = []
        self.p2_single_points = {}  # 存储单加密点

        for id, value in id_value_pairs:
            # 哈希标识符到点
            raw_point = self.hash_to_point(id)

            # 乘以k2: [k2] * H(id)
            encrypted_point = self.scalar_multiply_point(raw_point, self.private_scalar)
            encrypted_bytes = self.point_to_bytes(encrypted_point)
            enc_value = self.paillier_public_key.encrypt(value)

            # 存储配对关系
            self.paired_data.append((encrypted_bytes, enc_value))
            self.p2_single_points[encrypted_bytes] = encrypted_point

        # 随机打乱所有列表以增强隐私
        random.shuffle(self.double_encrypted)
        random.shuffle(self.paired_data)

        # 拆分为点和值列表
        p2_single_encrypted = [item[0] for item in self.paired_data]
        p2_encrypted_values = [item[1] for item in self.paired_data]

        return self.double_encrypted, p2_single_encrypted, p2_encrypted_values, self.paillier_public_key

    def protocol_p1_step3(self, p2_double_encrypted, p2_single_encrypted, encrypted_values, paillier_pk):
        """P1 - 步骤3: 计算交集和关联值之和"""
        # 计算k1的模逆元
        k1_inv = pow(self.private_scalar, -1, self.order)  # 使用模逆更安全

        # 部分解密: [k1^{-1}] * ([k2][k1]H(v_i)) = [k2]H(v_i)
        partial_decrypted = {}
        for enc_bytes in p2_double_encrypted:
            point = self.bytes_to_point(enc_bytes)
            # 乘以k1^{-1}: [k1_inv] * point
            decrypted_point = self.scalar_multiply_point(point, k1_inv)
            decrypted_bytes = self.point_to_bytes(decrypted_point)
            partial_decrypted[decrypted_bytes] = True

        # 查找交集并求和
        intersection_count = 0
        sum_ciphertext = paillier_pk.encrypt(0)

        # 遍历P2的点和值对
        for point_bytes, enc_value in zip(p2_single_encrypted, encrypted_values):
            if point_bytes in partial_decrypted:
                intersection_count += 1
                sum_ciphertext += enc_value

        return sum_ciphertext, intersection_count

    def protocol_p2_step4(self, sum_ciphertext):
        """P2 - 步骤4: 解密结果"""
        return self.paillier_private_key.decrypt(sum_ciphertext)


if __name__ == "__main__":
    print("=" * 50)
    print("私有交集求和协议实现 - 修复版")
    print("=" * 50)

    # 调试模式 - 小数据集
    print("\n调试模式: 小数据集测试")
    p1_identifiers = ["user1", "user2"]
    p2_data = [("user1", 100), ("user3", 200)]

    print("P1 标识符:", p1_identifiers)
    print("P2 数据:", p2_data)

    p1 = PrivateIntersectionSum()
    p2 = PrivateIntersectionSum()

    # 步骤1: P1 发送处理后的标识符
    p1_points = p1.protocol_p1_step1(p1_identifiers)
    print(f"P1 发送 {len(p1_points)} 个加密点")

    # 步骤2: P2 处理数据并返回
    p2_double, p2_single, p2_values, paillier_pk = p2.protocol_p2_step2(p1_points, p2_data)
    print(f"P2 返回: {len(p2_double)} 个双加密点, {len(p2_single)} 个单加密点")

    # 步骤3: P1 计算交集和值之和
    sum_ct, count = p1.protocol_p1_step3(p2_double, p2_single, p2_values, paillier_pk)
    print(f"P1 计算: 交集大小 = {count}, 加密和 = {sum_ct.ciphertext() if count > 0 else 'N/A'}")

    # 步骤4: P2 解密结果
    total = p2.protocol_p2_step4(sum_ct) if count > 0 else 0
    print(f"P2 解密: 值之和 = {total}")

    # 预期结果
    expected_intersection = set(p1_identifiers) & set(id for id, _ in p2_data)
    expected_sum = sum(value for id, value in p2_data if id in expected_intersection)

    print(f"预期: 交集大小 = {len(expected_intersection)}, 值之和 = {expected_sum}")
    print(f"状态: {'通过' if count == len(expected_intersection) and total == expected_sum else '失败'}")

    # 测试数据集1 - 基础测试
    print("\n测试用例1: 基础功能验证")
    p1_identifiers = ["202200460165@SDU.com", "202200460166@SDU.com", "202200460167@SDU.com", "202200460169@SDU.com"]
    p2_data = [
        ("202200460165@SDU.com", 100),
        ("202200460166@SDU.com", 200),
        ("202200460168@SDU.com", 300),
        ("202200460169@SDU.com", 400)
    ]

    p1 = PrivateIntersectionSum()
    p2 = PrivateIntersectionSum()

    p1_points = p1.protocol_p1_step1(p1_identifiers)
    p2_double, p2_single, p2_values, paillier_pk = p2.protocol_p2_step2(p1_points, p2_data)
    sum_ct, count = p1.protocol_p1_step3(p2_double, p2_single, p2_values, paillier_pk)
    total = p2.protocol_p2_step4(sum_ct) if count > 0 else 0

    # 预期结果计算
    expected_intersection = set(p1_identifiers) & set(id for id, _ in p2_data)
    expected_sum = sum(value for id, value in p2_data if id in expected_intersection)

    print(f"结果: 交集大小 = {count}, 关联值之和 = {total}")
    print(f"预期: 交集大小 = {len(expected_intersection)}, 关联值之和 = {expected_sum}")
    print(f"状态: {'通过' if count == len(expected_intersection) and total == expected_sum else '失败'}")

    # 测试数据集2 - 无交集
    print("\n测试用例2: 无交集")
    p1_identifiers = ["userA@test.com", "userB@test.com"]
    p2_data = [("userC@test.com", 500), ("userD@test.com", 600)]

    p1 = PrivateIntersectionSum()
    p2 = PrivateIntersectionSum()

    p1_points = p1.protocol_p1_step1(p1_identifiers)
    p2_double, p2_single, p2_values, paillier_pk = p2.protocol_p2_step2(p1_points, p2_data)
    sum_ct, count = p1.protocol_p1_step3(p2_double, p2_single, p2_values, paillier_pk)
    total = p2.protocol_p2_step4(sum_ct) if count > 0 else 0

    print(f"结果: 交集大小 = {count}, 关联值之和 = {total}")
    print(f"预期: 交集大小 = 0, 关联值之和 = 0")
    print(f"状态: {'通过' if count == 0 and total == 0 else '失败'}")

    # 测试数据集3 - 完全交集
    print("\n测试用例3: 完全交集")
    shared_ids = ["id1", "id2", "id3"]
    p1_identifiers = shared_ids
    p2_data = [(id, i * 100) for i, id in enumerate(shared_ids, 1)]

    p1 = PrivateIntersectionSum()
    p2 = PrivateIntersectionSum()

    p1_points = p1.protocol_p1_step1(p1_identifiers)
    p2_double, p2_single, p2_values, paillier_pk = p2.protocol_p2_step2(p1_points, p2_data)
    sum_ct, count = p1.protocol_p1_step3(p2_double, p2_single, p2_values, paillier_pk)
    total = p2.protocol_p2_step4(sum_ct) if count > 0 else 0

    expected_sum = sum(value for _, value in p2_data)
    print(f"结果: 交集大小 = {count}, 关联值之和 = {total}")
    print(f"预期: 交集大小 = {len(shared_ids)}, 关联值之和 = {expected_sum}")
    print(f"状态: {'通过' if count == len(shared_ids) and total == expected_sum else '失败'}")

    print("\n所有测试完成!")