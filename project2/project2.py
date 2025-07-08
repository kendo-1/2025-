import cv2
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.fftpack import dct, idct
import os
from tqdm import tqdm
from tabulate import tabulate

# 设置matplotlib使用支持中文的字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


# ===================== 辅助函数 =====================
def string_to_bits(s):
    """将字符串转换为二进制位序列"""
    bits = []
    for char in s:
        byte = format(ord(char), '08b')
        bits.extend([int(bit) for bit in byte])
    return bits


def bits_to_string(bits):
    """将二进制位序列转换回字符串"""
    chars = []
    for i in range(0, len(bits), 8):
        if i + 8 > len(bits):
            break
        byte = bits[i:i + 8]
        try:
            chars.append(chr(int(''.join(map(str, byte)), 2)))
        except:
            chars.append('?')  # 无效字符替换
    return ''.join(chars)


def calculate_ber(original_bits, extracted_bits):
    """计算位错误率(BER)"""
    min_len = min(len(original_bits), len(extracted_bits))
    if min_len == 0:
        return 1.0  # 全部错误
    errors = sum(1 for a, b in zip(original_bits[:min_len], extracted_bits[:min_len]) if a != b)
    return errors / min_len


def apply_attack(img, attack_type, **params):
    """应用指定的攻击到图像"""
    if attack_type == "rotation":
        angle = params.get('angle', 10)
        h, w = img.shape[:2]
        center = (w // 2, h // 2)
        M = cv2.getRotationMatrix2D(center, angle, 1.0)
        return cv2.warpAffine(img, M, (w, h))

    elif attack_type == "flip":
        direction = params.get('direction', 'horizontal')
        return cv2.flip(img, 1) if direction == 'horizontal' else cv2.flip(img, 0)

    elif attack_type == "crop":
        percent = params.get('percent', 0.1)
        h, w = img.shape[:2]
        crop_h = int(h * (1 - percent))
        crop_w = int(w * (1 - percent))
        return img[:crop_h, :crop_w]

    elif attack_type == "contrast":
        factor = params.get('factor', 1.5)
        img = img.astype(np.float32)
        img = (img - 127.5) * factor + 127.5
        return np.clip(img, 0, 255).astype(np.uint8)

    elif attack_type == "brightness":
        value = params.get('value', 50)
        hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
        h, s, v = cv2.split(hsv)
        v = cv2.add(v, value)
        v = np.clip(v, 0, 255)
        hsv = cv2.merge((h, s, v))
        return cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)

    elif attack_type == "jpeg":
        quality = params.get('quality', 50)
        result, encimg = cv2.imencode('.jpg', img, [int(cv2.IMWRITE_JPEG_QUALITY), quality])
        return cv2.imdecode(encimg, 1)

    elif attack_type == "blur":
        kernel_size = params.get('kernel_size', 5)
        return cv2.GaussianBlur(img, (kernel_size, kernel_size), 0)

    elif attack_type == "resize":
        scale = params.get('scale', 0.8)
        h, w = img.shape[:2]
        new_w, new_h = int(w * scale), int(h * scale)
        resized = cv2.resize(img, (new_w, new_h), interpolation=cv2.INTER_LINEAR)
        return cv2.resize(resized, (w, h), interpolation=cv2.INTER_LINEAR)

    elif attack_type == "translation":
        x = params.get('x', 20)
        y = params.get('y', 20)
        h, w = img.shape[:2]
        M = np.float32([[1, 0, x], [0, 1, y]])
        return cv2.warpAffine(img, M, (w, h))

    return img.copy()


def safe_string(s, max_len=30):
    """安全处理字符串，移除不可打印字符并截断"""
    # 移除非打印字符
    clean = ''.join(c for c in s if c.isprintable() or c == ' ')
    # 截断长度
    return clean[:max_len] + "..." if len(clean) > max_len else clean


# ===================== 水印算法 =====================
def embed_watermark_lsb(img, watermark_str):
    """LSB水印嵌入"""
    watermark_bits = string_to_bits(watermark_str)
    n_bits = len(watermark_bits)

    watermarked = img.copy()
    flat = watermarked.flatten()

    if n_bits > len(flat):
        raise ValueError(f"水印过长! 最大支持 {len(flat)} 位, 当前 {n_bits} 位")

    for i in range(n_bits):
        flat[i] = (flat[i] & 0xFE) | watermark_bits[i]

    return flat.reshape(watermarked.shape)


def extract_watermark_lsb(img, watermark_length):
    """LSB水印提取"""
    n_bits = watermark_length * 8
    flat = img.flatten()
    extracted_bits = [pixel & 1 for pixel in flat[:n_bits]]
    return bits_to_string(extracted_bits)


def embed_watermark_dct(img, watermark_str, strength=10):
    """DCT水印嵌入"""
    watermark_bits = string_to_bits(watermark_str)
    n_bits = len(watermark_bits)

    yuv = cv2.cvtColor(img, cv2.COLOR_BGR2YUV)
    Y = yuv[:, :, 0].astype(np.float32)

    h, w = Y.shape
    block_size = 8
    blocks_h = h // block_size
    blocks_w = w // block_size
    max_blocks = blocks_h * blocks_w

    if n_bits > max_blocks:
        raise ValueError(f"水印过长! 最大支持 {max_blocks} 位, 当前 {n_bits} 位")

    random.seed(42)
    block_indices = list(range(max_blocks))
    random.shuffle(block_indices)
    selected_blocks = block_indices[:n_bits]

    for bit_idx, block_idx in enumerate(selected_blocks):
        i = block_idx // blocks_w
        j = block_idx % blocks_w
        y = i * block_size
        x = j * block_size

        block = Y[y:y + block_size, x:x + block_size]
        if block.shape != (block_size, block_size):
            continue

        dct_block = dct(dct(block.T, norm='ortho').T, norm='ortho')

        coeff1 = dct_block[2, 3]
        coeff2 = dct_block[3, 2]

        if watermark_bits[bit_idx] == 1:
            if coeff1 <= coeff2:
                diff = abs(coeff2 - coeff1) + strength
                dct_block[2, 3] = coeff1 + diff
                dct_block[3, 2] = coeff2 - diff
        else:
            if coeff1 >= coeff2:
                diff = abs(coeff1 - coeff2) + strength
                dct_block[2, 3] = coeff1 - diff
                dct_block[3, 2] = coeff2 + diff

        idct_block = idct(idct(dct_block.T, norm='ortho').T, norm='ortho')
        Y[y:y + block_size, x:x + block_size] = idct_block

    yuv[:, :, 0] = np.clip(Y, 0, 255)
    return cv2.cvtColor(yuv, cv2.COLOR_YUV2BGR)


def extract_watermark_dct(img, watermark_length):
    """DCT水印提取"""
    n_bits = watermark_length * 8
    extracted_bits = []

    yuv = cv2.cvtColor(img, cv2.COLOR_BGR2YUV)
    Y = yuv[:, :, 0].astype(np.float32)

    h, w = Y.shape
    block_size = 8
    blocks_h = h // block_size
    blocks_w = w // block_size
    max_blocks = blocks_h * blocks_w

    random.seed(42)
    block_indices = list(range(max_blocks))
    random.shuffle(block_indices)
    selected_blocks = block_indices[:n_bits]

    for block_idx in selected_blocks:
        i = block_idx // blocks_w
        j = block_idx % blocks_w
        y = i * block_size
        x = j * block_size

        block = Y[y:y + block_size, x:x + block_size]
        if block.shape != (block_size, block_size):
            extracted_bits.append(0)
            continue

        dct_block = dct(dct(block.T, norm='ortho').T, norm='ortho')
        coeff1 = dct_block[2, 3]
        coeff2 = dct_block[3, 2]
        extracted_bits.append(1 if coeff1 > coeff2 else 0)

    return bits_to_string(extracted_bits)


# ===================== 鲁棒性测试 =====================
def robustness_test(watermarked_img, watermark_str, algorithm="dct"):
    """
    水印鲁棒性测试
    :param watermarked_img: 含水印图像
    :param watermark_str: 原始水印字符串
    :param algorithm: 使用的算法 ('lsb' 或 'dct')
    :return: 测试结果列表
    """
    # 定义攻击配置（减少噪声和锐化部分）
    attacks = [
        {"name": "无攻击", "type": "none", "params": {}},
        {"name": "水平翻转", "type": "flip", "params": {"direction": "horizontal"}},
        {"name": "垂直翻转", "type": "flip", "params": {"direction": "vertical"}},
        {"name": "旋转10度", "type": "rotation", "params": {"angle": 10}},
        {"name": "旋转30度", "type": "rotation", "params": {"angle": 30}},
        {"name": "平移(20,20)", "type": "translation", "params": {"x": 20, "y": 20}},
        {"name": "裁剪10%", "type": "crop", "params": {"percent": 0.1}},
        {"name": "裁剪30%", "type": "crop", "params": {"percent": 0.3}},
        {"name": "对比度增强(1.5x)", "type": "contrast", "params": {"factor": 1.5}},
        {"name": "对比度减弱(0.7x)", "type": "contrast", "params": {"factor": 0.7}},
        {"name": "亮度增加(50)", "type": "brightness", "params": {"value": 50}},
        {"name": "亮度减少(-50)", "type": "brightness", "params": {"value": -50}},
        {"name": "JPEG压缩(70)", "type": "jpeg", "params": {"quality": 70}},
        {"name": "JPEG压缩(30)", "type": "jpeg", "params": {"quality": 30}},
        {"name": "模糊(5x5)", "type": "blur", "params": {"kernel_size": 5}},
        {"name": "缩放(0.8x)", "type": "resize", "params": {"scale": 0.8}},
        {"name": "缩放(1.2x)", "type": "resize", "params": {"scale": 1.2}},
    ]

    original_bits = string_to_bits(watermark_str)
    results = []
    extract_func = extract_watermark_lsb if algorithm == "lsb" else extract_watermark_dct

    for attack in tqdm(attacks, desc=f"{algorithm.upper()}鲁棒性测试"):
        # 应用攻击
        if attack["type"] != "none":
            attacked_img = apply_attack(watermarked_img.copy(), attack["type"], **attack["params"])
        else:
            attacked_img = watermarked_img.copy()

        # 提取水印
        try:
            extracted_str = extract_func(attacked_img, len(watermark_str))
        except Exception as e:
            extracted_str = ""

        # 计算指标
        extracted_bits = string_to_bits(extracted_str) if extracted_str else []
        ber = calculate_ber(original_bits, extracted_bits)
        match = watermark_str == extracted_str

        # 安全处理提取的水印字符串
        safe_extracted = safe_string(extracted_str)

        # 保存结果
        results.append({
            "攻击类型": attack["name"],
            "提取的水印": safe_extracted,
            "完全匹配": "是" if match else "否",
            "BER": f"{ber:.4f}",
            "算法": algorithm.upper()
        })

    return results


# ===================== 可视化结果 =====================
def visualize_results(results, algorithm):
    """可视化测试结果"""
    attacks = [res["攻击类型"] for res in results]
    bers = [float(res["BER"]) for res in results]

    plt.figure(figsize=(12, 6))
    plt.bar(attacks, bers, color='skyblue')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('位错误率 (BER)')
    plt.title(f'{algorithm.upper()}水印鲁棒性测试')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f'{algorithm}_robustness_results.png', dpi=200)
    plt.show()

    # 打印表格结果
    print(f"\n{algorithm.upper()}水印鲁棒性测试结果:")
    print(tabulate(results, headers="keys", tablefmt="grid", stralign="left"))


# ===================== 主程序 =====================
def main():
    # 配置参数
    IMAGE_PATH = "host.jpg"
    WATERMARK = "kendo-202200460165"

    # 读取原始图像
    print("读取原始图像...")
    img = cv2.imread(IMAGE_PATH)
    if img is None:
        # 尝试使用其他路径
        IMAGE_PATH = os.path.join("images", "host.jpg")
        img = cv2.imread(IMAGE_PATH)
        if img is None:
            # 创建黑色图像作为后备
            img = np.zeros((512, 512, 3), dtype=np.uint8)
            print("警告: 无法加载图像，使用默认黑色图像")

    print(f"图像尺寸: {img.shape[1]}x{img.shape[0]}")
    print(f"水印内容: '{WATERMARK}' ({len(WATERMARK)} 字符)")

    # 嵌入两种水印
    print("\n嵌入LSB水印...")
    try:
        watermarked_lsb = embed_watermark_lsb(img, WATERMARK)
    except Exception as e:
        print(f"LSB嵌入失败: {str(e)}")
        watermarked_lsb = img.copy()

    print("嵌入DCT水印...")
    try:
        watermarked_dct = embed_watermark_dct(img, WATERMARK, strength=15)
    except Exception as e:
        print(f"DCT嵌入失败: {str(e)}")
        watermarked_dct = img.copy()

    # LSB水印鲁棒性测试
    print("\n开始LSB水印鲁棒性测试...")
    lsb_results = robustness_test(watermarked_lsb, WATERMARK, algorithm="lsb")
    visualize_results(lsb_results, "lsb")

    # DCT水印鲁棒性测试
    print("\n开始DCT水印鲁棒性测试...")
    dct_results = robustness_test(watermarked_dct, WATERMARK, algorithm="dct")
    visualize_results(dct_results, "dct")

    # 保存结果到文件 (使用UTF-8编码)
    with open("robustness_results.txt", "w", encoding="utf-8") as f:
        f.write("水印鲁棒性测试结果\n")
        f.write("=" * 50 + "\n")
        f.write(f"水印内容: {WATERMARK}\n")
        f.write(f"图像尺寸: {img.shape[1]}x{img.shape[0]}\n\n")

        f.write("LSB水印测试结果:\n")
        f.write(tabulate(lsb_results, headers="keys", tablefmt="grid"))
        f.write("\n\n")

        f.write("DCT水印测试结果:\n")
        f.write(tabulate(dct_results, headers="keys", tablefmt="grid"))

    print("\n所有测试完成! 结果已保存到robustness_results.txt")


if __name__ == "__main__":
    main()