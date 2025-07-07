#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// SM4官方SBox
static const uint8_t SM4_SBOX[256] = {
    0xD6, 0x90, 0xE9, 0xFE, 0xCC, 0xE1, 0x3D, 0xB7, 0x16, 0xB6, 0x14, 0xC2, 0x28, 0xFB, 0x2C, 0x05,
    0x2B, 0x67, 0x9A, 0x76, 0x2A, 0xBE, 0x04, 0xC3, 0xAA, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
    0x9C, 0x42, 0x50, 0xF4, 0x91, 0xEF, 0x98, 0x7A, 0x33, 0x54, 0x0B, 0x43, 0xED, 0xCF, 0xAC, 0x62,
    0xE4, 0xB3, 0x1C, 0xA9, 0xC9, 0x08, 0xE8, 0x95, 0x80, 0xDF, 0x94, 0xFA, 0x75, 0x8F, 0x3F, 0xA6,
    0x47, 0x07, 0xA7, 0xFC, 0xF3, 0x73, 0x17, 0xBA, 0x83, 0x59, 0x3C, 0x19, 0xE6, 0x85, 0x4F, 0xA8,
    0x68, 0x6B, 0x81, 0xB2, 0x71, 0x64, 0xDA, 0x8B, 0xF8, 0xEB, 0x0F, 0x4B, 0x70, 0x56, 0x9D, 0x35,
    0x1E, 0x24, 0x0E, 0x5E, 0x63, 0x58, 0xD1, 0xA2, 0x25, 0x22, 0x7C, 0x3B, 0x01, 0x21, 0x78, 0x87,
    0xD4, 0x00, 0x46, 0x57, 0x9F, 0xD3, 0x27, 0x52, 0x4C, 0x36, 0x02, 0xE7, 0xA0, 0xC4, 0xC8, 0x9E,
    0xEA, 0xBF, 0x8A, 0xD2, 0x40, 0xC7, 0x38, 0xB5, 0xA3, 0xF7, 0xF2, 0xCE, 0xF9, 0x61, 0x15, 0xA1,
    0xE0, 0xAE, 0x5D, 0xA4, 0x9B, 0x34, 0x1A, 0x55, 0xAD, 0x93, 0x32, 0x30, 0xF5, 0x8C, 0xB1, 0xE3,
    0x1D, 0xF6, 0xE2, 0x2E, 0x82, 0x66, 0xCA, 0x60, 0xC0, 0x29, 0x23, 0xAB, 0x0D, 0x53, 0x4E, 0x6F,
    0xD5, 0xDB, 0x37, 0x45, 0xDE, 0xFD, 0x8E, 0x2F, 0x03, 0xFF, 0x6A, 0x72, 0x6D, 0x6C, 0x5B, 0x51,
    0x8D, 0x1B, 0xAF, 0x92, 0xBB, 0xDD, 0xBC, 0x7F, 0x11, 0xD9, 0x5C, 0x41, 0x1F, 0x10, 0x5A, 0xD8,
    0x0A, 0xC1, 0x31, 0x88, 0xA5, 0xCD, 0x7B, 0xBD, 0x2D, 0x74, 0xD0, 0x12, 0xB8, 0xE5, 0xB4, 0xB0,
    0x89, 0x69, 0x97, 0x4A, 0x0C, 0x96, 0x77, 0x7E, 0x65, 0xB9, 0xF1, 0x09, 0xC5, 0x6E, 0xC6, 0x84,
    0x18, 0xF0, 0x7D, 0xEC, 0x3A, 0xDC, 0x4D, 0x20, 0x79, 0xEE, 0x5F, 0x3E, 0xD7, 0xCB, 0x39, 0x48
};

// 固定参数 
static const uint32_t FK[4] = {
    0xA3B1BAC6, 0x56AA3350, 0x677D9197, 0xB27022DC
};

static const uint32_t CK[32] = {
    0x00070E15, 0x1C232A31, 0x383F464D, 0x545B6269,
    0x70777E85, 0x8C939AA1, 0xA8AFB6BD, 0xC4CBD2D9,
    0xE0E7EEF5, 0xFC030A11, 0x181F262D, 0x343B4249,
    0x50575E65, 0x6C737A81, 0x888F969D, 0xA4ABB2B9,
    0xC0C7CED5, 0xDCE3EAF1, 0xF8FF060D, 0x141B2229,
    0x30373E45, 0x4C535A61, 0x686F767D, 0x848B9299,
    0xA0A7AEB5, 0xBCC3CAD1, 0xD8DFE6ED, 0xF4FB0209,
    0x10171E25, 0x2C333A41, 0x484F565D, 0x646B7279
};

// 查表法优化使用的4KB表
static uint32_t SM4_TTABLE[4][256];
static int sm4_ttable_initialized = 0;

// 32位循环左移
static inline uint32_t rotl32(uint32_t x, uint8_t n) {
    return (x << n) | (x >> (32 - n));
}

// 线性变换L
static inline uint32_t sm4_L_transform(uint32_t x) {
    return x ^ rotl32(x, 2) ^ rotl32(x, 10) ^ rotl32(x, 18) ^ rotl32(x, 24);
}

// 初始化T-Table 
void sm4_init_ttable() {
    if (sm4_ttable_initialized) return;

    for (int i = 0; i < 256; i++) {
        uint8_t s = SM4_SBOX[i];
        // 正确方法：独立计算每个字节位置的变换
        SM4_TTABLE[0][i] = sm4_L_transform((uint32_t)s << 24);
        SM4_TTABLE[1][i] = sm4_L_transform((uint32_t)s << 16);
        SM4_TTABLE[2][i] = sm4_L_transform((uint32_t)s << 8);
        SM4_TTABLE[3][i] = sm4_L_transform((uint32_t)s);
    }

    sm4_ttable_initialized = 1;
}

// 密钥扩展
void sm4_key_schedule(const uint8_t key[16], uint32_t rk[32]) {
    uint32_t k[36];

    for (int i = 0; i < 4; i++) {
        k[i] = ((uint32_t)key[4 * i] << 24) |
            ((uint32_t)key[4 * i + 1] << 16) |
            ((uint32_t)key[4 * i + 2] << 8) |
            key[4 * i + 3];
        k[i] ^= FK[i];
    }

    for (int i = 0; i < 32; i++) {
        uint32_t tmp = k[i + 1] ^ k[i + 2] ^ k[i + 3] ^ CK[i];

        tmp = (SM4_SBOX[(tmp >> 24) & 0xFF] << 24) |
            (SM4_SBOX[(tmp >> 16) & 0xFF] << 16) |
            (SM4_SBOX[(tmp >> 8) & 0xFF] << 8) |
            SM4_SBOX[tmp & 0xFF];

        tmp = tmp ^ rotl32(tmp, 13) ^ rotl32(tmp, 23);

        k[i + 4] = k[i] ^ tmp;
        rk[i] = k[i + 4];
    }
}

// 基础轮函数
static uint32_t sm4_round_func_basic(uint32_t x) {
    uint32_t b = ((uint32_t)SM4_SBOX[(x >> 24) & 0xFF] << 24) |
        ((uint32_t)SM4_SBOX[(x >> 16) & 0xFF] << 16) |
        ((uint32_t)SM4_SBOX[(x >> 8) & 0xFF] << 8) |
        SM4_SBOX[x & 0xFF];

    return sm4_L_transform(b);
}

// 查表法轮函数
static uint32_t sm4_round_func_ttable(uint32_t x) {
    return SM4_TTABLE[0][(x >> 24) & 0xFF] ^
        SM4_TTABLE[1][(x >> 16) & 0xFF] ^
        SM4_TTABLE[2][(x >> 8) & 0xFF] ^
        SM4_TTABLE[3][x & 0xFF];
}

// 通用加密函数
static void _sm4_crypt(const uint8_t input[16], uint8_t output[16],
    const uint32_t rk[32], int use_ttable) {
    uint32_t block[4];

    for (int i = 0; i < 4; i++) {
        block[i] = ((uint32_t)input[4 * i] << 24) |
            ((uint32_t)input[4 * i + 1] << 16) |
            ((uint32_t)input[4 * i + 2] << 8) |
            input[4 * i + 3];
    }

    for (int round = 0; round < 32; round++) {
        uint32_t tmp = block[1] ^ block[2] ^ block[3] ^ rk[round];

        tmp = use_ttable ?
            sm4_round_func_ttable(tmp) :
            sm4_round_func_basic(tmp);

        tmp ^= block[0];

        block[0] = block[1];
        block[1] = block[2];
        block[2] = block[3];
        block[3] = tmp;
    }

    uint32_t final_block[4] = { block[3], block[2], block[1], block[0] };

    for (int i = 0; i < 4; i++) {
        output[4 * i] = (final_block[i] >> 24) & 0xFF;
        output[4 * i + 1] = (final_block[i] >> 16) & 0xFF;
        output[4 * i + 2] = (final_block[i] >> 8) & 0xFF;
        output[4 * i + 3] = final_block[i] & 0xFF;
    }
}

// 加密函数
void sm4_encrypt_basic(const uint8_t input[16], uint8_t output[16], const uint32_t rk[32]) {
    _sm4_crypt(input, output, rk, 0);
}

void sm4_encrypt_ttable(const uint8_t input[16], uint8_t output[16], const uint32_t rk[32]) {
    if (!sm4_ttable_initialized) sm4_init_ttable();
    _sm4_crypt(input, output, rk, 1);
}

// 解密函数
void sm4_decrypt_ttable(const uint8_t input[16], uint8_t output[16], const uint32_t rk[32]) {
    if (!sm4_ttable_initialized) sm4_init_ttable();

    uint32_t rk_dec[32];
    for (int i = 0; i < 32; i++) {
        rk_dec[i] = rk[31 - i];
    }

    _sm4_crypt(input, output, rk_dec, 1);
}

// 性能测试函数
void benchmark_sm4() {
    uint8_t key[16] = { 0x01,0x23,0x45,0x67,0x89,0xAB,0xCD,0xEF,0xFE,0xDC,0xBA,0x98,0x76,0x54,0x32,0x10 };
    uint8_t plain[16] = { 0x01,0x23,0x45,0x67,0x89,0xAB,0xCD,0xEF,0xFE,0xDC,0xBA,0x98,0x76,0x54,0x32,0x10 };
    uint8_t cipher[16];
    uint32_t rk[32];

    sm4_key_schedule(key, rk);

    const int iterations = 1000000;
    clock_t start, end;

    // 基础加密性能测试
    start = clock();
    for (int i = 0; i < iterations; i++) {
        sm4_encrypt_basic(plain, cipher, rk);
    }
    end = clock();
    double basic_time = (double)(end - start) / CLOCKS_PER_SEC;

    // 查表法性能测试
    start = clock();
    for (int i = 0; i < iterations; i++) {
        sm4_encrypt_ttable(plain, cipher, rk);
    }
    end = clock();
    double ttable_time = (double)(end - start) / CLOCKS_PER_SEC;

    printf("性能测试 (%d次加密):\n", iterations);
    printf("基础方法: %.4f 秒 (%.2f MB/s)\n",
        basic_time,
        (iterations * 16) / (basic_time * 1000000));
    printf("查表方法: %.4f 秒 (%.2f MB/s)\n",
        ttable_time,
        (iterations * 16) / (ttable_time * 1000000));
    printf("速度提升: %.2f 倍\n", basic_time / ttable_time);
}

// 打印十六进制数据
void print_hex(const char* label, const uint8_t* data, size_t len) {
    printf("%s: ", label);
    for (size_t i = 0; i < len; i++) {
        printf("%02X", data[i]);
    }
    printf("\n");
}

int main() {
    uint8_t key[16] = {
        0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF,
        0xFE, 0xDC, 0xBA, 0x98, 0x76, 0x54, 0x32, 0x10
    };

    uint8_t plaintext[16] = {
        0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF,
        0xFE, 0xDC, 0xBA, 0x98, 0x76, 0x54, 0x32, 0x10
    };

    uint8_t ciphertext[16];
    uint8_t decrypted[16];
    uint32_t rk[32];

    sm4_key_schedule(key, rk);

    // 基础加密
    sm4_encrypt_basic(plaintext, ciphertext, rk);
    print_hex("基础加密结果", ciphertext, 16);

    // 查表法加密
    sm4_encrypt_ttable(plaintext, ciphertext, rk);
    print_hex("查表法加密结果", ciphertext, 16);

    // 解密验证
    sm4_decrypt_ttable(ciphertext, decrypted, rk);
    print_hex("解密结果", decrypted, 16);

    // 标准测试
    uint8_t expected[16] = {
        0x68, 0x1E, 0xDF, 0x34, 0xD2, 0x06, 0x96, 0x5E,
        0x86, 0xB3, 0xE9, 0x4F, 0x53, 0x6E, 0x42, 0x46
    };

    if (memcmp(ciphertext, expected, 16) == 0) {
        printf("加密验证通过\n");
    }
    else {
        printf("加密验证失败\n");
    }

    if (memcmp(decrypted, plaintext, 16) == 0) {
        printf("解密验证通过\n");
    }
    else {
        printf("解密验证失败\n");
    }

    // T-table
    benchmark_sm4();

    return 0;
}