#include <iostream>
#include <cstring>
#include <vector>
#include <chrono>
#include <immintrin.h>
#include <intrin.h>
#include <cstdlib>

// 编译器兼容性处理
#if defined(_MSC_VER)
#define bswap32(x) _byteswap_ulong(x)
#define ALIGN32 __declspec(align(32))
#define ALIGNED_MALLOC(size, align) _aligned_malloc(size, align)
#define ALIGNED_FREE(ptr) _aligned_free(ptr)
#elif defined(__GNUC__) || defined(__clang__)
#define bswap32(x) __builtin_bswap32(x)
#define ALIGN32 __attribute__((aligned(32)))
#define ALIGNED_MALLOC(size, align) aligned_alloc(align, size)
#define ALIGNED_FREE(ptr) free(ptr)
#else
#error "Unsupported compiler"
#endif

// 算法常量定义
ALIGN32 const uint32_t IV[8] = {
    0x7380166F, 0x4914B2B9, 0x172442D7, 0xDA8A0600,
    0xA96F30BC, 0x163138AA, 0xE38DEE4D, 0xB0FB0E4E
};
const uint32_t T0 = 0x79CC4519;
const uint32_t T1 = 0x7A879D8A;

// 核心算法宏定义
#define ROTL(x, n) (((x) << (n)) | ((x) >> (32 - (n))))
#define FF0(x, y, z) ((x) ^ (y) ^ (z))
#define GG0(x, y, z) ((x) ^ (y) ^ (z))
#define FF1(x, y, z) (((x) & (y)) | ((x) & (z)) | ((y) & (z)))
#define GG1(x, y, z) (((x) & (y)) | ((~(x)) & (z)))

inline uint32_t FF(int j, uint32_t x, uint32_t y, uint32_t z) {
    return j < 16 ? FF0(x, y, z) : FF1(x, y, z);
}

inline uint32_t GG(int j, uint32_t x, uint32_t y, uint32_t z) {
    return j < 16 ? GG0(x, y, z) : GG1(x, y, z);
}

// 置换函数
inline uint32_t P0(uint32_t x) { return x ^ ROTL(x, 9) ^ ROTL(x, 17); }
inline uint32_t P1(uint32_t x) { return x ^ ROTL(x, 15) ^ ROTL(x, 23); }

// AVX2辅助函数
inline __m256i mm256_rotl_epi32(__m256i x, int n) {
    return _mm256_or_si256(_mm256_slli_epi32(x, n), _mm256_srli_epi32(x, 32 - n));
}

inline __m256i mm256_p0(__m256i x) {
    __m256i rot9 = mm256_rotl_epi32(x, 9);
    __m256i rot17 = mm256_rotl_epi32(x, 17);
    return _mm256_xor_si256(x, _mm256_xor_si256(rot9, rot17));
}

// ==================== 基础版本实现 ====================
void message_expansion_basic(const uint32_t* block, uint32_t* W, uint32_t* W_prime) {
    for (int i = 0; i < 16; ++i) {
        W[i] = bswap32(block[i]);
    }

    for (int i = 16; i < 68; ++i) {
        W[i] = P1(W[i - 16] ^ W[i - 9] ^ ROTL(W[i - 3], 15))
            ^ ROTL(W[i - 13], 7) ^ W[i - 6];
    }

    for (int i = 0; i < 64; ++i) {
        W_prime[i] = W[i] ^ W[i + 4];
    }
}

void compress_basic(uint32_t* V, const uint32_t* W, const uint32_t* W_prime) {
    uint32_t A = V[0], B = V[1], C = V[2], D = V[3];
    uint32_t E = V[4], F = V[5], G = V[6], H = V[7];

    for (int j = 0; j < 64; ++j) {
        const uint32_t Tj = j < 16 ? T0 : T1;
        const uint32_t SS1 = ROTL(ROTL(A, 12) + E + ROTL(Tj, j % 32), 7);
        const uint32_t SS2 = SS1 ^ ROTL(A, 12);
        const uint32_t TT1 = FF(j, A, B, C) + D + SS2 + W_prime[j];
        const uint32_t TT2 = GG(j, E, F, G) + H + SS1 + W[j];

        D = C;
        C = ROTL(B, 9);
        B = A;
        A = TT1;
        H = G;
        G = ROTL(F, 19);
        F = E;
        E = P0(TT2);
    }

    V[0] ^= A; V[1] ^= B; V[2] ^= C; V[3] ^= D;
    V[4] ^= E; V[5] ^= F; V[6] ^= G; V[7] ^= H;
}

void sm3_hash_basic(const uint8_t* data, size_t len, uint8_t digest[32]) {
    uint32_t V[8];
    memcpy(V, IV, sizeof(IV));

    const uint64_t bit_len = static_cast<uint64_t>(len) * 8;
    const size_t pad_len = ((len + 1 + 8) + 63) & ~63;
    std::vector<uint8_t> padded(pad_len, 0);
    memcpy(padded.data(), data, len);
    padded[len] = 0x80;
    for (int i = 0; i < 8; ++i) {
        padded[pad_len - 8 + i] = static_cast<uint8_t>(bit_len >> (56 - i * 8));
    }

    for (size_t i = 0; i < pad_len; i += 64) {
        uint32_t W[68], W_prime[64];
        message_expansion_basic(reinterpret_cast<const uint32_t*>(padded.data() + i), W, W_prime);
        compress_basic(V, W, W_prime);
    }

    for (int i = 0; i < 8; ++i) {
        V[i] = bswap32(V[i]);
        memcpy(digest + i * 4, &V[i], 4);
    }
}

// ==================== 优化的AVX2版本 - 真正并行处理 ====================
// 转置函数：将8个连续块重组为AVX2友好格式
void transpose_8x8(uint32_t* in, __m256i* out) {
    for (int i = 0; i < 8; i++) {
        __m256i row = _mm256_set_epi32(
            in[i + 7 * 8], in[i + 6 * 8], in[i + 5 * 8], in[i + 4 * 8],
            in[i + 3 * 8], in[i + 2 * 8], in[i + 1 * 8], in[i + 0 * 8]
        );
        out[i] = row;
    }
}

// 消息扩展的AVX2优化版本
void message_expansion_avx2(const __m256i* block, __m256i* W) {
    // 前16个字直接加载
    for (int i = 0; i < 16; i++) {
        W[i] = block[i];
    }

    // 并行消息扩展 (16-67)
    for (int i = 16; i < 68; i++) {
        __m256i w16 = W[i - 16];
        __m256i w9 = W[i - 9];
        __m256i w3 = mm256_rotl_epi32(W[i - 3], 15);
        __m256i w13 = mm256_rotl_epi32(W[i - 13], 7);
        __m256i w6 = W[i - 6];

        // P1(x) = x ^ ROTL(x, 15) ^ ROTL(x, 23)
        __m256i p1_in = _mm256_xor_si256(_mm256_xor_si256(w16, w9), w3);
        __m256i p1 = _mm256_xor_si256(p1_in,
            _mm256_xor_si256(mm256_rotl_epi32(p1_in, 15),
                mm256_rotl_epi32(p1_in, 23)));

        // W[i] = P1(...) ^ ROTL(W[i-13], 7) ^ W[i-6]
        W[i] = _mm256_xor_si256(p1, _mm256_xor_si256(w13, w6));
    }
}

// 压缩函数的AVX2优化版本
void compress_avx2(__m256i* V, const __m256i* W) {
    // 初始化状态
    __m256i A = V[0], B = V[1], C = V[2], D = V[3];
    __m256i E = V[4], F = V[5], G = V[6], H = V[7];

    // 预计算轮常数
    const __m256i T0_vec = _mm256_set1_epi32(T0);
    const __m256i T1_vec = _mm256_set1_epi32(T1);

    // 预计算W_prime
    ALIGN32 __m256i W_prime[64];
    for (int i = 0; i < 64; i++) {
        W_prime[i] = _mm256_xor_si256(W[i], W[i + 4]);
    }

    for (int j = 0; j < 64; ++j) {
        // 选择常量T
        const __m256i Tj_vec = j < 16 ? T0_vec : T1_vec;
        const int rot_j = j % 32;

        // 计算SS1 = ROTL(ROTL(A, 12) + E + ROTL(Tj, j % 32), 7)
        __m256i rotA12 = mm256_rotl_epi32(A, 12);
        __m256i rotTj = mm256_rotl_epi32(Tj_vec, rot_j);
        __m256i SS1 = _mm256_add_epi32(rotA12, E);
        SS1 = _mm256_add_epi32(SS1, rotTj);
        SS1 = mm256_rotl_epi32(SS1, 7);

        // 计算SS2 = SS1 ^ rotA12
        __m256i SS2 = _mm256_xor_si256(SS1, rotA12);

        // 计算TT1 = FF(j, A, B, C) + D + SS2 + W_prime[j]
        __m256i TT1;
        if (j < 16) {
            TT1 = _mm256_xor_si256(A, _mm256_xor_si256(B, C)); // FF0
        }
        else {
            // FF1 = (A & B) | (A & C) | (B & C)
            __m256i AB = _mm256_and_si256(A, B);
            __m256i AC = _mm256_and_si256(A, C);
            __m256i BC = _mm256_and_si256(B, C);
            TT1 = _mm256_or_si256(_mm256_or_si256(AB, AC), BC);
        }
        TT1 = _mm256_add_epi32(TT1, D);
        TT1 = _mm256_add_epi32(TT1, SS2);
        TT1 = _mm256_add_epi32(TT1, W_prime[j]);

        // 计算TT2 = GG(j, E, F, G) + H + SS1 + W[j]
        __m256i TT2;
        if (j < 16) {
            TT2 = _mm256_xor_si256(E, _mm256_xor_si256(F, G)); // GG0
        }
        else {
            // GG1 = (E & F) | ((~E) & G)
            __m256i EF = _mm256_and_si256(E, F);
            __m256i notE = _mm256_andnot_si256(E, _mm256_set1_epi32(0xFFFFFFFF));
            __m256i notE_and_G = _mm256_and_si256(notE, G);
            TT2 = _mm256_or_si256(EF, notE_and_G);
        }
        TT2 = _mm256_add_epi32(TT2, H);
        TT2 = _mm256_add_epi32(TT2, SS1);
        TT2 = _mm256_add_epi32(TT2, W[j]);

        // 更新状态
        D = C;
        C = mm256_rotl_epi32(B, 9);
        B = A;
        A = TT1;
        H = G;
        G = mm256_rotl_epi32(F, 19);
        F = E;
        E = mm256_p0(TT2);
    }

    // 更新最终状态
    V[0] = _mm256_xor_si256(V[0], A);
    V[1] = _mm256_xor_si256(V[1], B);
    V[2] = _mm256_xor_si256(V[2], C);
    V[3] = _mm256_xor_si256(V[3], D);
    V[4] = _mm256_xor_si256(V[4], E);
    V[5] = _mm256_xor_si256(V[5], F);
    V[6] = _mm256_xor_si256(V[6], G);
    V[7] = _mm256_xor_si256(V[7], H);
}

// 处理完整块的主函数
void sm3_hash_avx2_parallel(const uint8_t* data, size_t len, uint8_t digest[32]) {
    // 初始化状态 - 每个通道独立
    ALIGN32 __m256i V[8];
    for (int i = 0; i < 8; i++) {
        V[i] = _mm256_setr_epi32(
            IV[i], IV[i], IV[i], IV[i],
            IV[i], IV[i], IV[i], IV[i]
        );
    }

    // 计算完整块数和剩余字节
    size_t total_blocks = len / 64;
    size_t remaining_bytes = len % 64;
    size_t processed = 0;

    // 处理完整块 (每次处理8个块)
    ALIGN32 uint32_t block_data[8 * 16] = { 0 };
    ALIGN32 __m256i W[68];

    for (size_t i = 0; i < total_blocks; i += 8) {
        // 一次处理8个块
        size_t blocks_this_round = (total_blocks - i) > 8 ? 8 : total_blocks - i;

        // 重组数据为AVX2友好格式
        for (size_t b = 0; b < blocks_this_round; b++) {
            const uint32_t* src = reinterpret_cast<const uint32_t*>(data + (i + b) * 64);
            for (int j = 0; j < 16; j++) {
                block_data[j * 8 + b] = bswap32(src[j]);
            }
        }

        // 填充剩余块（不足8个块时）
        for (size_t b = blocks_this_round; b < 8; b++) {
            for (int j = 0; j < 16; j++) {
                block_data[j * 8 + b] = 0;
            }
        }

        // 转置数据
        ALIGN32 __m256i transposed[16];
        transpose_8x8(block_data, transposed);

        // 消息扩展
        message_expansion_avx2(transposed, W);

        // 压缩
        compress_avx2(V, W);

        processed += blocks_this_round * 64;
    }

    // 处理剩余数据
    if (remaining_bytes > 0) {
        // 提取第一个通道的状态 (其他通道保持不变)
        uint32_t V_basic[8];
        for (int i = 0; i < 8; i++) {
            V_basic[i] = _mm256_extract_epi32(V[i], 0);
        }

        // 准备最后一个块
        ALIGN32 uint8_t last_block[64] = { 0 };
        memcpy(last_block, data + processed, remaining_bytes);
        last_block[remaining_bytes] = 0x80;

        // 添加长度信息
        uint64_t bit_len = static_cast<uint64_t>(len) * 8;
        for (int i = 0; i < 8; i++) {
            last_block[56 + i] = static_cast<uint8_t>(bit_len >> (56 - i * 8));
        }

        // 使用基础版本处理最后一个块
        uint32_t W_basic[68], W_prime_basic[64];
        message_expansion_basic(reinterpret_cast<const uint32_t*>(last_block),
            W_basic, W_prime_basic);
        compress_basic(V_basic, W_basic, W_prime_basic);

        // 仅更新第一个通道
        for (int i = 0; i < 8; i++) {
            V[i] = _mm256_insert_epi32(V[i], V_basic[i], 0);
        }
    }

    // 提取最终状态（仅第一个通道有效）
    uint32_t final_state[8];
    for (int i = 0; i < 8; i++) {
        final_state[i] = _mm256_extract_epi32(V[i], 0);
    }

    // 转换为大端序并输出
    for (int i = 0; i < 8; i++) {
        final_state[i] = bswap32(final_state[i]);
        memcpy(digest + i * 4, &final_state[i], 4);
    }
}

// ==================== 性能测试工具 ====================
double test_performance(const char* name, void (*hash_func)(const uint8_t*, size_t, uint8_t*),
    const std::vector<uint8_t>& data, int runs = 10) {
    uint8_t digest[32];

    // 预热缓存
    for (int i = 0; i < 3; i++) {
        hash_func(data.data(), data.size(), digest);
    }

    // 正式测试
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < runs; i++) {
        hash_func(data.data(), data.size(), digest);
    }
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(end - start).count();
    double speed = (data.size() * runs) / (elapsed * 1024 * 1024);

    std::cout << name << " 吞吐量: " << speed << " MB/s\n";
    return speed;
}

// CPU特性检测
bool cpu_supports_avx2() {
    int cpuInfo[4];
    __cpuid(cpuInfo, 1);
    if (!(cpuInfo[2] & (1 << 27))) return false; // 检查OSXSAVE

    __cpuidex(cpuInfo, 7, 0);
    return (cpuInfo[1] & (1 << 5)) != 0; // 检查AVX2
}

// 测试用例
int main() {
    // 标准测试向量（"abc"）
    const char* test_str = "abc";
    uint8_t digest_basic[32];
    uint8_t digest_avx2[32];
    const uint8_t expected[32] = {
        0x66, 0xC7, 0xF0, 0xF4, 0x62, 0xEE, 0xED, 0xD9,
        0xD1, 0xF2, 0xD4, 0x6B, 0xDC, 0x10, 0xE4, 0xE2,
        0x41, 0x67, 0xC4, 0x87, 0x5C, 0xF2, 0xF7, 0xA2,
        0x29, 0x7D, 0xA0, 0x2B, 0x8F, 0x4B, 0xA8, 0xE0
    };

    std::cout << "=== 正确性验证 ===\n";

    // 基础版测试
    sm3_hash_basic(reinterpret_cast<const uint8_t*>(test_str), 3, digest_basic);
    std::cout << "基础版 SM3(\"abc\"): ";
    for (int i = 0; i < 32; i++) printf("%02x", digest_basic[i]);
    if (memcmp(digest_basic, expected, 32) == 0) std::cout << " [通过]";
    else std::cout << " [失败]";
    std::cout << "\n";

    // AVX2版测试
    if (cpu_supports_avx2()) {
        sm3_hash_avx2_parallel(reinterpret_cast<const uint8_t*>(test_str), 3, digest_avx2);

        std::cout << "AVX2优化版 SM3(\"abc\"): ";
        for (int i = 0; i < 32; i++) printf("%02x", digest_avx2[i]);
        if (memcmp(digest_avx2, expected, 32) == 0) std::cout << " [通过]";
        else std::cout << " [失败]";
        std::cout << "\n";
    }
    else {
        std::cout << "当前CPU不支持AVX2，跳过AVX2版测试\n";
    }

    // 长消息测试 (1000字节)
    std::vector<uint8_t> long_data(1000, 'a');
    uint8_t digest_long_basic[32], digest_long_avx2[32];
    const uint8_t expected_long[32] = {
        0xf4, 0xbe, 0xdc, 0xa9, 0x73, 0x22, 0x7d, 0x45,
        0xc5, 0xb8, 0x22, 0x55, 0x1d, 0x2e, 0x76, 0x2d,
        0x4c, 0xfb, 0x0e, 0x9a, 0xf7, 0x0b, 0x24, 0x14,
        0x52, 0x54, 0x57, 0x27, 0xb5, 0xfb, 0x04, 0x6f
    };

    // 基础版长消息测试
    sm3_hash_basic(long_data.data(), long_data.size(), digest_long_basic);
    std::cout << "基础版长消息测试: ";
    if (memcmp(digest_long_basic, expected_long, 32) == 0) std::cout << " [通过]";
    else std::cout << " [失败]";
    std::cout << "\n";

    if (cpu_supports_avx2()) {
        sm3_hash_avx2_parallel(long_data.data(), long_data.size(), digest_long_avx2);
        std::cout << "AVX2版长消息测试: ";
        if (memcmp(digest_long_avx2, expected_long, 32) == 0) {
            std::cout << " [通过]\n";
        }
        else {
            std::cout << " [失败]\n";
            std::cout << "期望: ";
            for (int i = 0; i < 32; i++) printf("%02x", expected_long[i]);
            std::cout << "\n实际: ";
            for (int i = 0; i < 32; i++) printf("%02x", digest_long_avx2[i]);
            std::cout << "\n";
        }
    }

    // 性能测试（100MB数据）
    const size_t TEST_SIZE = 1024 * 1024 * 100;
    std::vector<uint8_t> bulk_data(TEST_SIZE, 0x61); // 填充'a'

    std::cout << "\n=== 性能测试 (100MB数据) ===\n";

    double basic_speed = test_performance("基础版", sm3_hash_basic, bulk_data, 5);
    double avx2_speed = 0;

    if (cpu_supports_avx2()) {
        avx2_speed = test_performance("AVX2并行优化版", sm3_hash_avx2_parallel, bulk_data, 10);
        if (basic_speed > 0) {
            std::cout << "加速比: " << avx2_speed / basic_speed << "x\n";
        }
    }
    else {
        std::cout << "当前CPU不支持AVX2，跳过AVX2性能测试\n";
    }

    return 0;
}