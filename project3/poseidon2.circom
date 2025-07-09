````
pragma circom 2.0.0;

// Poseidon2哈希函数实现 (t=2, d=5)
template Poseidon2() {
    signal input in_private[2];      
    signal output hash_output;       
    
    // 安全参数：t=2（状态大小），d=5（轮数）
    
    signal state[2];
    state[0] <== in_private[0];
    state[1] <== in_private[1];
    
    
    var round_constants = [
        // 5轮 * 2个元素 = 10个常数
        0x0e7d0d2c7e5c5d57, 0x0c850c0c0c0c0c0c,
        0x0d2d2d2d2d2d2d2d, 0x0a5a5a5a5a5a5a5a,
        0x0b6b6b6b6b6b6b6b, 0x0c7c7c7c7c7c7c7c,
        0x0d8d8d8d8d8d8d8d, 0x0e9e9e9e9e9e9e9e,
        0x0fafafafafafafaf, 0x0101010101010101
    ];
    
    // 预定义MDS矩阵
    var mds_matrix = [
        [2, 1],
        [1, 3]
    ];
    
    // 轮函数处理 (d=5轮)
    for (var r = 0; r < 5; r++) {
        // 1. 添加轮常数
        signal after_arc[2];
        after_arc[0] <== state[0] + round_constants[r*2];
        after_arc[1] <== state[1] + round_constants[r*2+1];
        
        // 2. S-box应用（非线性变换）
        signal after_sbox[2];
        
        // 完整轮：所有元素应用S-box (x^5)
        if (r < 2 || r >= 3) { // 第0,1,3,4轮
            // 计算x^5: (x^2)^2 * x
            signal sq0;
            sq0 <== after_arc[0] * after_arc[0];
            signal sq_sq0;
            sq_sq0 <== sq0 * sq0;
            after_sbox[0] <== sq_sq0 * after_arc[0];
            
            signal sq1;
            sq1 <== after_arc[1] * after_arc[1];
            signal sq_sq1;
            sq_sq1 <== sq1 * sq1;
            after_sbox[1] <== sq_sq1 * after_arc[1];
        } 
        
        else { // 第2轮
            // 第一个元素应用S-box
            signal sq0;
            sq0 <== after_arc[0] * after_arc[0];
            signal sq_sq0;
            sq_sq0 <== sq0 * sq0;
            after_sbox[0] <== sq_sq0 * after_arc[0];
            
            // 第二个元素保持不变
            after_sbox[1] <== after_arc[1];
        }
        
        // 3. MDS矩阵混合
        signal new_state[2];
        new_state[0] <== mds_matrix[0][0] * after_sbox[0] + mds_matrix[0][1] * after_sbox[1];
        new_state[1] <== mds_matrix[1][0] * after_sbox[0] + mds_matrix[1][1] * after_sbox[1];
        
       
        if (r < 4) {
            state[0] <== new_state[0];
            state[1] <== new_state[1];
        }
    }
    
    
    hash_output <== state[0];
}


template Main() {
    signal input in_private[2];      // 私有输入
    signal input hash_output;        
    
    
    component hasher = Poseidon2();
    
    // 连接输入
    hasher.in_private[0] <== in_private[0];
    hasher.in_private[1] <== in_private[1];
    
    
    hasher.hash_output === hash_output;
}


component main = Main();