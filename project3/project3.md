# Project3用circom实现poseidon2哈希算法的电路
要求： 
1) poseidon2哈希算法参数参考参考文档1的Table1，选用哈希函数的参数(n,t,d)=(256,2,5)
2）电路的公开输入用poseidon2哈希值，隐私输入为哈希原象，哈希算法的输入只考虑一个block即可。
3) 用Groth16算法生成证明

## 在 **Ubuntu 系统** 上安装 Circom 和 SnarkJS 的步骤。

```bash

sudo apt update

sudo apt install -y nodejs

cargo install --locked circom --version 2.1.7

sudo npm install -g snarkjs@0.7.7

git clone https://github.com/iden3/circomlib.git ~/circomlib

sudo apt install -y cmake
npm install -g wasm-pack

circom --version  
snarkjs --version 

```

![image-20250709112217728](C:\Users\28795\Desktop\project3\image-20250709112217728.png)

Circom及SnarkJS工作流程：

![image-20250709112926309](C:\Users\28795\Desktop\project3\image-20250709112926309.png)

* Circom：编写ZKP电路的DSL.

* Snarkjs：实现ZKP逻辑，实现Groth16.

* Node.js:Witness
首先编译电路，生成R1CS文件
```
circom circuits/poseidon2.circom --r1cs --wasm --sym 
```
采用和Circom电路相同逻辑的python函数来得到Poseidon2的哈希后的实际结果。

poseidon2_hash.py的结果如下：

![image-20250709144817110](C:\Users\28795\Desktop\project3\image-20250709144817110.png)

根据函数结果准备输入文件（input.json)

电路需要两个输入：私有输入 `in_private` 和公共输出 `hash_output`。`hash_output` 必须是 `in_private` 经过 Poseidon2 哈希后的实际结果。

**`input.json` 示例：**
JSON
```
{
    "in_private": ["10", "20"],
    "hash_output": "8576556543856606835948037659587238157455622082871502003842405520578322528367 "
  }

```
生成witness：
```
node generate_witness.js poseidon2.wasm input.json witness.wtns
```

最后生成证明（Groth16协议）
我们将使用 snarkjs 的 Groth16 协议来生成零知识证明。
生成并验证proof

```
snarkjs groth16 prove poseidon2_final.zkey witness.wtns proof.json public.json
snarkjs zkey export verificationkey poseidon2.r1cs poseidon2_0000.zkey
snarkjs groth16 verify verification_key.json public.json proof.json
```

![image-20250709164321174](C:\Users\28795\Desktop\project3\image-20250709164321174.png)

![image-20250709164342566](C:\Users\28795\Desktop\project3\image-20250709164342566.png)