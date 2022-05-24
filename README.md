# 基2、4、8 DIT&DIF FFT

## 1. 离散傅里叶变换（DFT）

给定一个离散的实数序列$x[n]$，我们可以用DFT得到一个离散的**频谱**(spectrum)$X[k]$，其中频谱第$k$个点计算公式为：

$$
F[k]=\sum_{n=0}^{N-1}x[n]e^{-j\dfrac{2\pi}{N}k}\quad(0\leq k\leq N-1)
$$

一般我们记$W_N=e^{-j\dfrac{2\pi}{N}}$，则DFT可以写成：

$$
F[k]=FFT(x[n])=\sum_{n=0}^{N-1}x[n]W_N^{n,k}\quad(0\leq k\leq N-1)
$$


### 1.1 频域抽取基2FFT

基2快速傅里叶变换作用在$N=2^m$的序列$ x[n]$上，快速傅里叶变换的核心思想为**分而治之**，即**分治法**，该思想的核心是将一个长度为$N$的问题，分级为两个长度为$\dfrac{N}{2}$的问题，首先进行如下变换：

$$
F[k]=\sum_{n=0}^{\dfrac{N}{2}-1}x[n]\times W_{N}^{k,n}+\sum_{n=0}^{\dfrac{N}{2}-1}x[n+\dfrac{N}{2}]\times W_{N}^{k,n+\dfrac{N}{2}}
$$

根据性质$W_N^{k,n+\dfrac{N}{2}}=W_N^{k,n}\times (-1)^k$，代入(3)式得：

$$
F[k]=\sum_{n=0}^{\dfrac{N}{2}-1}W_N^{k,n}(x[n]+(-1)^kx[n+\dfrac{N}{2}])
$$

分类讨论：

+ 代入k=2r，根据性质
 
$$
W_N^{2r,n}=W_\dfrac{N}{2}^{r,n}
$$

有：

$$
F[2r]=\sum_{n=0}^{\dfrac{N}{2}-1}W_\dfrac{N}{2}^{r,n}(x[n]+x[n+\dfrac{N}{2}])=FFT(x[n]+x[n+\dfrac{N}{2}])
$$

+ 代入k=2r+1，根据性质

$$
W_N^{2r+1,n}=W_\dfrac{N}{2}^{r,n}\times e^{-j\dfrac{2\pi}{N}n}
$$

有：

$$
F[2r+1]=\sum_{n=0}^{\dfrac{N}{2}-1}W_\dfrac{N}{2}^{r,n}\times e^{-j\dfrac{2\pi}{N}n}(x[n]-x[n+\dfrac{N}{2}])=FFT(W_N^n(x[n]-x[n+\dfrac{N}{2}]))
$$

![image-20220524094910020](https://github.com/SWang-FD/FFT-DIT-DIF-radix-2-4-8/tree/main/README.assets/202205240949093.png)

### 1.2 时域抽取基2FFT

![image-20220524095207783](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205240952815.png?token=ARMJFAPPEX65ZONTX3V3C23CRQ5II)

#### 1.2.1 蝶形运算实现DIT基2FFT

蝶形运算的公式如下，蝶形运算输入为$X_0$和$X_1$，输出为$Y_0$和$Y_1$系数为$W$：

$$
\begin{cases}
Y_0=X_0+W\times X_1\\
Y_1=X_0-W\times X_1
\end{cases}
$$

+ 2点FFT：序列为$X_2=\{x[0],x[1]\}$，下标2表示2点FFT

  $$
  F_2[0]=x[0]\times W_2^{0,0}+x[1]\times W_2^{0,1}=x[0]+x[1]=X_2[0]+X_2[1]\\
  F_2[1]=x[0]\times W_2^{1,0}+x[1]\times W_2^{1,1}=x[0]-x[1]=X_2[0]-X_2[1]
  $$

+ 4点FFT：序列为$X_4=\{x[0],x[1],x[2],x[3]\}$，将4点FFT分成两个2点FFT：

  $$
  X_{2,0}=\{X_4[0]+X_4[2],X_4[1]+X_4[3]\}\\
  X_{2,1}=\{W_4^0(X_4[0]-X_4[2]),W_4^1(X_4[1]-X_4[3])\}
  $$
  
  那么：
  
  $$
  F_4[0]=X_{2,0}[0]+X_{2,0}[1]\\
  F_4[1]=X_{2,0}[0]-X_{2,0}[1]\\
  F_4[2]=X_{2,1}[0]+X_{2,1}[1]\\
  F_4[3]=X_{2,1}[0]-X_{2,1}[1]
  $$

+ N点FFT，按时间抽取，**蝶形运算**：

  ![img](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205240953639.png?token=ARMJFAOFLZ4V6MAL2EPN5MDCRQ5MK)

我们把输入序列$x[n]$的序号用$D=log_2N$比特整数表示，然后，我们把第$i$个元素和第$reverseBit(i)$个元素调换位置。

$reverseBit(i)$就是把$D$ bit整数的二进制位**倒序排列**，得到一个新的值。例如，$6=110_2$，要跟他对调位置的是$011_2=3$。

对于第$i$次迭代，有$2^{D-i}$个蝶形单元，每个蝶形单元的大小为$2^i$。

#### 1.2.2 抽取规律

**DIT**：对于不同的Radix，第$i,\quad i=0,1,2...$次迭代，有$groupNum=radix^{D-1-i}$个蝶形计算组，每个蝶形计算组有$groupSize=radix^i$个蝶算单元。

**DIF**：对于不同的Radix，第$i,\quad i=0,1,2...$次迭代，有$groupNum=radix^i$个蝶形计算组，每个蝶形计算组有$groupSize=radix^{D-1-i}$个蝶算单元。

对于第$j,\quad j=0,1,2...$个蝶形计算组中的第$k,\quad k=0,1,2...$个蝶算单元，蝶算单元中第$m,\quad m=0,1,2...$个元素在（重新排序后的）输入数组中的坐标为

$$
index[m]=j\times groupSize\times radix+k+m\times groupSize
$$


## 2. 基4FFT

![image-20220524095440235](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205240954305.png?token=ARMJFAOCGBMAXQBJRYJOKITCRQ5R2)

![image-20220524101235560](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205241012615.png?token=ARMJFAISJ7CPV77ZKMFOXFDCRQ7VA)

![7241055-058837b165f2cf1b](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205241027871.png?token=ARMJFALU5PMJCMO672BYQ4DCRRBLO)

$reverseBit(i)$排序规则，每两个bit为一组，将各组位置倒序排列，对于长度为16的FFT，bitwidth = 4，例如：$6=0110_2$，将LST的两个bit $10_2$和MST的两个bit $01_2$交换位置，即$1001_2=9$。



## 3. 基8 FFT

![image-20220524102854416](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205241028477.png?token=ARMJFALOLXABVV3W5PK4OCLCRRBSG)

![image-20220524102920138](https://github.com/SWang-FD/Picture-for-Typora/main/img/202205241029184.png?token=ARMJFAKYJHUJ4RRVYFKOFTTCRRBT2)

$reverseBit(i)$排序规则，每三个bit为一组，将各组位置倒序排列，对于长度为64的FFT，bitwidth = 6，例如：$10=001010_2$，将LST的三个bit $010_2$和MST的三个bit $001_2$交换位置，即$010001_2=17$。
