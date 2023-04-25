#include "fft.h"

int fft(
    complexFloat *input, 
    complexFloat *output,
    bool IFFT
){
    complexFloat *butterfly = (complexFloat*) malloc(sizeof(complexFloat) * fftLength);
    
    // DIT需要重新排列输入顺序，作为FFT第一阶段的输入数据
    for(u32 i = 0; i < fftLength; i++){
        // 时域抽取FFT, 逆序输入, 顺序输出
        // 频域抽取FFT, 顺序输入, 逆序输出
        u32 iOrder;
        iOrder = (fftType == DIT) ? reverseBit(i) : i;
        if(IFFT){
            butterfly[iOrder].re = input[i].re;
            butterfly[iOrder].im = -input[i].im;
        }else{
            butterfly[iOrder] = input[i];
        }
    }
    
    complexFloat tmpButterfly[(int)(log(fftRadix)/log(2))][fftRadix];
    int indexButterfly[fftRadix];
    int indexWeight;

    // 共有fftStage个阶段
    // 每个阶段有groupNum个蝶形计算组
    // 每个蝶形计算组有groupSize个蝶形单元
    for(int i = 0; i < fftStage; i++){
        int groupNum, groupSize;
        groupNum = (fftType == DIT) ? pow(fftRadix, fftStage-1 - i) : pow(fftRadix, i);
        groupSize = (fftType == DIT) ? pow(fftRadix, i) : pow(fftRadix, fftStage-1 - i);
        for(int j = 0; j < groupNum; j++){
            for(int k = 0; k < groupSize; k++){
                for(int m = 0; m < fftRadix; m++){
                    indexButterfly[m] = j * groupSize * fftRadix + k + m * groupSize;
                }
                indexWeight = k * groupNum;
                if(fftRadix == 2){
                    if(fftType == DIT){
                        tmpButterfly[0][0] = complexMul(butterfly[indexButterfly[0]], getWeight(0,indexWeight));
                        tmpButterfly[0][1] = complexMul(butterfly[indexButterfly[1]], getWeight(1,indexWeight));
                        butterfly[indexButterfly[0]] = complexAdd(tmpButterfly[0][0], tmpButterfly[0][1]);
                        butterfly[indexButterfly[1]] = complexSub(tmpButterfly[0][0], tmpButterfly[0][1]);
                    }else if(fftType == DIF){
                        tmpButterfly[0][0] = complexAdd(butterfly[indexButterfly[0]], butterfly[indexButterfly[1]]);
                        tmpButterfly[0][1] = complexSub(butterfly[indexButterfly[0]], butterfly[indexButterfly[1]]);
                        butterfly[indexButterfly[0]] = complexMul(tmpButterfly[0][0], getWeight(0,indexWeight));
                        butterfly[indexButterfly[1]] = complexMul(tmpButterfly[0][1], getWeight(1,indexWeight));
                    }
                }else if(fftRadix == 4){
                    if(fftType == DIT){
                        tmpButterfly[0][0] = complexMul(butterfly[indexButterfly[0]], getWeight(0,indexWeight));
                        tmpButterfly[0][1] = complexMul(butterfly[indexButterfly[1]], getWeight(1,indexWeight));
                        tmpButterfly[0][2] = complexMul(butterfly[indexButterfly[2]], getWeight(2,indexWeight));
                        tmpButterfly[0][3] = complexMul(butterfly[indexButterfly[3]], getWeight(3,indexWeight));
                        tmpButterfly[1][0] = complexAdd(tmpButterfly[0][0], tmpButterfly[0][2]);
                        tmpButterfly[1][2] = complexSub(tmpButterfly[0][0], tmpButterfly[0][2]);
                        tmpButterfly[1][1] = complexAdd(tmpButterfly[0][1], tmpButterfly[0][3]);
                        tmpButterfly[1][3] = complexSub(tmpButterfly[0][1], tmpButterfly[0][3]);
                        tmpButterfly[1][3] = complexMul(tmpButterfly[1][3], (complexFloat){0,-1});
                        butterfly[indexButterfly[0]] = complexAdd(tmpButterfly[1][0], tmpButterfly[1][1]);
                        butterfly[indexButterfly[1]] = complexAdd(tmpButterfly[1][2], tmpButterfly[1][3]);
                        butterfly[indexButterfly[2]] = complexSub(tmpButterfly[1][0], tmpButterfly[1][1]);
                        butterfly[indexButterfly[3]] = complexSub(tmpButterfly[1][2], tmpButterfly[1][3]);
                    }else if(fftType == DIF){
                        tmpButterfly[0][0] = complexAdd(butterfly[indexButterfly[0]], butterfly[indexButterfly[2]]);
                        tmpButterfly[0][2] = complexSub(butterfly[indexButterfly[0]], butterfly[indexButterfly[2]]);
                        tmpButterfly[0][1] = complexAdd(butterfly[indexButterfly[1]], butterfly[indexButterfly[3]]);
                        tmpButterfly[0][3] = complexSub(butterfly[indexButterfly[1]], butterfly[indexButterfly[3]]);
                        tmpButterfly[0][3] = complexMul(tmpButterfly[0][3], (complexFloat){0,-1});
                        tmpButterfly[1][0] = complexAdd(tmpButterfly[0][0], tmpButterfly[0][1]);
                        tmpButterfly[1][1] = complexAdd(tmpButterfly[0][2], tmpButterfly[0][3]);
                        tmpButterfly[1][2] = complexSub(tmpButterfly[0][0], tmpButterfly[0][1]);
                        tmpButterfly[1][3] = complexSub(tmpButterfly[0][2], tmpButterfly[0][3]);
                        butterfly[indexButterfly[0]] = complexMul(tmpButterfly[1][0], getWeight(0,indexWeight));
                        butterfly[indexButterfly[1]] = complexMul(tmpButterfly[1][1], getWeight(1,indexWeight));
                        butterfly[indexButterfly[2]] = complexMul(tmpButterfly[1][2], getWeight(2,indexWeight));
                        butterfly[indexButterfly[3]] = complexMul(tmpButterfly[1][3], getWeight(3,indexWeight));
                    }
                }else if(fftRadix == 8){
                    if(fftType == DIT){
                        tmpButterfly[0][0] = complexMul(butterfly[indexButterfly[0]], getWeight(0,indexWeight));
                        tmpButterfly[0][1] = complexMul(butterfly[indexButterfly[1]], getWeight(1,indexWeight));
                        tmpButterfly[0][2] = complexMul(butterfly[indexButterfly[2]], getWeight(2,indexWeight));
                        tmpButterfly[0][3] = complexMul(butterfly[indexButterfly[3]], getWeight(3,indexWeight));
                        tmpButterfly[0][4] = complexMul(butterfly[indexButterfly[4]], getWeight(4,indexWeight));
                        tmpButterfly[0][5] = complexMul(butterfly[indexButterfly[5]], getWeight(5,indexWeight));
                        tmpButterfly[0][6] = complexMul(butterfly[indexButterfly[6]], getWeight(6,indexWeight));
                        tmpButterfly[0][7] = complexMul(butterfly[indexButterfly[7]], getWeight(7,indexWeight));
                        tmpButterfly[1][0] = complexAdd(tmpButterfly[0][0], tmpButterfly[0][4]);
                        tmpButterfly[1][4] = complexSub(tmpButterfly[0][0], tmpButterfly[0][4]);
                        tmpButterfly[1][2] = complexAdd(tmpButterfly[0][2], tmpButterfly[0][6]);
                        tmpButterfly[1][6] = complexSub(tmpButterfly[0][2], tmpButterfly[0][6]);
                        tmpButterfly[1][1] = complexAdd(tmpButterfly[0][1], tmpButterfly[0][5]);
                        tmpButterfly[1][5] = complexSub(tmpButterfly[0][1], tmpButterfly[0][5]);
                        tmpButterfly[1][3] = complexAdd(tmpButterfly[0][3], tmpButterfly[0][7]);
                        tmpButterfly[1][7] = complexSub(tmpButterfly[0][3], tmpButterfly[0][7]);
                        tmpButterfly[1][6] = complexMul(tmpButterfly[1][6], (complexFloat){0,-1});
                        tmpButterfly[1][7] = complexMul(tmpButterfly[1][7], (complexFloat){0,-1});
                        tmpButterfly[2][0] = complexAdd(tmpButterfly[1][0], tmpButterfly[1][2]);
                        tmpButterfly[2][2] = complexSub(tmpButterfly[1][0], tmpButterfly[1][2]);
                        tmpButterfly[2][4] = complexAdd(tmpButterfly[1][4], tmpButterfly[1][6]);
                        tmpButterfly[2][6] = complexSub(tmpButterfly[1][4], tmpButterfly[1][6]);
                        tmpButterfly[2][1] = complexAdd(tmpButterfly[1][1], tmpButterfly[1][3]);
                        tmpButterfly[2][3] = complexSub(tmpButterfly[1][1], tmpButterfly[1][3]);
                        tmpButterfly[2][5] = complexAdd(tmpButterfly[1][5], tmpButterfly[1][7]);
                        tmpButterfly[2][7] = complexSub(tmpButterfly[1][5], tmpButterfly[1][7]);
                        tmpButterfly[2][5] = complexMul(tmpButterfly[2][5], (complexFloat){cos(M_PI/4),-sin(M_PI/4)});
                        tmpButterfly[2][3] = complexMul(tmpButterfly[2][3], (complexFloat){0,-1});
                        tmpButterfly[2][7] = complexMul(tmpButterfly[2][7], (complexFloat){cos(3*M_PI/4),-sin(3*M_PI/4)});
                        butterfly[indexButterfly[0]] = complexAdd(tmpButterfly[2][0], tmpButterfly[2][1]);
                        butterfly[indexButterfly[1]] = complexAdd(tmpButterfly[2][4], tmpButterfly[2][5]);
                        butterfly[indexButterfly[2]] = complexAdd(tmpButterfly[2][2], tmpButterfly[2][3]);
                        butterfly[indexButterfly[3]] = complexAdd(tmpButterfly[2][6], tmpButterfly[2][7]);
                        butterfly[indexButterfly[4]] = complexSub(tmpButterfly[2][0], tmpButterfly[2][1]);
                        butterfly[indexButterfly[5]] = complexSub(tmpButterfly[2][4], tmpButterfly[2][5]);
                        butterfly[indexButterfly[6]] = complexSub(tmpButterfly[2][2], tmpButterfly[2][3]);
                        butterfly[indexButterfly[7]] = complexSub(tmpButterfly[2][6], tmpButterfly[2][7]);
                    }else if(fftType == DIF){
                        tmpButterfly[0][0] = complexAdd(butterfly[indexButterfly[0]], butterfly[indexButterfly[4]]);
                        tmpButterfly[0][4] = complexSub(butterfly[indexButterfly[0]], butterfly[indexButterfly[4]]);
                        tmpButterfly[0][2] = complexAdd(butterfly[indexButterfly[2]], butterfly[indexButterfly[6]]);
                        tmpButterfly[0][6] = complexSub(butterfly[indexButterfly[2]], butterfly[indexButterfly[6]]);
                        tmpButterfly[0][1] = complexAdd(butterfly[indexButterfly[1]], butterfly[indexButterfly[5]]);
                        tmpButterfly[0][5] = complexSub(butterfly[indexButterfly[1]], butterfly[indexButterfly[5]]);
                        tmpButterfly[0][3] = complexAdd(butterfly[indexButterfly[3]], butterfly[indexButterfly[7]]);
                        tmpButterfly[0][7] = complexSub(butterfly[indexButterfly[3]], butterfly[indexButterfly[7]]);
                        tmpButterfly[0][6] = complexMul(tmpButterfly[0][6], (complexFloat){0,-1});
                        tmpButterfly[0][7] = complexMul(tmpButterfly[0][7], (complexFloat){0,-1});
                        tmpButterfly[1][0] = complexAdd(tmpButterfly[0][0], tmpButterfly[0][2]);
                        tmpButterfly[1][2] = complexSub(tmpButterfly[0][0], tmpButterfly[0][2]);
                        tmpButterfly[1][4] = complexAdd(tmpButterfly[0][4], tmpButterfly[0][6]);
                        tmpButterfly[1][6] = complexSub(tmpButterfly[0][4], tmpButterfly[0][6]);
                        tmpButterfly[1][1] = complexAdd(tmpButterfly[0][1], tmpButterfly[0][3]);
                        tmpButterfly[1][3] = complexSub(tmpButterfly[0][1], tmpButterfly[0][3]);
                        tmpButterfly[1][5] = complexAdd(tmpButterfly[0][5], tmpButterfly[0][7]);
                        tmpButterfly[1][7] = complexSub(tmpButterfly[0][5], tmpButterfly[0][7]);
                        tmpButterfly[1][5] = complexMul(tmpButterfly[1][5], (complexFloat){cos(M_PI/4),-sin(M_PI/4)});
                        tmpButterfly[1][3] = complexMul(tmpButterfly[1][3], (complexFloat){0,-1});
                        tmpButterfly[1][7] = complexMul(tmpButterfly[1][7], (complexFloat){cos(3*M_PI/4),-sin(3*M_PI/4)});
                        tmpButterfly[2][0] = complexAdd(tmpButterfly[1][0], tmpButterfly[1][1]);
                        tmpButterfly[2][1] = complexAdd(tmpButterfly[1][4], tmpButterfly[1][5]);
                        tmpButterfly[2][2] = complexAdd(tmpButterfly[1][2], tmpButterfly[1][3]);
                        tmpButterfly[2][3] = complexAdd(tmpButterfly[1][6], tmpButterfly[1][7]);
                        tmpButterfly[2][4] = complexSub(tmpButterfly[1][0], tmpButterfly[1][1]);
                        tmpButterfly[2][5] = complexSub(tmpButterfly[1][4], tmpButterfly[1][5]);
                        tmpButterfly[2][6] = complexSub(tmpButterfly[1][2], tmpButterfly[1][3]);
                        tmpButterfly[2][7] = complexSub(tmpButterfly[1][6], tmpButterfly[1][7]);
                        butterfly[indexButterfly[0]] = complexMul(tmpButterfly[2][0], getWeight(0,indexWeight));
                        butterfly[indexButterfly[1]] = complexMul(tmpButterfly[2][1], getWeight(1,indexWeight));
                        butterfly[indexButterfly[2]] = complexMul(tmpButterfly[2][2], getWeight(2,indexWeight));
                        butterfly[indexButterfly[3]] = complexMul(tmpButterfly[2][3], getWeight(3,indexWeight));
                        butterfly[indexButterfly[4]] = complexMul(tmpButterfly[2][4], getWeight(4,indexWeight));
                        butterfly[indexButterfly[5]] = complexMul(tmpButterfly[2][5], getWeight(5,indexWeight));
                        butterfly[indexButterfly[6]] = complexMul(tmpButterfly[2][6], getWeight(6,indexWeight));
                        butterfly[indexButterfly[7]] = complexMul(tmpButterfly[2][7], getWeight(7,indexWeight));
                    }
                }else{
                    printf("Error: only support fftRadix = 2, 4, 8");
                    return 1;
                }
            }
        }
    }

    // DIF需要重新排列输出顺序
    for(u32 i = 0; i < fftLength; i++){
        // 时域抽取FFT, 逆序输入, 顺序输出
        u32 iOrder;
        iOrder = (fftType == DIF) ? reverseBit(i) : i;
        if(IFFT){
            output[iOrder].re = butterfly[i].re / fftLength;
    	    output[iOrder].im = -butterfly[i].im / fftLength;
        }else{
            output[iOrder] = butterfly[i];
        }
    }
    free(butterfly);
    butterfly = NULL;
    return 0;
}

static complexFloat getWeight(
    int iButterfly,
    int indexWeight
){
    complexFloat result;
    result.re = cos(2 * M_PI * iButterfly * indexWeight / fftLength);
    result.im = -sin(2 * M_PI * iButterfly * indexWeight / fftLength);
    return result;
}

static u32 reverseBit(
    u32 n
){
    u32 bitWidth = 0;
    while((1<<bitWidth) < fftLength) bitWidth++;
    if(fftRadix == 8){
        n = ((n & 0b111000111000111000111000) >> 3 ) | ((n & 0b000111000111000111000111) << 3 );
        n = ((n & 0b111111000000111111000000) >> 6 ) | ((n & 0b000000111111000000111111) << 6 );
        n = ((n & 0b111111111111000000000000) >> 12 ) | ((n & 0b000000000000111111111111) << 12 );
        n = n >> (24 - bitWidth);
    }else{
        if(fftRadix == 2) n = ((n & 0xAAAAAAAA) >> 1 ) | ((n & 0x55555555) << 1 );
        n = ((n & 0xCCCCCCCC) >> 2 ) | ((n & 0x33333333) << 2 );
        n = ((n & 0xF0F0F0F0) >> 4 ) | ((n & 0x0F0F0F0F) << 4 );
        n = ((n & 0xFF00FF00) >> 8 ) | ((n & 0x00FF00FF) << 8 );
        n = ((n & 0xFFFF0000) >> 16 ) | ((n & 0x0000FFFF) << 16 );
        n = n >> (32 - bitWidth);
    }
    return n;
}

static complexFloat complexAdd(
    complexFloat A,
    complexFloat B
){
    complexFloat result;
    result.re = A.re + B.re;
    result.im = A.im + B.im;
    return result;
}

static complexFloat complexSub(
    complexFloat A,
    complexFloat B
){
    complexFloat result;
    result.re = A.re - B.re;
    result.im = A.im - B.im;
    return result;
}

static complexFloat complexMul(
    complexFloat A,
    complexFloat B
){
    complexFloat result;
    result.re = A.re * B.re - A.im * B.im;
    result.im = A.im * B.re + A.re * B.im;
    return result;
}
