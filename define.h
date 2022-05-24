#ifndef __DEFINE_H__
#define __DEFINE_H__

#define fftRadix 8 // 基2, 4, 8FFT
#define fftStage 4 // FFT阶段数
#define fftLength pow(fftRadix, fftStage) // FFT点数, max = 2^24

#define DIT 0 // 时域抽取FFT
#define DIF 1 // 频域抽取FFT
#define fftType DIF


typedef unsigned int u32;
typedef unsigned short u16;
typedef char u8;

typedef struct _complexFloat{
    float re;
    float im;
}complexFloat;

#endif
