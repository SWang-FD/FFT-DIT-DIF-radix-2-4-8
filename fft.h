#ifndef __FFT_H__
#define __FFT_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "define.h"

static complexFloat complexAdd(
    complexFloat A,
    complexFloat B
);
static complexFloat complexSub(
    complexFloat A,
    complexFloat B
);
static complexFloat complexMul(
    complexFloat A,
    complexFloat B
);
static u32 reverseBit(
    u32 n
);
static complexFloat getWeight(
    int iButterfly,
    int indexWeight
);
extern int fft(
    complexFloat *input,
    complexFloat *output,
    bool IFFT
);

#endif
