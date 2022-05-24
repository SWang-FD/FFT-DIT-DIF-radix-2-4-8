#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "fft.h"

int main(void){
    FILE *fpsignal = fopen("0signal.txt","r");
    FILE *fpfftInput = fopen("1fftInput.txt","w");
    FILE *fpspectrum = fopen("2spectrum.txt","w");
    FILE *fpifftOutput = fopen("3ifftOutput.txt","w");
    if (fpsignal == NULL ||
        fpfftInput == NULL ||
        fpspectrum == NULL ||
        fpifftOutput == NULL){
        printf("Error: open file failed!");
        return 1;
    }

    float *signal = (float*) malloc(sizeof(float) * fftLength);
    complexFloat *fftInput = (complexFloat*) malloc(sizeof(complexFloat) * fftLength);
    complexFloat *spectrum = (complexFloat*) malloc(sizeof(complexFloat) * fftLength);
    complexFloat *ifftOutput = (complexFloat*) malloc(sizeof(complexFloat) * fftLength);

    // 读取输入信号数据
    for(int i = 0; i < fftLength; i++){
        fscanf(fpsignal,"%f", &signal[i]);
    }
    fclose(fpsignal);

    // FFT输入数据
    for(int i = 0; i < fftLength; i++){
        fftInput[i].re = signal[i];
        fftInput[i].im = 0;
        fprintf(fpfftInput, "%f\t%f\n", fftInput[i].re, fftInput[i].im);
    }
    fclose(fpfftInput);
    printf("read input signal data finished !\n");

    // FFT输出频谱
    fft(fftInput, spectrum, false);
    for(int i = 0; i < fftLength; i++){
        fprintf(fpspectrum, "%.5f \t%.5f\n", spectrum[i].re, spectrum[i].im);
    }
    fclose(fpspectrum);
    printf("write fft output spectrum finished !\n");

    // IFFT输出结果
    fft(spectrum, ifftOutput, true);
    for(int i = 0; i < fftLength; i++){
        fprintf(fpifftOutput, "%f\t%f\n", ifftOutput[i].re, ifftOutput[i].im);
    }
    fclose(fpifftOutput);
    printf("write ifft output data finished !\n");

    free(signal);
    free(fftInput);
    free(spectrum);
    free(ifftOutput);
    signal = NULL;
    fftInput = NULL;
    spectrum = NULL;
    ifftOutput = NULL;
    return 0;
}