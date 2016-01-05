#ifndef _UTIL_H
#define _UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void expand2d(float *uo, float *ui, int nzpad, int nxpad, int nz, int nx);

float *damp_make(int nb);

void damp2d_apply(float **uu, float *damp, int nz, int nx, int nb);

void wwin2d(float **uo, float **ui, int nzo, int nxo, int nb);

void wfld2d_inject(float **uo, float **ui, int nzo, int nxo, int nb);

bool approximatelyEqual(float a, float b, float epsilon);
bool essentiallyEqual(float a, float b, float epsilon);

#endif
