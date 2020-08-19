#include "mex.h"
#include <stdio.h>

void printHelloWorld() {
  printf("Hello, world!\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  printHelloWorld();
}
