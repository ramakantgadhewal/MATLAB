#include "mex.h"
#include <iostream>

using namespace std;

void printHello() {
  cout << "Hello, world!\n";
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  printHello();
}
