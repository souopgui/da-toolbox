#ifndef CPP_TOOLS_H
#define CPP_TOOLS_H

#include <algorithm>
#include <iostream>
  //#include <iomanip>
  //#include <vector>

using namespace std;

#define T_Rell double
#define T_Int int

void dvector_sort(T_Rell *rda_x, int id_nbelt);
void ivector_sort(T_Int *ida_x, int id_nbelt);

//Interface for Fortran call
extern "C"{
  void dvector_sort_(T_Rell *rda_x, int *id_nbelt);
  void ivector_sort_(T_Int *ida_x, int *id_nbelt);
}

#endif
