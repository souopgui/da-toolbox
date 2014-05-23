//OBSERVATION.CPP
//Defines C++ functions to be used for some purposes
#include "cpp_tools.h"

void dvector_sort_(T_Rell *rda_x, int *id_nbelt){
  dvector_sort(rda_x, *id_nbelt);
}

void dvector_sort(T_Rell *rda_x, int id_nbelt){
    //vector<T_Rell> val(rda_x, rda_x + id_nbelt);
    //std::sort(val.begin(), val.end());
    //std::copy(val.begin(), val.end(), rda_x) ;
  std::sort(rda_x, rda_x + id_nbelt);
}


void ivector_sort_(T_Int *ida_x, int *id_nbelt){
  ivector_sort(ida_x, *id_nbelt);
}

void ivector_sort(T_Int *ida_x, int id_nbelt){
  std::sort(ida_x, ida_x + id_nbelt);
}
