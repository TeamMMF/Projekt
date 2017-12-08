#include <cstdio>
#include "Dynamic.h"

int main(){

  // sample array. 
  int b[] = { 1, 3, 2, 4, 3, 5, 4, 6 };
  
  // should return 5
  printf("%d\n", lis( b, 8 ) );
}