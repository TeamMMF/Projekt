#include <iostream>
#include <chrono>
#include <unistd.h>
#include "Dynamic.h"

using namespace std;

int main(){

  // sample array. 
  int b[] = { 1, 14,15,13,1, 14,15,13,3, 2, 4,5,4,2,1,3,1, 14,15,13,1, 14,15,13,3, 2, 4,5,4,2,1,3, 5, 4,0,11,7,6,6,7,8,9,10,3, 2, 4,5,4,2,1,3, 5, 4,0,11,7,6,6,7,8,9,10, 5,1, 14,15,13,1, 14,15,13,3, 2, 4,5,4,2,1,3, 5, 4,0,11,7,6,6,7,8,9,10,3, 2, 4,5,4,2,1,3, 5, 4,0,11,7,6,6,7,8,9,10, 4,0,11,7,6,6,7,8,9,10,3, 2, 4,5,4,2,1,3, 5, 4,0,11,7,6,6,7,8,9,10 };


  
  // should return 5
  
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  // int result = lis( b,  140 );
  usleep(500);
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  cout << "Duration: " << chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count() << endl;
  // printf("%d\n", result );
}