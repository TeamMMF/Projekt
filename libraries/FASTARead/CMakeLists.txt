cmake_minimum_required(VERSION 3.2)

add_library(FASTARead FASTARead.cpp)
install (TARGETS FASTARead DESTINATION bin)
install (FILES FASTARead.h DESTINATION include)