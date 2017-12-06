cmake_minimum_required(VERSION 3.2)

add_library(lcskpp lcskpp.cpp)
install (TARGETS lcskpp DESTINATION bin)
install (FILES lcskpp.h DESTINATION include)