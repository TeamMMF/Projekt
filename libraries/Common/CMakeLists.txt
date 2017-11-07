add_library(Common Common.cpp)
install (TARGETS Common DESTINATION bin)
install (FILES Common.hpp DESTINATION include)