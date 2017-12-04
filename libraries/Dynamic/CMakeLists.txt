add_library(Dynamic Dynamic.cpp)
install (TARGETS Dynamic DESTINATION bin)
install (FILES Dynamic.hpp DESTINATION include)