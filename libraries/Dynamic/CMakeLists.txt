add_library(Dynamic Dynamic.cpp)
install (TARGETS Dynamic DESTINATION bin)
target_link_libraries(Dynamic Common)
install (FILES Dynamic.hpp DESTINATION include)