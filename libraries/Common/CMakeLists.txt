add_library(Common Common.cpp)
install (TARGETS Common DESTINATION bin)
target_link_libraries(Common CustomTypes)
install (FILES Common.hpp DESTINATION include)
