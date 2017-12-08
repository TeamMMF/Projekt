add_library(CustomTypes CustomTypes.cpp)
install (TARGETS CustomTypes DESTINATION bin)
install (FILES Common.cpp DESTINATION include)
