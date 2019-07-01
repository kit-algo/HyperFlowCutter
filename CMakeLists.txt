cmake_minimum_required(VERSION 3.5)

project(Hyper)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Wall -Wpedantic")
set(CMAKE_CXX_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native")

add_subdirectory(extern/tlx)

add_executable(OwnPaToH OwnPaToH.cpp io/hmetisreader.cpp util/random.cpp)
target_link_libraries(OwnPaToH tlx ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)

add_executable(HyperFlowCutter staggered_hyperflowcutter.cpp util/random.cpp io/hmetisreader.cpp util/logging_information.cpp)
target_link_libraries(HyperFlowCutter tlx ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)

add_executable(ReBaHFC RefineAndBalanceHyperFlowCutter.cpp util/random.cpp io/hmetisreader.cpp util/logging_information.cpp)
target_link_libraries(ReBaHFC tlx ${CMAKE_SOURCE_DIR}/extern/libpatoh.a)
