cmake_minimum_required(VERSION 3.9)
project(maxdp)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")
set(CMAKE_CXX_FLAGS_PROFILING "-pg ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE}) 
message(STATUS "CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_CXX_FLAGS_PROFILING ${CMAKE_CXX_FLAGS_PROFILING}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fconcepts -DIL_STD")

find_package(cplex REQUIRED)
find_package(fmt REQUIRED)
set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost 1.60.0 COMPONENTS program_options system filesystem REQUIRED)

include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS} $ENV{HOME}/libf77)

add_executable(maxdp main.cpp 2ge.cpp bal.cpp balvns.cpp cmdline.cpp color.cpp constructive.cpp coudert.cpp dynamicbalance.cpp dynamicdispersion.cpp ec.cpp lowerbound.cpp oscillate.cpp readinstance.cpp solution.cpp stats.cpp ub.cpp ubcommon.cpp ubi.cpp ubk.cpp ubrb.cpp ubs.cpp umdp.cpp umdpfernandez.cpp tabucol/PartialColAndTabuCol.cpp tabucol/initializeColoring.cpp tabucol/manipulateArrays.cpp tabucol/reactcol.cpp tabucol/tabu.cpp tabucol/util/Graph.cpp tabucol/util/inputGraph.cpp tabucol/util/manipulateArrays.cpp tabucol/util/solution.cpp tabucol/util/tabu.cpp tabucol/util/util.cpp btdsatur/bkdmain.cpp btdsatur/bktdsat.cpp btdsatur/colorrtns.cpp btdsatur/graph.cpp btdsatur/maxclique.cpp )
target_link_libraries(maxdp ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt dl pthread fmt::fmt)

add_executable(generateinstance generateinstance.cpp)
target_link_libraries(generateinstance ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt dl pthread fmt::fmt)
