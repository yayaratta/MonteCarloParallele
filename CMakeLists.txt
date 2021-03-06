cmake_minimum_required(VERSION 2.8)
project(mc-pricerParallele C CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

#Peut etre erreur
if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


#Detect MPI
set(MPI_LIBS)
find_package(MPI REQUIRED)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set (MPI_LIBS ${MPI_LIBS} ${MPI_CXX_LIBRARIES})
    add_definitions(-DOMPI_SKIP_MPICXX)
else (MPI_CXX_FOUND)
    message(STATUS "MPI not found.")
endif (MPI_CXX_FOUND)

# Detect Boost
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
message(STATUS "Boost include directory ${Boost_INCLUDE_DIR}")

# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Add math lib
set(LIBS ${LIBS} m)

add_executable(pricer src/parser.cpp src/pricerPacked.cpp src/Options/Option.cpp src/BlackScholesModel.cpp src/MonteCarlo.cpp src/Options/OptionBasket.cpp src/Options/AsianOption.cpp src/Options/PerformanceOption.cpp)
target_link_libraries(pricer ${MPI_LIBS} ${LIBS})
install(TARGETS pricer DESTINATION bin)

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()
