cmake_minimum_required(VERSION 3.10)
project(fieldline_tracing)

find_package(NetCDF REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/opt/homebrew/include)
link_directories(/opt/homebrew/lib)

#add_executable(example main.cpp)
add_executable(fieldline_tracing fieldline_tracing.cxx SplineInterpolation.cxx NetCDFLoader.cxx RK4.cxx)
#target_link_libraries(fieldline_tracing NetCDF::NetCDF)
target_link_libraries(fieldline_tracing netcdf)
