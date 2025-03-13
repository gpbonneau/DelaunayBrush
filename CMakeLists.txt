cmake_minimum_required(VERSION 3.16)
project(circumsphereDLL)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
set(Boost_NO_BOOST_CMAKE OFF)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)
cmake_policy(SET CMP0167 NEW)  # Use the new method for finding Boost


find_package(Boost REQUIRED)  # Make sure Boost is found

find_package(CGAL REQUIRED)

# Set Compiler for Unix system
if(UNIX)
   set(CMAKE_CXX_COMPILER g++)                      # Can switch to clang++ if prefered
   add_definitions(-g -O2 -std=c++14 -Wall -Wextra -Wfatal-errors -Wno-pragmas -Wno-unknown-warning-option) # Can adapt compiler flags if needed
   add_definitions(-Wno-sign-compare -Wno-type-limits) # Remove some warnings
endif()


# Set Compiler for Windows/Visual Studio
if(MSVC)
   set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT  ${executable_name} ) # default project (avoids AllBuild)
   set_target_properties( ${executable_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}$<0:> ) # default output in root dir
   set_target_properties( ${executable_name} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) # default debug execution in root dir
   
   # Avoids the warning /W3 overided by /W4 when using Ninja
   if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
   endif()

    add_definitions(/MP /wd4244 /wd4127 /wd4267 /wd4706 /wd4458 /wd4996 /wd26495 /openmp)   # Parallel build (/MP) + disable some warnings
    source_group(TREE ${CMAKE_SOURCE_DIR} FILES ${src_files})  #Allow to explore source directories as a tree in Visual Studio
endif()

# Add your project files
file(GLOB SRC_FILES *.cpp)
# add_executable(${PROJECT_NAME} ${SRC_FILES})
add_library(${PROJECT_NAME} SHARED ${SRC_FILES})

target_link_libraries(${PROJECT_NAME}  PUBLIC CGAL::CGAL)



