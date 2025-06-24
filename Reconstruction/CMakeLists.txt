cmake_minimum_required(VERSION 3.10)
project(TrackReconstruction)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree Matrix Physics MathCore Gui Graf3d)

# Include ROOT libraries and setup
include(${ROOT_USE_FILE})

# Set source and include directories
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Set source files
file(GLOB SOURCES 
    ${SOURCE_DIR}/*.cpp
)

# Verify sources were found
if(NOT SOURCES)
    message(FATAL_ERROR "No source files found in ${SOURCE_DIR}")
endif()

# Print found sources for debugging
message(STATUS "Found sources: ${SOURCES}")

# Set include directories
include_directories(
    ${INCLUDE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

# Create executable
add_executable(track_reconstruction ${SOURCES})

# Link ROOT libraries
target_link_libraries(track_reconstruction
    ${ROOT_LIBRARIES}
)

# Install target (optional)
install(TARGETS track_reconstruction DESTINATION bin)
