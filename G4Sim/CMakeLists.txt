cmake_minimum_required(VERSION 3.16)
project(MuonSpectrometerSim)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED gdml ui_all vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

# Find ROOT
find_package(ROOT REQUIRED CONFIG)

# Add ROOT flags
include(${ROOT_USE_FILE})

##include(${Geant4_USE_FILE})
##include_directories(include)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(MuonSpectrometerSim
	MuonSpectrometerSim.cc
	${sources}
	${headers})

target_include_directories(MuonSpectrometerSim PRIVATE
    ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(MuonSpectrometerSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
