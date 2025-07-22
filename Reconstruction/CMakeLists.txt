cmake_minimum_required(VERSION 3.16)
project(MuonTrackingAnalysis)


# ------------------------------------------------------------------------------
# Find ROOT
# ------------------------------------------------------------------------------
#find_package(ROOT REQUIRED CONFIG)
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist Minuit2 Geom Gdml)
include(${ROOT_USE_FILE})

# ------------------------------------------------------------------------------
# GENFIT Config
# ------------------------------------------------------------------------------
set(GENFITINCDIR /Users/ukose/sw/kits/GenFit-install/include)
set(GENFITLIBDIR /Users/ukose/sw/kits/GenFit-install/lib)

include_directories(${GENFITINCDIR})
link_directories(${GENFITLIBDIR})

add_executable(MuonTrackingAnalysis
    src/main.cc
    src/TrackAnalysis.cc
)

target_include_directories(MuonTrackingAnalysis PRIVATE
    ${ROOT_INCLUDE_DIRS}
    ${GENFITINCDIR}
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(MuonTrackingAnalysis
    ${ROOT_LIBRARIES}
    genfit2
)
