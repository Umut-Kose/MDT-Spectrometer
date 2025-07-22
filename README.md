# Studying magnet options for FASERCAL


Code to simulate and analyze muons in Magnet


# G4Sim

Consideration:

Each module consists of:

1 Magnet Plane:

Made of ARMCO steel, 100×100 cm², 3 cm thick

Two horizontal slits (50 cm long anf 1cm wide, 25 cm away from Top/Bottom and centered)

Magnetized with ±1.5 T uniform field (top: +1.5 T, center: –1.5 T, bottom: +1.5 T)

Aluminum coil strip (5 cm wide and 4 mm thick)

1 SciFi Tracker Plane (4 layers/plane) immediately after the magnet

Position resolution: typically 50 um

30 cm air gap between modules

There are 5 such modules in total.

cd build
cmake -DGeant4_DIR=$G4LIBDIR ..
make -j

./MuonSpectrometerSim macro/run_1000events.mac
getting scifi_hits.root output file


# Display

using TEve root class

compiling the code
make clean
make distclean
make

./eve_display scifi_hits.root ../G4Sim/Geometry/detector_geometry.gdml

# Analysis

still need to be developed

