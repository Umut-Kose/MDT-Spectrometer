# Studying MDT magnet options for FASERCAL


Code to simulate and analyze muons in Magnet


# G4Sim

# ATLAS Monitored Drift Chamber
* 3-Layer honeycomb MDT chamber
* 34 tubes per layer
* 102 tubes per station
* Y pitch within the layer 30mm (tubes touching)
* Z pitch between layers 25.98 mm
# All tubes in all stations oriented along the X axis
* Tube length= 860 mm
* Tube diameter = 30 mm 
    (0.4 mm thick Al pipe, inner radius 14.6 mm)
* Aluminum endplates = 10 mm
* Tracking station thickness 102 mm
* Spatial resolution 80 um
# Magnet same as Baby-Mind 
* muons bend in YZ plane
* Each magnet 400 mm thick iron
* Magnetic field only on the iron
* Bx = -1.5 T for |y| < 250 mm; +1.5 T for 250<|y| < 500 mm
* B = 0 outside magnet volume

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

