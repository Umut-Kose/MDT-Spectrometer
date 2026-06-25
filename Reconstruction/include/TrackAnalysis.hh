#ifndef TRACKANALYSIS_HH
#define TRACKANALYSIS_HH

#include <vector>
#include <memory>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TTree.h"
#include <TGeoManager.h>  
#include <tuple>

// Include complete GENFIT headers
#include <AbsBField.h>
#include <RKTrackRep.h>
#include <FieldManager.h>
#include <Track.h>
#include <MeasuredStateOnPlane.h>
#include <PlanarMeasurement.h>

class TrackAnalysis {
public:
    struct Hit {
        int stationID;
        int layerID;
        int trackID; // GEANT track ID
        int pdg; // PDG code of the particle
        TVector3 position; // mm
        TVector3 momentum; // MeV/c
    };

    struct Tracklet {
        int station;
        TVector3 position;      // Center position in station
        TVector3 dir;           // Unit direction vector
        int pdg;                // PDG code of the particle
        int trackID;            // GEANT track ID
        double truthMomentum;   // Momentum in GeV/c
    };

    struct Result {
        int eventID;
        int trackID;
        int pdg; // PDG code of the particle

        double truthMomentum;   // Momentum from truth information
        
        double sagittaMomentum; // Momentum from sagitta method
        double sagittaMomentum_trk; // Momentum from sagitta method for tracklets        
        
        double kasaMomentum; // Momentum from Kasa circle fit
        double kasaMomentum_trk;    // Momentum from Kasa circle fit for tracklets

        double taubinMomentum; // Momentum from Taubin circle fit
        double taubinMomentum_trk;   // Momentum from Taubin circle fit for tracklets

        double trackletDeflectionMomentum; // Momentum from tracklet deflection
        double trackletDeflectionMomentum_trk; // Momentum from tracklet deflection for track

        double kalmanMomentum;  // Momentum from Kalman filter        
        double kalmanMomentumErr; // Error on Kalman momentum
        double kalmanChi2;  // Chi2 from Kalman filter
        
        double genfitMomentum;  // Momentum from GENFIT
        double genfitChi2;  // Chi2 from GENFIT
        int genfitCharge;         // +1 or -1
        double genfitMomentumErr; // Error on GENFIT momentum

        // Add charge fields for each method
        int sagittaCharge = 0;           // Sign from sagitta method
        int sagittaCharge_trk = 0;       // Sign from sagitta method for tracklets
        int kasaCharge = 0;              // Sign from Kasa circle fit
        int kasaCharge_trk = 0;          // Sign from Kasa circle fit for tracklets
        int taubinCharge = 0;            // Sign from Taubin circle fit
        int taubinCharge_trk = 0;        // Sign from Taubin circle fit for tracklets
        int trackletDeflectionCharge = 0; // Sign from tracklet deflection
        int trackletDeflectionCharge_trk = 0; // Sign from tracklet deflection for tracklets
        int kalmanCharge = 0;            // Sign from Kalman filter
    
         // Resolution fields (momentum errors)
        double sagittaResolution = -1;
        double sagittaResolution_trk = -1;
        double kasaResolution = -1;
        double kasaResolution_trk = -1;
        double taubinResolution = -1;
        double taubinResolution_trk = -1;
        double trackletDeflectionMomentumErr = -1;
        double trackletDeflectionMomentumErr_trk = -1;

    };

struct CircleFitResult {
    double yc, zc;   // center in meters
    double R;        // radius in meters
    double sigma;    // RMS error (m)
    int ndf;         // N-3
    double p;       // momentum estimate (GeV/c)
    double p_err;   // error on momentum (GeV/c)
    int charge;     // +1 or -1
    bool ok;
};

struct State {
    double y;
    double slope;
    double p; // momentum (GeV)
    double z;
    double cov[3][3]; // Covariance matrix
};

    int verbose = 0; // 0 = silent, 1 = info, 2 = debug, etc.

    static constexpr double MagnetFieldStrength = 1.8; // Tesla changed from 1.8T
    const double IronThickness = 0.15; // meters, default iron thickness


    const double mm_to_cm = 0.1;
    const double cm_to_mm = 10.0;
    const double MeV_to_GeV = 0.001;
    const double GeV_to_MeV = 1000.0;

    double GetIronThickness() const;  // in meters
    double GetMagneticField() const;  // in Tesla

    //TrackAnalysis();
    ~TrackAnalysis();

    void ProcessFile(const char* filename);
    void WriteResults(const char* filename);

    TrackAnalysis(TGeoManager* geo = nullptr);
    void BuildGeometryMaps();  // Call in constructor
    double GetFieldX(double x, double y, double z) const;
    double GetFieldBX(double y) const;
    bool IsInMagnet(double z) const;
    bool IsMagnetExist(double z1,double z2) const;


private:
    class MagneticField : public genfit::AbsBField {
    public:
        MagneticField(const TrackAnalysis* ta) : ta_(ta) {}

        TVector3 get(const TVector3& position) const override {
            //std::cout << "[MAGFIELD] Query at y=" << position.Y() << std::endl;
            double y = position.Y();
            double B = ta_->MagnetFieldStrength; 
            if (std::abs(y) >= 25.0 && std::abs(y) <= 50.0) return TVector3(B, 0, 0);
            if (std::abs(y) < 25.0) return TVector3(-1*B, 0, 0);
            return TVector3(0, 0, 0);
        }
        
        genfit::AbsBField* clone() const { 
            return new MagneticField(*this); 
        }
    private:
        const TrackAnalysis* ta_;

    };

    class GeoField : public genfit::AbsBField {
    public:
        // Take a pointer to the parent analysis so it can access GetFieldX
        GeoField(const TrackAnalysis* ta) : ta_(ta) {}
        TVector3 get(const TVector3& position) const override 
        {
            double x_mm = position.X() * 10.0; // cm to mm
            double y_mm = position.Y() * 10.0; // cm to mm  
            double z_mm = position.Z() * 10.0; // cm to mm
            
            double bx = ta_->GetFieldX(x_mm, y_mm, z_mm); // in Tesla
            
            // ADD DEBUG OUTPUT HERE:
            //if (ta_->verbose >= 2) {
            //    std::cout << "[GEOFIELD] Query at (x=" << x_mm << ", y=" << y_mm 
            //            << ", z=" << z_mm << ") mm --> Bx=" << bx << " T" << std::endl;
            // }
            
            return TVector3(bx, 0, 0);
        }
    
        genfit::AbsBField* clone() const { return new GeoField(*this); }

    private:
        const TrackAnalysis* ta_;
    };

    TGeoManager* geo_; // Geometry pointer
    std::vector<std::pair<double,double>> magnet_z_ranges_; // zmin, zmax
    std::vector<double> scifi_layer_z_; // z position for each SciFi layer
    //genfit::AbsBField* magneticField_; // Not MagneticField*


    void ProcessEvent(int eventID, 
                        const std::vector<int>& trackID,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z,
                        const std::vector<double>& px,
                        const std::vector<double>& py,
                        const std::vector<double>& pz,
                        const std::vector<int>& layerID,
                        const std::vector<int>& stationID,
                        const std::vector<int>& pdg);

    void RunKalmanFilter(const std::vector<Hit>& hits,double seed_p);
    void RunKalmanFilterTracklets(const std::vector<Tracklet>& hits, double p_seed= 30.0);
    void RunGenfit(const std::vector<Hit>& hits);
    void RunGenfitTracklets(const std::vector<Tracklet>& tracklets, double p_seed= 30.0);
    void RunGenfitGlobalHits(const std::vector<Hit>& hits, double seed_p=30.0);
    double CalculateSagittaMomentum(const std::vector<Tracklet>& hits);
    double CalculateSagittaMomentum(const std::vector<Hit>& hits);
    void PrintOut(const std::vector<Hit>& hits);
    void PrintStationFits(const std::vector<Hit>& hits);
    //void FitStationTracklets(const std::vector<Hit>& hits);    
    
    
    void PlotHitsAndTrack(const std::vector<Hit>& hits,
                      const std::vector<TVector3>& fittedTrackPoints,
                      const std::string& title);
    void PlotMCViews(const std::vector<Hit>& hits, const std::string& basename);
    void PlotHitsAndTracklets(const std::vector<Hit>& hits, const std::vector<Tracklet>& tracklets, const std::string& title);
    std::vector<Tracklet> FitStationTracklets(const std::vector<Hit>& hits);
    void PlotTrackletsViews(const std::vector<Hit>& hits, const std::vector<Tracklet>& tracklets, const std::string& basename);

    void ReconstructGlobalTrack(const std::vector<Tracklet>& tracklets);
    void PlotGlobalTrackletFit(const std::vector<Tracklet>& tracklets, const std::string& basename);
    double MomentumFromSlopeChange(const std::vector<TVector3>& avgs,
            std::function<double(double, double, double)> GetField,
            const std::vector<std::pair<double, double>>& magnet_z_ranges);
    


    void PlotCircleFitYZ(const std::vector<Hit>& hits, double yc, double zc, double R, const std::string& title);
        double KinkMomentum_Bx(const std::vector<Hit>& hits, double bx_tesla, double iron_thickness_mm);        
    std::vector<Result> results_;
    genfit::AbsTrackRep* trackRepMinus_;
    genfit::AbsTrackRep* trackRepPlus_;

    //MagneticField* magneticField_;
    genfit::AbsBField* magneticField_; // Use AbsBField for GENFIT compatibility
    double EstimateMomentumKickMethod(
        const std::vector<Tracklet>& tracklets, 
        const std::vector<std::pair<int,int>>& iron_regions, 
        double B, double iron_thickness);
   // std::pair<double, double> 
    std::tuple<double, double, int> MomentumFromTrackletDeflection(const std::vector<Tracklet>& tracklets,double bx_field, double iron_thickness_m,bool use_3D_angle); 
   double GlobalDirectionDeflectionMomentum(const std::vector<Tracklet>& tracklets,double bx_field, double iron_thickness_m);      
  
    // Circle fit in (y, z) plane for Bx field using KAsa and Taubin methods
    CircleFitResult CircleFitYZ_Taubin(const std::vector<Tracklet>& tracklets, double bx_field = 0.0);
    CircleFitResult CircleFitYZ_Taubin(const std::vector<Hit>& hits, double bx_field = 0.0);
    CircleFitResult CircleFitYZ_Kasa(const std::vector<Tracklet>& tracklets, double bx_field);
    CircleFitResult CircleFitYZ_Kasa(const std::vector<Hit>& hits, double bx_field);


  std::vector<TVector3> SmoothedTrackletPositions(const std::vector<Tracklet>& tracklets);
void LinearFit(const std::vector<TVector3>& pts, double& slope, double& intercept, char coord);

};

#endif