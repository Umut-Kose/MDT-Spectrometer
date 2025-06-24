#ifndef MAGNETICFIELD_HH
#define MAGNETICFIELD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

class MagneticField : public G4MagneticField {
public:
    MagneticField() = default;
    ~MagneticField() = default;

    virtual void GetFieldValue(const G4double point[4], G4double* Bfield) const override;
};

#endif

