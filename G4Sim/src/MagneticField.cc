#include "MagneticField.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>


void MagneticField::GetFieldValue(const G4double point[4], G4double* Bfield) const {
    // point[0] = x, point[1] = y, point[2] = z

    // Assume magnet box extends: x ∈ [-500, 500] mm, y ∈ [-500, 500] mm, z ∈ [–25, +25] mm (centered)
    G4double x = point[0];
    G4double y = point[1];
    G4double z = point[2];

    // Default: no field
    Bfield[0] = Bfield[1] = Bfield[2] = 0.0;

    // Assume ±1.5 Tesla in steel, depending on y (top/bottom vs center)
    if (std::abs(y) >= 250. * mm && std::abs(y) <= 500. * mm) {
      // Top or Bottom: +1.5 T
      Bfield[0] = +1.5 * tesla;
    } else if (std::abs(y) < 250. * mm) {
      // Middle: –1.5 T
      Bfield[0] = -1.5 * tesla;
    }

    
}
/*
// previously used code
void MagneticField::GetFieldValue(const G4double point[4], G4double* Bfield) const {
    // point[0] = x, point[1] = y, point[2] = z

    // Assume magnet box extends: x ∈ [-500, 500] mm, y ∈ [-500, 500] mm, z ∈ [–25, +25] mm (centered)
    G4double x = point[0];
    G4double y = point[1];
    G4double z = point[2];

    // Default: no field
    Bfield[0] = Bfield[1] = Bfield[2] = 0.0;

    // Assume ±1.5 Tesla in steel, depending on y (top/bottom vs center)
    if (std::abs(y) >= 250. * mm && std::abs(y) <= 500. * mm) {
      // Top or Bottom: +1.5 T
      Bfield[0] = +1.5 * tesla;
    } else if (std::abs(y) < 250. * mm) {
      // Middle: –1.5 T
      Bfield[0] = -1.5 * tesla;
    }

    
}
*/