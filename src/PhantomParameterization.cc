//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "PhantomParameterization.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhantomParameterization
::PhantomParameterization(const G4ThreeVector& voxelSize,
                                    G4int nz,
                                    std::vector<G4Material*>& mat):
  G4VNestedParameterisation(),
  fdX(voxelSize.x()),fdY(voxelSize.y()),fdZ(voxelSize.z()),
  fNz(nz),fMat(mat)
{
  // Pre-calculate the positions of all the voxels.
  // X and Y positions are already defined in DetectorConstruction
  // by using replicated volume.
  // Z position is calculated here based on the Z copy number.
  fpZ.clear();
  G4double zp;
  for ( G4int iz = 0; iz < fNz; iz++){
    zp = (-fNz+1+2*iz)*fdZ;
    fpZ.push_back(zp);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhantomParameterization::~PhantomParameterization(){
  fpZ.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Material assignment to geometry.
//
G4Material* PhantomParameterization
::ComputeMaterial(G4VPhysicalVolume* currentVol, const G4int copyNo,
                  const G4VTouchable* parentTouch)
{
  if(parentTouch==0) return fMat[0]; // protection for initialization and
                                     // vis at idle state
  // Copy number of voxels. 
  // Copy number of X and Y are obtained from replication number.
  // Copy nymber of Z is the copy number of current voxel.
  G4int ix = parentTouch->GetReplicaNumber(0);
  G4int iy = parentTouch->GetReplicaNumber(1);
  G4int iz = copyNo;
	
  // If we were actually reading in a patient scan, we would use information
  // from that scan to assign the appropriate material to each voxel.
  // But for this simple tutorial example, we just hard code that:
  // center voxels are bone, colored red
  // voxels around those are water, colored blue
  // voxels around those are air, colored grey
  G4Material* mat = 0;
  G4VisAttributes* visAttributes = 0;
	
  if (iz > 1 && iz < 4 && iy > 3 && iy < 6 && ix > 3 && ix < 6) {
	mat = fMat[2];
    visAttributes = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  } else if (iz > 0 && iz < 5  && iy > 1 && iy < 8 && ix > 1 && ix < 8) {
    mat = fMat[1];
    visAttributes = new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.02));
  } else {
    mat = fMat[0];
    visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.02));
  }
  
  if (iz >= 10 && iz < 12 && iy >= 10 && iy < 13 && ix >= 10 && ix < 13) {
    // tumor material (red)
    mat = fMat[1];
    visAttributes = new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.1));
  } else {
    mat = fMat[0];
    visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.02));
  }

  visAttributes -> SetForceSolid(true);
  currentVol->GetLogicalVolume()->SetVisAttributes(visAttributes);

  return mat;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// This method is called by Geant4's material scanner so that it can
// know how much memory to set aside for the materials tables
//
G4int PhantomParameterization::GetNumberOfMaterials() const{
  return fMat.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// This method is called by Geant4's material scanner so that it can
// tell what materials to instantiate
//
G4Material* PhantomParameterization::GetMaterial(G4int i) const{
  return fMat[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Translate to the appropriate one of the pre-calculated positions.
//
void PhantomParameterization
::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const{
  G4ThreeVector position(0.,0.,fpZ[copyNo]);
  physVol->SetTranslation(position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// This example has no variation in the dimensions of the voxels.
//
void PhantomParameterization
::ComputeDimensions(G4Box& box, const G4int, const G4VPhysicalVolume* ) const{
  box.SetXHalfLength(fdX);
  box.SetYHalfLength(fdY);
  box.SetZHalfLength(fdZ);
}
