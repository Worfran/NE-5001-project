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
// $Id: PhantomHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file PhantomHit.hh
/// \brief Definition of the PhantomHit class

#ifndef PhantomHit_h
#define PhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Phantom hit
///
/// It records:
/// - the voxel ID
/// - the dose

// =============================================
// Exercise 2 Step 1:
// Create a hit class for the phantom
//    A hit is characterized by: an index (which phantom
//    voxel has fired, and the time of the fire
//    Complete this class as appropriate
//  The following methods should be preapred:
//    - Constructor with two parameters (index and time)
//    - Operator new and delete that use an allocator
//    - Setters and Getters for index and dose
//    - A Print method that dumps on G4cout the information
//      contained in the hit

class PhantomHit : public G4VHit
{
public:
    PhantomHit(G4int i,G4double dose);
    virtual ~PhantomHit() {}
    
    inline void *operator new(size_t);
    inline void operator delete(void*aHit);

    void Print();
    
    G4int GetID() const { return fId; }

    void SetDose(G4double val) { fDose = val; }
    G4double GetDose() const { return fDose; }
    
private:
    G4int fId;
    G4double fDose;
};

typedef G4THitsCollection<PhantomHit> PhantomHitsCollection;

extern G4ThreadLocal G4Allocator<PhantomHit>* PhantomHitAllocator;

inline void* PhantomHit::operator new(size_t)
{
    if (!PhantomHitAllocator)
        PhantomHitAllocator = new G4Allocator<PhantomHit>;
    return (void*)PhantomHitAllocator->MallocSingle();
}

inline void PhantomHit::operator delete(void*aHit)
{
    PhantomHitAllocator->FreeSingle((PhantomHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
