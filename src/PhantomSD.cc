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
// $Id: PhantomSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file PhantomSD.cc
/// \brief Implementation of the PhantomSD class

#include "PhantomSD.hh"
#include "PhantomHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4VSolid.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhantomSD::PhantomSD(G4String name, G4int nxVoxels, G4int nzVoxels)
: G4VSensitiveDetector(name),
fnxVoxels(nxVoxels), fnzVoxels(nzVoxels),
fHitsCollection(0), fHCID(-1)
{
    G4String HCname = "phantomColl";
    collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhantomSD::~PhantomSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhantomSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new PhantomHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PhantomSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    // =============================================
    // Exercise 2 Step 2:
    // Collect dose information from G4Step.
    // First: if there is no energy deposition in the voxel, the detector
    // did not trigger. Just return.
    // Second: Get copyNo and dose from Step. You should use PreStepPoint
    // Third: check if hits collection has already an entry for given voxel.
    //   If not create a new entry for this voxel and add it to the detector,.
	//   if there has already been some dose recorded in this voxel,
    //   just accumulate the new dose into the existing entry.

	G4double edep = step->GetTotalEnergyDeposit();
	if (edep==0.) return true;
	
	G4StepPoint* preStepPoint = step->GetPreStepPoint();

    G4TouchableHistory* touchable
      = (G4TouchableHistory*)(preStepPoint->GetTouchable());
	
	G4int iY = touchable->GetReplicaNumber(2);
	G4int iX = touchable->GetReplicaNumber(1);
	G4int iZ = touchable->GetReplicaNumber(0);

	G4int copyNo = iY * fnxVoxels * fnzVoxels + iX * fnzVoxels + iZ;

	G4cout << "iY: " << iY << ", iX:" << iX << ", iZ: " << iZ << ", copyNo: " << copyNo << G4endl;
	
	G4double density = preStepPoint->GetMaterial()->GetDensity();
	G4double dose = edep / ( density * touchable->GetSolid()->GetCubicVolume() );
	
    // check if this voxel already has a hit
    G4int ix = -1;
    for (G4int i=0;i<fHitsCollection->entries();i++)
    {
        if ((*fHitsCollection)[i]->GetID()==copyNo)
        {
            ix = i;
            break;
        }
    }

    if (ix>=0)
        // if it has, then accumulate from previous dose
    {
		dose += (*fHitsCollection)[ix]->GetDose();
        (*fHitsCollection)[ix]->SetDose(dose);
    } else
        // if not, create a new hit and set it to the collection
    {
        PhantomHit* hit = new PhantomHit(copyNo,dose);
        fHitsCollection->insert(hit);
    }    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
