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
// $Id: DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "PhantomParameterization.hh"
#include "PhantomSD.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(0),
  fHodoscope1Logical(0),
  fWirePlane1Logical(0),
  fVisAttributes(),
  fArmAngle(0.*deg), fArmRotation(0), fSecondArmPhys(0)

{
    fArmRotation = new G4RotationMatrix();
    fArmRotation->rotateY(fArmAngle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fArmRotation;
    delete fMessenger;
    
    for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Construct materials
    ConstructMaterials();
    G4Material* air = G4Material::GetMaterial("G4_AIR");
    G4Material* argonGas = G4Material::GetMaterial("G4_Ar");
    G4Material* scintillator 
      = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4bool checkOverlaps = true;

    
    // geometries --------------------------------------------------------------
    // experimental hall (world volume)
    G4VSolid* worldSolid 
      = new G4Box("worldBox",1.*m,1.*m,1.*m);
    G4LogicalVolume* worldLogical
      = new G4LogicalVolume(worldSolid,air,"worldLogical");
    G4VPhysicalVolume* worldPhysical
      = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                          false,0,checkOverlaps);
    
    // first arm
    G4VSolid* firstArmSolid 
      = new G4Box("firstArmBox",0.5*m,0.5*m,0.5*m);
    G4LogicalVolume* firstArmLogical
      = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,-30*cm),firstArmLogical,
                      "firstArmPhysical",worldLogical,
                      false,0,checkOverlaps);
    
    // second arm
    G4VSolid* secondArmSolid 
      = new G4Box("secondArmBox",10.*cm,10.*cm,10*cm); // TO DO: adjust size as needed currently generating overlaps and outside range
    G4LogicalVolume* secondArmLogical
      = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical");
    G4double x = -10.*cm * std::sin(fArmAngle);
    G4double z = 10.*cm * std::cos(fArmAngle);
    fSecondArmPhys
      = new G4PVPlacement(fArmRotation,G4ThreeVector(x,0.,z),secondArmLogical,
                          "fSecondArmPhys",worldLogical,
                          false,0,checkOverlaps);
    
    
    // =============================================
    // Target absorber in second arm
    // =============================================
    G4Material* absorberMaterial = G4Material::GetMaterial("ADIPOSE_TISSUE_ICRP");
    
    // Define the phantom size and voxelization parameters
    G4ThreeVector phantomSize = G4ThreeVector(8.*cm, 8.*cm, 8.*cm); // 8 cm cube
    G4int nVoxelsX = 20; // Number of voxels along X (8 cm / 0.4 cm)
    G4int nVoxelsY = 20; // Number of voxels along Y (8 cm / 0.4 cm)
    G4int nVoxelsZ = 20; // Number of voxels along Z (8 cm / 0.4 cm)

    // Calculate voxel size
    G4ThreeVector voxelSize;
    voxelSize.setX(phantomSize.x() / nVoxelsX);
    voxelSize.setY(phantomSize.y() / nVoxelsY);
    voxelSize.setZ(phantomSize.z() / nVoxelsZ);

    // Create the overall phantom volume
    G4VSolid* phantomSolid = new G4Box("PhantomBox", phantomSize.x() / 2., phantomSize.y() / 2., phantomSize.z() / 2.);
    G4LogicalVolume* fPhantomLogical = new G4LogicalVolume(phantomSolid, absorberMaterial, "PhantomLogical");
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fPhantomLogical, "PhantomPhysical", worldLogical, false, 0, checkOverlaps);

    // Slice the phantom into Y slices using the Replica Volume technique
    G4String yRepName("RepY");
    G4VSolid* solYRep = new G4Box(yRepName, phantomSize.x() / 2., voxelSize.y() / 2., phantomSize.z() / 2.);
    G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep, absorberMaterial, yRepName);
    new G4PVReplica(yRepName, logYRep, fPhantomLogical, kYAxis, nVoxelsY, voxelSize.y());

    // Further slice these along X using the Replica Volume technique
    G4String xRepName("RepX");
    G4VSolid* solXRep = new G4Box(xRepName, voxelSize.x() / 2., voxelSize.y() / 2., phantomSize.z() / 2.);
    G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep, absorberMaterial, xRepName);
    new G4PVReplica(xRepName, logXRep, logYRep, kXAxis, nVoxelsX, voxelSize.x());

    // Define a single voxel box
    G4String zVoxName("Voxel");
    G4VSolid* solVoxel = new G4Box(zVoxName, voxelSize.x() / 2., voxelSize.y() / 2., voxelSize.z() / 2.);
    G4LogicalVolume* fVoxelLogical = new G4LogicalVolume(solVoxel, absorberMaterial, zVoxName);

    // Create a vector containing the materials for the phantom
    std::vector<G4Material*> phantomMat(3);
    phantomMat[0] = absorberMaterial;
    phantomMat[1] = absorberMaterial;
    phantomMat[2] = absorberMaterial;

    // Define the parameterization for the phantom
    PhantomParameterization* paramPhantom = new PhantomParameterization(voxelSize / 2., nVoxelsZ, phantomMat);

    // Apply the parameterization to the Z slices
    new G4PVParameterised("PhantomVoxels", fVoxelLogical, logXRep, kUndefined, nVoxelsZ, paramPhantom);
    
    
    // visualization attributes ------------------------------------------------
    
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    visAttributes->SetVisibility(false);
    worldLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    visAttributes->SetVisibility(false);
    firstArmLogical->SetVisAttributes(visAttributes);
    secondArmLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.8888,0.8888,0.0));
    visAttributes->SetVisibility(false);
    fVisAttributes.push_back(visAttributes);
    
    // return the world physical volume ----------------------------------------
    
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // sensitive detector -----------------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;
 
  G4VSensitiveDetector* phantomSD = new PhantomSD(SDname="/phantom", 20, 20);
  SDman->AddNewDetector(phantomSD);
  fVoxelLogical->SetSensitiveDetector(phantomSD);
 
 
  // Register the field and its manager for deleting
  G4AutoDelete::Register(fMagneticField);
  G4AutoDelete::Register(fFieldMgr); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    // Air 
    nistManager->FindOrBuildMaterial("G4_AIR");
  
    // Argon gas
    nistManager->FindOrBuildMaterial("G4_Ar");

    // Scintillator
    // (PolyVinylToluene, C_9H_10)
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    // CsI
    // =============================================
    // Exercise 1a
    //    Create material Cesium Iodide.
    //  Chemical formula CsI.
    // Some properties:
    //    Cs : A=55 Zeff=132.9*g/mol
    //     I : A=53 , Zeff=126.9*g/mol
    // Density of christal of CsI is rho=4.51*g/cm^3
    // =============================================

    G4Element* el_Cs = new G4Element("Cesium","Cs",55.,132.9*g/mole);
    G4Element* el_I  = new G4Element("Iodine","I",53.,126.9*g/mole);
    G4Material* mat_CsI = new G4Material("CsI", 4.51*g/cm3,2);
    mat_CsI->AddElement(el_Cs,1);
    mat_CsI->AddElement(el_I,1);
    
    // Lead
    // =============================================
    // Exercise 1b
    //    Create material Lead from Nist Manager.
    // Note that it is actually a mixture of several isotopes.
    // If you want to build it by hand starting from isotopes you
    // can:
    //    G4Isotope* iso1= new G4Isotope( ... )
    //    ....
    //    G4Element* elem = new G4Element("name","symbol",nisotopes);
    //    elem->AddIsotope( iso1 );
    //    elem->AddIsotope( iso2 );
    //    ...
    // =============================================
    nistManager->FindOrBuildMaterial("G4_Pb");

    // ==========================================================
    // Custom material from user-provided composition (mass fractions)
    // Composition (Atomic number : fraction by weight):
    //  1 : 0.119477
    //  6 : 0.637240
    //  7 : 0.007970
    //  8 : 0.232333
    // 11 : 0.000500
    // 12 : 0.000020
    // 15 : 0.000160
    // 16 : 0.000730
    // 17 : 0.001190
    // 19 : 0.000320
    // 20 : 0.000020
    // 26 : 0.000020
    // 30 : 0.000020
    // ==========================================================
    G4Element* el_H  = nistManager->FindOrBuildElement(1);  
    G4Element* el_C  = nistManager->FindOrBuildElement(6);  
    G4Element* el_N  = nistManager->FindOrBuildElement(7);  
    G4Element* el_O  = nistManager->FindOrBuildElement(8);  
    G4Element* el_Na = nistManager->FindOrBuildElement(11); 
    G4Element* el_Mg = nistManager->FindOrBuildElement(12); 
    G4Element* el_P  = nistManager->FindOrBuildElement(15); 
    G4Element* el_S  = nistManager->FindOrBuildElement(16); 
    G4Element* el_Cl = nistManager->FindOrBuildElement(17); 
    G4Element* el_K  = nistManager->FindOrBuildElement(19); 
    G4Element* el_Ca = nistManager->FindOrBuildElement(20); 
    G4Element* el_Fe = nistManager->FindOrBuildElement(26); 
    G4Element* el_Zn = nistManager->FindOrBuildElement(30); 

    // ADIPOSE TISSUE (ICRP)
    // Density (g/cm3) = 9.20000E-01
    G4double density = 9.2e-01 * g/cm3; 
    G4Material* ICRP = new G4Material("ADIPOSE_TISSUE_ICRP", density, 13);

    // Add elements with mass fractions (fractions by weight sum should be ~1.0)
    ICRP->AddElement(el_H,  0.119477);
    ICRP->AddElement(el_C,  0.637240);
    ICRP->AddElement(el_N,  0.007970);
    ICRP->AddElement(el_O,  0.232333);
    ICRP->AddElement(el_Na, 0.000500);
    ICRP->AddElement(el_Mg, 0.000020);
    ICRP->AddElement(el_P,  0.000160);
    ICRP->AddElement(el_S,  0.000730);
    ICRP->AddElement(el_Cl, 0.001190);
    ICRP->AddElement(el_K,  0.000320);
    ICRP->AddElement(el_Ca, 0.000020);
    ICRP->AddElement(el_Fe, 0.000020);
    ICRP->AddElement(el_Zn, 0.000020);
    
    
    // Important: Never use a real gas with 0 density as materials.
    //            Physics processes requires density>0 and may give wrong
    //            results if they encounter 0 density material.
    //            Instead use one of the following methods if you need
    //            "vacuum" (e.g. a beam pipe internal)
    // Vacuum "Galactic"
    // nistManager->FindOrBuildMaterial("G4_Galactic");

    // Vacuum "Air with low density"
    // G4Material* air = G4Material::GetMaterial("G4_AIR");
    // G4double density = 1.0e-5*air->GetDensity();
    // nistManager
    //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
