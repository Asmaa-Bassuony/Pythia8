// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.
/*
#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << mult;
  return 0;
}


*/

//=======================================================================================================================================================



// File: hist.cc
// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.
// Copyright (C) 2014 Torbjorn Sjostrand

// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TCanvas.h"
#include "TFile.h"
using namespace Pythia8;

int main(int argc, char* argv[]) {

 int nListJets = 5;

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
 // pythia.readString("HardQCD:all = on");
  pythia.readString("SoftQCD:all = on");
  //pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Beams:eCM = 14000.");
  pythia.init();    //(2212 /* p */, 2212 /* p */, 13000 /* TeV */)

//===========================================================================================

// Common parameters for the two jet finders.
  double etaMax   = 7;
 // double etaMaxC  = 1.5;
 // double etaMaxF  = 6.6;
  double radius   = 0.7;
  double pTjetMin = 10.;
  // Exclude neutrinos (and other invisible) from study.
  int    nSel     = 2;
  // Range and granularity of CellJet jet finder.
  int    nEta     = 80;
  int    nPhi     = 64;
  int    nEvent   = 100000000;
  
   //===========================================================================================================//
 //============================== Tracker Region  =============================================================//
//===========================================================================================================//
// Set up SlowJet jet finder, with anti-kT clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);
  //SlowJet slowJet( -1, radius, pTjetMin, etaMaxC, nSel, 1);
 // SlowJet slowJet( -1, radius, pTjetMin, etaMaxF, nSel, 1);

  // Set up CellJet jet finder.
  CellJet cellJet( etaMax, nEta, nPhi, nSel);
 //  CellJet cellJet( etaMaxC, nEta, nPhi, nSel);
 // CellJet cellJet( etaMaxF, nEta, nPhi, nSel);
 
  
  // Create file on which histogram(s) can be saved.
   

// ================  Book histogram ======================================================================

   TH1F *mult      = new TH1F("mult2_5","charged multiplicity", 800, -0.5, 799.5);
   TH1F *mult11    = new TH1F("mult2_5111","charged multiplicity", 800, -0.5, 799.5);
   TH1F *mult05    = new TH1F("mult_5","charged multiplicity", 800, -0.5, 799.5);
   TH1F *mult1     = new TH1F("mult_1","charged multiplicity", 800, -0.5, 799.5);
   TH1F *mult1_5   = new TH1F("mult1_5","charged multiplicity", 800, -0.5, 799.5);
   TH1F *mult2     = new TH1F("mult2_0","charged multiplicity", 800, -0.5, 799.5);
   
  TH1F *eta1        = new TH1F("eta2_5","eta_All", 120, -3.0, 3.0);
  TH1F *eta         = new TH1F("eta_ALL","eta_All", 140, -7.0, 7.0);

//NEW work 12/2021
  TH1F *eta_1       = new TH1F("eta_1","eta_All", 250, -2.5, 2.5);
  TH1F *eta_2       = new TH1F("eta_2","eta_All", 250, -2.5, 2.5);
  TH1F *eta_3       = new TH1F("eta_3","eta_All", 250, -2.5, 2.5);
  TH1F *eta_4       = new TH1F("eta_4","eta_All", 250, -2.5, 2.5);
  TH1F *eta_5       = new TH1F("eta_5","eta_All", 250, -2.5, 2.5);
  TH1F *pT_05       = new TH1F("pT_0_5","pT_All", 2000, -0.5, 199.5);
  TH1F *pT_1        = new TH1F("pT_1_0","pT_All", 2000, -0.5, 199.5);
  TH1F *pT_15       = new TH1F("pT_1_5","pT_All", 2000, -0.5, 199.5);
  TH1F *pT_2        = new TH1F("pT_2_0","pT_All", 2000, -0.5, 199.5);
  TH1F *pT_25       = new TH1F("pT_2_5","pT_All", 2000, -0.5, 199.5);



  // Histograms. Note similarity in names, even when the two jet finders
  // do not calculate identically the same property (pT vs. ET, y vs. eta).
TH1F *nJetsSI     = new TH1F("Nu_Jet_Slow_inside","number of jets, SlowJet", 50, -0.5, 49.5);
TH1F *nJetsCI     = new TH1F("Nu_Jet_Cell_inside","number of jets, CellJet", 50, -0.5, 49.5);

TH1F *Mass_jetS   = new TH1F("invariant_mass_Slow", "invariant_massof this four-vector", 10000, -0.5,9999.5);
TH1F *Mass_jetC   = new TH1F("invariant_mass_Cell", "invariant_massof this four-vector", 10000, -0.5,9999.5);

TH1F *nJetsS     = new TH1F("Nu_Jet_Slow","number of jets, SlowJet", 50, -0.5, 49.5);
TH1F *nJetsC     = new TH1F("Nu_Jet_Cell","number of jets, CellJet", 50, -0.5, 49.5);
TH1F *nJetsD     = new TH1F("Nu_Jet_Cell_Slow","number of jets, CellJet - SlowJet", 45, -22.5, 22.5);
TH1F *JetMult_S  = new TH1F("Jet_mult_Slow", "Jet_mult_Slow", 500, -0.5, 499.5);
TH1F *JetMult_C  = new TH1F("Jet_mult_Cell", "Jet_mult_Cell", 500, -0.5, 499.5);
TH1F *eTjetsS    = new TH1F("eT_Slow","pT for jets, SlowJet", 100, 0., 500.);
TH1F *eTjetsC    = new TH1F("eT_Cell","eT for jets, CellJet", 100, 0., 500.);
TH1F *etaJetsS   = new TH1F("y_Slow","y for jets, SlowJet", 100, -5., 5.);
TH1F *etaJetsC   = new TH1F("eta_Cell","eta for jets, CellJet", 100, -5., 5.);
TH1F *phiJetsS   = new TH1F("phi_Slow","phi for jets, SlowJwt", 100, -M_PI, M_PI);
TH1F *phiJetsC   = new TH1F("Phi_Cell","phi for jets, CellJet", 100, -M_PI, M_PI);
TH1F *distJetsS  = new TH1F("R_Slow","R distance between jets, SlowJet", 100, 0., 10.);
TH1F *distJetsC  = new TH1F("R_Cell","R distance between jets, CellJet", 100, 0., 10.);
TH1F *eTdiffS    = new TH1F("pT_Diff_Slow","pT difference, SlowJet", 100, -100., 400.);
TH1F *eTdiffC    = new TH1F("et_diss_Slow","eT difference, CellJet", 100, -100., 400.);
  

 

  // Book histogram.
    TH1F *mult_Cen     = new TH1F("mult_Cen","charged multiplicity", 800, -0.5, 799.5);
    TH1F *eta_Cen      = new TH1F("eta_Cen","eta_All", 120, -3.0, 3.0);
    TH1F *pT_Cen       = new TH1F("pT-Cen","pT_All", 200, -0.5, 199.5);

  // Histograms. Note similarity in names, even when the two jet finders
  // do not calculate identically the same property (pT vs. ET, y vs. eta).

    TH1F *nJetsS_CenI     = new TH1F("Nu_Jet_Slow_Central_inside","number of jets, SlowJet", 50, -0.5, 49.5);
    TH1F *nJetsC_CenI     = new TH1F("Nu_Jet_Cell_Central_inside","number of jets, CellJet", 50, -0.5, 49.5);
    
    TH1F *nJetsS_Cen     = new TH1F("Nu_Jet_Slow_Central","number of jets, SlowJet", 50, -0.5, 49.5);
    TH1F *nJetsC_Cen     = new TH1F("Nu_Jet_Cell_Central","number of jets, CellJet", 50, -0.5, 49.5);
    
    TH1F *JetMult_SC     = new TH1F("Jet_mult_Slow_Cen", "Jet_mult_Slow", 500, 0.5, 499.5);
    TH1F *JetMult_CC     = new TH1F("Jet_mult_Cell_Cen", "Jet_mult_Cell", 500, 0.5, 499.5);
    TH1F *Mass_jetSC     = new TH1F("invariant_mass_Slow_Cen", "invariant_massof this four-vector", 10000, -0.5,9999.5);
    TH1F *Mass_jetCC     = new TH1F("invariant_mass__Cell_Cen", "invariant_massof this four-vector", 10000, -0.5,9999.5);
    
    TH1F *nJetsD_Cen     = new TH1F("Nu_Jet_Cell_Slow_Central","number of jets, CellJet - SlowJet", 45, -22.5, 22.5);
    TH1F *eTjetsS_Cen    = new TH1F("eT_Slow_Central","pT for jets, SlowJet", 100, 0., 500.);
    TH1F *eTjetsC_Cen    = new TH1F("eT_Cell_Central","eT for jets, CellJet", 100, 0., 500.);
    TH1F *etaJetsS_Cen   = new TH1F("y_Slow_Central","y for jets, SlowJet", 100, -5., 5.);
    TH1F *etaJetsC_Cen   = new TH1F("eta_Cell_Central","eta for jets, CellJet", 100, -5., 5.);
    TH1F *phiJetsS_Cen   = new TH1F("phi_Slow_Central","phi for jets, SlowJwt", 100, -M_PI, M_PI);
    TH1F *phiJetsC_Cen   = new TH1F("Phi_Cell_Central","phi for jets, CellJet", 100, -M_PI, M_PI);
    TH1F *distJetsS_Cen  = new TH1F("R_Slow_Central","R distance between jets, SlowJet", 100, 0., 10.);
    TH1F *distJetsC_Cen  = new TH1F("R_Cell_Central","R distance between jets, CellJet", 100, 0., 10.);
    TH1F *eTdiffS_Cen    = new TH1F("pT_Diff_Slow_Central","pT difference, SlowJet", 100, -100., 400.);
    TH1F *eTdiffC_Cen    = new TH1F("et_diss_Slow_Central","eT difference, CellJet", 100, -100., 400.);

   

  // Book histogram.
  TH1F *multF     = new TH1F("mult_Forward","charged multiplicity", 800, -0.5, 799.5);
  TH1F *etaF      = new TH1F("eta_All_Forward","eta_All", 120, 5.0, 7.0);
  TH1F *pTF       = new TH1F("pT-All_Forward","pT_All", 200, -0.5, 199.5);
 

  // Histograms. Note similarity in names, even when the two jet finders
  // do not calculate identically the same property (pT vs. ET, y vs. eta).
TH1F *nJetsSF     = new TH1F("Nu_Jet_Slow_Foward","number of jets, SlowJet", 50, -0.5, 49.5);
TH1F *nJetsCF     = new TH1F("Nu_Jet_Cell_Foward","number of jets, CellJet", 50, -0.5, 49.5);

TH1F *nJetsSFI     = new TH1F("Nu_Jet_Slow_Foward_inside","number of jets, SlowJet", 50, -0.5, 49.5);
TH1F *nJetsCFI     = new TH1F("Nu_Jet_Cell_Foward_inside","number of jets, CellJet", 50, -0.5, 49.5);

TH1F *JetMult_SF   = new TH1F("Jet_mult_Slow_Forward", "Jet_mult_Slow", 500, 0.5, 499.5);
TH1F *JetMult_CF   = new TH1F("Jet_mult_Cell_Forward", "Jet_mult_Cell", 500, 0.5, 499.5);
TH1F *Mass_jetSF   = new TH1F("invariant_mass_Slow_Forward", "invariant_massof this four-vector", 10000, -0.5,9999.5);
TH1F *Mass_jetCF   = new TH1F("invariant_mass_Cell_Forward", "invariant_massof this four-vector", 10000, -0.5,9999.5);

TH1F *nJetsDF     = new TH1F("Nu_Jet_Cell_Slow_Foward","number of jets, CellJet - SlowJet", 45, -22.5, 22.5);
TH1F *eTjetsSF    = new TH1F("eT_Slow_Foward","pT for jets, SlowJet", 100, 0., 500.);
TH1F *eTjetsCF    = new TH1F("eT_Cell_Foward","eT for jets, CellJet", 100, 0., 500.);
TH1F *etaJetsSF   = new TH1F("y_Slow_Foward","y for jets, SlowJet", 100, -5., 5.);
TH1F *etaJetsCF   = new TH1F("eta_Cell_Foward","eta for jets, CellJet", 100, -5., 5.);
TH1F *phiJetsSF   = new TH1F("phi_Slow_Foward","phi for jets, SlowJwt", 100, -M_PI, M_PI);
TH1F *phiJetsCF   = new TH1F("Phi_Cell_Foward","phi for jets, CellJet", 100, -M_PI, M_PI);
TH1F *distJetsSF  = new TH1F("R_Slow_Foward","R distance between jets, SlowJet", 100, 0., 10.);
TH1F *distJetsCF  = new TH1F("R_Cell_Foward","R distance between jets, CellJet", 100, 0., 10.);
TH1F *eTdiffSF    = new TH1F("pT_Diff_Slow_Foward","pT difference, SlowJet", 100, -100., 400.);
TH1F *eTdiffCF    = new TH1F("et_diss_Slow_Foward","eT difference, CellJet", 100, -100., 400.);

TH1F *mult_eta5   = new TH1F("mult_eta5","charged multiplicity", 800, -0.5, 799.5);
TH1F *mult_eta1   = new TH1F("mult_eta1","charged multiplicity", 800, -0.5, 799.5);
TH1F *mult_eta15  = new TH1F("mult_eta15","charged multiplicity", 800, -0.5, 799.5);
TH1F *mult_eta2   = new TH1F("mult_eta2","charged multiplicity", 800, -0.5, 799.5);
TH1F *mult_eta25  = new TH1F("mult_eta25","charged multiplicity", 800, -0.5, 799.5);
    
//======   Begin event loop. Generate event; skip if generation aborted ================

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles.
  //  int count     = 0;
  //  int nCharged  = 0;                 // for |eta| < 2.5
  //  int nChargedC = 0;                // for |eta| > 3.0 &&|eta| < 5.0
  //  int nChargedF = 0;                // for |eta| > 5.2.0 &&|eta| < 6.6
    
    int counter  = 0;                   // for |eta| < 0.5
    int counter1 = 0;                  // for |eta| < 1.0
    int counter2 = 0;                  // for |eta| < 1.5
    int counter3 = 0;                  // for |eta| < 2.0
    int counter4 = 0;                  // for |eta| < 2.5

    int counteta5   = 0;
    int counteta1   = 0;
    int counteta15  = 0;
    int counteta2   = 0;
    int counteta25  = 0;
    
    for (int i = 0; i < pythia.event.size(); ++i){
        
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()&&pythia.event[i].pT()> 0.1){
          
  

        if (abs(pythia.event[i].eta()) < 0.5  && pythia.event[i].pT()>= 0.4) counter++;
        if (abs(pythia.event[i].eta()) < 1.0  && pythia.event[i].pT()>= 0.4) counter1++; 
        if (abs(pythia.event[i].eta()) < 1.5  && pythia.event[i].pT()>= 0.4) counter2++;
        if (abs(pythia.event[i].eta()) < 2.0  && pythia.event[i].pT()>= 0.4) counter3++;
        if (abs(pythia.event[i].eta()) < 2.5  && pythia.event[i].pT()>= 0.4) counter4++; 

        if (abs(pythia.event[i].eta()) < 0.5  ) 
		{
		 	counteta5++;		
			pT_05 ->Fill(pythia.event[i].pT());
			}
        if (abs(pythia.event[i].eta()) < 1.0  ) 
		{
		 	counteta1++;
			pT_1  ->Fill(pythia.event[i].pT()); 
			}
        if (abs(pythia.event[i].eta()) < 1.5  )
		{
		 	counteta15++;
		 	pT_15 ->Fill(pythia.event[i].pT());
			}

        if (abs(pythia.event[i].eta()) < 2.0  )
		{
		 	counteta2++;
 			pT_2  ->Fill(pythia.event[i].pT());
			}
        if (abs(pythia.event[i].eta()) <= 2.4 &&pythia.event[i].pT()>= 0.1  ) 
		{
		 	counteta25++;
			pT_25 ->Fill(pythia.event[i].pT());   	 
        		}
 	if (abs(pythia.event[i].eta()) < 2.5 && pythia.event[i].pT()> 0.1) eta_1  -> Fill(pythia.event[i].eta());
	if (abs(pythia.event[i].eta()) < 2.5 && pythia.event[i].pT()> 0.2) eta_2  -> Fill(pythia.event[i].eta());
	if (abs(pythia.event[i].eta()) < 2.5 && pythia.event[i].pT()> 0.3) eta_3  -> Fill(pythia.event[i].eta());
	if (abs(pythia.event[i].eta()) < 2.5 && pythia.event[i].pT()> 0.4) eta_4  -> Fill(pythia.event[i].eta());
	if (abs(pythia.event[i].eta()) < 2.5 && pythia.event[i].pT()> 0.5) eta_5  -> Fill(pythia.event[i].eta());
        
 /*   if ( abs(pythia.event[i].eta())>3. && pythia.event[i].eta() < 5.)
        {
          nChargedC++;
          eta      -> Fill(pythia.event[i].eta());
          eta_Cen  -> Fill(pythia.event[i].eta());
          pT_Cen   -> Fill(pythia.event[i].pT());
        }
        
    if ( abs(pythia.event[i].eta())>5.2&& abs(pythia.event[i].eta())<6.6)
        {
           nChargedF++;
           eta         -> Fill(pythia.event[i].eta());
           etaF        -> Fill(pythia.event[i].eta());
           pTF         -> Fill(pythia.event[i].pT());
                   }*/
        
     }
    }
    // Fill charged multiplicity in histogram. End event loop.
  //  mult11       -> Fill( count );
   // mult_Cen     -> Fill( nChargedC );
 //   multF        -> Fill( nChargedF );   
    
    //Multiplicity in different eta ranges
   
  if (counter > 1)
             mult         -> Fill( counter4 );
   if (counter > 3)
            mult05       -> Fill (counter);
  if (counter > 3) 
            mult1        -> Fill (counter1);
  if (counter > 3) 
            mult1_5      -> Fill (counter2);
  if (counter > 3) 
            mult2        -> Fill (counter3);

if (counteta5  > 0)	mult_eta5   -> Fill(counteta5);
if (counteta1  > 0)	mult_eta1   -> Fill(counteta1);
if (counteta15 > 0)	mult_eta15  -> Fill(counteta15);
if (counteta2  > 0)	mult_eta2   -> Fill(counteta2);
if (counteta25 > 0)	mult_eta25  -> Fill(counteta25);
  }
  
 //=====================================================================//
 //                            Jet Work                                //
//=====================================================================//
   
   /*
   for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    
    
    int slow_count = 0;
    int cell_count = 0;

    int slow_countC =0;
    int cell_countC =0;
    
    int slow_countF = 0;
    int cell_countF = 0;
    
// Analyze Slowet jet properties. List first few.
    slowJet. analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Fill SlowJet inclusive jet distributions.
    
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
    
  //          jet study in   central region ========================= tracker region      
        
        if (abs(slowJet.y(i))< 2.5){ 
            eTjetsS   -> Fill( slowJet.pT(i) );
            etaJetsS  -> Fill( slowJet.y(i) );
            phiJetsS  -> Fill( slowJet.phi(i) );
            JetMult_S   -> Fill( slowJet.multiplicity(i));
            Mass_jetS-> Fill( slowJet.m(i));
            slow_count++; 
            nJetsS-> Fill( slow_count );
         }
 
  //          jet study in   forward region ========================= Hadron forward region

       if (abs(slowJet.y(i))> 3. && abs(slowJet.y(i)) < 5) {
        slow_countC++;
       eTjetsS_Cen   -> Fill( slowJet.pT(i) );
       etaJetsS_Cen  -> Fill( slowJet.y(i) );
       phiJetsS_Cen  -> Fill( slowJet.phi(i) );
       nJetsS_Cen    -> Fill( slow_countC );
       Mass_jetSC    -> Fill( slowJet.m(i));
       JetMult_SC    -> Fill( slowJet.multiplicity(i));
                }
 
 //          jet study in very forward region ========================= castor region
 
         if (abs(slowJet.y(i)) < 5.2 && abs(slowJet.y(i) ) > 6.6){     
                eTjetsSF   -> Fill( slowJet.pT(i) );
                etaJetsSF  -> Fill( slowJet.y(i) );
                phiJetsSF  -> Fill( slowJet.phi(i) );
                slow_countF++;  
               nJetsSFI    -> Fill( slow_countF );
               Mass_jetSF  -> Fill( slowJet.m(i));
              JetMult_SF   -> Fill( slowJet.multiplicity(i));
            
            }
        
    }
    
    nJetsS-> Fill( slow_count );
    
    // Fill SlowJet distance between jets.
    for (int i = 0; i < slowJet.sizeJet() - 1; ++i)
    for (int j = i +1; j < slowJet.sizeJet(); ++j) {
        if (abs(slowJet.y(i))< 2.5 && abs (cellJet.etaWeighted(j)<2.5)){  
            double dEta = slowJet.y(i) - slowJet.y(j);
            double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
            if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
            double dR = sqrt( pow2(dEta) + pow2(dPhi) );
             distJetsS-> Fill( dR );
            }
         
       if (abs(slowJet.y(i)) > 3. && abs(slowJet.y(i)) < 5. && 
           abs(slowJet.y(j))> 3. && abs(slowJet.y(j))< 5. ) {
            double dEta = slowJet.y(i) - slowJet.y(j);
            double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
            if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
            double dR = sqrt( pow2(dEta) + pow2(dPhi) );
            distJetsS_Cen -> Fill( dR );
           }
           
       if (abs(slowJet.y(i))>5.2&& abs(slowJet.y(i) )<6.6 &&
             abs(slowJet.y(j))>5.2&& abs(slowJet.y(j) )<6.6) {
                double dEta = slowJet.y(i) - slowJet.y(j);
                double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
            
                if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
                double dR = sqrt( pow2(dEta) + pow2(dPhi) );
                distJetsSF-> Fill( dR );
                }
    }

    
    
    // Fill SlowJet pT-difference between jets (to check ordering of list).
    for (int i = 1; i < slowJet.sizeJet(); ++i){
         if (abs(slowJet.y(i))< 2.5) {
            eTdiffS-> Fill( slowJet.pT(i-1) - slowJet.pT(i) );
         }
    
        if (abs(slowJet.y(i)) > 3. && abs(slowJet.y(i)) < 5. ){
            eTdiffS_Cen  -> Fill( slowJet.pT(i-1) - slowJet.pT(i) );     
        }
        
        if (abs(slowJet.y(i))>5.2&& abs(slowJet.y(i) )<6.6 &&
            abs(cellJet.etaWeighted(i))>5.2&& abs(cellJet.etaWeighted(i) )<6.6){
            
            eTdiffSF-> Fill( slowJet.pT(i-1) - slowJet.pT(i) );
            }
        }    
   
   
// Analyze CellJet jet properties. List first few.
   cellJet. analyze( pythia.event, pTjetMin, radius );
    if (iEvent < nListJets) cellJet.list();
    
    // Fill CellJet inclusive jet distributions.
   
    nJetsC-> Fill( cellJet.size() );
    for (int i = 0; i < cellJet.size(); ++i) {
         if (abs(cellJet.etaWeighted(i))< 2.5){ 
            eTjetsC   -> Fill( cellJet.eT(i) );  
            etaJetsC  -> Fill( cellJet.etaWeighted(i) );
            phiJetsC  -> Fill( cellJet.phiWeighted(i) );
            JetMult_C -> Fill( cellJet.multiplicity(i));
            Mass_jetC -> Fill( slowJet.m(i));
                     }
        if ( abs(cellJet.etaWeighted(i))  > 3. && abs(cellJet.etaWeighted(i)) < 5.) {
            eTjetsC_Cen   -> Fill( cellJet.eT(i) );
            etaJetsC_Cen  -> Fill( cellJet.etaWeighted(i) );
            phiJetsC_Cen  -> Fill( cellJet.phiWeighted(i) );
            cell_countC++;
            nJetsC_Cen   -> Fill( cell_countC );
            Mass_jetCC   -> Fill( cellJet.m(i));
            JetMult_CC   -> Fill( cellJet.multiplicity(i));
            }             
      
      if (abs(cellJet.etaWeighted(i))>5.2&& abs(cellJet.etaWeighted(i) )<6.6) {
         
            eTjetsCF     -> Fill( cellJet.eT(i) );
            etaJetsCF    -> Fill( cellJet.etaWeighted(i) );
            phiJetsCF    -> Fill( cellJet.phiWeighted(i) );
            cell_countF++;
            Mass_jetCF   -> Fill( cellJet.m(i));
            nJetsCI      -> Fill( cell_countF );
            JetMult_CF   -> Fill( cellJet.multiplicity(i));
            
            }         
               
              
        }

    // Fill CellJet distance between jets.
    for (int i = 0; i < cellJet.size() - 1; ++i)
    for (int j = i +1; j < cellJet.size(); ++j) {
        if (abs(cellJet.etaWeighted(i))< 2.5 && abs(cellJet.etaWeighted(j))< 2.5){ 
      double dEta = cellJet.etaWeighted(i)- cellJet.etaWeighted(j);
      double dPhi = abs( cellJet.phiWeighted(i)- cellJet.phiWeighted(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsC-> Fill( dR );
        }
    
    if (abs(cellJet.etaWeighted(i))> 3. && abs(cellJet.etaWeighted(i)) < 5. && 
        abs(cellJet.etaWeighted(j))> 3. && abs(cellJet.etaWeighted(j)) < 5.) {
      double dEta = cellJet.etaWeighted(i) - cellJet.etaWeighted(j);
      double dPhi = abs( cellJet.phiWeighted(i) - cellJet.phiWeighted(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsC_Cen -> Fill( dR );
        }
   
       if (abs(cellJet.etaWeighted(i))>5.2&& abs(cellJet.etaWeighted(i) )<6.6 &&
            abs(cellJet.etaWeighted(j))>5.2&& abs(cellJet.etaWeighted(j) )<6.6) {
            double dEta = cellJet.etaWeighted(i) - cellJet.etaWeighted(j);
            double dPhi = abs( cellJet.phiWeighted(i)  - cellJet.phiWeighted(j) );
            if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
            double dR = sqrt( pow2(dEta) + pow2(dPhi) );
            distJetsCF-> Fill( dR );
        }  
       
    }
    // Fill CellJet ET-difference between jets (to check ordering of list).
    for (int i = 1; i < cellJet.size(); ++i){
     if (abs(cellJet.etaWeighted(i))< 2.5 && abs(slowJet.y(i))< 2.5){

      eTdiffC-> Fill( cellJet.eT(i-1) - cellJet.eT(i) );
      
     }
    
    if ( abs(cellJet.etaWeighted(i)) > 3. && abs(cellJet.etaWeighted(i)) < 5.){
      eTdiffC_Cen      -> Fill( cellJet.eT(i-1) - cellJet.eT(i) );
          }
 
   if(abs(cellJet.etaWeighted(i))>5.2 && abs(cellJet.etaWeighted(i) )<6.6)
       eTdiffCF-> Fill( cellJet.eT(i-1) - cellJet.eT(i) );
   
   // Compare number of jets for the two finders.
    nJetsD-> Fill( cellJet.size() - slowJet.sizeJet() );
    
    
}*/
//===================================================================================================
 
  // Statistics on event generation.
  pythia.stat();

  /* Show histogram. Possibility to close it.
  mult->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  */
  // Save histogram on file and close file.

//=============== Tracker region ====================================================================

TFile* f = new TFile("Asmaa_20TeV.root", "RECREATE");

mult         -> Write();
mult11        -> Write();
mult05       -> Write();
mult1        -> Write();
mult1_5      -> Write();
mult2        -> Write();
multF        -> Write();
etaF         -> Write();
pTF          -> Write();
eta          -> Write();
eta1         -> Write();

mult_Cen     -> Write();
eta_Cen      -> Write();
pT_Cen       -> Write();

pT_05         -> Write();
pT_1          -> Write();
pT_15         -> Write();
pT_2          -> Write();
pT_25         -> Write();
eta_1         -> Write();
eta_2         -> Write();
eta_3         -> Write();
eta_4         -> Write();
eta_5         -> Write();


mult_eta5     -> Write();
mult_eta1     -> Write();
mult_eta15    -> Write();
mult_eta2     -> Write();
mult_eta25    -> Write();
 

/*
nJetsS       -> Write();
nJetsC       -> Write();
nJetsD       -> Write();
eTjetsS      -> Write();
eTjetsC      -> Write();
etaJetsS     -> Write();
etaJetsC     -> Write();
phiJetsS     -> Write();
phiJetsC     -> Write();
distJetsS    -> Write();
distJetsC    -> Write();
eTdiffS      -> Write();
eTdiffC      -> Write();

//=============== Central region ====================================================================

nJetsS_Cen       -> Write();
nJetsC_Cen       -> Write();
nJetsD_Cen       -> Write();
eTjetsS_Cen      -> Write();
eTjetsC_Cen      -> Write();
etaJetsS_Cen     -> Write();
etaJetsC_Cen     -> Write();
phiJetsS_Cen     -> Write();
phiJetsC_Cen     -> Write();
distJetsS_Cen    -> Write();
distJetsC_Cen    -> Write();
eTdiffS_Cen      -> Write();
eTdiffC_Cen      -> Write();

//============== Forward region =====================================================================



nJetsSF       -> Write();
nJetsCF       -> Write();
nJetsDF       -> Write();
eTjetsSF      -> Write();
eTjetsCF      -> Write();
etaJetsSF     -> Write();
etaJetsCF     -> Write();
phiJetsSF     -> Write();
phiJetsCF     -> Write();
distJetsSF    -> Write();
distJetsCF    -> Write();
eTdiffSF      -> Write();
eTdiffCF      -> Write();

JetMult_S     -> Write(); 
JetMult_C     -> Write();

JetMult_SC     -> Write(); 
JetMult_CC     -> Write();

JetMult_SF     -> Write(); 
JetMult_CF     -> Write();

Mass_jetS      -> Write(); 
Mass_jetC      -> Write(); 
Mass_jetSC     -> Write(); 
Mass_jetCC     -> Write(); 
Mass_jetSF     -> Write(); 
Mass_jetCF     -> Write(); 
*/
   f->Close();
//  delete outFile;

  // Done.
  return 0;
}





