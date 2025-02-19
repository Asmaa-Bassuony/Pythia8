# Pythia8
The provided code includes sections for using the PYTHIA event generator and ROOT histogramming libraries to study charged particle multiplicities and jet properties in high-energy physics simulations. However, due to the incomplete truncation at the end, a few details might be missing.

# Code Overview
Event Generation with PYTHIA:

Initializes the PYTHIA generator for simulations of proton-proton collisions at a specified center-of-mass energies 
The event loop generates particles and collects kinematic information for analysis.
Histogramming with ROOT:

Various histograms (TH1F) are defined to store information about particle multiplicities, transverse momentum (pT), pseudorapidity (eta), jet properties, and other event-level features.
Histograms are grouped for analyzing central (e.g., mult_Cen, eta_Cen) and forward (e.g., multF, etaF) regions of the detector.
Jet Reconstruction:

Implements two jet reconstruction algorithms:
SlowJet (anti-kT algorithm): Clusters particles to form jets with a specified radius and transverse momentum threshold.
CellJet: Grid-based clustering for a fast estimation of jet properties.
Event Analysis:

Event data such as charged multiplicity, jet count, invariant mass, and kinematic properties are filled into histograms.
Separate histograms track differences between SlowJet and CellJet results, providing insights into jet-finder performance.
