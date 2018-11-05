# dune_foils_simulation_solidanglemethod
Low energy event light simulations in DUNE.

Code capable of simulating the induced scintillation of low energy events (specifically supernova neutrinos, argon-39 beta decay and radon-222 alpha decay) occuring within the active volume of DUNE including TPB-coated wavelength-shifting reflector foils on the cathode.

## Using the code
The code generates .root files which contain both an "event_tree" and a "data_tree", where the former contains information on each event entry (e.g. a supernova event) and the latter essentially contains an entry for each and every DETECTED photon from ALL events simulated (=> data_tree >> event_tree).

* Before you start, you'll need to set up ROOT. Do this by entering "source setupRoot.sh" into the terminal
* Then, compile the code by typing "make -B", to recompile everything.
  * As many of the configurable parameters are in the 'analyze_light_histo.h' header file, if you 	change a parameter in this file, it is best to recompile everything in the project.
* The Makefile generates an executable that can be run with "./analyze_light_histo" (or whatever you change the name to). If you happen to be missing the data file, a segmentation violation will occur. Before the crash readout, you will find that the requested file could not be found. Change your path, and it should then run fine.

The code creates a root file with an event tree and data trees for all,vuv and vis photons. The event_tree has data on an event-by-event basis, and data_tree has the information based on DETECTED photons from ALL events.

## IMPORTANT NOTES:
1.  Use the "analyze_light_histo.h" file to change:
	  - The type of event
	  - The number of events
	  - The position of the events (make sure only one Boolean is true)
	  - The foil configuration (full coverage, cathode-only coverage or no coverage)
	  - The general properties of an event type, e.g. electron scintillation yield (most likely will not need to be changed)

2. All of the main functionailty used to generate the .root files is found in "analyze_light_histo.cc". This is well commented and I would advise going through it slowly to understand how the trees are created/filled.

3. solid_angle_functions.h and solid_angle_functions.cc contain the functions with determine the number of hits on each optical channel using the solid angle instead + corrections for Rayleigh scattering instead of an optical library. Preliminary. 

4. timingparam.h and timingparam.cc contain the functions which determine the timing of the photons with the updated timing parameterisations + discretisation for efficiency. Visible light timings preliminary.

5. utility_functions.cc contains the energy spectrum functions.
