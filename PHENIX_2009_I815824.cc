// -*- C++ -*-
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet
{
/// **********************************************************************************
///   plots over events in the 10-30% and 50-80% centrality range, but that data is seperate
///   Triggers particles are defined as Pi0s with 13 GeV < pT < 20 GeV or
///   gamma particles with 8 GeV < pT < 20 GeV. Both are |eta| < 1.0
///   Associated partices are all charged hadrons with pT > 1.2 GeV and less than trigger pT


	class PHENIX_2009_I815824 : public HeavyIonAnalysis
	{
	// Variables pertaining to just this class
	// I moved the private section above the public to introduce
	//   the variables you will see in init and later
	private:
		const int N_CENT_TYPES = 3; // number of centrality bins that you need to worry about, this case three, 0-20, 20-40, 40-60
		const int CENT_TYPE_EDGES[N_CENT_TYPES][3] = { {0,20},{20,40},{40,60} }; // 0-20, 20-40, 40-60
		
		
		//Histo1DPtr _h1dPhi[N_CENT_TYPES]; // delta phi, split by centrality
		scatter2DPtr _yield[24];
		scatter2DPtr _IAA[4];
		scatter2DPtr _IAAz[4];

		// Number of trigger particles for each centrality over
		//  all events. Needed in finalize() for normalization
		unsigned long long int nTrigger[N_CENT_TYPES];

	public:
		PHENIX_2009_I815824():HeavyIonAnalysis("PHENIX_2009_I815824") {}

		/// Book histograms and initialise projections before the run
		///   Every histogram you will use must be booked, and every projection
		///   you will use must be declared
		void init()
		{
			//**** Select centrality method ****
			// The 50 here is the number of events used JUST to determine future events'
			//   centrality. These events will show a centrality of -1.0, which is not
			//   physical. You can adjust this number for testing, but never have it lower
			//   than the number of events you are testing against.
			//   The test file ampt.hepmc has 100 events
			addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod");

			//**** Trigger particle set ****
			// Many of you will have "trigger" particles. This is code to trigger off of a
			//   specific type or particle in a certain range. There will ALWAYS be an "abseta"
			//   cut on all of your particles. It may differ from your trigger and associated
			//   particles. This example uses Pi0 and gamma of differet energies as triggers.
			const int pidGAMMA = 22; // The PID for Gamma
			// Find other pid codes at: http://home.fnal.gov/~mrenna/lutp0613man2/node44.html
			// I apply the following cuts for trigger particles:
			//  |eta| < 1.0  for all particles
			//  or  A gamma with pT between 5 and 7, 7 and 9, 9 and 12, 12 and 15 GeV
			Cuts cutTrigger = (Cuts::pid == pidGAMMA && Cuts::pt > 5.0 * GeV && Cuts::pt < 7.0 * GeV)) || (Cuts::pid == pidGAMMA && Cuts::pt > 7.0 * GeV && Cuts::pt < 9.0 * GeV)) || (Cuts::pid == pidGAMMA && Cuts::pt > 9.0 * GeV && Cuts::pt < 12.0 * GeV)) || (Cuts::pid == pidGAMMA && Cuts::pt > 12.0 * GeV && Cuts::pt < 15.0 * GeV))
			FinalState fs(cutTrigger);
			declare(fs, "partTrigger");

			//**** Associated particle set ****
			// Charged Particles with cuts
			// Max pT for an associated particle will be limited by the pT of the trigger
			//  but that is based on the specific trigger. So that check will be done later
			//  inside a loop in analysis
			Cuts cutAssoc = Cuts::abseta < 1.0 && Cuts::pt > 1.2 * GeV && Cuts::pt < 20.0 * GeV;
			ChargedFinalState cfs(cutAssoc);
			declare(cfs, "partAssoc");

			//booking scatter plots
			//yield scatter plots
			_yield[0] = bookScatter2D(1, 1, 1);
			_yield[1] = bookScatter2D(1, 1, 2);
			_yield[2] = bookScatter2D(2, 1, 1);
			_yield[3] = bookScatter2D(2, 1, 2);
			_yield[4] = bookScatter2D(3, 1, 1);
			_yield[5] = bookScatter2D(3, 1, 2);
			_yield[6] = bookScatter2D(4, 1, 1);
			_yield[7] = bookScatter2D(4, 1, 2);
			_yield[8] = bookScatter2D(5, 1, 1);
			_yield[9] = bookScatter2D(5, 1, 2);
			_yield[10] = bookScatter2D(6, 1, 1);
			_yield[11] = bookScatter2D(6, 1, 2);
			_yield[12] = bookScatter2D(7, 1, 1);
			_yield[13] = bookScatter2D(7, 1, 2);
			_yield[14] = bookScatter2D(8, 1, 1);
			_yield[15] = bookScatter2D(8, 1, 2);
			_yield[16] = bookScatter2D(9, 1, 1);
			_yield[17] = bookScatter2D(9, 1, 2);
			_yield[18] = bookScatter2D(10, 1, 1);
			_yield[19] = bookScatter2D(10, 1, 2);
			_yield[20] = bookScatter2D(11, 1, 1);
			_yield[21] = bookScatter2D(11, 1, 2);
			_yield[22] = bookScatter2D(12, 1, 1);
			_yield[23] = bookScatter2D(12, 1, 2);

			//IAA scatter plots
			_IAA[0] = bookScatter2D(13, 1, 1);
			_IAA[1] = bookScatter2D(14, 1, 1);
			_IAA[2] = bookScatter2D(15, 1, 1);
			_IAA[3] = bookScatter2D(16, 1, 1);

			//IAAz scatter plots
			_IAAz[0] = bookScatter2D(17, 1, 1);
			_IAAz[1] = bookScatter2D(18, 1, 1);
			_IAAz[2] = bookScatter2D(19, 1, 1);
			_IAAz[3] = bookScatter2D(20, 1, 1);

			//**** Initialize counters ****
			// Set number of events in counters to zero
			// These are used later for normalizing histograms
			//  in finalize
			for (int i = 0; i < N_CENT_TYPES; ++i) nTrigger[i] = 0;
		}

		/// Per event calculations
		///  Do your per event calculations and
		///  fill your histograms here
		void analyze(const Event& event)
		{
			// Get the centrality for each event
			const double c = centrality(event, "IPMethod");

			// The first 50 events (number from addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod") )
			//  will give a centrality outside of 0-100, specifically -1.0
			if ((c < 0.) || (c > 100.)) vetoEvent;

			/// commented out the code below since I think the above code is cleaner, but left it incase I or someone else decided to use it instead in the future.
			/// ***
			// Alternatively, you may want to set the acceptance to a narrower range based on what your paper plots
			// For my fake paper there are two ranges. You can use these as indices for an array of histograms instead.
			//int centralityIndex = -1; // -1 - not in plot range
									 // 0 - between 10-30
									 // 1 - between 50-80
			// Find which centrality range of interest this falls into
			//for (int i = 0; i < N_CENT_TYPES; ++i) {
			//	if (c > CENT_TYPE_EDGES[i][0] && c <= CENT_TYPE_EDGES[i][1]) {
			//		centralityIndex = i;
			//		break; // since theres no centrality overlap, break the loop when you find the correct index
			//	}
			//}
			// If not a centrality of interest, stop processing the event
			//if (centralityIndex == -1) vetoEvent;



			// In the background RIVET is handling the the declared projections
			//  but you need to access them inside analyze to use tCFShem
			const FinalState & fs = apply<FinalState>(event,"partTrigger");
			const ChargedFinalState & cfs = apply<ChargedFinalState>(event,"partAssoc");

			// Now that you have the projection objects, you need to get the particles from them.
			Particles tracksTrigger = fs.particlesByPt();
			Particles tracksAssoc = cfs.particlesByPt();

			// Increase appropriate counter
			nTrigger[centralityIndex] += tracksTrigger.size();

			// foreach goes through each Particle in tracksTrigger, and references it in partTrigger.
			//  Its a simplification of a for loop where you would have had to find the number of entries
			//  and loop over them.
			// I am assuming I am calculating delta-phi between trigger particles and charged hadrons
			double deltaPhi;
			foreach (const Particle& partTrigger, tracksTrigger) {
				// Loop over all associated particles
				foreach (const Particle& partAssoc, tracksAssoc) {
					// Only include associated particles with pT less than the trigger
					if (partAssoc.pt() < partTrigger.pt()) {
						deltaPhi = partAssoc.phi() - partTrigger.phi();
					
						//   so make sure deltaPhi falls in this range
						// M_PI is part of the <math.h> header. You need to use the define: #define _USE_MATH_DEFINES
						//   See header includes at top of file for example
						while (deltaPhi < 0) deltaPhi += 2 * M_PI;

						// Fill histogram with delta phi information
						//  In this case use a 1 as the second parameter, the weight.
						//  With event weighting factored in, that number would be variable
						//  and nTrigger would need to factor in weighting.
						_h1dPhi[centralityIndex]->fill(deltaPhi,1);
					}
				}
			}
		}

		// After all events are processed this is called.
		//  Here I need to normalize my histogram by dividing by
		//  the number of triggers.
		void finalize()
		{
			for (int i = 0; i <= 24; ++i){
				// calculate the inverse of the number of triggers
				double scaleFactor = 1./(double)nTrigger[i];

				// scale by that factor
				_yield[i]->scaleW(scaleFactor);
				// Histo1Dptr are just pointers to YODA objects, and their documentation
				// is in the YODA docs.
				// e.g.: https://yoda.hepforge.org/doxy/classYODA_1_1Histo1D.html
			}
			for (int i = 0; i <= 4; ++i) {
				double scaleFactor = 1. / (double)nTrigger[i];
				_IAA[i]->scaleW(scaleFactor);
			}
			for (int i = 0; i <= 4; ++i) {
				double scaleFactor = 1. / (double)nTrigger[i];
				_IAAz[i]->scaleW(scaleFactor);
			}
		}
	};


	DECLARE_RIVET_PLUGIN(PHENIX_2009_I815824);
}
