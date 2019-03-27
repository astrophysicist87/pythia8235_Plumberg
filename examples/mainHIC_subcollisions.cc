// main113.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This test program will generate Pb-Pb collisions at
// sqrt(S_NN)=2.76TeV using the Angantyr model for Heavy Ion
// collisions. The analysis will divide the event in centrality
// classes using the same observable as was used for p-Pb in the ATLAS
// analysis in arXiv:1508.00848 [hep-ex] (see main112.cc). The
// centrality classes are same as in the ALICE analysis in
// arXiv:1012.1657 [nucl-ex] although the actual observable used is
// not the same. Histograms of multiplicity distributions are measured
// for each centrality percentile.

// Note that heavy ion collisions are computationally quite CPU
// intensive and generating a single event will take around a second
// on a reasonable desktop. To get reasonable statistics, this program
// will take a couple of hours to run.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;
using namespace std;

int main(int argc, char *argv[])
{

	Pythia pythia;

	// Can easily add more later
	std::unordered_map<string, string> particle_IDs
		= {
			{ "p"  , "2212"       },
			{ "d"  , "1000010020" },
			{ "t"  , "1000010030" },
			{ "He3", "1000020030" },
			{ "He4", "1000020040" },
			{ "Li" , "1000030060" },
			{ "C"  , "1000060120" },
			{ "O"  , "1000080160" },
			{ "Cu" , "1000290630" },
			{ "Xe" , "1000541290" },
			{ "Au" , "1000791970" },
			{ "Pb" , "1000822080" },
			{ "U"  , "1000922380" }
		  };

	// Turn on tracking of space-time information
	//pythia.readString("Fragmentation:setVertices = on");
	pythia.readString("PartonVertex:setVertex = on");
	
	if (argc != 6)
	{
		cerr << "Incorrect number of arguments!" << endl;
		cerr << "Usage: ./mainHIC [Projectile nucleus] [Target nucleus] [Beam energy in GeV]"
				<< " [Number of events] [Results directory]" << endl;
		exit(8);
	}

	// Setup the beams.
	pythia.readString("Beams:idA = " + particle_IDs[string(argv[1])]);
	pythia.readString("Beams:idB = " + particle_IDs[string(argv[2])]);
	pythia.readString("Beams:eCM = " + string(argv[3]));
	pythia.readString("SoftQCD:all = on");
	//pythia.readString("HardQCD:all = on");
	pythia.readString("Beams:frameType = 1");

	cout << "Set " << "Beams:idA = " + particle_IDs[string(argv[1])]
		<< " and " << "Beams:idB = " + particle_IDs[string(argv[2])]
		<< " and " << "Beams:eCM = " + string(argv[3]) << endl;

	// Initialize the Angantyr model to fit the total and semi-inclusive
	// cross sections in Pythia within some tolerance.
	pythia.readString("HeavyIon:SigFitErr = "
						"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
	// These parameters are typically suitable for sqrt(S_NN)=5 TeV
	pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
	// A simple genetic algorithm is run for 20 generations to fit the
	// parameters.
	pythia.readString("HeavyIon:SigFitNGen = 20");

	// Initialise Pythia.
	pythia.init();

	//const int total_number_of_events = 100000;
	const int total_number_of_events = atoi(argv[4]);

	string systemSpecs = string(argv[1]) + string(argv[2]) + "_" + string(argv[3]) + "GeV_Nev" + string(argv[4]) + ".dat";

	string path = string(argv[5]);

	ofstream out( ( path + systemSpecs ).c_str());

	// Loop over events.
	for ( int iEvent = 0; iEvent < total_number_of_events; ++iEvent )
	{
		// Generate the next event
		if ( !pythia.next() )
			continue;

		for (int i = 0; i < pythia.event.size(); ++i)
		{
			Particle & p = pythia.event[i];
			int p_status = abs( p.status() );
			int p_id = abs( p.id() );
			if ( p_status == 21 or p_status == 31 or p_status == 23 or p_status == 33
					or p_id == 21 or p_id == 31 or p_id == 23 or p_id == 33 )
			{
		 			out << iEvent << "   " << i
		 				<< "   " << p.id()
		 				<< "   " << p.e()
		 				<< "   " << p.px()
		 				<< "   " << p.py()
		 				<< "   " << p.pz()
		 				<< "   " << p.tProd()
		 				<< "   " << p.xProd()
		 				<< "   " << p.yProd()
		 				<< "   " << p.zProd()
		 				<< endl;

			}

		}

	}

	out.close();

	// And we're done!
	return 0;
}
