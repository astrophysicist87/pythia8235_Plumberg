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
	pythia.readString("Fragmentation:setVertices = on");
	//pythia.readString("PartonVertex:setVertex = on");

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
	const int max_events_per_file = 10000;
	int current_file_index = 0;
	string file_index_string = "";
	if (total_number_of_events > max_events_per_file)
	{
		file_index_string = "_" + to_string(current_file_index);
	}

	//string systemSpecs = string(argv[1]) + string(argv[2]) + "_" + string(argv[3]) + "GeV_"
	//						+ "C" + string(argv[4]) + "_" + string(argv[5]) + "_Nev" + string(argv[4]);
	string systemSpecs = string(argv[1]) + string(argv[2]) + "_" + string(argv[3]) + "GeV_Nev" + string(argv[4]);

	//string path = "./results/";
	//string path = "/scratch/blixen/plumberg/results/";
	string path = string(argv[5]) + "/";
	ostringstream filename_stream, mult_fn_stream;
	filename_stream //<< path
					<< systemSpecs
					<< file_index_string
					<< ".dat";
	ofstream outmain( ( path + filename_stream.str()).c_str());
	mult_fn_stream //<< path
					<< systemSpecs << "_mult"
					<< ".dat";
	ofstream outMultiplicities( (path + mult_fn_stream.str()).c_str());

	//cerr << "Pythia(): printing to "
	//		<< path + systemSpecs + "_S_x_p_filenames.dat"
	//		<< " and " << path + systemSpecs + "_total_N_filename.dat" << endl;

	ofstream outfilenames(path + systemSpecs + "_S_x_p_filenames.dat");
	outfilenames << filename_stream.str() << endl;

	ofstream outmult_filenames(path + systemSpecs + "_total_N_filename.dat");
	outmult_filenames << mult_fn_stream.str() << endl;
	outmult_filenames.close();

	int count = 0;

	// Loop over events.
	for ( int iEvent = 0; iEvent < total_number_of_events; ++iEvent )
	{
		if ( !pythia.next() )
			continue;

		int event_multiplicity = 0;
		int pion_multiplicity = 0;

		for (int i = 0; i < pythia.event.size(); ++i)
		{
			Particle & p = pythia.event[i];
			if ( p.isFinal() and p.isHadron() )
			{
				//count all final hadrons in multiplicity
				event_multiplicity++;

				// p.status() < 91 means p is not decay particle
		 		if ( p.id() == 211 )	// i.e., is pi^+
				{
					// i.e., only do it once
					if (count < 1)
					{
						ofstream out(path + "HBT_particle.dat");
						out << "name = " << pythia.event[i].name()
							<< "\t\t# Particle name" << endl
							<< "monval = " << pythia.event[i].id()
							<< "\t\t# Monte-Carlo number" << endl
							<< "mass = " << pythia.event[i].m()
							<< "\t\t# mass" << endl
							<< "charge = " << pythia.event[i].charge()
							<< "\t\t# charge" << endl
							<< "spinType = " << pythia.event[i].spinType()
							<< "\t\t# spin type" << endl
							<< "chargeType = " << pythia.event[i].chargeType()
							<< "\t\t# charge type" << endl;
						out.close();
						++count;
					}

		 			outmain
		 				<< iEvent << "   " << i
		 				<< "   " << p.e()
		 				<< "   " << p.px()
		 				<< "   " << p.py()
		 				<< "   " << p.pz()
		 				<< "   " << p.tProd()
		 				<< "   " << p.xProd()
		 				<< "   " << p.yProd()
		 				<< "   " << p.zProd()
		 				<< endl;

					pion_multiplicity++;
				}
			}
		}

		// If too many events for single file, set-up new file here
		if ( (iEvent + 1) % max_events_per_file == 0
			and total_number_of_events > max_events_per_file
			and iEvent + 1 < total_number_of_events )
		{
			outmain.close();
			file_index_string = "_" + to_string(++current_file_index);
			filename_stream.str("");
			filename_stream //<< path
				<< systemSpecs
				<< file_index_string
				<< ".dat";
			outmain.open(( path + filename_stream.str()).c_str());
			outfilenames << filename_stream.str() << endl;

		}

		outMultiplicities << iEvent << "   "
					<< event_multiplicity << "   "
					<< pion_multiplicity << endl;
	}

	outmain.close();
	outMultiplicities.close();
	outfilenames.close(); 

	// And we're done!
	return 0;
}
