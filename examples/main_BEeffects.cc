// main_BEeffects.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;
using namespace std;

vector<int> get_centrality_limits(
			const double centrality_class_lower_limit,
			const double centrality_class_upper_limit,
			const int n_events_to_use, Pythia & pythia );

void print_particle_record(
		int iEvent, vector<Particle> & particles_to_output,
		ofstream & record_stream );

int main(int argc, char *argv[])
{
	// Check number of command-line arguments.
	if (argc != 8)
	{
		cerr << "Incorrect number of arguments!" << endl;
		cerr << "Usage: ./main_BEeffects [Projectile nucleus] [Target nucleus] [Beam energy in GeV]"
				<< " [Number of events] [Results directory]"
				<< " [Lower centrality %] [Upper centrality %]" << endl;
		exit(8);
	}

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

	// turn on and set Bose-Einstein effects
	pythia.readString("HadronLevel:BoseEinstein = on");
	//pythia.readString("BoseEinstein:QRef = " + string(argv[6]));
	pythia.readString("BoseEinstein:widthSep = 1.0");

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

	ofstream outfilenames(path + systemSpecs + "_S_x_p_filenames.dat");
	outfilenames << filename_stream.str() << endl;

	ofstream outmult_filenames(path + systemSpecs + "_total_N_filename.dat");
	outmult_filenames << mult_fn_stream.str() << endl;
	outmult_filenames.close();

	int count = 0;

	// Estimate centrality class limits
	const int n_events_to_use = 10000;
	const double centrality_class_lower_limit = atof( argv[6] );
	const double centrality_class_upper_limit = atof( argv[7] );

	cout << "Read in these centrality limits: "
			<< centrality_class_lower_limit << " to "
			<< centrality_class_upper_limit << endl;

	vector<int> centrality_limits
		= get_centrality_limits(
			centrality_class_lower_limit,
			centrality_class_upper_limit,
			n_events_to_use, pythia );

	const int multiplicity_lower_limit = centrality_limits[0];
	const int multiplicity_upper_limit = centrality_limits[1];

	cout << "Accepted multiplicity range: "
			<< multiplicity_lower_limit << " to "
			<< multiplicity_upper_limit << endl;

	// Loop over events.
	int iEvent = 0;
	do
	{
		if ( !pythia.next() )
			continue;

		int event_multiplicity = 0;
		int pion_multiplicity = 0;

		vector<Particle> particles_to_output;

		for (int i = 0; i < pythia.event.size(); ++i)
		{
			Particle & p = pythia.event[i];
			if ( p.isFinal() and p.isHadron() )
			{
				//count all final hadrons in multiplicity
				event_multiplicity++;

		 		if ( p.id() == 211 )	// i.e., is pi^+
				{
					// i.e., only do it once
					if ( count < 1 )
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

					if ( thermal_only and p.status() != 99 )	// only works for mom.-space modifications
						continue;

					particles_to_output.push_back( p );

					pion_multiplicity++;
				}
			}
		}

		bool event_in_chosen_centrality_class
			= ( event_multiplicity >= multiplicity_lower_limit)
				and ( event_multiplicity <= multiplicity_upper_limit);

		// only save this event if N^{ch} is correct range
		if ( not event_in_chosen_centrality_class )
			continue;

		print_particle_record( iEvent, particles_to_output, outmain );

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

		++iEvent;

	} while ( iEvent < total_number_of_events );

	outmain.close();
	outMultiplicities.close();
	outfilenames.close(); 

	// And we're done!
	return 0;
}









bool decreasing (int i,int j) { return (i>j); }


vector<int> get_centrality_limits(
			const double centrality_class_lower_limit,
			const double centrality_class_upper_limit,
			const int n_events_to_use, Pythia & pythia )
{
	int iEvent = 0;
	vector<int> event_multiplicities;

	vector<int> results(2);

	// if we're just doing all the events, no matter what
	if ( centrality_class_lower_limit < 1.e-6
			and 100.0 - centrality_class_upper_limit < 1.e-6 )
	{
		cout << "main_BEeffects(): setting default centralities" << endl;
		results[0] = 0;
		results[1] = 1e+9;
	}
	else
	{
		do
		{
			if ( !pythia.next() )
				continue;

			int event_multiplicity = 0;

			for (int i = 0; i < pythia.event.size(); ++i)
			{
				Particle & p = pythia.event[i];
				if ( p.isFinal() and p.isHadron() )
				{
					//count all final hadrons in multiplicity
					//to estimate centrality
					event_multiplicity++;
				}
			}

			event_multiplicities.push_back( event_multiplicity );

			++iEvent;

		} while ( iEvent < n_events_to_use );

		sort( event_multiplicities.begin(), event_multiplicities.end(), decreasing );

		int lower_index = max( 0, (int)floor( 0.01*centrality_class_lower_limit*n_events_to_use+0.5 ) );
		int upper_index = min( n_events_to_use - 1, (int)floor( 0.01*centrality_class_upper_limit*n_events_to_use+0.5 ) );
		results[0] = ( 100.0 - centrality_class_upper_limit < 1.e-6 ) ?
						0 : 
						event_multiplicities[upper_index];	//smaller multiplicity limit first
		results[1] = ( centrality_class_lower_limit < 1.e-6 ) ?
						1e+9 : 
						event_multiplicities[lower_index];	//larger multiplicity limit second
	}

	return ( results );
}




void print_particle_record(
		int iEvent, vector<Particle> & particles_to_output,
		ofstream & record_stream )
{
	for (int i = 0; i < (int)particles_to_output.size(); ++i)
	{
		Particle & p = particles_to_output[i];

		record_stream
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
	}

	return;
}

// End of file
