// mainCP.cc is not a part of the PYTHIA event generator.
// Copyright (C) 2018 Christopher Plumberg.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is program generates pion distributions for a given
// number of events.  Output can be passed on to HBT event
// generator

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Pythia8/Pythia.h"

using namespace Pythia8;
using namespace std;

int main()
{
	// Generator. Process selection. LHC initialization.
	Pythia pythia;
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("HardQCD:all = on");
	pythia.readString("PhaseSpace:pTHatMin = 20.");

	// Set tracking of space-time information
	bool space_time_tracking = true;
	string space_time_tracking_string = "";
	if (space_time_tracking)
	{
		pythia.readString("Fragmentation:setVertices = on");
		space_time_tracking_string = "_wSTtracking";
	}

	// Set modeling of Bose-Einstein effects
	bool Bose_Einstein_included = true;
	string Bose_Einstein_included_string = "";
	if (Bose_Einstein_included)
	{
		pythia.readString("HadronLevel:BoseEinstein = on");
		Bose_Einstein_included_string = "_wBEeffects";
	}

	// Initialize event generator
	pythia.init();

	const int total_number_of_events = 100000;
	const int max_events_per_file = 10000;
	int current_file_index = 0;
	string file_index_string = "";
	if (total_number_of_events > max_events_per_file)
	{
		file_index_string = "_" + to_string(current_file_index);
	}

	string path = "./results/";
	ostringstream filename_stream;
	filename_stream << path << "pp_13TeV"
					<< space_time_tracking_string
					<< Bose_Einstein_included_string
					<< file_index_string
					<< ".dat";
	ofstream outmain(filename_stream.str().c_str());

	int count = 0;

	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < total_number_of_events; ++iEvent)
	{
		if (!pythia.next())
			continue;

		// Loop over all particles in event
		for (int i = 0; i < pythia.event.size(); ++i)
			if (pythia.event[i].isFinal() && pythia.event[i].isHadron())
			{
				if ( pythia.event[i].id() == 211 )	// i.e., is pi^+
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
						<< "   " << pythia.event[i].e()
						<< "   " << pythia.event[i].px()
						<< "   " << pythia.event[i].py()
						<< "   " << pythia.event[i].pz()
						<< "   " << pythia.event[i].tProd()
						<< "   " << pythia.event[i].xProd()
						<< "   " << pythia.event[i].yProd()
						<< "   " << pythia.event[i].zProd()
						<< endl;
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
			filename_stream << path << "pA_13TeV"
				<< space_time_tracking_string
				<< Bose_Einstein_included_string
				<< file_index_string
				<< ".dat";
			outmain.open(filename_stream.str().c_str());
		}

	}

	outmain.close();

	return 0;
}
