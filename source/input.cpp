/*
 *	DESCRIPTION:	Read input data for a simulation from the file 'settings_file.txt'
 *
 */


#include "input.hpp"
#include <fstream>
#include <sstream>
#include <cassert>



void settingsfile :: read_settings_file ()
{
	std::ifstream infile("settings_file.txt");


	std::string BCRchoice;
	std::string BCLchoice;
	std::string RS_purechoice;
	std::string RS_mixedchoice;
	std::string ICchoice;
	std::string GFMchoice;
	std::string simchoice;
	std::string eos1choice;
	std::string eos2choice;
	std::string FSchoice;



	std::string line;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string inputname;

		iss >> inputname;

		if (inputname == "length") iss >> length;

		if (inputname == "numGC") iss >> numGC;

		if (inputname == "lsnumGC") iss >> lsnumGC;

		if (inputname == "lslength") iss >> lslength;

		if (inputname == "fluid1_gamma") iss >> fluid1_gamma;

		if (inputname == "fluid1_B") iss >> fluid1_B;

		if (inputname == "fluid2_gamma") iss >> fluid2_gamma;

		if (inputname == "fluid2_B") iss >> fluid2_B;

		if (inputname == "CFL") iss >> CFL;
		
		if (inputname == "RS_pure")
		{
			iss >> RS_purechoice;
			if (RS_purechoice == "HLLC_idealgas") RS_pure = HLLC_idealgas;
			else if (RS_purechoice == "exact_idealgas") RS_pure = exact_idealgas;
			else assert(!"Invalid pure RS");
		}

		if (inputname == "RS_mixed")
		{
			iss >> RS_mixedchoice;
			if (RS_mixedchoice == "M_HLLC") RS_mixed = M_HLLC;
			else if (RS_mixedchoice == "exact_idealgas") RS_mixed = exact_idealgas;
			else assert(!"Invalid mixed RS");
		}
	
		if (inputname == "FS")
		{
			iss >> FSchoice;
			if (FSchoice == "Godunov") FS = Godunov;
			if (FSchoice == "MUSCL_FS") FS = MUSCL_FS;
		}

		if (inputname == "GFM")
		{
			iss >> GFMchoice;
			if (GFMchoice == "Original") GFM = Original;
			else if (GFMchoice == "Isobaricfix") GFM = Isobaricfix;
			else if (GFMchoice == "Real") GFM = Real;
			else assert(!"Invalid GFM");
		}

		if (inputname == "IC")
		{
			iss >> ICchoice;
			if (ICchoice == "TTC1") IC = TTC1;
			else if (ICchoice == "TTC2") IC = TTC2;
			else if (ICchoice == "TTC3") IC = TTC3;
			else if (ICchoice == "TTC4") IC = TTC4;
			else if (ICchoice == "TTC5") IC = TTC5;
			else assert(!"Invalid IC");
		}

		if (inputname == "eos1")
		{
			iss >> eos1choice;
			if (eos1choice == "ideal") eos1 = ideal;
		}

		if (inputname == "eos2")
		{
			iss >> eos2choice;
			if (eos2choice == "ideal") eos2 = ideal;
		}

		if (inputname == "BC_L")
		{
			iss >> BCLchoice;
			if (BCLchoice == "reflective") BC_L = reflective;
			if (BCLchoice == "transmissive") BC_L = transmissive;
			if (BCLchoice == "periodic") BC_L = periodic;
		}
		
		if (inputname == "BC_R")
		{
			iss >> BCRchoice;
			if (BCRchoice == "reflective") BC_R = reflective;
			if (BCRchoice == "transmissive") BC_R = transmissive;
			if (BCRchoice == "periodic") BC_R = periodic;
		}

		if (inputname == "sim")
		{
			iss >> simchoice;
			if (simchoice == "onefluid") sim = onefluid;
			if (simchoice == "twofluid") sim = twofluid;
		}

		if (inputname == "outputpath") iss >> outputpath;

		

	}

	basename = outputpath + GFMchoice + "_" + FSchoice + "_" + ICchoice + "_" 
		+ eos1choice + "-" + eos2choice + "_" + std::to_string(length) + "_";

	if (simchoice == "onefluid") basename = outputpath + "onefluid_" + FSchoice + "_" + ICchoice + "_" + eos1choice + "_" + std::to_string(length) + "_";
	
	infile.close();
}
