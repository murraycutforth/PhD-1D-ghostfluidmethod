#ifndef INPUT_H
#define INPUT_H


#include <string>


struct settingsfile {

	int length;
	double x0;
	double dx;
	int numGC;
	int lsnumGC;
	int lslength;
	double lsdx;
	double fluid1_gamma;
	double fluid1_B;
	double fluid2_gamma;
	double fluid2_B;
	double T;
	double CFL;
	std::string RS_pure;
	std::string RS_mixed;
	std::string FS;
	std::string GFM;
	std::string IC;
	std::string eos1;
	std::string eos2;
	std::string BC_L;
	std::string BC_R;
	std::string sim;
	std::string outputpath;
	std::string basename;

	void read_settings_file ();
};

#endif
