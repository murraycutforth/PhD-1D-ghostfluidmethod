#ifndef INPUT_H
#define INPUT_H


#include <string>


enum RS_type {HLLC_idealgas, M_HLLC, exact_idealgas};

enum FS_type {Godunov, MUSCL_FS};

enum GFM_type {Original, Isobaricfix, Real};

enum IC_type {TTC1, TTC2, TTC3, TTC4, TTC5};

enum eos_type {ideal, tait};

enum BC_type {reflective, transmissive, periodic, nothing};

enum sim_type {onefluid, twofluid};



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

	RS_type RS_pure;
	RS_type RS_mixed;
	FS_type FS;
	GFM_type GFM;
	IC_type IC;
	eos_type eos1;
	eos_type eos2;
	BC_type BC_L;
	BC_type BC_R;

	sim_type sim;

	std::string outputpath;
	std::string basename;

	void read_settings_file ();
};

#endif
