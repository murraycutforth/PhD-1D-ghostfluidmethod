#ifndef DATA_STORAGE
#define DATA_STORAGE




class eos_base;



#include <blitz/array.h>
#include <memory>
#include <string>






enum RS_type {HLLC, PVRS};

enum FS_type {Godunov};

enum GFM_type {Original};

enum IC_type {TC1, TC2};

enum ls_IC_type {T1};

enum eos_type {ideal, tammann};

enum BC_type {reflective, transmissive, nothing};



struct settingsfile {

	int length;
	double x0;
	double dx;
	int numGC;

	double fluid1_gamma;
	double fluid1_A;
	double fluid2_gamma;
	double fluid2_A;

	RS_type RS;
	FS_type FS;
	GFM_type GFM;
	IC_type IC;
};




struct arrayinfo {

	int length;
	double x0;
	double dx;
	int numGC;
	BC_type leftBC;
	BC_type rightBC;

	double cellcentre_coord (int i);
	int cellindex (double x);
};


bool operator==(const arrayinfo& rhs, const arrayinfo& lhs);




class levelset_array;




class twofluid_array {

	public:

	blitz::Array<double,2> fluid1;
	blitz::Array<double,2> fluid2;
	arrayinfo array;
	std::shared_ptr<eos_base> eos1;
	std::shared_ptr<eos_base> eos2;


	twofluid_array (arrayinfo array, std::shared_ptr<eos_base> eos1, std::shared_ptr<eos_base> eos2);
	

	void apply_BCs ();

	void output_to_file (std::string name);

	void output_realfluid_to_file (std::string name, levelset_array& ls);

};





class onefluid_array {

	public:

	blitz::Array<double,2> fluid;
	arrayinfo array;
	std::shared_ptr<eos_base> eos;


	onefluid_array (arrayinfo array, std::shared_ptr<eos_base> eos);

	onefluid_array (blitz::Array<double,2>& fluid, arrayinfo array, std::shared_ptr<eos_base> eos);


	void apply_BCs();

	void output_to_file (std::string name);

	blitz::Array<double,1> total_conserved_quantities ();
};




class levelset_array {

	public:

	blitz::Array<double,1> phi;
	arrayinfo array;


	
	levelset_array (arrayinfo array);

	double operator() (double x);

	double linear_interpolation (double x);

	void apply_BCs();

	void output_to_file (std::string name);
};



#endif
