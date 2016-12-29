#ifndef DATA_STORAGE
#define DATA_STORAGE


#include "input.hpp"
#include <blitz/array.h>
#include <memory>
#include <string>


class eos_base;



struct arrayinfo {

	int length;
	double x0;
	double dx;
	int numGC;
	std::string leftBC;
	std::string rightBC;

	double cellcentre_coord (int i);
	int cellindex (double x);
};

bool operator==(const arrayinfo& rhs, const arrayinfo& lhs);



class fluid_state_array {

	public:
	
	arrayinfo array;
	blitz::Array<double,2> CV;
	std::shared_ptr<eos_base> eos;
	

	fluid_state_array (arrayinfo array, std::shared_ptr<eos_base> eos);

	fluid_state_array (const fluid_state_array& other);

	fluid_state_array ();
	

	fluid_state_array copy ();

	double get_u (int i);
	
	void apply_BCs ();

	void output_to_file (std::string name);

	blitz::Array<double,1> total_conserved_quantities ();
};



class levelset_array {

	public:

	arrayinfo array;
	blitz::Array<double,1> phi;

	
	levelset_array (arrayinfo array);
	
	levelset_array (const levelset_array& other);

	levelset_array ();

	
	levelset_array copy ();

	double linear_interpolation (double x);

	double normal (double x);

	void advection_step (double dt, blitz::Array<double,1> vfield);

	void apply_BCs();

	void output_to_file (std::string name);
};



#endif
