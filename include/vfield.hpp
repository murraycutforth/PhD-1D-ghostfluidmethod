#ifndef VFIELD_H
#define VFIELD_H



#include "data_storage.hpp"
#include "flow_solver.hpp"
#include "ghost_fluid_method.hpp"



#include <memory>




class vfield_base {

	public:

	double t;

	
	vfield_base ();

	virtual double get_u (double x) =0;
};





class vfield_test1 : public vfield_base {

	public:


	vfield_test1 ();

	double get_u (double x);
};




class vfield_oldstate : public vfield_base {

	public:

	twofluid_array& states;
	levelset_array& ls;


	vfield_oldstate (twofluid_array& oldstates, levelset_array& oldls);

	double get_u (double x);
};



class vfield_starstate : public vfield_base {

	public:

	std::shared_ptr<flow_solver_base> FS;
	levelset_array& ls;
	arrayinfo statearray;


	vfield_starstate (	std::shared_ptr<flow_solver_base> FS,
				levelset_array& ls,
				arrayinfo statearray);

	double get_u (double x);
};






class vfield_mixedRPsolution : public vfield_base {

	public:

	std::shared_ptr<ghost_fluid_method_base> GFM;
	levelset_array& ls;


	vfield_mixedRPsolution (std::shared_ptr<ghost_fluid_method_base> GFM, levelset_array& ls);

	double get_u (double x);
};




#endif
