#ifndef GHOST_FLUID_METHOD
#define GHOST_FLUID_METHOD



#include "data_storage.hpp"
#include "riemann_solver.hpp"


#include <memory>



class ghost_fluid_method_base {

	public:

	blitz::Array<double,1> extension_interface_velocity;

	
	ghost_fluid_method_base (arrayinfo array);

	virtual void set_ghost_cells (	twofluid_array& states, 
					levelset_array& ls, 
					std::shared_ptr<riemann_solver_base> rs) =0;
};



class original_GFM : public ghost_fluid_method_base {

	public:

	original_GFM (arrayinfo array);

	void set_ghost_cells (	twofluid_array& states, 
				levelset_array& ls,
				std::shared_ptr<riemann_solver_base> rs);

};






class isobaric_fix_GFM : public ghost_fluid_method_base {

	public:

	isobaric_fix_GFM (arrayinfo array);

	void set_ghost_cells (	twofluid_array& states, 
				levelset_array& ls,
				std::shared_ptr<riemann_solver_base> rs);

};



class rGFM : public ghost_fluid_method_base {

	public:

	rGFM (arrayinfo array);

	void set_ghost_cells (	twofluid_array& states,
				levelset_array& ls,
				std::shared_ptr<riemann_solver_base> rs);

};



#endif
