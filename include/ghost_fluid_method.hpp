#ifndef GHOST_FLUID_METHOD
#define GHOST_FLUID_METHOD


#include "data_storage.hpp"
#include "riemann_solver.hpp"
#include <memory>


class GFM_base {

	public:

	blitz::Array<double,1> extension_interface_velocity;

	
	GFM_base (arrayinfo array);

	virtual void set_ghost_cells (	
	
		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls, 
		std::shared_ptr<multimat_RS_base> RS
	) =0;
};


class Original_GFM : public GFM_base {

	public:

	Original_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls, 
		std::shared_ptr<multimat_RS_base> RS
	);
};


class R_GFM : public GFM_base {

	public:

	R_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls, 
		std::shared_ptr<multimat_RS_base> RS
	);
};
	

#endif
