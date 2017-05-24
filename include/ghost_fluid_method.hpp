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
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	) =0;
};


class Original_GFM : public GFM_base {

	public:

	Original_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls,  
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	);
};

class newR_GFM : public GFM_base {

	public:

	newR_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls,  
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	);
};



class R_GFM : public GFM_base {

	public:

	R_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls,  
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	);
};


class M_GFM : public GFM_base {

	public:

	M_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls,  
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	);
};


class P_GFM : public GFM_base {

	public:

	P_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls,  
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	);
};


class newmethod1_GFM : public GFM_base {

	public:

	newmethod1_GFM (arrayinfo array);

	void set_ghost_cells (

		fluid_state_array& state1,
		fluid_state_array& state2,
		levelset_array& ls,  
		levelset_array& ls_prev,
		std::shared_ptr<multimat_RS_base> RS,
		double dt
	);
};
	

#endif
