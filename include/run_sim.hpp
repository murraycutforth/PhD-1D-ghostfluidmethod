#ifndef RUNSIM_H
#define RUNSIM_H


#include "input.hpp"
#include "eos.hpp"
#include "flow_solver.hpp"
#include "riemann_solver.hpp"
#include "ghost_fluid_method.hpp"
#include "data_storage.hpp"


class sim_base {

	public:

	virtual void run_sim (settingsfile SF) =0;
};



class onefluid_sim : public sim_base {

	public:

	void run_sim (settingsfile SF);

	double compute_dt (double CFL, fluid_state_array& state, double T, double t);
};



class twofluid_sim : public sim_base {
	
	private:
	
	std::shared_ptr<eos_base> eos1;
	std::shared_ptr<eos_base> eos2;
	std::shared_ptr<singlefluid_RS_base> RS_pure;
	std::shared_ptr<multimat_RS_base> RS_mixed;
	std::shared_ptr<flow_solver_base> FS;
	std::shared_ptr<GFM_base> GFM;
	fluid_state_array statearr1;
	fluid_state_array statearr2;
	levelset_array ls;
	levelset_array prev_ls;
	blitz::Array<double,1> FL1;
	blitz::Array<double,1> FR1;
	blitz::Array<double,1> FL2;
	blitz::Array<double,1> FR2;
	blitz::Array<double,1> U0;
	blitz::Array<double,1> Ut;
	

	public:
	
	twofluid_sim ();

	void run_sim (settingsfile SF);

	double compute_dt (double CFL, double T, double t, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);

	void output_endoftimestep (int numsteps, settingsfile& SF, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);
	
	void output_endofsimulation (int numsteps, settingsfile& SF, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);

	void output_realfluidonly (std::string name, settingsfile& SF, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);
};

#endif
