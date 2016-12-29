#ifndef RUNSIM_H
#define RUNSIM_H


#include "input.hpp"
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

	public:

	void run_sim (settingsfile SF);

	double compute_dt (double CFL, double T, double t, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);

	void output_endoftimestep (int numsteps, settingsfile& SF, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);
	
	void output_endofsimulation (int numsteps, settingsfile& SF, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);

	void output_realfluidonly (std::string name, settingsfile& SF, fluid_state_array& state1, fluid_state_array& state2, levelset_array& ls);
};

#endif
