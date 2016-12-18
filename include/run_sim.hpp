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



//class twofluid_sim : public sim_base {
//
//	public:
//
//	void run_sim (settingsfile SF);
//};

#endif
