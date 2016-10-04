#ifndef RUNSIM_H
#define RUNSIM_H



#include "data_storage.hpp"



class sim_base {

	public:

	virtual void run_sim (settingsfile SF) =0;
};



class serial_onefluid_sim : public sim_base {

	public:


	void run_sim (settingsfile SF);
};



class serial_twofluid_sim : public sim_base {

	public:


	void run_sim (settingsfile SF);
};





#endif
