#ifndef INITIAL_CONDITIONS
#define INITIAL_CONDITIONS



#include "data_storage.hpp"




void initialise_onefluid (onefluid_array& state, IC_type IC);



void initialise_twofluid (twofluid_array& states, IC_type IC);



void initialise_levelset (levelset_array& ls, ls_IC_type IC);



#endif
