#ifndef LEVELSET_ADVECTION
#define LEVELSET_ADVECTION



#include "data_storage.hpp"
#include "vfield.hpp"




#include <memory>





class levelset_advection_base {

	public:

	virtual void advection_operator (	levelset_array& ls,
						levelset_array& newls,
						std::shared_ptr<vfield_base> vfield,
						double dt) =0;
};





class levelset_advection_firstorderup : public levelset_advection_base {

	public:

	void advection_operator (	levelset_array& ls,
					levelset_array& newls,
					std::shared_ptr<vfield_base> vfield,
					double dt);
};




#endif
