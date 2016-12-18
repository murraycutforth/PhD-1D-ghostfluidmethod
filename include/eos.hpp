#ifndef EOS
#define EOS

#include "data_storage.hpp"
#include <blitz/array.h>
#include <memory>


double specific_ie_cv (blitz::Array<double,1> state);

bool is_state_physical (blitz::Array<double,1> state);


class eos_base {

	public:

	eos_base ();


	virtual double a (blitz::Array<double,1> state) =0;

	virtual double p (blitz::Array<double,1> state) =0;

	virtual double E (blitz::Array<double,1> primitives) =0;

	virtual double E (double rho, double u, double p) =0;

	virtual double rho_pS (double p, double S) =0;

	virtual double S_prho (double p, double rho) =0;

	virtual double specific_ie_prim (blitz::Array<double,1> primitives) =0;

	virtual double specific_ie_prim (double rho, double u, double p) =0;

	virtual double get_gamma () =0;

	virtual double get_Tau (blitz::Array<double,1> state) =0;

	virtual double get_Psi (blitz::Array<double,1> state) =0;

	virtual double postshock_density (double P_star, double P_K, double rho_K) =0;

	virtual double postrarefaction_density (double P_star, double P_K, double rho_K) =0;

	virtual eos_type get_eos_type () =0;
};





class eos_idealgas : public eos_base {

	public:

	double gamma;

	eos_idealgas (double gamma);


	double a (blitz::Array<double,1> state);

	double p (blitz::Array<double,1> state);

	double E (blitz::Array<double,1> primitives);

	double E (double rho, double u, double p);	

	double rho_pS (double p, double S);

	double S_prho (double p, double rho);

	double specific_ie_prim (blitz::Array<double,1> primitives);

	double specific_ie_prim (double rho, double u, double p);
	
	double get_Tau (blitz::Array<double,1> state);

	double get_Psi (blitz::Array<double,1> state);

	double get_gamma ();
	
	double postshock_density (double P_star, double P_K, double rho_K);

	double postrarefaction_density (double P_star, double P_K, double rho_K);

	eos_type get_eos_type ();
};




class eos_tait : public eos_base {

	public:

	double gamma;
	double B;

	eos_tait (double gamma, double B);


	double a (blitz::Array<double,1> state);

	double p (blitz::Array<double,1> state);

	double E (blitz::Array<double,1> primitives);

	double E (double rho, double u, double p);	

	double rho_pS (double p, double S);

	double S_prho (double p, double rho);

	double specific_ie_prim (blitz::Array<double,1> primitives);

	double specific_ie_prim (double rho, double u, double p);
	
	double get_Tau (blitz::Array<double,1> state);

	double get_Psi (blitz::Array<double,1> state);

	double get_gamma ();

	double postshock_density (double P_star, double P_K, double rho_K);

	double postrarefaction_density (double P_star, double P_K, double rho_K);

	eos_type get_eos_type ();
};


#endif
