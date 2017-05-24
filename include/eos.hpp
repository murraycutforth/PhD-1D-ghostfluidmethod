#ifndef EOS
#define EOS

#include "data_storage.hpp"
#include <blitz/array.h>
#include <memory>
#include <string>


class eos_base {

	public:

	eos_base ();


	// General usage EOS functions

	virtual std::string get_eos_type () =0;

	virtual double a (blitz::Array<double,1> state) =0;
	
	virtual double a (double rho, double p) =0;

	virtual double p (blitz::Array<double,1> state) =0;

	virtual double E (blitz::Array<double,1> primitives) =0;

	virtual double E (double rho, double u, double p) =0;

	virtual double specific_ie_prim (blitz::Array<double,1> primitives) =0;

	virtual double specific_ie_prim (double rho, double u, double p) =0;

	virtual double get_gamma () =0;
	
	virtual double get_Pinf () =0;


	// Functions used by the original ghost fluid method

	virtual double rho_constant_entropy (double p_old, double rho_old, double p_new) =0;

	virtual double rho (double p, double S) =0;

	virtual double S (double p, double rho) =0;


	// Functions used by the M-HLLC solver

	virtual double get_Tau (blitz::Array<double,1> state) =0;

	virtual double get_Psi (blitz::Array<double,1> state) =0;

	virtual double postshock_density (double P_star, double P_K, double rho_K) =0;

	virtual double postrarefaction_density (double P_star, double P_K, double rho_K) =0;


};





class eos_idealgas : public eos_base {

	public:

	double gamma;

	eos_idealgas (double gamma);


	std::string get_eos_type ();

	double get_gamma ();
	
	double get_Pinf ();

	double a (blitz::Array<double,1> state);
	
	double a (double rho, double p);

	double p (blitz::Array<double,1> state);

	double E (blitz::Array<double,1> primitives);

	double E (double rho, double u, double p);	
	
	double specific_ie_prim (blitz::Array<double,1> primitives);

	double specific_ie_prim (double rho, double u, double p);


	double rho_constant_entropy (double p_old, double rho_old, double p_new);

	double rho (double p, double S);

	double S (double p, double rho);



	double get_Tau (blitz::Array<double,1> state);

	double get_Psi (blitz::Array<double,1> state);

	double postshock_density (double P_star, double P_K, double rho_K);

	double postrarefaction_density (double P_star, double P_K, double rho_K);

	
};




class eos_stiffenedgas : public eos_base {

	public:

	double gamma;
	
	double Pinf;

	eos_stiffenedgas (double gamma, double Pinf);


	std::string get_eos_type ();

	double get_gamma ();
	
	double get_Pinf ();

	double a (blitz::Array<double,1> state);
	
	double a (double rho, double p);

	double p (blitz::Array<double,1> state);

	double E (blitz::Array<double,1> primitives);

	double E (double rho, double u, double p);	
	
	double specific_ie_prim (blitz::Array<double,1> primitives);

	double specific_ie_prim (double rho, double u, double p);


	double rho_constant_entropy (double p_old, double rho_old, double p_new);

	double rho (double p, double S);

	double S (double p, double rho);



	double get_Tau (blitz::Array<double,1> state);

	double get_Psi (blitz::Array<double,1> state);

	double postshock_density (double P_star, double P_K, double rho_K);

	double postrarefaction_density (double P_star, double P_K, double rho_K);

	
};



#endif
