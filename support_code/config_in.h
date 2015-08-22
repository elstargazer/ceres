/*
 * config_in.h
 *
 *  Created on: Aug 17, 2015
 *      Author: antonermakov
 */


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <libconfig.h++>

#include "local_math.h"

using namespace std;
using namespace libconfig;

namespace Step22 {
namespace system_parameters {

// Mesh file name
string mesh_filename;

// No idea what there parameters are
double eq_r;
double polar_r;

// Body parameters
double* depths;
unsigned int sizeof_depths;
double crust_thickness;
double r_core_eq;
double r_core_polar;
double mantle_rho;
double core_rho;
double period;
double omegasquared;

// Rheology parameters
double eta_ceiling;
double eta_floor;
double pressure_scale;
double q;
double ice_G;
double rock_G;
bool cylindrical;
bool continue_plastic_iterations;

// plasticity variables
bool plasticity_on;
unsigned int max_plastic_iterations;
double smoothing_radius;

// viscoelasticity variables
unsigned int initial_elastic_iterations;
double elastic_time;
double viscous_time;
double initial_disp_target;
double final_disp_target;
double current_time_interval;

//mesh refinement variables
unsigned int global_refinement;
unsigned int small_r_refinement;
unsigned int crustal_refinement;
unsigned int surface_refinement;

//solver variables
int iteration_coefficient;
double tolerance_coefficient;

//time step variables
double present_time;
double present_timestep;
double total_viscous_steps;
}

class config_in
{
public:
	config_in(char*);
};

config_in::config_in(char* filename)
{

	// This example reads the configuration file 'example.cfg' and displays
	// some of its contents.

	  Config cfg;

	  // Read the file. If there is an error, report it and exit.
	  try
	  {
	    cfg.readFile(filename);
	  }
	  catch(const FileIOException &fioex)
	  {
	    std::cerr << "I/O error while reading file." << std::endl;

	  }
	  catch(const ParseException &pex)
	  {
	    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
	              << " - " << pex.getError() << std::endl;

	  }

	  // Get mesh name.
	  try
	  {
        string msh = cfg.lookup("mesh_filename");
	    system_parameters::mesh_filename = msh;
	    cout << "mesh filename: " << system_parameters::mesh_filename << endl << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
	    cerr << "No 'mesh_filename' setting in configuration file." << endl;
	  }

	  // get radii

	  const Setting& root = cfg.getRoot();

	  try
	  {
	    const Setting& radii = root["radii"];
	    radii.lookupValue("eq_r", system_parameters::eq_r);
	    radii.lookupValue("polar_r", system_parameters::polar_r);

	    cout << "eq_r = " << system_parameters::eq_r << endl;
	    cout << "polar_r = " << system_parameters::polar_r << endl;

	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the get radii block" << endl;
	  }

	  // get body parameters
	  try
	  {

		const Setting& set_depths = cfg.lookup("body_parameters.depths");

		unsigned int ndepths = set_depths.getLength();
		system_parameters::sizeof_depths = ndepths;

        cout << "Number of depth = " << ndepths << endl;
		system_parameters::depths = new double[ndepths];

		double d;
		for(int i=0; i<ndepths; i++)
		{
			d = set_depths[i];
			system_parameters::depths[i] = d;
		    cout << "depth[" << i << "] = " << system_parameters::depths[i] << endl;
		}

	    const Setting& body_parameters = root["body_parameters"];
	    body_parameters.lookupValue("crust_thickness", system_parameters::crust_thickness);
	    body_parameters.lookupValue("r_core_eq", system_parameters::r_core_eq);
	    body_parameters.lookupValue("r_core_polar", system_parameters::r_core_polar);
	    body_parameters.lookupValue("mantle_rho", system_parameters::mantle_rho);
	    body_parameters.lookupValue("core_rho", system_parameters::core_rho);
	    body_parameters.lookupValue("period", system_parameters::period);
	    system_parameters::omegasquared = pow(TWOPI / 3600 / system_parameters::period, 2.0);

	    cout << "crust_thickness = " << system_parameters::crust_thickness << endl;
	    cout << "r_core_eq = " << system_parameters::r_core_eq << endl;
	    cout << "r_core_polar = " << system_parameters::r_core_polar << endl;
	    cout << "mantle_rho = " << system_parameters::mantle_rho << endl;
	    cout << "core_rho = " << system_parameters::core_rho << endl;
	    cout << "period = " << system_parameters::period << endl;
	    cout << "omegasquared = " << system_parameters::omegasquared << endl;

	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the body parameters block" << endl;

	  }

	// Rheology parameters
	  try
	  {

	    const Setting& rheology_parameters = root["rheology_parameters"];
	    rheology_parameters.lookupValue("eta_ceiling", system_parameters::eta_ceiling);
	    rheology_parameters.lookupValue("eta_floor", system_parameters::eta_floor);
	    rheology_parameters.lookupValue("pressure_scale", system_parameters::pressure_scale);
	    rheology_parameters.lookupValue("q", system_parameters::q);
	    rheology_parameters.lookupValue("ice_G", system_parameters::ice_G);
	    rheology_parameters.lookupValue("rock_G", system_parameters::rock_G);
	    rheology_parameters.lookupValue("cylindrical", system_parameters::cylindrical);
	    rheology_parameters.lookupValue("continue_plastic_iterations", system_parameters::continue_plastic_iterations);

	    cout << "eta_ceiling = " << system_parameters::eta_ceiling << endl;
	    cout << "eta_floor = " << system_parameters::eta_floor << endl;
	    cout << "pressure_scale = " << system_parameters::pressure_scale << endl;
	    cout << "q = " << system_parameters::q << endl;
	    cout << "ice_G = " << system_parameters::ice_G << endl;
	    cout << "rock_G = " << system_parameters::rock_G << endl;
	    cout << "cylindrical = " << system_parameters::cylindrical << endl;
	    cout << "continue_plastic_iterations = " << system_parameters::continue_plastic_iterations << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the rheology parameters block" << endl;
	  }

	  // Plasticity parameters
	  try
	  {

	    const Setting& plasticity_parameters = root["plasticity_parameters"];
	    plasticity_parameters.lookupValue("plasticity_on", system_parameters::plasticity_on);
	    plasticity_parameters.lookupValue("max_plastic_iterations", system_parameters::max_plastic_iterations);
	    plasticity_parameters.lookupValue("smoothing_radius", system_parameters::smoothing_radius);

	    cout << "plasticity_on = " << system_parameters::plasticity_on << endl;
	    cout << "max_plastic_iterations = " << system_parameters::max_plastic_iterations << endl;
	    cout << "smoothing_radius = " << system_parameters::smoothing_radius << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the plasticity parameters block" << endl;
	  }

	// Viscoelasticity parameters

	  try
	  {

	    const Setting& viscoelasticity_parameters = root["viscoelasticity_parameters"];
	    viscoelasticity_parameters.lookupValue("initial_elastic_iterations", system_parameters::initial_elastic_iterations);
	    viscoelasticity_parameters.lookupValue("elastic_time", system_parameters::elastic_time);
	    viscoelasticity_parameters.lookupValue("viscous_time", system_parameters::viscous_time);
	    viscoelasticity_parameters.lookupValue("initial_disp_target", system_parameters::initial_disp_target);
	    viscoelasticity_parameters.lookupValue("final_disp_target", system_parameters::final_disp_target);
	    viscoelasticity_parameters.lookupValue("current_time_interval", system_parameters::current_time_interval);

	    system_parameters::viscous_time *= SECSINYEAR;


	    cout << "initial_elastic_iterations = " << system_parameters::initial_elastic_iterations << endl;
	    cout << "elastic_time = " << system_parameters::elastic_time << endl;
	    cout << "viscous_time = " << system_parameters::viscous_time << endl;
	    cout << "initial_disp_target = " << system_parameters::initial_disp_target << endl;
	    cout << "final_disp_target = " << system_parameters::final_disp_target << endl;
	    cout << "current_time_interval = " << system_parameters::current_time_interval << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the viscoelasticity parameters block" << endl;
	  }

	// Mesh refinement parameters
	  try
	  {

	    const Setting& mesh_refinement_parameters = root["mesh_refinement_parameters"];
	    mesh_refinement_parameters.lookupValue("global_refinement", system_parameters::global_refinement);
	    mesh_refinement_parameters.lookupValue("small_r_refinement", system_parameters::small_r_refinement);
	    mesh_refinement_parameters.lookupValue("crustal_refinement", system_parameters::crustal_refinement);
	    mesh_refinement_parameters.lookupValue("surface_refinement", system_parameters::surface_refinement);

	    cout << "global_refinement = " << system_parameters::global_refinement << endl;
	    cout << "small_r_refinement = " << system_parameters::small_r_refinement << endl;
	    cout << "crustal_refinement = " << system_parameters::crustal_refinement << endl;
	    cout << "surface_refinement = " << system_parameters::surface_refinement << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		 cerr << "We've got a problem in the mesh refinement parameters block" << endl;
	  }

	  // Solver parameters
	  try
	  {
	    const Setting& solve_parameters = root["solve_parameters"];
	    solve_parameters.lookupValue("iteration_coefficient", system_parameters::iteration_coefficient);
	    solve_parameters.lookupValue("tolerance_coefficient", system_parameters::tolerance_coefficient);

	    cout << "iteration_coefficient = " << system_parameters::iteration_coefficient << endl;
	    cout << "tolerance_coefficient = " << system_parameters::tolerance_coefficient << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the solver parameters block" << endl;
	  }

	  // Time step parameters
	  try
	  {
	    const Setting& time_step_parameters = root["time_step_parameters"];
	    time_step_parameters.lookupValue("present_time", system_parameters::present_time);
	    time_step_parameters.lookupValue("present_timestep", system_parameters::present_timestep);
	    time_step_parameters.lookupValue("total_viscous_steps", system_parameters::total_viscous_steps);

	    cout << "present_time = " << system_parameters::present_time << endl;
	    cout << "present_timestep = " << system_parameters::present_timestep << endl;
	    cout << "total_viscous_steps = " << system_parameters::total_viscous_steps << endl;
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the time step parameters block" << endl;
	  }
}
}





