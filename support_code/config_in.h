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
string output_folder;

// No idea what there parameters are
double eq_r;
double polar_r;

// Body parameters
double* depths;
double transition_zone;
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
double contaminant_effect;
double cmb_contrast;
double pressure_scale;
double q;
double ice_k;
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
unsigned int present_timestep;
unsigned int total_viscous_steps;
}

class config_in
{
public:
	config_in(char*);
	
private:
	void write_config();
};

void config_in::write_config()
{
	std::ostringstream config_parameters;
	config_parameters << system_parameters::output_folder << "/run_parameters.txt";
	std::ofstream fout_config(config_parameters.str().c_str());
	fout_config << "mesh filename: " << system_parameters::mesh_filename << endl << endl;
    fout_config << "eq_r = " << system_parameters::eq_r << endl;
    fout_config << "polar_r = " << system_parameters::polar_r << endl;
    fout_config << "transition_zone = " << system_parameters::transition_zone << endl;
    fout_config << "crust_thickness = " << system_parameters::crust_thickness << endl;
    fout_config << "r_core_eq = " << system_parameters::r_core_eq << endl;
    fout_config << "r_core_polar = " << system_parameters::r_core_polar << endl;
    fout_config << "mantle_rho = " << system_parameters::mantle_rho << endl;
    fout_config << "core_rho = " << system_parameters::core_rho << endl;
    fout_config << "period = " << system_parameters::period << endl;
    fout_config << "omegasquared = " << system_parameters::omegasquared << endl;
    fout_config << "eta_ceiling = " << system_parameters::eta_ceiling << endl;
    fout_config << "eta_floor = " << system_parameters::eta_floor << endl;
	fout_config << "contaminant_effect = " << system_parameters::contaminant_effect << endl;
    fout_config << "pressure_scale = " << system_parameters::pressure_scale << endl;
    fout_config << "cmb_contrast = " << system_parameters::cmb_contrast << endl;
    fout_config << "q = " << system_parameters::q << endl;
	fout_config << "ice_k = " << system_parameters::ice_k << endl;
    fout_config << "ice_G = " << system_parameters::ice_G << endl;
    fout_config << "rock_G = " << system_parameters::rock_G << endl;
    fout_config << "cylindrical = " << system_parameters::cylindrical << endl;
    fout_config << "continue_plastic_iterations = " << system_parameters::continue_plastic_iterations << endl;
    fout_config << "plasticity_on = " << system_parameters::plasticity_on << endl;
    fout_config << "max_plastic_iterations = " << system_parameters::max_plastic_iterations << endl;
    fout_config << "smoothing_radius = " << system_parameters::smoothing_radius << endl;
    fout_config << "initial_elastic_iterations = " << system_parameters::initial_elastic_iterations << endl;
    fout_config << "elastic_time = " << system_parameters::elastic_time << endl;
    fout_config << "viscous_time = " << system_parameters::viscous_time << endl;
    fout_config << "initial_disp_target = " << system_parameters::initial_disp_target << endl;
    fout_config << "final_disp_target = " << system_parameters::final_disp_target << endl;
    fout_config << "current_time_interval = " << system_parameters::current_time_interval << endl;
    fout_config << "global_refinement = " << system_parameters::global_refinement << endl;
    fout_config << "small_r_refinement = " << system_parameters::small_r_refinement << endl;
    fout_config << "crustal_refinement = " << system_parameters::crustal_refinement << endl;
    fout_config << "surface_refinement = " << system_parameters::surface_refinement << endl;
	fout_config << "iteration_coefficient = " << system_parameters::iteration_coefficient << endl;
	fout_config << "tolerance_coefficient = " << system_parameters::tolerance_coefficient << endl;
	fout_config << "present_time = " << system_parameters::present_time << endl;
	fout_config << "present_timestep = " << system_parameters::present_timestep << endl;
	fout_config << "total_viscous_steps = " << system_parameters::total_viscous_steps << endl;
	
	fout_config.close();
}

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
        string out = cfg.lookup("output_folder");

        system_parameters::output_folder = out;
	    system_parameters::mesh_filename = msh;

        string output = cfg.lookup("output_folder");
	    system_parameters::output_folder = output;

	  }
	  catch(const SettingNotFoundException &nfex)
	  {
	    cerr << "No 'mesh_filename' or 'output_folder' setting in configuration file." << endl;
	  }

	  // get radii

	  const Setting& root = cfg.getRoot();

	  try
	  {
	    const Setting& radii = root["radii"];
	    radii.lookupValue("eq_r", system_parameters::eq_r);
	    radii.lookupValue("polar_r", system_parameters::polar_r);
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

//        cout << "Number of depth = " << ndepths << endl;
		system_parameters::depths = new double[ndepths];

		double d;
		for(unsigned int i=0; i<ndepths; i++)
		{
			d = set_depths[i];
			system_parameters::depths[i] = d;
//		    cout << "depth[" << i << "] = " << system_parameters::depths[i] << endl;
		}

	    const Setting& body_parameters = root["body_parameters"];
	    body_parameters.lookupValue("crust_thickness", system_parameters::crust_thickness);
		body_parameters.lookupValue("transition_zone", system_parameters::transition_zone);
	    body_parameters.lookupValue("r_core_eq", system_parameters::r_core_eq);
	    body_parameters.lookupValue("r_core_polar", system_parameters::r_core_polar);
	    body_parameters.lookupValue("mantle_rho", system_parameters::mantle_rho);
	    body_parameters.lookupValue("core_rho", system_parameters::core_rho);
	    body_parameters.lookupValue("period", system_parameters::period);
	    system_parameters::omegasquared = pow(TWOPI / 3600 / system_parameters::period, 2.0);
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
		rheology_parameters.lookupValue("contaminant_effect", system_parameters::contaminant_effect);
	    rheology_parameters.lookupValue("pressure_scale", system_parameters::pressure_scale);
		rheology_parameters.lookupValue("cmb_contrast", system_parameters::cmb_contrast);
	    rheology_parameters.lookupValue("q", system_parameters::q);
		rheology_parameters.lookupValue("ice_k", system_parameters::ice_k);
	    rheology_parameters.lookupValue("ice_G", system_parameters::ice_G);
	    rheology_parameters.lookupValue("rock_G", system_parameters::rock_G);
	    rheology_parameters.lookupValue("cylindrical", system_parameters::cylindrical);
	    rheology_parameters.lookupValue("continue_plastic_iterations", system_parameters::continue_plastic_iterations);
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
	  }
	  catch(const SettingNotFoundException &nfex)
	  {
		  cerr << "We've got a problem in the time step parameters block" << endl;
	  }
	  write_config();
}
}





