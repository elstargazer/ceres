/*
Viscoelastoplastic relaxation of Ceres
Author: Roger R. Fu
Adapted from Fu et al. 2014 Icarus 240, 133-145 starting Oct. 19, 2014
 */

/*
Summary of output files:

One per run:

- initial_mesh.eps					:  Visualization of initially imported mesh
- physical_times.txt				:  Columns are (1) step number corresponding to other files, (2) physical times at the time when each calculation is run in sec, (3) number of the final plasticity iteration in each timestep.  Written in do_elastic_steps() for elastic steps and do_flow_step() for viscous steps

One per timestep:

- timeXX_elastic_displacements.txt	:  Vtk-readable file with columns (1) x, (2) y, (3) u_x, (4) u_y, (5) P.  Written in output_results() function, which is run immediately after solve().
- timeXX_baseviscosities.txt		:  Columns (1) cell x, (2) cell y, (3) base viscosity in Pa s.  Written in solution_stresses().
- timeXX_surface.txt				:  Surface (defined as where P=0 boundary condition is applied) vertices at the beginning of timestep, except for the final timestep.  Written in write_vertices() function, which is called immediately after setup_dofs() except for the final iteration, when it is called after move_mesh()

One per plasticity step:

- timeXX_flowYY.txt					:  Same as timeXX_elastic_displacements.txt above
- timeXX_principalstressesYY.txt	:  Columns with sigma1 and sigma3 at each cell.  Same order as timeXX_baseviscosities.txt.  Written in solution_stresses().
- timeXX_stresstensorYY.txt			:  Columns with components 11, 22, 33, and 13 of stress tensor at each cell.  Written in solution_stresses().
- timeXX_failurelocations00.txt		:  Gives x,y coordinates of all cells where failure occurred.  Written in solution_stresses().
- timeXX_viscositiesregYY.txt		:  Gives smoothed and regularized (i.e., floor and ceiling-filtered) effective viscosities.  Written at end of solution_stresses().

*/

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>

// Then we need to include the header file
// for the sparse direct solver UMFPACK:
#include <deal.II/lac/sparse_direct.h>

// This includes the libary for the
// incomplete LU factorization that will
// be used as a preconditioner in 3D:
#include <deal.II/lac/sparse_ilu.h>
#include "support_code/ellipsoid_grav.h"
#include <fstream>
#include <sstream>
#include <time.h>
#include <armadillo>

#include "support_code/ellipsoid_fit.h"


// As in all programs, the namespace dealii
// is included:
namespace Step22 {
namespace system_parameters {
double eq_r = 0;
double polar_r = 0;

double depths[] = {62000, 58000, -10};
double crust_thickness = depths[0];
double r_core_eq = 0;
double r_core_polar = 0;
double mantle_rho = 1465.7;
double core_rho = 2491.1;
double period = 1000;
double omegasquared = 2 * 3.141592653589793 / 3600 / period * 2
		* 3.141592653589793 / 3600 / period;
double eta_ceiling = 1e17 * 1e5;
double eta_floor = eta_ceiling / 1e5;
double pressure_scale = eta_ceiling / 1e3;
double q = 2; // heat flux in mW/m2
double ice_G = 9.33e9;//Bland et al. 2013
double rock_G = 40e9;
bool cylindrical = true;
bool continue_plastic_iterations = true;
char mesh_filename[] = "meshes/mesh_test_def_quad_2.inp";

//plasticity variables
bool plasticity_on = false;
unsigned int max_plastic_iterations = 50;
double smoothing_radius = 10000;

//viscoelasticity variables
unsigned int initial_elastic_iterations = 1;
double elastic_time = 1;
double viscous_time = 3e3 * 3.1557e7;
double initial_disp_target = 6000;
double final_disp_target = 300;
double current_time_interval = 0;

//mesh variables
unsigned int global_refinement = 0;
unsigned int small_r_refinement = 0;
unsigned int crustal_refinement = 0;
unsigned int surface_refinement = 0;

//solver variables
int iteration_coefficient = 3000;
double tolerance_coefficient = 1e-10;

//time step variables
double present_time = 0;
double present_timestep = 0;
double total_viscous_steps = 2;
}
using namespace dealii;
using namespace arma;

template<int dim>
struct InnerPreconditioner;

template<>
struct InnerPreconditioner<2> {
	typedef SparseDirectUMFPACK type;
};

template<>
struct InnerPreconditioner<3> {
	typedef SparseILU<double> type;
};

// Auxiliary functions

template<int dim>
class AuxFunctions {
public:
	Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2> > &grad_u);
};

template<int dim>
Tensor<2, 2> AuxFunctions<dim>::get_rotation_matrix(
		const std::vector<Tensor<1, 2> > &grad_u) {
	const double curl = (grad_u[1][0] - grad_u[0][1]);

	const double angle = std::atan(curl);

	const double t[2][2] = { { cos(angle), sin(angle) }, { -sin(angle), cos(
			angle) } };
	return Tensor<2, 2>(t);
}

// Class for remembering material state/properties at each quadrature point

template<int dim>
struct PointHistory {
	SymmetricTensor<2, dim> old_stress;
	double old_phiphi_stress;
	double first_eta;
	double new_eta;
	double G;
};

// Primary class of this problem

template<int dim>
class StokesProblem {
public:
	StokesProblem(const unsigned int degree);
	void run();

private:
	void setup_dofs();
	void assemble_system();
	void solve();
	void output_results() const;
	void refine_mesh();
	void solution_stesses();
	std::vector<double> flow_law(double r, double z);
	double get_log_local_viscosity(double &r, double &z);
	double get_local_G(double &r, double &z);

	void setup_initial_mesh();
	void do_elastic_steps();
	void do_flow_step();
	void update_time_interval();
	void initialize_eta_and_G();
	void move_mesh();
	void write_vertices();
	void setup_quadrature_point_history();
	void update_quadrature_point_history();

	const unsigned int degree;

	Triangulation<dim> triangulation;
	const MappingQ1<dim> mapping;
	FESystem<dim> fe;
	DoFHandler<dim> dof_handler;
	unsigned int n_u = 0, n_p = 0;
	unsigned int plastic_iteration = 0;

	QGauss<dim> quadrature_formula;
	Vector<double> node_viscosities;
	Vector<double> quad_viscosities;
	Vector<double> cell_viscosities;
	std::vector<PointHistory<dim> > quadrature_point_history;

	ConstraintMatrix constraints;

	BlockSparsityPattern sparsity_pattern;
	BlockSparseMatrix<double> system_matrix;

	BlockVector<double> solution;
	BlockVector<double> system_rhs;

	std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;

	ellipsoid_fit<dim>   ellipsoid;
};

// Class for boundary conditions and rhs

template<int dim>
class BoundaryValuesP: public Function<dim> {
public:
	BoundaryValuesP() :
			Function<dim>(dim + 1) {
	}

	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

	virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template<int dim>
double BoundaryValuesP<dim>::value(const Point<dim> &p,
		const unsigned int component) const {
	Assert(component < this->n_components,
			ExcIndexRange (component, 0, this->n_components));

	Assert(p[0] >= -10, ExcLowerRange (p[0], 0)); //POSSIBLY FUDGED- LOOK INTO

	return 0;
}

template<int dim>
void BoundaryValuesP<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values) const {
	for (unsigned int c = 0; c < this->n_components; ++c)
		values(c) = BoundaryValuesP<dim>::value(p, c);
}

template<int dim>
class RightHandSide: public Function<dim> {
public:
	RightHandSide () : Function<dim>(dim+1) {}

	virtual double value(const Point<dim> &p, const unsigned int component,
			A_Grav_namespace::AnalyticGravity<dim> *aGrav) const;

	virtual void vector_value(const Point<dim> &p, Vector<double> &value,
			A_Grav_namespace::AnalyticGravity<dim> *aGrav) const;

	virtual void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &values,
			A_Grav_namespace::AnalyticGravity<dim> *aGrav) const;

};

template<int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
		const unsigned int component,
		A_Grav_namespace::AnalyticGravity<dim> *aGrav) const {

	std::vector<double> temp_vector(2);
	aGrav->get_gravity(p, temp_vector);

	if (component == 0) {
		return temp_vector[0] + system_parameters::omegasquared * p[0];	// * 1.2805;
	} else {
		if (component == 1)
			return temp_vector[1];
		else
			return 0;
	}
}

template<int dim>
void RightHandSide<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values,
		A_Grav_namespace::AnalyticGravity<dim> *aGrav) const {
	for (unsigned int c = 0; c < this->n_components; ++c)
		values(c) = RightHandSide<dim>::value(p, c, aGrav);
}

template<int dim>
void RightHandSide<dim>::vector_value_list(
		const std::vector<Point<dim> > &points,
		std::vector<Vector<double> > &values,
		A_Grav_namespace::AnalyticGravity<dim> *aGrav) const {
	// check whether component is in
	// the valid range is up to the
	// derived class
	Assert(values.size() == points.size(),
			ExcDimensionMismatch(values.size(), points.size()));

	for (unsigned int i = 0; i < points.size(); ++i)
		this->vector_value(points[i], values[i], aGrav);
}

// Class for linear solvers and preconditioners

template<class Matrix, class Preconditioner>
class InverseMatrix: public Subscriptor {
public:
	InverseMatrix(const Matrix &m, const Preconditioner &preconditioner);

	void vmult(Vector<double> &dst, const Vector<double> &src) const;

private:
	const SmartPointer<const Matrix> matrix;
	const SmartPointer<const Preconditioner> preconditioner;
};

template<class Matrix, class Preconditioner>
InverseMatrix<Matrix, Preconditioner>::InverseMatrix(const Matrix &m,
		const Preconditioner &preconditioner) :
		matrix(&m), preconditioner(&preconditioner) {
}

template<class Matrix, class Preconditioner>
void InverseMatrix<Matrix, Preconditioner>::vmult(Vector<double> &dst,
		const Vector<double> &src) const {
	SolverControl solver_control(1000 * src.size(), 1e-9 * src.l2_norm());

	SolverCG<> cg(solver_control);

	dst = 0;

	cg.solve(*matrix, dst, src, *preconditioner);
}

// Class for the SchurComplement

template<class Preconditioner>
class SchurComplement: public Subscriptor {
public:
	SchurComplement(const BlockSparseMatrix<double> &system_matrix,
			const InverseMatrix<SparseMatrix<double>, Preconditioner> &A_inverse);

	void vmult(Vector<double> &dst, const Vector<double> &src) const;

private:
	const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
	const SmartPointer<const InverseMatrix<SparseMatrix<double>, Preconditioner> > A_inverse;

	mutable Vector<double> tmp1, tmp2;
};

template<class Preconditioner>
SchurComplement<Preconditioner>::SchurComplement(
		const BlockSparseMatrix<double> &system_matrix,
		const InverseMatrix<SparseMatrix<double>, Preconditioner> &A_inverse) :
		system_matrix(&system_matrix), A_inverse(&A_inverse), tmp1(
				system_matrix.block(0, 0).m()), tmp2(
				system_matrix.block(0, 0).m()) {
}

template<class Preconditioner>
void SchurComplement<Preconditioner>::vmult(Vector<double> &dst,
		const Vector<double> &src) const {
	system_matrix->block(0, 1).vmult(tmp1, src);
	A_inverse->vmult(tmp2, tmp1);
	system_matrix->block(1, 0).vmult(dst, tmp2);
}

// StokesProblem::StokesProblem

template<int dim>
StokesProblem<dim>::StokesProblem(const unsigned int degree) :
		degree(degree),
		mapping(),
		triangulation(Triangulation<dim>::maximum_smoothing),
		fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1),
		dof_handler(triangulation),
		quadrature_formula(degree + 2),
		ellipsoid(&triangulation)
				{}

// Set up dofs

template<int dim>
void StokesProblem<dim>::setup_dofs() {
	A_preconditioner.reset();
	system_matrix.clear();

	dof_handler.distribute_dofs(fe);
	DoFRenumbering::Cuthill_McKee(dof_handler);

	std::vector<unsigned int> block_component(dim + 1, 0);
	block_component[dim] = 1;
	DoFRenumbering::component_wise(dof_handler, block_component);

//========================================Apply Boundary Conditions=====================================
	{
		constraints.clear();
		std::vector<bool> component_maskP(dim + 1, false);
		component_maskP[dim] = true;
		DoFTools::make_hanging_node_constraints(dof_handler, constraints);
		VectorTools::interpolate_boundary_values(dof_handler, 1,
				BoundaryValuesP<dim>(), constraints, component_maskP);
	}
	{
		std::set<unsigned char> no_normal_flux_boundaries;
		no_normal_flux_boundaries.insert(2);
		VectorTools::compute_no_normal_flux_constraints(dof_handler, 0,
				no_normal_flux_boundaries, constraints);
	}

	constraints.close();

	std::vector<unsigned int> dofs_per_block(2);
	DoFTools::count_dofs_per_block(dof_handler, dofs_per_block,
			block_component);
	n_u = dofs_per_block[0];
	n_p = dofs_per_block[1];

	std::cout << "   Number of active cells: " << triangulation.n_active_cells()
			<< std::endl << "   Number of degrees of freedom: "
			<< dof_handler.n_dofs() << " (" << n_u << '+' << n_p << ')'
			<< std::endl;

	{
		BlockCompressedSimpleSparsityPattern csp(2, 2);

		csp.block(0, 0).reinit(n_u, n_u);
		csp.block(1, 0).reinit(n_p, n_u);
		csp.block(0, 1).reinit(n_u, n_p);
		csp.block(1, 1).reinit(n_p, n_p);

		csp.collect_sizes();

		DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
		sparsity_pattern.copy_from(csp);
	}

	system_matrix.reinit(sparsity_pattern);

	solution.reinit(2);
	solution.block(0).reinit(n_u);
	solution.block(1).reinit(n_p);
	solution.collect_sizes();

	system_rhs.reinit(2);
	system_rhs.block(0).reinit(n_u);
	system_rhs.block(1).reinit(n_p);
	system_rhs.collect_sizes();
}

// Viscosity and Shear modulus functions

template<int dim>
std::vector<double> StokesProblem<dim>::flow_law(double r, double z)
{
	double lat = std::atan(z / r);
	double Tsurf = 178.045 * std::sqrt( std::sqrt( std::sin(3.14159265 / 2 - lat) ));
	double Tcmb = Tsurf * std::exp(3.0722e-6 * system_parameters::crust_thickness);
	// Grain Boundary Sliding
	double eta_surf = 0.00792447 * std::exp(5893.67 / Tsurf);
	double eta_cmb = 0.00792447 * std::exp(5893.67 / Tcmb);
	//usually, these are the silicate, cmb ice, and surface ice viscosities
	double eta_kinks[] = {system_parameters::eta_ceiling, eta_cmb, eta_surf};
//	double eta_kinks[] = {1e22, 1e22, 1e22};
	std::vector<double> etas(sizeof(system_parameters::depths) / sizeof(double *) * 2);
	for(unsigned int i=0; i < (sizeof(system_parameters::depths) / sizeof(double *)); i++)
	{
		etas[2*i] = system_parameters::depths[i];
		etas[2*i + 1] = eta_kinks[i];
	}

	return etas;
}

template<int dim>
double StokesProblem<dim>::get_log_local_viscosity(double &r, double &z) {
	double ecc = system_parameters::eq_r / system_parameters::polar_r;
	double Rminusr = system_parameters::eq_r - system_parameters::polar_r;
	double approx_a = std::sqrt(r * r + z * z * ecc * ecc);
	double approx_b = approx_a / ecc;
	double group1 = r * r + z * z - Rminusr * Rminusr;

	if (approx_b
			< (system_parameters::polar_r - system_parameters::crust_thickness)) {
		return system_parameters::eta_ceiling;
	} else {
		double a0 = approx_a;
		double error = 10000;
		// While loop finds the a axis of the "isodepth" ellipse for which the input point is on the surface.
		// An "isodepth" ellipse is defined as one whose axes a,b are related to the global axes A, B by: A-a = B-b
		while (error >= 10) {
			double a02 = a0 * a0;
			double a03 = a0 * a02;
			double a04 = a0 * a03;
			double fofa = a04 - (2 * Rminusr * a03) - (group1 * a02)
					+ (2 * r * r * Rminusr * a0) - (r * r * Rminusr * Rminusr);
			double fprimeofa = 4 * a03 - (6 * Rminusr * a02) - (2 * group1 * a0)
					+ (2 * r * r * Rminusr);
			double deltaa = -fofa / fprimeofa;
			a0 += deltaa;
			error = std::abs(deltaa);
		}
		double local_depth = system_parameters::eq_r - a0;
		if (local_depth < 0)
			local_depth = 0;

		std::vector<double> viscosity_function = flow_law(r, z);

		unsigned int n_visc_kinks = viscosity_function.size() / 2;

		//find the correct interval to do the interpolation in
		int n_minus_one = -1;
		for (unsigned int n = 1; n <= n_visc_kinks; n++) {
			unsigned int ndeep = 2 * n - 2;
			unsigned int nshallow = 2 * n;
			if (local_depth <= viscosity_function[ndeep] && local_depth >= viscosity_function[nshallow])
				n_minus_one = ndeep;
		}

		//find the viscosity interpolation
		if (n_minus_one == -1)
			return system_parameters::eta_ceiling;
		else {
			double visc_exponent =
					(viscosity_function[n_minus_one]
							- local_depth)
							/ (viscosity_function[n_minus_one]
									- viscosity_function[n_minus_one + 2]);
			double visc_base = viscosity_function[n_minus_one + 3]
					/ viscosity_function[n_minus_one + 1];
			// This is the true viscosity given the thermal profile
			double true_eta = viscosity_function[n_minus_one + 1] * std::pow(visc_base, visc_exponent);

			if(true_eta > system_parameters::eta_ceiling)
				return system_parameters::eta_ceiling;
			else
				if(true_eta < system_parameters::eta_floor)
					return system_parameters::eta_floor;
				else
					return true_eta;
		}
	}
}

template<int dim>
double StokesProblem<dim>::get_local_G(double &r, double &z) {
	//generates constant shear moduli in crust and mantle with a sharp discontinuity
	double a = system_parameters::eq_r - system_parameters::crust_thickness;
	double b = system_parameters::polar_r - system_parameters::crust_thickness;
	double expected_z = b * std::sqrt(1 - (r * r / a / a));

	if (z <= expected_z)
		return system_parameters::rock_G;
	else
		return system_parameters::ice_G;
}

// Initialize the eta and G parts of the quadrature_point_history object

template<int dim>
void StokesProblem<dim>::initialize_eta_and_G() {
	FEValues<dim> fe_values(fe, quadrature_formula, update_quadrature_points);

	const unsigned int n_q_points = quadrature_formula.size();

	for (typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(); cell != dof_handler.end(); ++cell) {
		PointHistory<dim> *local_quadrature_points_history =
				reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
		Assert(
				local_quadrature_points_history >= &quadrature_point_history.front(),
				ExcInternalError());
		Assert(
				local_quadrature_points_history < &quadrature_point_history.back(),
				ExcInternalError());
		fe_values.reinit(cell);

		for (unsigned int q = 0; q < n_q_points; ++q) {
			double r_value = fe_values.quadrature_point(q)[0];
			double z_value = fe_values.quadrature_point(q)[1];
			//defines local viscosity
			double local_viscosity = 0;
			local_viscosity = StokesProblem<dim>::get_log_local_viscosity(
						r_value, z_value);

			local_quadrature_points_history[q].first_eta = local_viscosity;
			local_quadrature_points_history[q].new_eta = local_viscosity;

			//defines local shear modulus
			double local_G = 0;
			local_G = StokesProblem<dim>::get_local_G(r_value, z_value);

			local_quadrature_points_history[q].G = local_G;

			//initializes the phi-phi stress
			local_quadrature_points_history[q].old_phiphi_stress = 0;
		}
	}
}

//====================== ASSEMBLE THE SYSTEM ======================

template<int dim>
void StokesProblem<dim>::assemble_system() {
	system_matrix = 0;
	system_rhs = 0;

	FEValues<dim> fe_values(fe, quadrature_formula,
			update_values | update_quadrature_points | update_JxW_values
					| update_gradients);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;

	const unsigned int n_q_points = quadrature_formula.size();

	FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> local_rhs(dofs_per_cell);

	std::vector<unsigned int> local_dof_indices(dofs_per_cell);

	// runs the gravity script function
	const RightHandSide<dim> right_hand_side;

	A_Grav_namespace::AnalyticGravity<dim> * aGrav =
			new A_Grav_namespace::AnalyticGravity<dim>;
	std::vector<double> grav_parameters;
	grav_parameters.push_back(system_parameters::eq_r);
	grav_parameters.push_back(system_parameters::polar_r);
	grav_parameters.push_back(system_parameters::r_core_eq);
	grav_parameters.push_back(system_parameters::r_core_polar);
	grav_parameters.push_back(system_parameters::mantle_rho);
	grav_parameters.push_back(system_parameters::core_rho);
	aGrav->setup_vars(grav_parameters);

	std::vector<Vector<double> > rhs_values(n_q_points,
			Vector<double>(dim + 1));

	const FEValuesExtractors::Vector velocities(0);
	const FEValuesExtractors::Scalar pressure(dim);

	std::vector<SymmetricTensor<2, dim> > phi_grads_u(dofs_per_cell);
	std::vector<double> div_phi_u(dofs_per_cell);
	std::vector<Tensor<1, dim> > phi_u(dofs_per_cell);
	std::vector<double> phi_p(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), first_cell = dof_handler.begin_active(),
			endc = dof_handler.end();

	for (; cell != endc; ++cell) {
		PointHistory<dim> *local_quadrature_points_history =
				reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
		Assert(
				local_quadrature_points_history >= &quadrature_point_history.front(),
				ExcInternalError());
		Assert(
				local_quadrature_points_history < &quadrature_point_history.back(),
				ExcInternalError());

		//initializes the rhs vector to the correct g values
		fe_values.reinit(cell);
		right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
				rhs_values, aGrav);

		std::vector<Vector<double> > new_viscosities(quadrature_formula.size(), Vector<double>(dim + 1));

		if (plastic_iteration != 0)
		{
			// Finds the cell viscosities from the previous plasticity cycle
			BlockVector<double> new_viscosities_block;
			new_viscosities_block.reinit(2);
			new_viscosities_block.block(0).reinit(n_u);
			new_viscosities_block.block(1).reinit(n_p);
			new_viscosities_block.collect_sizes();

			if ((n_u + n_p) != node_viscosities.size())
			{
				std::cout << "Node_viscosities vector is not the same length as the solution block vector.";
				std::exit(1);
			}

			for (unsigned int i = 0; i < node_viscosities.size(); i++)
			{
				if (i < n_u)
					new_viscosities_block.block(0)(i) = node_viscosities(i);
				else
					new_viscosities_block.block(1)(i-n_u) = node_viscosities(i);
			}

//			if (cell == first_cell)
//			{
//				std::cout << "block 0 and 1 sizes: " << n_u << " " << n_p << endl;
//				std::cout << "node_viscosities vector size: " << dof_handler.n_dofs() << endl;
//				std::cout << "number of vertices: " << triangulation.n_vertices() << endl;
//			}

			new_viscosities.resize(quadrature_formula.size());
			fe_values.get_function_values(new_viscosities_block, new_viscosities);
		}

		// Finds vertices where the radius is zero DIM
		bool is_singular = false;
		unsigned int singular_vertex_id = 0;
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
			if (cell->face(f)->center()[0] == 0) {
				is_singular = true;
				singular_vertex_id = f;
			}
		}

		if (is_singular == false || system_parameters::cylindrical == false) {
			local_matrix = 0;
			local_rhs = 0;

			// ===== outputs the local gravity
			std::vector<Point<dim> > quad_points_list(n_q_points);
			quad_points_list = fe_values.get_quadrature_points();

			if (plastic_iteration
					== (system_parameters::max_plastic_iterations - 1)) {
				if (cell != first_cell) {
					std::ofstream fout("gravity_field.txt", std::ios::app);
					fout << quad_points_list[0] << " " << rhs_values[0];
					fout.close();
				} else {
					std::ofstream fout("gravity_field.txt");
					fout << quad_points_list[0] << " " << rhs_values[0];
					fout.close();
				}
			}

			for (unsigned int q = 0; q < n_q_points; ++q) {
				const SymmetricTensor<2, dim> &old_stress =
						local_quadrature_points_history[q].old_stress;
				double &local_old_phiphi_stress =
						local_quadrature_points_history[q].old_phiphi_stress;
				double r_value = fe_values.quadrature_point(q)[0];
				double z_value = fe_values.quadrature_point(q)[1];
				double local_density = 0;
				double expected_core_z = system_parameters::r_core_polar
						* std::sqrt(
								1
										- r_value * r_value
												/ system_parameters::r_core_eq
												/ system_parameters::r_core_eq);

				if ((z_value - expected_core_z) < 1e-10)
					local_density = system_parameters::core_rho;
				else
					local_density = system_parameters::mantle_rho;

				//defines local viscosities
				double local_viscosity = 0;
				if (plastic_iteration == 0)
					local_viscosity = local_quadrature_points_history[q].first_eta;
				else
					local_viscosity = new_viscosities[q][0];

				// Define the local viscoelastic constants
				double local_eta_ve = 2
						/ ((1 / local_viscosity)
								+ (1 / local_quadrature_points_history[q].G
										/ system_parameters::current_time_interval));
				double local_chi_ve = 1
						/ (1
								+ (local_quadrature_points_history[q].G
										* system_parameters::current_time_interval
										/ local_viscosity));

				for (unsigned int k = 0; k < dofs_per_cell; ++k) {
					phi_grads_u[k] = fe_values[velocities].symmetric_gradient(k,
							q);
					div_phi_u[k] = (fe_values[velocities].divergence(k, q));
					phi_u[k] = (fe_values[velocities].value(k, q));
					if (system_parameters::cylindrical == true) {
						div_phi_u[k] *= (r_value);
						div_phi_u[k] += (phi_u[k][0]);
					}
					phi_p[k] = fe_values[pressure].value(k, q);
				}

				for (unsigned int i = 0; i < dofs_per_cell; ++i) {
					for (unsigned int j = 0; j <= i; ++j) {
						if (system_parameters::cylindrical == true) {
							local_matrix(i, j) += (phi_grads_u[i]
									* phi_grads_u[j] * 2 * local_eta_ve
									* r_value
									+ 2 * phi_u[i][0] * phi_u[j][0]
											* local_eta_ve / r_value
									- div_phi_u[i] * phi_p[j]
											* system_parameters::pressure_scale
									- phi_p[i] * div_phi_u[j]
											* system_parameters::pressure_scale
									+ phi_p[i] * phi_p[j] * r_value
											* system_parameters::pressure_scale)
									* fe_values.JxW(q);
						} else {
							local_matrix(i, j) += (phi_grads_u[i]
									* phi_grads_u[j] * 2 * local_eta_ve
									- div_phi_u[i] * phi_p[j]
											* system_parameters::pressure_scale
									- phi_p[i] * div_phi_u[j]
											* system_parameters::pressure_scale
									+ phi_p[i] * phi_p[j]) * fe_values.JxW(q);
						}
					}
					if (system_parameters::cylindrical == true) {
						const unsigned int component_i =
								fe.system_to_component_index(i).first;
						local_rhs(i) += (fe_values.shape_value(i, q)
								* rhs_values[q](component_i) * r_value
								* local_density
								- local_chi_ve * phi_grads_u[i] * old_stress
										* r_value
								- local_chi_ve * phi_u[i][0]
										* local_old_phiphi_stress)
								* fe_values.JxW(q);
					} else {
						const unsigned int component_i =
								fe.system_to_component_index(i).first;
						local_rhs(i) += fe_values.shape_value(i, q)
								* rhs_values[q](component_i) * fe_values.JxW(q)
								* local_density;
					}
				}
			}
		} // end of non-singular
		else {
			local_matrix = 0;
			local_rhs = 0;

			// ===== outputs the local gravity
			std::vector<Point<dim> > quad_points_list(n_q_points);
			quad_points_list = fe_values.get_quadrature_points();

			for (unsigned int q = 0; q < n_q_points; ++q) {
				const SymmetricTensor<2, dim> &old_stress =
						local_quadrature_points_history[q].old_stress;
				double &local_old_phiphi_stress =
						local_quadrature_points_history[q].old_phiphi_stress;
				double r_value = fe_values.quadrature_point(q)[0];
				double z_value = fe_values.quadrature_point(q)[1];
				double local_density = 0;
				double expected_core_z = system_parameters::r_core_polar
						* std::sqrt(
								1
										- r_value * r_value
												/ system_parameters::r_core_eq
												/ system_parameters::r_core_eq);

//					std::ofstream fout("all_quad_points.txt",std::ios::app);
//					fout <<  r_value << " " << z_value << "\n";
//					fout.close();

				if ((z_value - expected_core_z) < 1e-10)
					local_density = system_parameters::core_rho;
				else
					local_density = system_parameters::mantle_rho;

				//defines local viscosities
				double local_viscosity = 0;
				if (plastic_iteration == 0) {
					local_viscosity =
							local_quadrature_points_history[q].first_eta;
				} else
					local_viscosity = new_viscosities[q][0];

				// Define the local viscoelastic constants
				double local_eta_ve = 2
						/ ((1 / local_viscosity)
								+ (1 / local_quadrature_points_history[q].G
										/ system_parameters::current_time_interval));
				double local_chi_ve = 1
						/ (1
								+ (local_quadrature_points_history[q].G
										* system_parameters::current_time_interval
										/ local_viscosity));

				for (unsigned int k = 0; k < dofs_per_cell; ++k) {
					phi_grads_u[k] = fe_values[velocities].symmetric_gradient(k,
							q);
					div_phi_u[k] = (fe_values[velocities].divergence(k, q));
					phi_u[k] = (fe_values[velocities].value(k, q));
					if (system_parameters::cylindrical == true) {
						div_phi_u[k] *= (r_value);
						div_phi_u[k] += (phi_u[k][0]);
					}
					phi_p[k] = fe_values[pressure].value(k, q);
				}

				for (unsigned int i = 0; i < dofs_per_cell; ++i) {
					for (unsigned int j = 0; j <= i; ++j) {
						if (system_parameters::cylindrical == true) {
							local_matrix(i, j) += (phi_grads_u[i]
									* phi_grads_u[j] * 2 * local_eta_ve
									* r_value
									+ 2 * phi_u[i][0] * phi_u[j][0]
											* local_eta_ve / r_value
									- div_phi_u[i] * phi_p[j]
											* system_parameters::pressure_scale
									- phi_p[i] * div_phi_u[j]
											* system_parameters::pressure_scale
									+ phi_p[i] * phi_p[j] * r_value
											* system_parameters::pressure_scale)
									* fe_values.JxW(q);
						} else {
							local_matrix(i, j) += (phi_grads_u[i]
									* phi_grads_u[j] * 2 * local_eta_ve
									- div_phi_u[i] * phi_p[j]
											* system_parameters::pressure_scale
									- phi_p[i] * div_phi_u[j]
											* system_parameters::pressure_scale
									+ phi_p[i] * phi_p[j]) * fe_values.JxW(q);
						}
					}
					if (system_parameters::cylindrical == true) {
						const unsigned int component_i =
								fe.system_to_component_index(i).first;
						local_rhs(i) += (fe_values.shape_value(i, q)
								* rhs_values[q](component_i) * r_value
								* local_density
								- local_chi_ve * phi_grads_u[i] * old_stress
										* r_value
								- local_chi_ve * phi_u[i][0]
										* local_old_phiphi_stress)
								* fe_values.JxW(q);
					} else {
						const unsigned int component_i =
								fe.system_to_component_index(i).first;
						local_rhs(i) += fe_values.shape_value(i, q)
								* rhs_values[q](component_i) * fe_values.JxW(q)
								* local_density;
					}
				}
			}
		} // end of singular

		for (unsigned int i = 0; i < dofs_per_cell; ++i)
			for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
				local_matrix(i, j) = local_matrix(j, i);

		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(local_matrix, local_rhs,
				local_dof_indices, system_matrix, system_rhs);
	}

	std::cout << "   Computing preconditioner..." << std::endl << std::flush;

	A_preconditioner = std_cxx1x::shared_ptr<
			typename InnerPreconditioner<dim>::type>(
			new typename InnerPreconditioner<dim>::type());
	A_preconditioner->initialize(system_matrix.block(0, 0),
			typename InnerPreconditioner<dim>::type::AdditionalData());

	delete aGrav;

}

//====================== SOLVER ======================

template<int dim>
void StokesProblem<dim>::solve() {
	const InverseMatrix<SparseMatrix<double>,
			typename InnerPreconditioner<dim>::type> A_inverse(
			system_matrix.block(0, 0), *A_preconditioner);
	Vector<double> tmp(solution.block(0).size());

	{
		Vector<double> schur_rhs(solution.block(1).size());
		A_inverse.vmult(tmp, system_rhs.block(0));
		system_matrix.block(1, 0).vmult(schur_rhs, tmp);
		schur_rhs -= system_rhs.block(1);

		SchurComplement<typename InnerPreconditioner<dim>::type> schur_complement(
				system_matrix, A_inverse);

		int n_iterations = system_parameters::iteration_coefficient
				* solution.block(1).size();
		double tolerance_goal = system_parameters::tolerance_coefficient
				* schur_rhs.l2_norm();

		SolverControl solver_control(n_iterations, tolerance_goal);
		SolverCG<> cg(solver_control);

		std::cout << "\nMax iterations and tolerance are:  " << n_iterations
				<< " and " << tolerance_goal << std::endl;

		SparseILU<double> preconditioner;
		preconditioner.initialize(system_matrix.block(1, 1),
				SparseILU<double>::AdditionalData());

		InverseMatrix<SparseMatrix<double>, SparseILU<double> > m_inverse(
				system_matrix.block(1, 1), preconditioner);

		cg.solve(schur_complement, solution.block(1), schur_rhs, m_inverse);

		constraints.distribute(solution);


		std::cout << "  " << solver_control.last_step()
				<< " outer CG Schur complement iterations for pressure"
				<< std::endl;
	}

	{
		system_matrix.block(0, 1).vmult(tmp, solution.block(1));
		tmp *= -1;
		tmp += system_rhs.block(0);

		A_inverse.vmult(solution.block(0), tmp);
		constraints.distribute(solution);
		solution.block(1) *= (system_parameters::pressure_scale);
	}
}

//====================== OUTPUT RESULTS ======================
template<int dim>
void StokesProblem<dim>::output_results() const {
	std::vector < std::string > solution_names(dim, "velocity");
	solution_names.push_back("pressure");

	std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(
			dim, DataComponentInterpretation::component_is_part_of_vector);
	data_component_interpretation.push_back(
			DataComponentInterpretation::component_is_scalar);

	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names,
			DataOut<dim>::type_dof_data, data_component_interpretation);
	data_out.build_patches();

	std::ostringstream filename;
	if (system_parameters::present_timestep < system_parameters::initial_elastic_iterations)
	{
		filename << "time"
						<< Utilities::int_to_string(system_parameters::present_timestep, 2)
						<< "_elastic_displacements" << ".txt";
	}
	else
	{
		filename << "time"
				<< Utilities::int_to_string(system_parameters::present_timestep, 2)
				<< "_flow" << Utilities::int_to_string(plastic_iteration, 2) << ".txt";
	}


	std::ofstream output(filename.str().c_str());
	data_out.write_gnuplot(output);
}

//====================== FIND AND WRITE TO FILE THE STRESS TENSOR; IMPLEMENT PLASTICITY ======================

template<int dim>
void StokesProblem<dim>::solution_stesses() {
	//note most of this section only works with dim=2

	//name the output text files
	std::ostringstream stress_output;
	stress_output << "time"
			<< Utilities::int_to_string(system_parameters::present_timestep, 2)
			<< "_principalstresses" << Utilities::int_to_string(plastic_iteration, 2)
			<< ".txt";
	std::ofstream fout_snew(stress_output.str().c_str());
	fout_snew.close();

	std::ostringstream stresstensor_output;
	stresstensor_output << "time"
			<< Utilities::int_to_string(system_parameters::present_timestep, 2)
			<< "_stresstensor" << Utilities::int_to_string(plastic_iteration, 2)
			<< ".txt";
	std::ofstream fout_sfull(stresstensor_output.str().c_str());
	fout_sfull.close();

	std::ostringstream failed_cells_output;
	failed_cells_output << "time"
			<< Utilities::int_to_string(system_parameters::present_timestep, 2)
			<< "_failurelocations" << Utilities::int_to_string(plastic_iteration, 2)
			<< ".txt";
	std::ofstream fout_failed_cells(failed_cells_output.str().c_str());
	fout_failed_cells.close();

	std::ostringstream plastic_eta_output;
	plastic_eta_output << "time"
			<< Utilities::int_to_string(system_parameters::present_timestep, 2)
			<< "_viscositiesreg" << Utilities::int_to_string(plastic_iteration, 2)
			<< ".txt";
	std::ofstream fout_vrnew(plastic_eta_output.str().c_str());
	fout_vrnew.close();

	std::ostringstream initial_eta_output;
	if (plastic_iteration == 0)
	{
		initial_eta_output << "time"
			<< Utilities::int_to_string(system_parameters::present_timestep, 2)
			<< "_baseviscosities.txt";
		std::ofstream fout_baseeta(initial_eta_output.str().c_str());
		fout_baseeta.close();
	}

	std::cout << "Running stress calculations for plasticity iteration "
			<< plastic_iteration << "...\n\n";

	//This makes the set of points at which the stress tensor is calculated
	std::vector<Point<dim> > points_list(0);
	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	//This loop gets the gradients of the velocity field and saves it in the tensor_gradient_? objects DIM
	for (; cell != endc; ++cell) {
		points_list.push_back(cell->center());
	}
	// Make the FEValues object to evaluate values and derivatives at quadrature points
	FEValues<dim> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points | update_JxW_values);

	// Make the object that will hold the velocities and velocity gradients at the quadrature points
	std::vector < std::vector<Tensor<1, dim> >> velocity_grads(quadrature_formula.size(),
			std::vector < Tensor<1, dim> > (dim + 1));
	std::vector<Vector<double> > velocities(quadrature_formula.size(),
			Vector<double>(dim + 1));

	// Write the solution flow velocity and derivative for each cell
	std::vector<Vector<double> > vector_values(0);
	std::vector < std::vector<Tensor<1, dim> > > gradient_values(0);
	for (typename DoFHandler<dim>::active_cell_iterator cell =
		dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
	{
		fe_values.reinit(cell);
		fe_values.get_function_gradients(solution, velocity_grads);
		fe_values.get_function_values(solution, velocities);
		Vector<double> current_cell_velocity(dim+1);
		std::vector<Tensor<1, dim>> current_cell_grads(dim+1);
		double cell_area = 0;

		for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
		{
			cell_area += fe_values.JxW(q);
			velocities[q] *= fe_values.JxW(q);
			current_cell_velocity += velocities[q];
			for (unsigned int i = 0; i < (dim+1); i++)
			{
				velocity_grads[q][i] *= fe_values.JxW(q);
				current_cell_grads[i] += velocity_grads[q][i];
			}
		}
		current_cell_velocity /= cell_area;
		for (unsigned int i = 0; i < (dim+1); i++)
			current_cell_grads[i] /= cell_area;

		vector_values.push_back(current_cell_velocity);
		gradient_values.push_back(current_cell_grads);
	}

	//tracks where failure occurred
	std::vector<int> fail_ID;
	unsigned int total_fails = 0;

	//loop across all the cells to find and adjust eta of failing cells
	if (cell_viscosities.size() == 0)
		cell_viscosities.reinit(triangulation.n_active_cells());
	for (unsigned int i = 0; i < triangulation.n_active_cells(); i++) {
		double current_cell_viscosity = 0;

		//fill viscosities vector
		if (plastic_iteration == 0) {
			double local_viscosity =
					StokesProblem<dim>::get_log_local_viscosity(points_list[i][0], points_list[i][1]);
			current_cell_viscosity = local_viscosity;
		} else
		{
			current_cell_viscosity = cell_viscosities.operator()(i);
		}

		//find local pressure
		double cell_p = vector_values[i].operator()(2);
		//find stresses tensor
		//makes non-diagonalized local matrix A
		double sigma13 = 0.5
				* (gradient_values[i][0][1] + gradient_values[i][1][0]);
		mat A;
		A << gradient_values[i][0][0] << 0 << sigma13 << endr
		  << 0 << vector_values[i].operator()(0) / points_list[i].operator()(0)<< 0 << endr
		  << sigma13 << 0 << gradient_values[i][1][1] << endr;
		vec P;
		P << cell_p << cell_p << cell_p;
		mat Pmat = diagmat(P);
		mat B;
		B = (2 * current_cell_viscosity * A) - Pmat;

		//finds principal stresses
		vec eigval;
		mat eigvec;
		eig_sym(eigval, eigvec, B);
		double sigma1 = -min(eigval);
		double sigma3 = -max(eigval);

		// Writes text files for principal stresses, full stress tensor, base viscosities
		std::ofstream fout_snew(stress_output.str().c_str(), std::ios::app);
		fout_snew << " " << sigma1 << " " << sigma3 << "\n";
		fout_snew.close();

		std::ofstream fout_sfull(stresstensor_output.str().c_str(), std::ios::app);
		fout_sfull << A(0,0) << " " << A(1,1) << " " << A(2,2) << " " << A(0,2) << "\n";
		fout_sfull.close();

		if (plastic_iteration == 0)
		{
			std::ofstream fout_baseeta(initial_eta_output.str().c_str(), std::ios::app);
			fout_baseeta << points_list[i]<< " " << current_cell_viscosity << "\n";
			fout_baseeta.close();
		}

		// Finds adjusted effective viscosity
		double cell_effective_viscosity = 0;
		if (system_parameters::plasticity_on)
		{
			if (sigma1 >= 5 * sigma3) // this guarantees that viscosities only go down, never up
				{
					double reductionfactor = 1;
					if (sigma3 < 0)
						reductionfactor = 100;
					else
						reductionfactor = 1.9 * sigma1 / 5 / sigma3;
						cell_effective_viscosity = current_cell_viscosity / reductionfactor;
						fail_ID.push_back(1);
						total_fails++;

						std::ofstream fout_failed_cells(failed_cells_output.str().c_str(), std::ios::app);
						fout_failed_cells << points_list[i] << "\n";
						fout_failed_cells.close();
				} else {
					cell_effective_viscosity = current_cell_viscosity;
					fail_ID.push_back(0);
				}
		}
		else
		{
			cell_effective_viscosity = current_cell_viscosity;
			fail_ID.push_back(0);
		}


		// Prevents viscosities from dropping below the floor necessary for numerical stability
		if (current_cell_viscosity < system_parameters::eta_floor)
		{
			current_cell_viscosity = system_parameters::eta_floor;
			cell_viscosities.operator()(i) = cell_effective_viscosity;
		}
		else
			cell_viscosities.operator()(i) = cell_effective_viscosity;
	}

	std::cout << "   Number of failing cells: " << total_fails << "\n";
	if (total_fails <= 20)
		system_parameters::continue_plastic_iterations = false;

	node_viscosities.reinit(dof_handler.n_dofs());

	//maps effective viscosity field from cells to vertices to make the FEFieldFunction object for the smoothing
	DoFTools::distribute_cell_to_dof_vector(dof_handler, cell_viscosities,
			node_viscosities);
	Functions::FEFieldFunction<dim> new_viscosities1(dof_handler,
			node_viscosities);

	//Smoothes the new effective viscosity map
	for (unsigned int n = 0; n < triangulation.n_active_cells(); n++) {
		if (fail_ID[n] == 1) {
			double smooth_interval = 5000; //does work
			double n_radii = std::floor(
					system_parameters::smoothing_radius / smooth_interval);
			unsigned int n_angles = 8;
			std::vector<Point<dim> > smooth_points(0);
			smooth_points.push_back(points_list[n]);
			double x_sm = points_list[n][0];
			double y_sm = points_list[n][1];

			//get all the smoothing points for a certain cell
			for (unsigned int i = 1; i <= n_radii; i++) {
				double current_r = smooth_interval * i;
				double angles[16] = { 1, 0, .707, .707, 0, 1, -.707, .707, -1,
						0, -.707, -.707, 0, -1, .707, -.707 };
				for (unsigned int j = 0; j < n_angles; j++) {
					Point<dim> sm_point(x_sm + current_r * angles[2 * j],
							y_sm + current_r * angles[2 * j + 1]);
					if (sm_point[0] >= 0 && sm_point[1] >= 0) {
						double a = system_parameters::eq_r - 10000; //make sure none of the spots are inside the ellipse but outside of cells
						double b = system_parameters::polar_r - 10000;
						double expected_y_sq = b * b
								* (1 - (sm_point[0] * sm_point[0] / a / a));
						if ((sm_point[1] * sm_point[1]) <= expected_y_sq)
							smooth_points.push_back(sm_point);
					}
				}
			}

			std::vector<double> viscosity_at_smooth_points( smooth_points.size());
			new_viscosities1.value_list(smooth_points, viscosity_at_smooth_points);
			double new_cell_viscosity = 0;
			double sum_ln_visc = 0;
			int n_viscosities = 0;

			for (unsigned int k = 0; k < viscosity_at_smooth_points.size(); k++) {
				if (viscosity_at_smooth_points[k] > 0) //some of the viscosities very near the surface are negative for some reason
						{
					sum_ln_visc += std::log(viscosity_at_smooth_points[k]);
					n_viscosities++;
				}
			}
			new_cell_viscosity = std::exp(sum_ln_visc / n_viscosities);
			cell_viscosities.operator()(n) = new_cell_viscosity;
		}
	}

	//Writes the plasticity-corrected, smoothed effective viscosities to the node_viscosities vector
	DoFTools::distribute_cell_to_dof_vector(dof_handler, cell_viscosities,
			node_viscosities);

	for (unsigned int i = 0; i < triangulation.n_active_cells(); i++) {
		std::ofstream fout_vrnew(plastic_eta_output.str().c_str(), std::ios::app);
		fout_vrnew << " " << cell_viscosities[i] << "\n";
		fout_vrnew.close();
	}
}

//====================== SAVE STRESS TENSOR AT QUADRATURE POINTS ======================

template<int dim>
void StokesProblem<dim>::update_quadrature_point_history() {
	std::cout << "   Updating stress field...";

	FEValues<dim> fe_values(fe, quadrature_formula,
			update_values | update_gradients | update_quadrature_points);

	// Make the object that will hold the velocity gradients
	std::vector < std::vector<Tensor<1, dim> >> velocity_grads(quadrature_formula.size(),
			std::vector < Tensor<1, dim> > (dim + 1));
	std::vector<Vector<double> > velocities(quadrature_formula.size(),
			Vector<double>(dim + 1));

//		std::cout << std::vector<Tensor<1,dim> >(dim) << "\n";

	for (typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(); cell != dof_handler.end(); ++cell) {
		PointHistory<dim> *local_quadrature_points_history =
				reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
		Assert(
				local_quadrature_points_history >= &quadrature_point_history.front(),
				ExcInternalError());
		Assert(
				local_quadrature_points_history < &quadrature_point_history.back(),
				ExcInternalError());

		fe_values.reinit(cell);
		fe_values.get_function_gradients(solution, velocity_grads);
		fe_values.get_function_values(solution, velocities);

		for (unsigned int q = 0; q < quadrature_formula.size(); ++q) {
			// Define the local viscoelastic constants
			double local_eta_ve = 2
					/ ((1 / local_quadrature_points_history[q].new_eta)
							+ (1 / local_quadrature_points_history[q].G
									/ system_parameters::current_time_interval));
			double local_chi_ve =
					1
							/ (1
									+ (local_quadrature_points_history[q].G
											* system_parameters::current_time_interval
											/ local_quadrature_points_history[q].new_eta));

			// Compute new stress at each quadrature point
			SymmetricTensor<2, dim> new_stress;
			for (unsigned int i = 0; i < dim; ++i)
				new_stress[i][i] =
						local_eta_ve * velocity_grads[q][i][i]
								+ local_chi_ve
										* local_quadrature_points_history[q].old_stress[i][i];

			for (unsigned int i = 0; i < dim; ++i)
				for (unsigned int j = i + 1; j < dim; ++j)
					new_stress[i][j] =
							local_eta_ve
									* (velocity_grads[q][i][j]
											+ velocity_grads[q][j][i]) / 2
									+ local_chi_ve
											* local_quadrature_points_history[q].old_stress[i][j];

			// Rotate new stress
			AuxFunctions<dim> rotation_object;
			const Tensor<2, dim> rotation = rotation_object.get_rotation_matrix(
					velocity_grads[q]);
			const SymmetricTensor<2, dim> rotated_new_stress = symmetrize(
					transpose(rotation)
							* static_cast<Tensor<2, dim> >(new_stress)
							* rotation);
			local_quadrature_points_history[q].old_stress = rotated_new_stress;

//				std::cout << velocities[q] <<"\n";
			// For axisymmetric case, make the phi-phi element of stress tensor
			local_quadrature_points_history[q].old_phiphi_stress =
					(2 * local_eta_ve * velocities[q](0)
							/ fe_values.quadrature_point(q)[0]
							+ local_chi_ve
									* local_quadrature_points_history[q].old_phiphi_stress);
		}

	}
}

//====================== REDEFINE THE TIME INTERVAL FOR THE VISCOUS STEPS ======================
template<int dim>
void StokesProblem<dim>::update_time_interval()
{
	double move_goal_per_step = system_parameters::initial_disp_target -
			((system_parameters::initial_disp_target - system_parameters::final_disp_target) /
					system_parameters::total_viscous_steps *
					(system_parameters::present_timestep - system_parameters::initial_elastic_iterations));
	double max_velocity = 0;
	for(unsigned int i=0; i<solution.n_blocks()-1; i++)
		for(unsigned int j=0; j<solution.block(i).size(); j++)
			if(std::abs(solution.block(i)(j)) > max_velocity)
				max_velocity = std::abs(solution.block(i)(j));
	// NOTE: It is possible for this time interval to be very different from that used in the FE calculation.
	system_parameters::current_time_interval = move_goal_per_step / max_velocity;
	std::cout << "\n   Viscous time for moving mesh: " << system_parameters::current_time_interval << " s";
}

//====================== MOVE MESH ======================

template<int dim>
void StokesProblem<dim>::move_mesh() {
	std::cout << "\n" << "   Moving mesh...";


	std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
	for (typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
		for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
			if (vertex_touched[cell->vertex_index(v)] == false) {
				vertex_touched[cell->vertex_index(v)] = true;

				Point<dim> vertex_displacement;
				for (unsigned int d = 0; d < dim; ++d)
					vertex_displacement[d] = solution(
							cell->vertex_dof_index(v, d));

				cell->vertex(v) += vertex_displacement
						* system_parameters::current_time_interval;
			}

	std::vector<double> ellipse_axes(0);
	// compute fit to boundary 1
	ellipsoid.compute_fit(ellipse_axes, 1);

	std::cout << endl;
	std::cout << "New dimensions for best-fit ellipse: ";
	for(unsigned int j = 0; j < ellipse_axes.size(); j++)
		std::cout << ellipse_axes[j] << " ";

	system_parameters::eq_r = ellipse_axes[0];
	system_parameters::polar_r = ellipse_axes[1];
	system_parameters::r_core_eq = system_parameters::eq_r - system_parameters::crust_thickness;
	system_parameters::r_core_polar = system_parameters::polar_r - system_parameters::crust_thickness;

}

//====================== WRITE VERTICES TO FILE ======================

template<int dim>
void StokesProblem<dim>::write_vertices() {
	std::ostringstream vertices_output;
	vertices_output << "time"
		<< Utilities::int_to_string(system_parameters::present_timestep, 2)
		<< "_surface.txt";
	std::ofstream fout_final_vertices(vertices_output.str().c_str());
	fout_final_vertices.close();

	// Figure out if the vertex is on the boundary of the domain
	std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
	for (typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(); cell != triangulation.end(); ++cell)
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
		{
			unsigned char boundary_ids = cell->face(f)->boundary_indicator();
			if(boundary_ids == 1)
				{
				for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
					if (vertex_touched[cell->face(f)->vertex_index(v)] == false)
						{
							vertex_touched[cell->face(f)->vertex_index(v)] = true;
							std::ofstream fout_final_vertices(vertices_output.str().c_str(), std::ios::app);
							fout_final_vertices << cell->face(f)->vertex(v) << "\n";
							fout_final_vertices.close();
						}
				}
		}
}

//====================== SETUP INITIAL MESH ======================

template<int dim>
void StokesProblem<dim>::setup_initial_mesh() {
	GridIn<dim> grid_in;
	grid_in.attach_triangulation(triangulation);
	std::ifstream mesh_stream(system_parameters::mesh_filename,
			std::ifstream::in);
	grid_in.read_ucd(mesh_stream);
	std::ofstream out_eps ("initial_mesh.eps");
	GridOut grid_out;
	grid_out.write_eps (triangulation, out_eps);

//boundary indicator 1 is free surface, indicator 2 is inside
	double zero_tolerance = 1e-3;
	for (typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(); cell != triangulation.end(); ++cell)
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
			if (cell->face(f)->at_boundary()) {
				cell->face(f)->set_all_boundary_indicators(0);
				if (cell->face(f)->center()[0] > zero_tolerance) {
					if (cell->face(f)->center()[1] > zero_tolerance) {
						cell->face(f)->set_all_boundary_indicators(1);
						if (cell->face(f)->center()[0] <= 1) {
							cell->face(f)->set_all_boundary_indicators(3);
						}
					}
				}

			}
			if (cell->face(f)->center()[0] < zero_tolerance)
				cell->face(f)->set_all_boundary_indicators(2);
			if (cell->face(f)->center()[1] < zero_tolerance)
				cell->face(f)->set_all_boundary_indicators(2);
		}

	triangulation.refine_global(system_parameters::global_refinement);

//refines region near r=0
	if (system_parameters::small_r_refinement != 0) {
		for (unsigned int step = 0;
				step < system_parameters::small_r_refinement; ++step) {
			//		std::cout << "iteration " << step + 1 << "\n";
			typename dealii::Triangulation<dim>::active_cell_iterator cell =
					triangulation.begin_active(), endc = triangulation.end();

			for (; cell != endc; ++cell)
				for (unsigned int v = 0;
						v < GeometryInfo<dim>::vertices_per_cell; ++v) {
					Point<dim> current_vertex = cell->vertex(v);

					const double x_coord = current_vertex.operator()(0);

					if (std::fabs(x_coord) < 1e-10) {
						cell->set_refine_flag();
						break;
					}

				}
			triangulation.execute_coarsening_and_refinement();
		}
	}

//refines crustal region
	if (system_parameters::crustal_refinement != 0) {
		double a = system_parameters::eq_r - system_parameters::crust_thickness;
		double b = system_parameters::polar_r
				- system_parameters::crust_thickness;

		for (unsigned int step = 0;
				step < system_parameters::crustal_refinement; ++step) {
			//		std::cout << "iteration " << step + 1 << "\n";
			typename dealii::Triangulation<dim>::active_cell_iterator cell =
					triangulation.begin_active(), endc = triangulation.end();
			for (; cell != endc; ++cell)
				for (unsigned int v = 0;
						v < GeometryInfo<dim>::vertices_per_cell; ++v) {
					Point<dim> current_vertex = cell->vertex(v);

					const double x_coord = current_vertex.operator()(0);
					const double y_coord = current_vertex.operator()(1);
					double expected_z = -1;

					if ((x_coord - a) < -1e-10)
						expected_z = b
								* std::sqrt(1 - (x_coord * x_coord / a / a));

					if (y_coord >= expected_z) {
						cell->set_refine_flag();
						break;
					}
				}
			triangulation.execute_coarsening_and_refinement();
		}
	}

	//refines surface region
	if (system_parameters::surface_refinement != 0) {
		for (unsigned int step = 0;
				step < system_parameters::surface_refinement; ++step) {
			//		std::cout << "iteration " << step + 1 << "\n";
			typename dealii::Triangulation<dim>::active_cell_iterator cell =
					triangulation.begin_active(), endc = triangulation.end();
			for (; cell != endc; ++cell)
				for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell;
						++f) {
					if (cell->face(f)->at_boundary()) {
						if (cell->face(f)->center()[0] != 0) {
							if (cell->face(f)->center()[1] != 0)
								if (cell->face(f)->center()[0] > 1) {
									cell->set_refine_flag();
									break;
								}
						}
					}
				}
			triangulation.execute_coarsening_and_refinement();
		}
	}

	std::vector<double> ellipse_axes(0);
	// compute fit to boundary = 1
	ellipsoid.compute_fit(ellipse_axes, 1);

	std::cout << "The initial best-fit ellipse has dimensions: ";
	for(unsigned int j = 0; j < ellipse_axes.size(); j++)
		std::cout << ellipse_axes[j] << " ";
	std::cout << endl;

	system_parameters::eq_r = ellipse_axes[0];
	system_parameters::polar_r = ellipse_axes[1];
	system_parameters::r_core_eq = system_parameters::eq_r - system_parameters::crust_thickness;
	system_parameters::r_core_polar = system_parameters::polar_r - system_parameters::crust_thickness;

	write_vertices();
}

//====================== REFINE MESH ======================

template<int dim>
void StokesProblem<dim>::refine_mesh() {
	Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

	std::vector<bool> component_mask(dim + 1, false);
	component_mask[dim] = true;
	KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(degree + 1),
			typename FunctionMap<dim>::type(), solution,
			estimated_error_per_cell, component_mask);

	GridRefinement::refine_and_coarsen_fixed_number(triangulation,
			estimated_error_per_cell, 0.3, 0.0);
	triangulation.execute_coarsening_and_refinement();
}

//====================== SET UP THE DATA STRUCTURES TO REMEMBER STRESS FIELD ======================
template<int dim>
void StokesProblem<dim>::setup_quadrature_point_history() {
	unsigned int our_cells = 0;
	for (typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(); cell != triangulation.end(); ++cell)
		++our_cells;

	triangulation.clear_user_data();

	{
		std::vector<PointHistory<dim> > tmp;
		tmp.swap(quadrature_point_history);
	}
	quadrature_point_history.resize(our_cells * quadrature_formula.size());

	unsigned int history_index = 0;
	for (typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(); cell != triangulation.end(); ++cell) {
		cell->set_user_pointer(&quadrature_point_history[history_index]);
		history_index += quadrature_formula.size();
	}

	Assert(history_index == quadrature_point_history.size(), ExcInternalError());
}

//====================== DOES ELASTIC STEPS ======================
template<int dim>
void StokesProblem<dim>::do_elastic_steps()
{
	plastic_iteration = 0;
	unsigned int elastic_iteration = 0;

	// Writes files with the physical time and the highest number of plasticity iterations
	std::ostringstream times_filename;
	times_filename << "physical_times.txt";
	std::ofstream fout_times(times_filename.str().c_str());
	fout_times.close();

	while (elastic_iteration < system_parameters::initial_elastic_iterations)
	{

		std::ofstream fout_times(times_filename.str().c_str(), std::ios::app);
		fout_times << system_parameters::present_timestep << " "
									<< system_parameters::present_time << " 0" << "\n";
		fout_times.close();
		std::cout << "\n\nElastic iteration " << elastic_iteration
							<< "\n";
		setup_dofs();
		write_vertices();

		if (system_parameters::present_timestep == 0)
			initialize_eta_and_G();

		system_parameters::current_time_interval =
				system_parameters::elastic_time; //This is the time interval needed in assembling the problem

		std::cout << "   Assembling..." << std::endl << std::flush;
		assemble_system();

		std::cout << "   Solving..." << std::flush;
		solve();

		output_results();
		update_quadrature_point_history();

	//				std::cout << std::endl << "\a";
		elastic_iteration++;
		system_parameters::present_timestep++;
		system_parameters::present_time = system_parameters::present_time + system_parameters::current_time_interval;
		move_mesh();
	}
}

//====================== DO A SINGLE VISCOPLASTIC TIMESTEP ======================
template<int dim>
void StokesProblem<dim>::do_flow_step() {
	plastic_iteration = 0;
	while (plastic_iteration < system_parameters::max_plastic_iterations) {
		if (system_parameters::continue_plastic_iterations == true) {
			std::cout << "Plasticity iteration " << plastic_iteration << "\n";
			setup_dofs();
			write_vertices();

			std::cout << "   Assembling..." << std::endl << std::flush;
			assemble_system();

			std::cout << "   Solving..." << std::flush;
			solve();

			output_results();
			solution_stesses();

			if (system_parameters::continue_plastic_iterations == false) {
				// Writes the current timestep, physical time, and final plasticity_iteration
				std::ostringstream times_filename;
				times_filename << "physical_times.txt";
				std::ofstream fout_times(times_filename.str().c_str(), std::ios::app);
				fout_times << system_parameters::present_timestep << " "
							<< system_parameters::present_time << " " <<  plastic_iteration << "\n";
				fout_times.close();
				break;
			}
			plastic_iteration++;
			std::cout << "0 " << system_parameters::continue_plastic_iterations << endl;
//			std::cout << std::endl << "\a";
		}
	}
}

//====================== RUN ======================

template<int dim>
void StokesProblem<dim>::run() {

	// Sets up mesh and applies elastic displacement
	setup_initial_mesh();
	setup_quadrature_point_history();
	do_elastic_steps();

	// Computes viscous timesteps
	system_parameters::current_time_interval = system_parameters::viscous_time;
	unsigned int VEPstep = 0;
	while (system_parameters::present_timestep
			< (system_parameters::initial_elastic_iterations
					+ system_parameters::total_viscous_steps)) {
		if (system_parameters::continue_plastic_iterations == false)
			system_parameters::continue_plastic_iterations = true;
		std::cout << "\n\nViscoelastoplastic iteration " << VEPstep << "\n\n";
		if (system_parameters::present_timestep == 0)
			initialize_eta_and_G();

		// Computes plasticity
		do_flow_step();
		update_quadrature_point_history();
		update_time_interval();
		move_mesh();
		system_parameters::present_timestep++;
		system_parameters::present_time = system_parameters::present_time + system_parameters::current_time_interval;
		VEPstep++;
	}

	// Write the moved vertices time for the last viscous step
	write_vertices();
	std::ostringstream times_filename;
	times_filename << "physical_times.txt";
	std::ofstream fout_times(times_filename.str().c_str(), std::ios::app);
	fout_times << system_parameters::present_timestep << " "
					<< system_parameters::present_time << " 0\n";
	fout_times.close();
}

}

//====================== MAIN ======================

int main() {
	try {
		using namespace dealii;
		using namespace Step22;

		std::clock_t t1;
		std::clock_t t2;
		t1 = std::clock();

		deallog.depth_console(0);

		StokesProblem<2> flow_problem(1);
		flow_problem.run();

		t2 = std::clock();
		float diff (((float)t2 - (float)t1) / (float)CLOCKS_PER_SEC);
		std::cout  << "\n Program run in: " << diff << " seconds" << endl;
	} catch (std::exception &exc) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Exception on processing: " << std::endl << exc.what()
				<< std::endl << "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;

		return 1;
	} catch (...) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl << "Aborting!"
				<< std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	}

	return 0;
}
