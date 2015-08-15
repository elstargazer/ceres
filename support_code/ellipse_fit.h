/*
 * fe_test.cpp
 *
 *  Created on: Jul 24, 2015
 *      Author: antonermakov
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>
#include <sstream>

using namespace dealii;

template <int dim>
class FE_test
{
public:
	FE_test ();
  void run ();
  void output_results();
  void cp_tr(const Triangulation<dim,dim> &tr_in);
  void fit_ellipsoid(std::vector<double> &ell, unsigned int step, bool full_mesh);

private:

  void                 generate_ellipsoid(double* ell, int nrefine);
  Triangulation<dim,dim>   triangulation;
  FE_Q<dim>            fe;
  DoFHandler<dim>      dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;
};



template <int dim>
void FE_test<dim>::fit_ellipsoid(std::vector<double> &ell, unsigned int step, bool full_mesh)
{

	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
	typename Triangulation<dim>::active_cell_iterator endc = triangulation.end();

	FullMatrix<double> A(triangulation.n_vertices(),dim);
	Vector<double>     x(dim);
    Vector<double>     b(triangulation.n_vertices());

	std::vector<bool> vertex_touched (triangulation.n_vertices(),
		                                     false);
	unsigned int j = 0;
	std::vector<unsigned int> ind_bnry_row;
	std::vector<unsigned int> ind_bnry_col;

	// initial indicator
	bool indicator[3] = {true, true, true};

	// set tolerance, how close to the x = 0 or y = 0 or z = 0 planes
	double eps = 1e-3;

	// assemble the sensitivity matrix and r.h.s.
	for (; cell!=endc; ++cell)
        for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
	{
			unsigned char boundary_ids = cell->face(f)->boundary_indicator();
			if(boundary_ids == 1)
			{
				for (unsigned int v=1; v<GeometryInfo<dim>::vertices_per_face; ++v)
					if (vertex_touched[cell->face(f)->vertex_index(v)] == false)
					{
						vertex_touched[cell->face(f)->vertex_index(v)] = true;
						
						for (unsigned int i=0; i<dim; ++i)
						{
					    	// sensitivity matrix entry
					    	A(j,i) = pow(cell->face(f)->vertex(v)[i],2);
					    	// r.h.s. entry
					    	b[j] = 1.0;
					    	// if mesh if not full: set the indicator
					    	if (!full_mesh)
					    		indicator[i] = fabs(cell->face(f)->vertex(v)[i]) > eps;
						}
						if (!full_mesh) // check if all indicators are true, i.e. vertex is not at the cut
					       	if (std::all_of(std::begin(indicator),std::end(indicator),[](bool i){return i;}))
								ind_bnry_row.push_back(j);
						j++;
					}
			}
	}
    // maxtrix A'*A and vector A'*b;  A'*A*x = A'*b -- normal system of equations
	FullMatrix<double> AtA(dim,dim);
	Vector<double>     Atb(dim);


	 if (!full_mesh) // quadrant/octant mesh, extract only necessary values
	 {
	     FullMatrix<double> A_out(ind_bnry_row.size(),dim);
	     Vector<double>     b_out(ind_bnry_row.size());

	     for(unsigned int i=0;i<dim;i++)
		     ind_bnry_col.push_back(i);

	     for(unsigned int i=0;i<ind_bnry_row.size();i++)
		     b_out(i) = 1;

	     A_out.extract_submatrix_from(A, ind_bnry_row, ind_bnry_col);
	     A_out.Tmmult(AtA,A_out,true);
	     A_out.Tvmult(Atb,b_out,true);
	 }
	 else // full mesh, just go with what you already have
	 {
		 A.Tmmult(AtA,A,true);
		 A.Tvmult(Atb,b,true);
	 }

     // solve normal system of equations
	 SolverControl           solver_control (1000, 1e-12);
	 SolverCG<>              solver (solver_control);
	 solver.solve (AtA, x, Atb, PreconditionIdentity());

	 // find ellipsoidal axes
	 for(unsigned int i=0; i<dim; i++)
         ell.push_back(sqrt(1.0/x[i]));

}

template<int dim>
void FE_test<dim>::cp_tr(const Triangulation<dim,dim> &tr_in)
{
	triangulation.copy_triangulation(tr_in);
}

template <int dim>
FE_test<dim>::FE_test () :
  fe (1),
  dof_handler (triangulation)
{}

template <int dim>
void FE_test<dim>::generate_ellipsoid(double* ell, int nrefine)
{
	  const double radius = 1.0;

      const Point<dim> center = (dim == 2 ?
                                      Point<dim>(0.0,0.0) :
                                      Point<dim>(0.0,0.0,0.0));

	  GridGenerator::hyper_ball (triangulation, center, radius);

	  static const HyperBallBoundary<dim> boundary(center,radius);

	  triangulation.set_boundary(0, boundary);
	  triangulation.refine_global (nrefine);

	  // deforming mesh to create an ellipsoid

	  std::vector<bool> vertex_touched (triangulation.n_vertices(),
	                                     false);

	  for (typename DoFHandler<dim>::active_cell_iterator
	        cell = dof_handler.begin_active ();
	        cell != dof_handler.end(); ++cell)
	     for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	       if (vertex_touched[cell->vertex_index(v)] == false)
	         {
	           vertex_touched[cell->vertex_index(v)] = true;
	           for (unsigned int d=0; d<dim; ++d)
	        	   cell->vertex(v)[d] *= ell[d];
	         }
}

template <int dim>
void FE_test<dim>::output_results()
{
	  std::vector<bool> vertex_touched (triangulation.n_vertices(),
	                                     false);

	  // print all boundary vertices to a file
	  FILE * pFile;
	  pFile = fopen ("outer_points.txt","w");
	
	  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
	  typename Triangulation<dim>::active_cell_iterator endc = triangulation.end();

	  for (; cell!=endc; ++cell)
	          for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
	              if (cell->face(f)->at_boundary())
	            	  for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v)
	            		  if (vertex_touched[cell->face(f)->vertex_index(v)] == false)
	            		  {
	            			  vertex_touched[cell->face(f)->vertex_index(v)]  = true;
	            			  for (unsigned int i=0; i<dim; ++i)
	            		          fprintf(pFile,"%f ", cell->face(f)->vertex(v)[i]);
	            			  fprintf(pFile,"\n");
	            		  }

	  fclose (pFile);

	   // refining mesh
//	  for (unsigned int cycle=0; cycle<1; ++cycle)
//	      {
//	        std::cout << "Cycle " << cycle << ':' << std::endl;
//	        if (cycle != 0)
//	          triangulation.refine_global (1);
//	        std::cout << "   Number of active cells: "
//	                  << triangulation.n_active_cells()
//	                  << std::endl
//	                  << "   Total number of cells: "
//	                  << triangulation.n_cells()
//	                  << std::endl;
//	      }

	  std::ofstream out_eps ("mesh_test_def.eps");
	  std::ofstream out_vtk ("mesh_test_def.vtk");
//	  std::ofstream out_ucd ("mesh_test_def.inp");
	  GridOut grid_out;
	  grid_out.write_eps (triangulation, out_eps);
	  grid_out.write_vtk (triangulation, out_vtk);
//	  grid_out.write_ucd (triangulation, out_ucd);
}

template <int dim>
void FE_test<dim>::run ()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  Assert (((dim==3) | (dim == 2)), ExcInternalError())
  std::ifstream input_file("mesh_test_def.inp");
  grid_in.read_ucd (input_file);

  double ell[dim];
  fit_ellipsoid(ell, true);

  for(unsigned int i=0; i<dim; i++)
	  printf("axis[%i] = %f\n",i,ell[i]);

//	double ell[] =  {1.0, 1.0, 1.0};
//	generate_ellipsoid(ell, 6);
  output_results();

}