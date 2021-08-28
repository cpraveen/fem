#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>


namespace Step26
{
   using namespace dealii;
   
   
   template<int dim>
   class HeatEquation
   {
   public:
      HeatEquation();
      void run();
      
   private:
      void setup_system();
      void solve_time_step();
      void output_results() const;
      
      Triangulation<dim>   triangulation;
      FE_Q<dim>            fe;
      DoFHandler<dim>      dof_handler;
      SparsityPattern      sparsity_pattern;
      SparseMatrix<double> mass_matrix;
      SparseMatrix<double> laplace_matrix;
      SparseMatrix<double> system_matrix;
      
      Vector<double>       solution;
      Vector<double>       old_solution;
      Vector<double>       system_rhs;
      
      double               time;
      double               time_step;
      unsigned int         timestep_number;
      
      const double         theta;
   };
   
   
   
   
   template<int dim>
   class RightHandSide : public Function<dim>
   {
   public:
      RightHandSide ()
      :
      Function<dim>(),
      period (0.2)
      {}
      
      virtual double value (const Point<dim> &p,
                            const unsigned int component = 0) const;
      
   private:
      const double period;
   };
   
   
   
   template<int dim>
   double RightHandSide<dim>::value (const Point<dim> &p,
                                     const unsigned int component) const
   {
      Assert (component == 0, ExcInternalError());
      Assert (dim == 2, ExcNotImplemented());
      
      const double time = this->get_time();
      const double point_within_period = (time/period - std::floor(time/period));
      
      if ((point_within_period >= 0.0) && (point_within_period <= 0.2))
      {
         if ((p[0] > 0.5) && (p[1] > -0.5))
            return 1;
         else
            return 0;
      }
      else if ((point_within_period >= 0.5) && (point_within_period <= 0.7))
      {
         if ((p[0] > -0.5) && (p[1] > 0.5))
            return 1;
         else
            return 0;
      }
      else
         return 0;
   }
   
   
   
   template<int dim>
   class BoundaryValues : public Function<dim>
   {
   public:
      virtual double value (const Point<dim>  &p,
                            const unsigned int component = 0) const;
   };
   
   
   
   template<int dim>
   double BoundaryValues<dim>::value (const Point<dim> &/*p*/,
                                      const unsigned int component) const
   {
      Assert(component == 0, ExcInternalError());
      return 0;
   }
   
   
   
   template<int dim>
   HeatEquation<dim>::HeatEquation ()
   :
   fe(1),
   dof_handler(triangulation),
   time_step(1.0 / 500),
   theta(0.5)
   {}
   
   
   
   template<int dim>
   void HeatEquation<dim>::setup_system()
   {
      dof_handler.distribute_dofs(fe);
      
      std::cout << std::endl
      << "==========================================="
      << std::endl
      << "Number of active cells: " << triangulation.n_active_cells()
      << std::endl
      << "Number of degrees of freedom: " << dof_handler.n_dofs()
      << std::endl
      << std::endl;
      
      DynamicSparsityPattern dsp(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler,
                                      dsp);
      sparsity_pattern.copy_from(dsp);
      
      mass_matrix.reinit(sparsity_pattern);
      laplace_matrix.reinit(sparsity_pattern);
      system_matrix.reinit(sparsity_pattern);
      
      MatrixCreator::create_mass_matrix(dof_handler,
                                        QGauss<dim>(fe.degree+1),
                                        mass_matrix);
      MatrixCreator::create_laplace_matrix(dof_handler,
                                           QGauss<dim>(fe.degree+1),
                                           laplace_matrix);
      
      solution.reinit(dof_handler.n_dofs());
      old_solution.reinit(dof_handler.n_dofs());
      system_rhs.reinit(dof_handler.n_dofs());
   }
   
   
   template<int dim>
   void HeatEquation<dim>::solve_time_step()
   {
      SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
      SolverCG<> cg(solver_control);
      
      PreconditionSSOR<> preconditioner;
      preconditioner.initialize(system_matrix, 1.0);
      
      cg.solve(system_matrix, solution, system_rhs,
               preconditioner);
      
      std::cout << "     " << solver_control.last_step()
      << " CG iterations." << std::endl;
   }
   
   
   
   template<int dim>
   void HeatEquation<dim>::output_results() const
   {
      DataOut<dim> data_out;
      
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "U");
      
      data_out.build_patches();
      
      const std::string filename = ("solution-"
                                    + Utilities::int_to_string(timestep_number, 3) +
                                    ".vtk");
      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
   }
   
   template<int dim>
   void HeatEquation<dim>::run()
   {
      const unsigned int initial_global_refinement = 4;
      GridGenerator::hyper_L (triangulation);
      triangulation.refine_global (initial_global_refinement);
      
      setup_system();
      
      Vector<double> tmp;
      Vector<double> forcing_terms;
      
      tmp.reinit (solution.size());
      forcing_terms.reinit (solution.size());
      
      // Initial condition is zero
      VectorTools::interpolate(dof_handler,
                               ZeroFunction<dim>(),
                               old_solution);
      solution = old_solution;
      
      timestep_number = 0;
      time            = 0;
      
      output_results();
      
      while (time <= 0.5)
      {
         time += time_step;
         ++timestep_number;
         
         std::cout << "Time step " << timestep_number << " at t=" << time
         << std::endl;
         
         mass_matrix.vmult(system_rhs, old_solution);
         laplace_matrix.vmult(tmp, old_solution);
         system_rhs.add(-(1 - theta) * time_step, tmp);
         
         RightHandSide<dim> rhs_function;
         rhs_function.set_time(time);
         VectorTools::create_right_hand_side(dof_handler,
                                             QGauss<dim>(fe.degree+1),
                                             rhs_function,
                                             tmp);
         forcing_terms = tmp;
         forcing_terms *= time_step * theta;
         
         rhs_function.set_time(time - time_step);
         VectorTools::create_right_hand_side(dof_handler,
                                             QGauss<dim>(fe.degree+1),
                                             rhs_function,
                                             tmp);
         
         forcing_terms.add(time_step * (1 - theta), tmp);
         
         system_rhs += forcing_terms;
         
         system_matrix.copy_from(mass_matrix);
         system_matrix.add(theta * time_step, laplace_matrix);
         
         {
            BoundaryValues<dim> boundary_values_function;
            boundary_values_function.set_time(time);
            
            std::map<types::global_dof_index, double> boundary_values;
            VectorTools::interpolate_boundary_values(dof_handler,
                                                     0,
                                                     boundary_values_function,
                                                     boundary_values);
            
            MatrixTools::apply_boundary_values(boundary_values,
                                               system_matrix,
                                               solution,
                                               system_rhs);
         }
         
         solve_time_step();
         output_results();
         old_solution = solution;
      }
   }
}


int main()
{
   try
   {
      using namespace dealii;
      using namespace Step26;
      
      HeatEquation<2> heat_equation_solver;
      heat_equation_solver.run();
      
   }
   catch (std::exception &exc)
   {
      std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
      << std::endl << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
      
      return 1;
   }
   catch (...)
   {
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
