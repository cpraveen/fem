#include "dg.h"
#include "problem.h"

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int
main(int argc, char** argv)
{
   // Number of threads to use
   unsigned int n_threads = 1;
   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, n_threads);

   Problem<2> problem;

   ParameterHandler ph;
   declare_parameters(ph);
   if(argc < 2)
   {
      std::cout << "Specify input parameter file\n";
      std::cout << "It should contain following parameters.\n\n";
      ph.print_parameters(std::cout, ParameterHandler::Text);
      return 0;
   }
   ph.parse_input(argv[1]);
   if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      ph.print_parameters(std::cout, ParameterHandler::Text);

   Parameter param;
   param.final_time = problem.get_final_time(); // override this in input file
   parse_parameters(ph, param);

   Quadrature<1> quadrature_1d;
   if(param.basis == "gl")
   {
      quadrature_1d = QGauss<1>(param.degree+1);
   }
   else if(param.basis == "gll")
   {
      quadrature_1d = QGaussLobatto<1>(param.degree+1);
   }
   else
   {
      AssertThrow(false, ExcMessage("Unknown points"));
   }

   DGSystem<2> solver(param, problem, quadrature_1d);
   solver.run();

   return 0;
}
