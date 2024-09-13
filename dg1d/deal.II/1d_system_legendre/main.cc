#include "dg.h"
#include "problem_data.h"

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int
main(int argc, char** argv)
{
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
   ph.print_parameters(std::cout, ParameterHandler::Text);

   Parameter param;
   parse_parameters(ph, param);

   QGauss<1> cell_quadrature(param.degree+1);

   ScalarProblem<1> problem(param, cell_quadrature);
   problem.run();

   return 0;
}
