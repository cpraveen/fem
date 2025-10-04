#include "dg.h"
#include "problem.h"

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
      ph.print_parameters(std::cout, ParameterHandler::PRM);
      return 0;
   }
   ph.parse_input(argv[1]);
   ph.print_parameters(std::cout, ParameterHandler::PRM);

   Problem<2> problem;
   Parameter param;
   param.final_time = problem.get_final_time(); // override this in input file
   parse_parameters(ph, param);

   DGSystem<2> solver(param, problem);
   solver.run();

   return 0;
}
