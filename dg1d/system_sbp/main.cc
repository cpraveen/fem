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
   param.final_time = Problem::final_time; // will override if given in input
   parse_parameters(ph, param);

   DGSystem<1> solver(param);
   solver.run();

   return 0;
}
