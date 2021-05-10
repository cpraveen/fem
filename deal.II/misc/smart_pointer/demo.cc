#include <deal.II/base/smartpointer.h>
#include <deal.II/fe/fe_q.h>

using namespace dealii;

int main()
{
   SmartPointer<FiniteElement<2>> fe1;   // SmartPointer
   //FiniteElement<2>* fe1;              // normal pointer
   {
      FE_Q<2> fe(2);
      fe1 = &fe;
      std::cout << "Degree = " << fe1->degree << std::endl;
   }
   // fe has been destroyed, next line should give an error for SmartPointer
   // If you dont use SmartPointer, you may not get any error message,
   // some random thing may happen, result is not guaranteed to be correct.
   // It will be very difficult to catch such errors because code will run
   // but may give wrong answer.
   std::cout << "Degree = " << fe1->degree << std::endl;
   return 0;
}
