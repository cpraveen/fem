#include <deal.II/base/observer_pointer.h>
#include <deal.II/fe/fe_q.h>

using namespace dealii;

int main()
{
   FiniteElement<2>* fe1;                   // normal pointer
   ObserverPointer<FiniteElement<2>> fe2;   // ObserverPointer
   {
      FE_Q<2> fe(2);
      fe1 = &fe;
      fe2 = &fe;
      std::cout << "Normal   pointer: Degree = " << fe1->degree << std::endl;
      std::cout << "Observer pointer: Degree = " << fe2->degree << std::endl;
   }
   // fe has been destroyed, next line should give an error for ObserverPointer
   // If you dont use ObserverPointer, you may not get any error message,
   // some random thing may happen, result is not guaranteed to be correct.
   // It will be very difficult to catch such errors because code will run
   // but may give wrong answer.
   std::cout << "Normal   pointer: Degree = " << fe1->degree << std::endl;
   std::cout << "Observer pointer: Degree = " << fe2->degree << std::endl;
   return 0;
}
