using namespace dealii;

enum class TestCase {sine, hat, trihat, exp};

std::map<std::string, TestCase> TestCaseList{{"sine",   TestCase::sine}, 
                                             {"hat",    TestCase::hat},
                                             {"trihat", TestCase::trihat},
                                             {"exp",    TestCase::exp}};
//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition(TestCase test_case)
      :
      Function<dim>(),
      test_case(test_case)
   {
      if(test_case == TestCase::sine)
      {
         xmin = -1.0;
         xmax = +1.0;
      }
      else if(test_case == TestCase::hat)
      {
         xmin = -1.0;
         xmax = +1.0;
      }
      else if(test_case == TestCase::trihat)
      {
         xmin = -1.0;
         xmax = +1.0;
      }
      else if(test_case == TestCase::exp)
      {
         xmin = -20.0;
         xmax = +20.0;
      }
      else
      {
         std::cout << "Unknown test case\n";
      }
   }

   double value(const Point<dim>& p,
                const unsigned int component = 0) const override;
   double xmin, xmax;

private:
   const TestCase test_case;
};

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
template<int dim>
double
InitialCondition<dim>::value(const Point<dim>& p,
                             const unsigned int /* component */) const
{
   double x = p[0];
   double value = 0;

   // test case: sine
   if(test_case == TestCase::sine)
   {
      value = -std::sin(M_PI* x);
   }
   else if(test_case == TestCase::hat)
   {
      if(std::fabs(x) < 0.25)
         value = 1.0;
      else
         value = 0.0;
   }
   else if(test_case == TestCase::trihat)
   {
      while(x >  1.0) x = x - 2.0;
      while(x < -1.0) x = x + 2.0;
      if(std::fabs(x) > 0.5)
         value = 0.0;
      else if(x < 0.0)
         value = 1.0 + 2.0 * x;
      else
         value = 1.0 - 2.0 * x;
   }
   else if(test_case == TestCase::exp)
   {
      value = exp(-10.0 * x* x);
   }
   else
   {
      AssertThrow(false, ExcMessage("Unknown test case"));
   }

   return value;
}

//------------------------------------------------------------------------------
// Exact solution
//------------------------------------------------------------------------------
template <int dim>
class Solution : public Function<dim>
{
public:
   Solution(TestCase test_case, double final_time)
      :
      Function<dim>(),
      test_case(test_case),
      final_time(final_time)
   {}

   double value(const Point<dim>&    p,
                const unsigned int  component = 0) const override;
   Tensor<1, dim> gradient(const Point<dim>&    p,
                           const unsigned int  component = 0) const override;
private:
   const TestCase test_case;
   const double final_time;
};

//------------------------------------------------------------------------------
// Exact solution works correctly only for periodic case
//------------------------------------------------------------------------------
template<int dim>
double
Solution<dim>::value(const Point<dim>&    p,
                     const unsigned int) const
{
   Point<dim> pp(p);
   pp[0] -= final_time;
   InitialCondition<dim> initial_condition(test_case);
   return initial_condition.value(pp);
}

//------------------------------------------------------------------------------
// Exact solution works correctly only for periodic case
//------------------------------------------------------------------------------
template<int dim>
Tensor<1, dim>
Solution<dim>::gradient(const Point<dim>&    p,
                        const unsigned int) const
{
   Tensor<1, dim> values;
   double x = p[0] - final_time;

   // test case: sine
   if(test_case == TestCase::sine)
   {
      values[0] = -M_PI * std::cos(M_PI* x);
   }
   else if(test_case == TestCase::hat)
   {
      values[0] = 0;
   }
   else if(test_case == TestCase::trihat)
   {
      while(x >  1.0) x = x - 2.0;
      while(x < -1.0) x = x + 2.0;
      if(std::fabs(x) > 0.5)
         values[0] = 0.0;
      else if(x < 0.0)
         values[0] = 2.0;
      else
         values[0] = -2.0;
   }
   else if(test_case == TestCase::exp)
   {
      values[0] = -20.0 * x * exp(-10.0 * x* x);
   }
   else
   {
      AssertThrow(false, ExcMessage("Unknown test case"));
   }

   return values;
}
