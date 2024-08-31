using namespace dealii;

enum class TestCase {sin4pi, exp};

std::map<std::string, TestCase> TestCaseList{{"sin4pi", TestCase::sin4pi}, 
                                             {"exp",    TestCase::exp}};

//------------------------------------------------------------------------------
TestCase get_test_case(std::string value)
{
   auto search = TestCaseList.find(value);
   if (search != TestCaseList.end())
      return search->second;
   else
   {
      std::cout << "Available test cases\n";
      for (const auto &v : TestCaseList)
         std::cout << v.first << std::endl;
      AssertThrow(false, ExcMessage("Unknown test case"));
   }
}

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
      if(test_case == TestCase::sin4pi)
      {
         xmin = -1.0;
         xmax = +1.0;
      }
      else if(test_case == TestCase::exp)
      {
         xmin = -1.0;
         xmax = +1.0;
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
   if(test_case == TestCase::sin4pi)
   {
      value = std::sin(4.0 * M_PI * x);
   }
   else if(test_case == TestCase::exp)
   {
      value = exp(-50.0 * x * x);
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
// Exact solution not known
//------------------------------------------------------------------------------
template<int dim>
double
Solution<dim>::value(const Point<dim>&    p,
                     const unsigned int) const
{
   return 0.0;
}

//------------------------------------------------------------------------------
// Exact solution not known
//------------------------------------------------------------------------------
template<int dim>
Tensor<1, dim>
Solution<dim>::gradient(const Point<dim>&    p,
                        const unsigned int) const
{
   Tensor<1, dim> values;
   values[0] = 0.0;
   return values;
}
