using namespace dealii;

enum class TestCase {sine, hat, exact1};
std::map<std::string, TestCase> TestCaseList{{"sine",   TestCase::sine},
                                             {"hat",    TestCase::hat},
                                             {"exact1", TestCase::exact1}};

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
         xmin = 0.0;
         xmax = 1.0;
      }
      else if(test_case == TestCase::hat)
      {
         xmin = -1.0;
         xmax = +1.0;
      }
      else if(test_case == TestCase::exact1)
      {
         xmin = 0.0;
         xmax = 2.0 * M_PI;
      }
      else
      {
         AssertThrow(false, ExcMessage("Unknown test case"));
      }
   }

   double value(const Point<dim>& p,
                const unsigned int component = 0) const override;
   double xmin, xmax;

private:
   const TestCase test_case;
};

// Initial condition
template<int dim>
double
InitialCondition<dim>::value(const Point<dim>& p,
                             const unsigned int /* component */) const
{
   double x = p[0];
   double value = 0;

   if(test_case == TestCase::sine)
   {
      value = 1.5 + std::sin(2.0 * M_PI* x);
   }
   else if(test_case == TestCase::hat)
   {
      if(std::fabs(x) < 0.25)
         value = 1.0;
      else
         value = 0.0;
   }
   else if(test_case == TestCase::exact1)
   {
      value = 0.2 * std::sin(x);
   }
   else
   {
      AssertThrow(false, ExcMessage("Unknown test case"));
   }

   return value;
}

//------------------------------------------------------------------------------
// solve for u from
//       u = u0(x - t*u)
// where u0(x) = 0.2 * sin(x) using Newton method
//------------------------------------------------------------------------------
double
solution_exact1(const double x, const double t)
{
   const double EPS = 1.0e-14;
   const int ITMAX = 100;
   double u = 0.0;
   double f = u - 0.2 * std::sin(x - t* u);
   double df = 1.0 + 0.2 * t * std::cos(x - t* u);
   int iter = 0;
   while(std::fabs(f) > EPS && iter < ITMAX)
   {
      u = u - f / df;
      f = u - 0.2 * std::sin(x - t* u);
      df = 1.0 + 0.2 * t * std::cos(x - t* u);
      ++iter;
   }
   if(std::fabs(f) > EPS)
   {
      std::cout << "Newton did not converge\n";
      std::cout << "x,t,u,f = " << x << " " << t << " " << u << " "
                << f << std::endl;
      exit(0);
   }
   return u;
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

// Exact solution works correctly only for periodic case
template<int dim>
double
Solution<dim>::value(const Point<dim>&    p,
                     const unsigned int) const
{
   double x = p[0];
   double t = final_time;
   if(test_case == TestCase::exact1)
   {
      // u = u0(x - t*u)
      return solution_exact1(x, t);
   }
   else
   {
      return 0.0;
   }
}


// Exact solution works correctly only for periodic case
template<int dim>
Tensor<1, dim>
Solution<dim>::gradient(const Point<dim>&    p,
                        const unsigned int) const
{
   double x = p[0];
   double t = final_time;
   Tensor<1, dim> values;
   if(test_case == TestCase::exact1)
   {
      // u = u0(x - t*u)
      double u = solution_exact1(x, t);
      // u_x = u0'(x-t*u) * (1 - t * u_x)
      double du = 0.2 * std::cos(x - t* u);
      values[0] = du / (1.0 + t* du);
   }
   else
   {
      values[0] = 0;
   }
   return values;
}
