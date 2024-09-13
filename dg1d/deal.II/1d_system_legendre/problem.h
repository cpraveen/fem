namespace Problem
{
   extern double xmin, xmax;
   extern bool periodic;
   void initial_value(const Point<1>& p, 
                      Vector<double>& u);
   void boundary_value(const int       id, 
                       const double    t,
                       Vector<double>& ul,
                       Vector<double>& ur);
}
