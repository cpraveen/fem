#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <vector>
#include <iostream>

template <class T>
class Vector
{
   public:
      Vector (unsigned int nrow);
      Vector (const Vector<T>& v);
      ~Vector (){};
      unsigned int size () const;
      T& operator()(unsigned int i);
      T  operator()(unsigned int i) const;
      Vector<T>& operator=  (const T scalar);
      Vector<T>& operator*= (const T scalar);
      Vector<T>& operator-= (const Vector<T>& v);
      Vector<T>& operator+= (const Vector<T>& v);
      void       add(const T          scalar, 
                     const Vector<T>& v);
      friend std::ostream& operator<< (std::ostream&    os, 
                                       const Vector<T>& v)
      {
         for(unsigned int i=0; i<v.size(); ++i)
            os << v(i) << std::endl;
         return os;
      }

   private:
      unsigned int   nrow;
      std::vector<T> val;
};

#endif
