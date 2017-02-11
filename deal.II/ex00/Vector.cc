#include <iostream>
#include <cassert>

#include "Vector.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
Vector<T>::Vector (unsigned int nrow)
   :
      nrow (nrow)
{
   assert (nrow > 0);
   val.resize (nrow);
}

//-----------------------------------------------------------------------------
// Copy Constructor
//-----------------------------------------------------------------------------
template <class T>
Vector<T>::Vector (const Vector<T>& v)
{
   nrow = v.size ();
   val.resize (nrow);
   for(unsigned int i=0; i<nrow; ++i)
      val[i] = v(i);
}

//-----------------------------------------------------------------------------
// Return size of vector
//-----------------------------------------------------------------------------
template <class T>
unsigned int Vector<T>::size () const
{
   return nrow;
}

//-----------------------------------------------------------------------------
// Return i'th value
//-----------------------------------------------------------------------------
template <class T>
T Vector<T>::operator() (unsigned int i) const
{
   return val[i];
}

//-----------------------------------------------------------------------------
// Return reference to i'th value
//-----------------------------------------------------------------------------
template <class T>
T& Vector<T>::operator() (unsigned int i)
{
   return val[i];
}

//-----------------------------------------------------------------------------
// Set vector to a constant value
//-----------------------------------------------------------------------------
template <class T>
Vector<T>& Vector<T>::operator= (const T scalar)
{
   for(unsigned int i=0; i<nrow; ++i)
      val[i] = scalar;
   return *this;
}

//-----------------------------------------------------------------------------
// this = scalar * this
//-----------------------------------------------------------------------------
template <class T>
Vector<T>& Vector<T>::operator*= (const T scalar)
{
   for(unsigned int i=0; i<nrow; ++i)
      val[i] *= scalar;
   return *this;
}

//-----------------------------------------------------------------------------
// this = this - v
//-----------------------------------------------------------------------------
template <class T>
Vector<T>& Vector<T>::operator-= (const Vector<T>& v)
{
   assert (nrow == v.size());
   for(unsigned int i=0; i<nrow; ++i)
      val[i] -= v(i);
   return *this;
}

//-----------------------------------------------------------------------------
// this = this + v
//-----------------------------------------------------------------------------
template <class T>
Vector<T>& Vector<T>::operator+= (const Vector<T>& v)
{
   assert (nrow == v.size());
   for(unsigned int i=0; i<nrow; ++i)
      val[i] += v(i);
   return *this;
}

//-----------------------------------------------------------------------------
// this = this + scalar * v
//-----------------------------------------------------------------------------
template <class T>
void Vector<T>::add(const T          scalar, 
                    const Vector<T>& v)
{
   assert (nrow == v.size());
   for(unsigned int i=0; i<nrow; ++i)
      val[i] += scalar * v(i);
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class Vector<int>;
template class Vector<float>;
template class Vector<double>;
