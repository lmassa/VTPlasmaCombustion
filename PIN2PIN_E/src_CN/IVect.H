#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IVECT_H_
#define _IVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "SPACE.H"
#include "Misc.H"
#include "REAL.H"
#include "Vector.H"
#include "SPMD.H"
#include "NamespaceHeader.H"


/// A int vector in n-dimensional space
/**
  The class IVect is an implementation of a int vector in a
  n-dimensional space.
  IVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

class IVect: public Vector<int>
{
public: 

  /**
     \name Constructors and Accessors
  */
  /*@{*/

 
  IVect ();
 ///
  /**
     Construct a IVect whose components are zero.
  */
  IVect (unsigned int isize);

  IVect (unsigned int isize, unsigned int jone, const int& s);

  IVect (unsigned int i1, unsigned int i2);


  IVect (const Vector<int>& rhs);


  IVect (const std::vector<int>& invec);



  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this IVect is equivalent to argument IVect.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator== (const IVect& p) const;

  ///
  /**
     Returns true if this IVect is different from argument IVect.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!= (const IVect& p) const;

  ///
  /**
     Returns true if this IVect is less than argument IVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an IVect to be neither greater than, less than, nor equal to
     another.
  */
  bool operator< (const IVect& p) const;

  ///
  /**
     Returns true if this IVect is less than or equal to argument
     IVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an IVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  bool operator<= (const IVect& p) const;

  ///
  /**
     Returns true if this IVect is greater than argument IVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an IVect to be neither greater than, less than, nor equal
     to another.
  */
  bool operator> (const IVect& p) const;

  ///
  /**
     Returns true if this IVect is greater than or equal to argument
     IVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an IVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */

  bool operator>= (const IVect& p) const;

  /*@}*/

  /**
     \name Arithmetic Operators
  */
  /*@{*/

  ///
  /**
     Modifies this IVect by addition of a scalar to each component.
  */
  IVect& operator+= (int s);

  ///
  /**
     Returns a IVect that is this IVect with a scalar s added
     to each component.
  */
  IVect operator+ (int s) const;

  ///
  /**
     Modifies this IVect by component-wise addition by argument.
  */
  IVect& operator+= (const IVect& p);

  ///
  /**
     Modifies this IVect by subtraction of a scalar from each
     component.
  */
  IVect& operator-= (int s);

  IVect& operator = (int s);

  ///
  /**
     Modifies this IVect by component-wise subtraction by argument.
  */
  IVect& operator-= (const IVect& p);

  ///
  /**
     Returns a IVect that is this IVect with a scalar s subtracted
     from each component.
  */
  IVect operator- (int s) const;

  ///
  /**
     Modifies this IVect by multiplying each component by a scalar.
  */
  IVect& operator*= (int s);

  ///
  /**

  */
  int dotProduct(const IVect& a_rhs) const;
  IVect matProduct(const IVect& a_B, const int& a_nrows, const int& a_ncols, const int& a_ncolA) const;
  IVect matProduct(const IVect& a_B,const int& a_n) const;

  ///
  /**
     Modifies this IVect by component-wise multiplication by argument.
  */
  IVect& operator*= (const IVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this IVect and argument.
//XXX  */
//XXX  IVect operator* (const IVect& p) const;

  ///
  /**
     Returns a IVect that is this IVect with each component
     multiplied by a scalar.
  */
  IVect operator* (int s) const;

  ///
  /**
     Modifies this IVect by dividing each component by a scalar.
  */
  IVect& operator/= (int s);

  ///
  /**
     Modifies this IVect by component-wise division by argument.
  */
  IVect& operator/= (const IVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise quotient of this IVect by argument.
//XXX  */
//XXX  IVect operator/ (const IVect& p) const;

  ///
  /**
     Returns a IVect that is this IVect with each component
     divided by a scalar.
  */
  IVect operator/ (int s) const;

  ///
  /**
     Modifies this IVect by multiplying each component by a scalar.
  */
  IVect& scale (int s);

  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this IVect by taking component-wise min with IVect
     argument.
  */
  IVect& min (const IVect& p);

  ///
  /**
     Returns the IVect that is the component-wise minimum of two
     argument IVects.
  */
  friend inline IVect min (const IVect& p1,
			   const IVect& p2);

  ///
  /**
     Modifies this IVect by taking component-wise max with IVect
     argument.
  */
  IVect& max (const IVect& p);

  ///
  /**
     Returns the IVect that is the component-wise maximum of two
     argument IVects.
  */
  friend inline IVect max (const IVect& p1,
                              const IVect& p2);

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  IVect operator+ () const;

  ///
  /**
     Unary minus -- negates all components of this IVect.
  */
  IVect operator- () const;

  ///
  /**
     Sum of all components of this IVect.
  */
  int sum () const;

  ///
  /**
     sqrt(sum squares)
  */
  int vectorLength() const;

  ///
  /**
     Product of all components of this IVect.
  */
  int product () const;

  ///
  /**
     Component with the minimum value of this IVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int minDir(const bool& a_doAbs) const;

  ///
  /**
     Component with the maximum value of this IVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int maxDir(const bool& a_doAbs) const;

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Only for sending stuff to Fortran
   */
  const int* dataPtr() const;

  ///
  /**
     Only for sending stuff to Fortran
   */
  int* dataPtr() ;

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/


  static IVect Identity(const int& a_nrows, const int& a_ncols);
  static IVect Identity(const int& a_n);
  static IVect Identity(const int& a_nrows, const int& a_ncols, const int& s);

  /*@}*/

  /**
     \name Arithmetic friend functions
  */
  /*@{*/

  ///
  /**
     Returns a IVect that is a IVect <i>p</i> with
     a scalar <i>s</i> added to each component.
  */
  friend IVect operator+ (int            s,
			  const IVect&    p);

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend IVect operator- (int            s,
			  const IVect&    p);

  ///
  /**
     Returns a IVect that is a IVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend IVect operator* (int            s,
                         const IVect&     p);
  ///
  /**
     Returns a IVect that is a IVect <i>p</i> with each component
     divided by a scalar <i>s</i>.
  */
  friend IVect operator/ (int            s,
                         const IVect&     p);

  ///
  /**
     Returns component-wise sum of IVects <i>s</i> and <i>p</i>.
  */
  friend IVect operator+ (const IVect& s,
                          const IVect& p);

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend IVect operator- (const IVect& s,
                          const IVect& p);

  ///
  /**
     Returns component-wise product of <i>s</i> and <i>p</i>.
  */
  friend IVect operator* (const IVect& s,
                          const IVect& p);
  ///
  /**
     Returns component-wise quotient <i>p / s</i>.
  */
  friend IVect operator/ (const IVect& s,
                             const IVect& p);


  ///
  /**
     Returns a IVect obtained by multiplying each of the components
     of the given IVect by a scalar.
  */
  friend inline IVect scale (const IVect& p,
                                int            s);

  friend std::ostream& operator<< (std::ostream& ostr,
                                   const IVect& p);

  static size_t io_offset;

protected:


};

#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"

#include "NamespaceVar.H"
// IVect spcializations of linearization
template < >
int linearSize(const CH_XDIR::IVect& vindex);


#include "BaseNamespaceFooter.H"

#include "NamespaceHeader.H"


inline
IVect&
IVect::operator-= (int s)
{
  for(int i=0; i<v.size(); i++) v[i] -= s;
  return *this;
}


inline
IVect&
IVect::operator = (int s)
{
  for(int i=0; i<v.size(); i++) v[i] = s;
  return *this;
}

inline
IVect&
IVect::operator-= (const IVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] -= p.v[i];
  return *this;
}

inline
IVect
IVect::operator+ () const
{
  return IVect(*this);
}

inline
IVect
IVect::operator- () const
{
  IVect vr=v;
  for(int i=0; i<vr.size(); i++) vr[i] = -vr[i];
  return vr;
}

inline
IVect&
IVect::scale (int s)
{
  for(int i=0; i<v.size(); i++) v[i] *= s;
  return *this;
}

inline
int
IVect::sum () const
{
  int s=0;
  for(int i=0; i<v.size(); i++) s += v[i];

  return s;
}

inline
int
IVect::vectorLength () const
{
  int len = 0;
  for (int idir = 0; idir < v.size(); idir++)
    {
      len = len + v[idir]*v[idir];
    }
  len = sqrt(len);

  return len;
}

inline
int
IVect::product () const
{
  int s=1;
  for(int i=0; i<v.size(); i++) s *= v[i];
  return s;
}

inline
IVect
scale (const IVect& p,
       int         s)
{
  Vector<int> vr=p.v;
  for(int i=0; i<vr.size(); i++) vr[i] *= s;
  return IVect(vr);
}

inline
bool
IVect::operator< (const IVect& p) const
{
  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] < p.v[i]);
      if(!retval) break;
    }
  return retval;
}

inline
bool
IVect::operator<= (const IVect& p) const
{
  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] <= p.v[i]);
      if(!retval) break;
    }
  return retval;
}


inline
bool
IVect::operator> (const IVect& p) const
{
  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] > p.v[i]);
      if(!retval) break;
    }
  return retval;
}

inline
bool
IVect::operator>= (const IVect& p) const
{
  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] >= p.v[i]);
      if(!retval) break;
    }
  return retval;
}

inline
IVect&
IVect::min (const IVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] = Min(v[i], p.v[i]);
  return *this;
}

inline
IVect&
IVect::max (const IVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] = Max(v[i], p.v[i]);
  return *this;
}

inline
IVect
min (const IVect& p1,
     const IVect& p2)
{
  IVect p(p1);
  return p.min(p2);
}

inline
IVect
max (const IVect& p1,
     const IVect& p2)
{
  IVect p(p1);
  return p.max(p2);
}



#include "NamespaceFooter.H"
#endif
