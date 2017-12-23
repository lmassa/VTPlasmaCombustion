#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RVECT_H_
#define _RVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "SPACE.H"
#include "Misc.H"
#include "REAL.H"
#include "Vector.H"
#include "SPMD.H"
#include "IVect.H"
#include "NamespaceHeader.H"


/// A Real vector in n-dimensional space
/**
  The class RVect is an implementation of a Real vector in a
  n-dimensional space.
  RVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

class RVect: public Vector<Real>
{
public: 

  /**
     \name Constructors and Accessors
  */
  /*@{*/

 
  RVect ();
 ///
  /**
     Construct a RVect whose components are zero.
  */
  RVect (unsigned int isize);

  RVect (unsigned int isize, unsigned int jone);

  RVect (unsigned int isize, unsigned int jone, const Real& s);

  RVect (unsigned int isize, const Real& s);

  RVect (const Real& s1, const Real& s2, unsigned int isize);


  RVect (const Vector<Real>& rhs);
  RVect (const IVect& rhs);


  RVect (const std::vector<Real>& invec);

  RVect (unsigned int a_nbr, unsigned int a_nrowb, const RVect& a_blkMat, const RVect& a_bCoe);


  struct index_greater : std::unary_function<Real, bool>
  {
    index_greater(const Real &idx) : idx_(idx) {}
    bool operator()(const Real &arg) const { return arg > idx_; }
    const Real &idx_;
  };

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this RVect is equivalent to argument RVect.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator== (const RVect& p) const;

  ///
  /**
     Returns true if this RVect is different from argument RVect.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!= (const RVect& p) const;

  ///
  /**
     Returns true if this RVect is less than argument RVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an RVect to be neither greater than, less than, nor equal to
     another.
  */
  bool operator< (const RVect& p) const;

  ///
  /**
     Returns true if this RVect is less than or equal to argument
     RVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  bool operator<= (const RVect& p) const;

  ///
  /**
     Returns true if this RVect is greater than argument RVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an RVect to be neither greater than, less than, nor equal
     to another.
  */
  bool operator> (const RVect& p) const;

  ///
  /**
     Returns true if this RVect is greater than or equal to argument
     RVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */

  bool operator>= (const RVect& p) const;

  /*@}*/

  /**
     \name Arithmetic Operators
  */
  /*@{*/

  ///
  /**
     Modifies this RVect by addition of a scalar to each component.
  */
  RVect& operator+= (Real s);

  ///
  /**
     Returns a RVect that is this RVect with a scalar s added
     to each component.
  */
  RVect operator+ (Real s) const;

  ///
  /**
     Modifies this RVect by component-wise addition by argument.
  */
  RVect& operator+= (const RVect& p);

  ///
  /**
     Modifies this RVect by subtraction of a scalar from each
     component.
  */
  RVect& operator-= (Real s);

  RVect& operator = (Real s);
  RVect& operator = (RVect s);
  RVect& operator = (IVect s);

  ///
  /**
     Modifies this RVect by component-wise subtraction by argument.
  */
  RVect& operator-= (const RVect& p);

  ///
  /**
     Returns a RVect that is this RVect with a scalar s subtracted
     from each component.
  */
  RVect operator- (Real s) const;

  ///
  /**
     Modifies this RVect by multiplying each component by a scalar.
  */
  RVect& operator*= (Real s);

  ///
  /**

  */
  Real dotProduct(const RVect& a_rhs) const;
  RVect matProduct(const RVect& a_B, const int& a_nrows, const int& a_ncols, const int& a_ncolA) const;
  RVect matProduct(const RVect& a_B,const int& a_n) const;

  RVect pwr(const Real& a_s) const;
  //
//note matmat is the same as matProduct but returning *this
  RVect& matMat(const RVect& a_B, const int& a_nrows, const int& a_ncols, const int& a_ncolA);
  RVect& matMat(const RVect& a_B,const int& a_n);

  ///
  /**
     Modifies this RVect by component-wise multiplication by argument.
  */
  RVect& operator*= (const RVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this RVect and argument.
//XXX  */
//XXX  RVect operator* (const RVect& p) const;

  ///
  /**
     Returns a RVect that is this RVect with each component
     multiplied by a scalar.
  */
  RVect operator* (Real s) const;

  ///
  /**
     Modifies this RVect by dividing each component by a scalar.
  */
  RVect& operator/= (Real s);

  ///
  /**
     Modifies this RVect by component-wise division by argument.
  */
  RVect& operator/= (const RVect& p);

//XXX  ///
//XXX  /**
//XXX     Returns component-wise quotient of this RVect by argument.
//XXX  */
//XXX  RVect operator/ (const RVect& p) const;

  ///
  /**
     Returns a RVect that is this RVect with each component
     divided by a scalar.
  */
  RVect operator/ (Real s) const;

  ///
  /**
     Modifies this RVect by multiplying each component by a scalar.
  */
  RVect& scale (Real s);

  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this RVect by taking component-wise min with RVect
     argument.
  */
  RVect& min (const RVect& p);

  ///
  /**
     Returns the RVect that is the component-wise minimum of two
     argument RVects.
  */
  friend inline RVect min (const RVect& p1,
			   const RVect& p2);

  ///
  /**
     Modifies this RVect by taking component-wise max with RVect
     argument.
  */
  RVect& max (const RVect& p);

  ///
  /**
     Returns the RVect that is the component-wise maximum of two
     argument RVects.
  */
  friend inline RVect max (const RVect& p1,
                              const RVect& p2);

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  RVect operator+ () const;

  ///
  /**
     Unary minus -- negates all components of this RVect.
  */
  RVect operator- () const;

  ///
  /**
     Sum of all components of this RVect.
  */
  Real sum () const;

  ///
  /**
     sqrt(sum squares)
  */
  Real vectorLength() const;

  ///
  /**
     Product of all components of this RVect.
  */
  Real product () const;

  ///
  /**
     Component with the minimum value of this RVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int minDir(const bool& a_doAbs) const;

  ///
  /**
     Component with the maximum value of this RVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int maxDir(const bool& a_doAbs) const;

  int minloc() const;
  Real maxAbs() const;
  Real min() const;
  Real max() const;

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Only for sending stuff to Fortran
   */
  const Real* dataPtr() const;

  ///
  /**
     Only for sending stuff to Fortran
   */
  Real* dataPtr() ;

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/


  static RVect Identity(const int& a_nrows, const int& a_ncols);
  static RVect Identity(const int& a_n);
  static RVect Identity(const int& a_nrows, const int& a_ncols, const Real& s);

  int findLE(const Real a_r) const;

  RVect& setRow(const int& a_irow, const RVect& a_B);

  /*@}*/

  /**
     \name Arithmetic friend functions
  */
  /*@{*/

  ///
  /**
     Returns a RVect that is a RVect <i>p</i> with
     a scalar <i>s</i> added to each component.
  */
  friend RVect operator+ (Real            s,
			  const RVect&    p);

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend RVect operator- (Real            s,
			  const RVect&    p);

  ///
  /**
     Returns a RVect that is a RVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend RVect operator* (Real            s,
                         const RVect&     p);
  ///
  /**
     Returns a RVect that is a RVect <i>p</i> with each component
     divided by a scalar <i>s</i>.
  */
  friend RVect operator/ (Real            s,
                         const RVect&     p);

  ///
  /**
     Returns component-wise sum of RVects <i>s</i> and <i>p</i>.
  */
  friend RVect operator+ (const RVect& s,
                          const RVect& p);

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend RVect operator- (const RVect& s,
                          const RVect& p);

  ///
  /**
     Returns component-wise product of <i>s</i> and <i>p</i>.
  */
  friend RVect operator* (const RVect& s,
                          const RVect& p);
  ///
  /**
     Returns component-wise quotient <i>p / s</i>.
  */
  friend RVect operator/ (const RVect& s,
                             const RVect& p);


  ///
  /**
     Returns a RVect obtained by multiplying each of the components
     of the given RVect by a scalar.
  */
  friend inline RVect scale (const RVect& p,
                                Real            s);

  friend std::ostream& operator<< (std::ostream& ostr,
                                   const RVect& p);

  static size_t io_offset;

protected:


};

#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"

#include "NamespaceVar.H"
// RVect spcializations of linearization
template < >
int linearSize(const CH_XDIR::RVect& vindex);


#include "BaseNamespaceFooter.H"

#include "NamespaceHeader.H"


inline
RVect&
RVect::operator-= (Real s)
{
  for(int i=0; i<v.size(); i++) v[i] -= s;
  return *this;
}


inline
RVect&
RVect::operator = (Real s)
{
  for(int i=0; i<v.size(); i++) v[i] = s;
  return *this;
}


inline
RVect&
RVect::operator = (RVect s)
{
  s.swap (*this); // Non-throwing swap
  return *this;
}


inline
RVect&
RVect::operator = (IVect s)
{
  RVect tmp(s);
  tmp.swap (*this); // Non-throwing swap
  return *this;
}

inline
RVect&
RVect::operator-= (const RVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] -= p.v[i];
  return *this;
}

inline
RVect
RVect::operator+ () const
{
  return RVect(*this);
}

inline
RVect
RVect::operator- () const
{
  RVect vr=v;
  for(int i=0; i<vr.size(); i++) vr[i] = -vr[i];
  return vr;
}

inline
RVect&
RVect::scale (Real s)
{
  for(int i=0; i<v.size(); i++) v[i] *= s;
  return *this;
}

inline
Real
RVect::sum () const
{
  Real s=0;
  for(int i=0; i<v.size(); i++) s += v[i];

  return s;
}

inline
Real
RVect::vectorLength () const
{
  Real len = 0;
  for (int idir = 0; idir < v.size(); idir++)
    {
      len = len + v[idir]*v[idir];
    }
  len = sqrt(len);

  return len;
}

inline
Real
RVect::product () const
{
  Real s=1;
  for(int i=0; i<v.size(); i++) s *= v[i];
  return s;
}

inline
RVect
scale (const RVect& p,
       Real         s)
{
  Vector<Real> vr=p.v;
  for(int i=0; i<vr.size(); i++) vr[i] *= s;
  return RVect(vr);
}

inline
bool
RVect::operator< (const RVect& p) const
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
RVect::operator<= (const RVect& p) const
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
RVect::operator> (const RVect& p) const
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
RVect::operator>= (const RVect& p) const
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
RVect&
RVect::min (const RVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] = Min(v[i], p.v[i]);
  return *this;
}

inline
RVect&
RVect::max (const RVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] = Max(v[i], p.v[i]);
  return *this;
}

inline
RVect
min (const RVect& p1,
     const RVect& p2)
{
  RVect p(p1);
  return p.min(p2);
}

inline
RVect
max (const RVect& p1,
     const RVect& p2)
{
  RVect p(p1);
  return p.max(p2);
}



#include "NamespaceFooter.H"
#endif
