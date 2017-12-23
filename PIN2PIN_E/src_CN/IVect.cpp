#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SPACE.H"
#include "IVect.H"
#include "Misc.H"
using std::ostream;

#include "NamespaceHeader.H"
#define CHOFFSET(object, member) (int)((char*)&(object.member) - (char*)&object)
using std::ws;


const int*
IVect::dataPtr() const
{
  return v.data();
}

int*
IVect::dataPtr()
{
  //not sure this works for this version of C++
  return v.data();
}



IVect::IVect ():Vector<int>()
{

}

IVect::IVect (unsigned int isize):Vector<int>(isize, 0)
{

}


IVect::IVect (unsigned int isize, unsigned int jone, const int& s)
{
  v.resize(isize, 0);
  v[jone]=s;
}


IVect::IVect (unsigned int i1, unsigned int i2)
{
  int isize = i2-i1+1;
  v.resize(isize);
  for(int i=i1; i<=i2;i++) v[i-i1]=i;
}

IVect::IVect (const Vector<int>& rhs)
{
  v = rhs.constStdVector();
}


IVect::IVect (const std::vector<int>& invec)
{
  v = invec;
}





int IVect::dotProduct(const IVect& a_rhs) const
{

  int s=0;
  for(int i=0; i<v.size(); i++) s += v[i]*a_rhs.v[i];
  return s;
}


IVect IVect::matProduct(const IVect& a_B, const int& a_nrows, const int& a_ncols, const int& a_ncolA) const
{

  IVect Cmat ( a_nrows*a_ncols );
  int C;
  for(int i=0; i<a_nrows; i++) 
    for(int j=0; j<a_ncols; j++)
      {
	C =  v[a_ncolA*i]*a_B[j];
	for(int k=1; k<a_ncolA; k++) 
	  C =  C + v[a_ncolA*i + k]*a_B[a_ncols*k + j];
	Cmat[a_ncols*i+j] =  C;
      }
  return Cmat;
}

//overload for square matrices
IVect IVect::matProduct(const IVect& a_B, const int& a_n) const
{

  IVect Cmat ( a_n*a_n );
  int C;
  for(int i=0; i<a_n; i++) 
    for(int j=0; j<a_n; j++)
      {
	C =  v[a_n*i]*a_B[j];
	for(int k=1; k<a_n; k++) 
	  C =  C + v[a_n*i + k]*a_B[a_n*k + j];
	Cmat[a_n*i+j] =  C;
      }
  return Cmat;
}

IVect IVect::Identity(const int& a_nrows, const int& a_ncols)
{
  IVect Cmat ( a_nrows*a_ncols );
  for(int i=0; i<a_nrows; i++) Cmat[a_ncols*i+i] =  1.0;
  return Cmat;
}

IVect IVect::Identity(const int& a_n)
{
  IVect Cmat ( a_n*a_n );
  for(int i=0; i<a_n; i++) Cmat[a_n*i+i] =  1.0;
  return Cmat;
}
IVect IVect::Identity(const int& a_nrows, const int& a_ncols, const int& s)
{
  IVect Cmat ( a_nrows*a_ncols );
  for(int i=0; i<a_nrows; i++) Cmat[a_ncols*i+i] =  s;
  return Cmat;
}




bool
IVect::operator== (const IVect& p) const
{

  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] == p.v[i]);
      if(!retval) break;
    }
  return retval;
}

bool
IVect::operator!= (const IVect& p) const
{
  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] != p.v[i]);
      if(!retval) break;
    }
  return retval;
}

IVect&
IVect::operator+= (int s)
{
  for(int i=0; i<v.size(); i++) v[i] += s;
  return *this;
}

IVect&
IVect::operator+= (const IVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] += p[i];
  return *this;
}

IVect&
IVect::operator*= (int s)
{
  for(int i=0; i<v.size(); i++) v[i] *= s;
  return *this;
}

IVect&
IVect::operator*= (const IVect &p)
{
  for(int i=0; i<v.size(); i++) v[i] *= p[i];
  return *this;
}

//XXXIVect
//XXXIVect::operator* (const IVect& p) const
//XXX{
//XXX  IVect v(D_DECL6(vect[0]*p[0], vect[1]*p[1], vect[2]*p[2]));
//XXX  return v;
//XXX}

IVect
IVect::operator* (int s) const
{
  IVect w=v;
  for(int i=0; i<w.size(); i++) w[i] *= s;
  return w;
}

IVect
IVect::operator- (int s) const
{
  IVect w=v;
  for(int i=0; i<w.size(); i++) w[i] -= s;
  return w;
}


IVect
IVect::operator+ (int s) const
{
  IVect w=v;
  for(int i=0; i<w.size(); i++) w[i] += s;
  return w;
}

IVect&
IVect::operator/= (int s)
{
  for(int i=0; i<v.size(); i++) v[i] /= s;
  return *this;
}

IVect&
IVect::operator/= (const IVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] /= p[i];
  return *this;
}

//XXXIVect
//XXXIVect::operator/ (const IVect & p) const
//XXX{
//XXX  IVect result( D_DECL6( vect[0] / p[0], vect[1] / p[1], vect[2] / p[2] ) );
//XXX  return result ;
//XXX}

IVect
IVect::operator/ (int s) const
{
  IVect w=v;
  for(int i=0; i<w.size(); i++) w[i] /= s;
  return w;
}


int
IVect::minDir(const bool& a_doAbs) const
{
  int mDir = 0;
  for (int idir=0; idir<size(); idir++)
    {
      if (a_doAbs)
        {
          if (Abs(v[idir]) < Abs(v[mDir]))
            {
              mDir = idir;
            }
        }
      else
        {
          if (v[idir] < v[mDir])
            {
              mDir = idir;
            }
        }
    }
  return mDir;
}

int
IVect::maxDir(const bool& a_doAbs) const
{
  int mDir = 0;
  for (int idir=0; idir<size(); idir++)
    {
      if (a_doAbs)
        {
          if (Abs(v[idir]) > Abs(v[mDir]))
            {
              mDir = idir;
            }
        }
      else
        {
          if (v[idir] > v[mDir])
            {
              mDir = idir;
            }
        }
    }
  return mDir;
}


IVect
operator/ (int            s,
           const IVect& p)
{
  IVect w=p;
  for(int i=0; i<w.size(); i++) w[i] = s/w[i];
  return w;
}
IVect
operator+ (int            s,
           const IVect& p)
{
  IVect w=p;
  for(int i=0; i<w.size(); i++) w[i] += s;
  return w;
}

IVect
operator- (int            s,
           const IVect& p)
{
  IVect w=p;
  for(int i=0; i<w.size(); i++) w[i] = s-w[i];
  return w;
}

IVect
operator* (int            s,
           const IVect& p)
{
  IVect w=p;
  for(int i=0; i<w.size(); i++) w[i] *= s;
  return w;
}

IVect
operator/ (const IVect& s,
           const IVect& p)
{
  IVect w=s;
  for(int i=0; i<w.size(); i++) w[i] /= p[i];
  return w;
}

IVect
operator+ (const IVect& s,
           const IVect& p)
{
  IVect w=s;
  for(int i=0; i<w.size(); i++) w[i] += p[i];
  return w;
}

IVect
operator- (const IVect& s,
           const IVect& p)
{
  IVect w=s;
  for(int i=0; i<w.size(); i++) w[i] -= p[i];
  return w;
}

IVect
operator* (const IVect& s,
           const IVect& p)
{
  IVect w=s;
  for(int i=0; i<w.size(); i++) w[i] *= p[i];
  return w;
}

std::ostream&
operator<< (std::ostream& ostr, const IVect& p)
{
  if(p.size() <= 0) return ostr;

  ostr << "(" ;
  for(int i=0; i<p.size()-1; i++) ostr << p[i] << ", ";
  ostr <<  p[p.size()-1] <<  ")";
  // if (ostr.fail()) MayDay::Error("operator<<(ostream&,IVect&) failed");
  return ostr;
}

#include "NamespaceFooter.H"
#include "UsingNamespace.H"
#include "BaseNamespaceHeader.H"

template < >
int linearSize(const IVect& vindex)
{
  return sizeof(IVect);
}

#include "BaseNamespaceFooter.H"
