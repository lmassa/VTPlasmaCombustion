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
#include "RVect.H"
#include "Misc.H"
using std::ostream;

#include "NamespaceHeader.H"
#define CHOFFSET(object, member) (int)((char*)&(object.member) - (char*)&object)
using std::ws;

RVect anRVect;
size_t RVect::io_offset = CHOFFSET(anRVect, v);


const Real*
RVect::dataPtr() const
{
  return v.data();
}

Real*
RVect::dataPtr()
{
  //not sure this works for this version of C++
  return v.data();
}



RVect::RVect ():Vector<Real>()
{

}

RVect::RVect (unsigned int isize):Vector<Real>(isize, 0e0)
{

}


RVect::RVect (unsigned int isize, unsigned int jone)
{
  v.resize(isize, 0);
  v[jone]=1e0;
}

RVect::RVect (unsigned int isize, unsigned int jone, const Real& s)
{
  v.resize(isize, 0);
  v[jone]=s;
}


RVect::RVect (const Real& d1, const Real& d2, unsigned int isize)
{
  v.resize(isize);
  if(isize == 1)
    {
      v[0]=d1;
      return;
    }
  Real delta = (d2-d1)/(isize-1);
  for(int i=0;i<isize;i++) v[i]=d1+i*delta;
}


RVect::RVect (unsigned int isize, const Real& s)
{
  v.resize(isize, s);
}

RVect::RVect (const Vector<Real>& rhs)
{
  v = rhs.constStdVector();
}


RVect::RVect (const IVect& rhs)
{
  
  v.resize(rhs.size(),0e0);
  for (int i=0;i<rhs.size();i++)v[i]=rhs[i];
}


RVect::RVect (const std::vector<Real>& invec)
{
  v = invec;
}


RVect::RVect (unsigned int a_nbr, unsigned int a_blkRowSize, const RVect& a_blkMat, const RVect& a_bCoe)
{
  int nblocks = a_bCoe.size(); //total number of blocks (row X col)
  int nbc = nblocks/a_nbr;
  int bsize = a_blkMat.size();
  int isize = a_nbr*nbc*bsize;
  int blkColSize = bsize/a_blkRowSize;
  int ncolc = blkColSize*nbc;
  int nrowc = a_blkRowSize*a_nbr;
  v.resize(isize, 0);
  for(int ib=0;ib<a_nbr;ib++)
    {
      for(int ir=0;ir<a_blkRowSize;ir++)
	{
	  for(int jb=0;jb<nbc;jb++)
	    {
	      for(int jc=0;jc<blkColSize;jc++)
		{
		  int irowc = ir + ib*a_blkRowSize;
		  int icolc = jc + jb*blkColSize;
		  v[irowc*ncolc + icolc] = a_blkMat[ir*blkColSize + jc] * a_bCoe[ib*a_nbr + jb];
		}
	    }
	}
    }
}





Real RVect::dotProduct(const RVect& a_rhs) const
{

  Real s=0;
  for(int i=0; i<v.size(); i++) s += v[i]*a_rhs.v[i];
  return s;
}


RVect RVect::matProduct(const RVect& a_B, const int& a_nrows, const int& a_ncols, const int& a_ncolA) const
{

  RVect Cmat ( a_nrows*a_ncols );
  Real C;
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


RVect RVect::pwr(const Real& a_s) const
{
  RVect C(v.size());
  for(int i=0; i<v.size(); i++) C[i] = pow(v[i],a_s);
  return C;
}

//overload for square matrices
RVect RVect::matProduct(const RVect& a_B, const int& a_n) const
{

  RVect Cmat ( a_n*a_n );
  Real C;
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


//note matmat is the same as matProduct but return *this

RVect& RVect::matMat(const RVect& a_B, const int& a_nrows, const int& a_ncols, const int& a_ncolA)
{

  RVect Cmat ( a_nrows*a_ncols );
  Real C;
  for(int i=0; i<a_nrows; i++) 
    for(int j=0; j<a_ncols; j++)
      {
	C =  v[a_ncolA*i]*a_B[j];
	for(int k=1; k<a_ncolA; k++) 
	  C =  C + v[a_ncolA*i + k]*a_B[a_ncols*k + j];
	Cmat[a_ncols*i+j] =  C;
      }
  v= Cmat.v;
  return *this;
}

//overload for square matrices
RVect& RVect::matMat(const RVect& a_B, const int& a_n)
{

  RVect Cmat ( a_n*a_n );
  Real C;
  for(int i=0; i<a_n; i++) 
    for(int j=0; j<a_n; j++)
      {
	C =  v[a_n*i]*a_B[j];
	for(int k=1; k<a_n; k++) 
	  C =  C + v[a_n*i + k]*a_B[a_n*k + j];
	Cmat[a_n*i+j] =  C;
      }
  v= Cmat.v;
  return *this;
}

RVect RVect::Identity(const int& a_nrows, const int& a_ncols)
{
  RVect Cmat ( a_nrows*a_ncols );
  for(int i=0; i<a_nrows; i++) Cmat[a_ncols*i+i] =  1.0;
  return Cmat;
}

RVect RVect::Identity(const int& a_n)
{
  RVect Cmat ( a_n*a_n );
  for(int i=0; i<a_n; i++) Cmat[a_n*i+i] =  1.0;
  return Cmat;
}
RVect RVect::Identity(const int& a_nrows, const int& a_ncols, const Real& s)
{
  RVect Cmat ( a_nrows*a_ncols );
  for(int i=0; i<a_nrows; i++) Cmat[a_ncols*i+i] =  s;
  return Cmat;
}


//overload for square matrices
RVect& RVect::setRow(const int& a_irow, const RVect& a_B)
{
  int bSize = a_B.size();
  int mysize = v.size();
  
  int ncols = bSize;
  int nrows = mysize/bSize;

  for(int i=0; i<ncols; i++) 
    v[a_irow*ncols + i] = a_B[i];

  return *this;
}


int
RVect::findLE(const Real a_r) const
{
  if(size() == 0) return -1;
  vector<Real>& w = *(const_cast<vector<Real>* >(&v));
  std::vector<Real>::iterator it = std::find_if (w.begin(), w.end(), index_greater(a_r));
  if(it>=w.end()) return (w.size()-1);
  else if (it == w.begin()) return (0);
  else
    return (it-w.begin()-1);
}


bool
RVect::operator== (const RVect& p) const
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
RVect::operator!= (const RVect& p) const
{
  bool retval = true;
  for(int i=0; i<v.size(); i++) 
    {
      retval = retval && (v[i] != p.v[i]);
      if(!retval) break;
    }
  return retval;
}

RVect&
RVect::operator+= (Real s)
{
  for(int i=0; i<v.size(); i++) v[i] += s;
  return *this;
}

RVect&
RVect::operator+= (const RVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] += p[i];
  return *this;
}

RVect&
RVect::operator*= (Real s)
{
  for(int i=0; i<v.size(); i++) v[i] *= s;
  return *this;
}

RVect&
RVect::operator*= (const RVect &p)
{
  for(int i=0; i<v.size(); i++) v[i] *= p[i];
  return *this;
}

//XXXRVect
//XXXRVect::operator* (const RVect& p) const
//XXX{
//XXX  RVect v(D_DECL6(vect[0]*p[0], vect[1]*p[1], vect[2]*p[2]));
//XXX  return v;
//XXX}

RVect
RVect::operator* (Real s) const
{
  RVect w=v;
  for(int i=0; i<w.size(); i++) w[i] *= s;
  return w;
}

RVect
RVect::operator- (Real s) const
{
  RVect w=v;
  for(int i=0; i<w.size(); i++) w[i] -= s;
  return w;
}


RVect
RVect::operator+ (Real s) const
{
  RVect w=v;
  for(int i=0; i<w.size(); i++) w[i] += s;
  return w;
}

RVect&
RVect::operator/= (Real s)
{
  for(int i=0; i<v.size(); i++) v[i] /= s;
  return *this;
}

RVect&
RVect::operator/= (const RVect& p)
{
  for(int i=0; i<v.size(); i++) v[i] /= p[i];
  return *this;
}

//XXXRVect
//XXXRVect::operator/ (const RVect & p) const
//XXX{
//XXX  RVect result( D_DECL6( vect[0] / p[0], vect[1] / p[1], vect[2] / p[2] ) );
//XXX  return result ;
//XXX}

RVect
RVect::operator/ (Real s) const
{
  RVect w=v;
  for(int i=0; i<w.size(); i++) w[i] /= s;
  return w;
}


int
RVect::minDir(const bool& a_doAbs) const
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
RVect::maxDir(const bool& a_doAbs) const
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


Real
RVect::maxAbs() const
{

  vector<Real>::const_iterator mx, mn;  
  mx = max_element(v.begin(), v.end());
  mn = min_element(v.begin(), v.end());
  return(Max(Abs(*mx),Abs(*mn)));
}


int
RVect::minloc() const
{
  vector<Real>::const_iterator mn;  
  mn = min_element(v.begin(), v.end());
  return(mn-v.begin());
}

Real
RVect::min() const
{
  if (size() == 0) return(0);
  vector<Real>::const_iterator mn;  
  mn = min_element(v.begin(), v.end());
  return(*mn);
}

Real RVect::max() const
{
  if (size() == 0) return(0);
  vector<Real>::const_iterator mx;  
  mx = max_element(v.begin(), v.end());
  return(*mx);
}


RVect
operator/ (Real            s,
           const RVect& p)
{
  RVect w=p;
  for(int i=0; i<w.size(); i++) w[i] = s/w[i];
  return w;
}
RVect
operator+ (Real            s,
           const RVect& p)
{
  RVect w=p;
  for(int i=0; i<w.size(); i++) w[i] += s;
  return w;
}

RVect
operator- (Real            s,
           const RVect& p)
{
  RVect w=p;
  for(int i=0; i<w.size(); i++) w[i] = s-w[i];
  return w;
}

RVect
operator* (Real            s,
           const RVect& p)
{
  RVect w=p;
  for(int i=0; i<w.size(); i++) w[i] *= s;
  return w;
}

RVect
operator/ (const RVect& s,
           const RVect& p)
{
  RVect w=s;
  for(int i=0; i<w.size(); i++) w[i] /= p[i];
  return w;
}

RVect
operator+ (const RVect& s,
           const RVect& p)
{
  RVect w=s;
  for(int i=0; i<w.size(); i++) w[i] += p[i];
  return w;
}

RVect
operator- (const RVect& s,
           const RVect& p)
{
  RVect w=s;
  for(int i=0; i<w.size(); i++) w[i] -= p[i];
  return w;
}

RVect
operator* (const RVect& s,
           const RVect& p)
{
  RVect w=s;
  for(int i=0; i<w.size(); i++) w[i] *= p[i];
  return w;
}

std::ostream&
operator<< (std::ostream& ostr, const RVect& p)
{
  if(p.size() <= 0) return ostr;

  if(false)
    {
      ostr << "(" ;
      for(int i=0; i<p.size()-1; i++) ostr << p[i] << ", ";
      ostr <<  p[p.size()-1] <<  ")";
    }
  else
    {
      for(int i=0; i<p.size()-1; i++) ostr << p[i] << ", ";
      ostr <<  p[p.size()-1];
    }
  // if (ostr.fail()) MayDay::Error("operator<<(ostream&,RVect&) failed");
  return ostr;
}

#include "NamespaceFooter.H"
#include "UsingNamespace.H"
#include "BaseNamespaceHeader.H"

template < >
int linearSize(const RVect& vindex)
{
  return sizeof(RVect);
}

#include "BaseNamespaceFooter.H"
