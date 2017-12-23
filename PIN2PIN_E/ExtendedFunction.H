#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTENDEDFUNCTION_H_
#define _EXTENDEDFUNCTION_H_

#include "ParmParse.H"
#include "parstream.H"

//Physics classes
#include "LevelDataOps.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "Box.H"
#include "RVect.H"
#include "IVect.H"
#include "BaseEBFaceFAB.H"
#include "BaseEBCellFAB.H"
#include "EBLGIntegrator.H" //for CNUM

//addtional classes to test PETSC
#include <petscmat.h>
#include <cstdarg>

#include "Tuple.H"
#include <deque>

#include "NamespaceHeader.H"

enum EFdata_t {FACE, CELL, FLUX, UNDETERMINED};
template<typename T, size_t N>
T * end(T (&ra)[N]) {
    return ra + N;
}
template<typename T, size_t N>
int numel(T (&ra)[N]) {
    return N;
}

struct differential
{
  deque<Mat> grad;
  vector<string> dependent;
  vector<int> colSize;
  vector<EFdata_t> colType;
  vector<bool> isMatCreated;
  vector<BaseEBCellFAB<int>* > ColsfromCELLtoIb;
  vector<BaseEBFaceFAB<int>* > ColsfromFACEtoIb;
  vector<IVect> nnzs;

  string name;
  int rowSize;
  EFdata_t rowType;
  BaseEBCellFAB<int>* fromCELLtoIb;
  BaseEBFaceFAB<int>* fromFACEtoIb;
  differential() : rowType(UNDETERMINED), fromCELLtoIb(NULL), fromFACEtoIb(NULL) {};
  int size ()
  {
    return grad.size();
  }
};

//The Extended Function Class
template<class T>
class ExtendedFunction
{
public:
  //constructors
  ExtendedFunction():m_val(NULL), m_diff(NULL),m_diffAllocated(false),m_valAllocated(false),m_blockSize(CNUM)
  {
  }
  ExtendedFunction(T* a_val):m_val(a_val), m_diff(NULL),m_diffAllocated(false),m_valAllocated(false),m_blockSize(CNUM)
  {
  }
  ExtendedFunction(differential* a_diff):m_val(NULL), m_diff(a_diff),m_diffAllocated(false),m_valAllocated(false),m_blockSize(CNUM)
  {
  }
  ExtendedFunction(T* a_val, differential* a_diff):m_val(a_val), m_diff(a_diff),m_diffAllocated(false),m_valAllocated(false),m_blockSize(CNUM)
  {
  }
  ExtendedFunction(T* a_val, const string& a_name, const vector <string>& a_names):m_val(a_val),m_blockSize(CNUM),m_valAllocated(false)
  {
    m_diff = new differential;
    m_diffAllocated = true;
    m_diff->name = a_name;
    m_diff->dependent = a_names;
    int nele = a_names.size();
    m_diff->grad.resize(nele);
    m_diff->isMatCreated.resize(nele,false);
  }

  ExtendedFunction(T* a_val, const string& a_name, const char* cnames[],const int& a_N):m_val(a_val),m_blockSize(CNUM),m_valAllocated(false)
  {
    std::vector<std::string> vnames(cnames, cnames + a_N);
    m_diff = new differential;
    m_diffAllocated = true;
    m_diff->name = a_name;
    m_diff->dependent = vnames;
    int nele = vnames.size();
    m_diff->grad.resize(nele);
    m_diff->isMatCreated.resize(nele,false);
  }
  ExtendedFunction(T* a_val, int a_newComps, const string& a_name, const char* cnames[],const int& a_N):m_val(a_val),m_blockSize(CNUM),m_valAllocated(false)
  {
    std::vector<std::string> vnames(cnames, cnames + a_N);
    m_diff = new differential;
    m_diffAllocated = true;
    m_diff->name = a_name;
    m_diff->dependent = vnames;
    int nele = vnames.size();
    m_diff->grad.resize(nele);
    m_diff->isMatCreated.resize(nele,false);

    m_FAB.resize(nele,NULL);
    for (int k=0;k<m_FAB.size();k++)
      {
	m_FAB = new T( &a_val, a_newComps);
      }
    m_FABAllocated = true;
  }
  //particularized constructors for different DataTypes
  //for the EBCellFAB
  ExtendedFunction(const EBISBox& a_ebisBox, const Box& a_region, int a_nComp, const string& a_name):m_val(NULL), m_diff(NULL),m_diffAllocated(false),m_valAllocated(false),m_blockSize(CNUM)
  {
    setDefaultvalues; //yet to be written
    m_val = new T(a_ebisBox,a_region,a_nComp);
    m_valAllocated=true;
    m_diff = new differential;
    m_diffAllocated = true;
    m_diff->name = a_name;
    
  }
  ~ExtendedFunction()
  {
    if(m_diffAllocated && m_diff!= NULL) 
      {
	delete m_diff;
	m_diff=NULL;
      }
    if(m_FABAllocated) 
      {
	for (int k=0;k<m_FAB.size();k++)
	  {
	    if(m_FAB[k] != NULL)
	      {
		delete m_FAB[k];
		m_FAB[k]=NULL;
	      }
	  }
      }
    if(m_valAllocated && m_val!= NULL) 
      {
	delete m_val;
	m_val=NULL;
      }
  }

  //some get functions
  T& getVal();
  BaseEBCellFAB<int> & getCELLmap();
  BaseEBFaceFAB<int> & getFACEmap();
  BaseEBFaceFAB<int> & getFACEmap(const string& a_str);
  Mat& getDeriv(const int& a_i);
  Mat& getDeriv(const string& a_str);
  int getCol(const string& a_str);
  differential& getDifferential();
  differential* getDiffPtr();
  string& getDependentName(const int& a_i);
  BaseFab<Real>& getSingleValuedFAB()
  {
    return(getVal().getSingleValuedFAB());
  }
  BaseFab<Real>& getSingleValuedFAB(const int& a_n); // expanded it down after the FAB utilities

  BaseEBCellFAB<int>& getCellMap(){return *(m_diff->fromCELLtoIb);};
  IVect& getNNZ(const string& a_str);
  int& getCellNum(){return *(m_diff->rowSize);};

  const Real& operator() (const VolIndex& a_ndin, int a_nVarLoc) const;
  Real& operator() (const VolIndex& a_ndin,int  a_nVarLoc);
  const Real& operator() (const VolIndex& a_ndin, int a_nVarLoc, int a_nDeriv) const;
  Real& operator() (const VolIndex& a_ndin,int  a_nVarLoc, int a_nDeriv);

  int size ()
  {
    if(m_diff == NULL) 
      return 0;
    else			 
      return m_diff->size();
  };


  string& name()
  {
    if(m_diff == NULL) 
      return string();
    else			 
      return m_diff->name;
  };

  bool isMatCreated(const int& a_i) const
  {
    return m_diff->isMatCreated[a_i];
  };


  void toggleCreated(const int& a_i)
  {
    m_diff->isMatCreated[a_i] = true;
  };



  void toggleCreated(const string& a_str);

  void allocateBlockMatrices(const int&a_nblk, const int&a_nrows, const int&a_ncols, const int&a_NZBperRow)
  {
    
  for(int k=0;k<size();k++)
    {
      if(!isMatCreated(k))
	{
	  MatCreateSeqBAIJ(PETSC_COMM_SELF,a_nblk, a_nrows, a_ncols, a_NZBperRow,NULL, &getDeriv(k));
	  toggleCreated(k);
	}
	
    }
  };

  void allocateBlockMatrices(const int&a_NZBperRow)
  {
    
  for(int k=0;k<size();k++)
    {
      if(isMatCreated(k))
	{
	  MatDestroy(&getDeriv(k));
	  m_diff->isMatCreated[k]=false;
	}
      MatCreateSeqBAIJ(PETSC_COMM_SELF,m_blockSize, m_numVofs, m_diff->colSize[k], a_NZBperRow,NULL, &getDeriv(k));
      toggleCreated(k);
	
    }
  };

  void allocateBlockMatrices()
  {
    
  for(int k=0;k<size();k++)
    {
      if(isMatCreated(k))
	{
	  MatDestroy(&getDeriv(k));
	  m_diff->isMatCreated[k]=false;
	}
      MatCreateSeqBAIJ(PETSC_COMM_SELF,m_blockSize, m_numVofs, m_diff->colSize[k], PETSC_DEFAULT,m_diff->nnzs[k].dataPtr(), &getDeriv(k));
      toggleCreated(k);
	
    }
  };

  void assembleMatrices(const MatAssemblyType type=MAT_FINAL_ASSEMBLY)
  {
    
  for(int k=0;k<size();k++)
    {
      PetscBool isAssembled;
      MatAssembled(getDeriv(k),&isAssembled);
      if(!isAssembled)
	{
	  int ierr = MatAssemblyBegin(getDeriv(k),type);//CHKERRQ(ierr);
	  ierr = MatAssemblyEnd      (getDeriv(k),type);//CHKERRQ(ierr);
	}
    }
  }

  //FAB utilities
  void setupFABS()
  {
    int nele = size();
    int nComps = m_blockSize*m_blockSize;
    m_FAB.resize(nele,NULL);
    for (int k=0;k<m_FAB.size();k++)
      {
	m_FAB[k] = new T( *m_val, nComps);
	m_FAB[k]->setVal(0e0);
      }
    m_FABAllocated = true;
  }
  T& getFAB(const string& a_str);
  T& getFAB(const int& a_n) {return *m_FAB[a_n];};
  void convertFaceFABtoMAT();
  void convertIVFABtoMAT(); //case in which <T> is an IVFAB
  void convertCellFABtoMAT();//case in which <T> is a CellFAB


  //mapping utilities
  int getIb(const VolIndex a_vof)
  {
    int ival = (*m_diff->fromCELLtoIb)(a_vof,0);
    return(ival);
  }
  int getIb(const FaceIndex a_face)
  {
    int ival = (*m_diff->fromFACEtoIb)(a_face,0);
    return(ival);
  }

  //operation routines
  void registerEF(const int num, ...);
  void chainEF(const int num, ...);
  void setVal(const Real& a_val);
  ExtendedFunction<T>&   copy( const ExtendedFunction<T>& source);
  ExtendedFunction<T>&   incr( const ExtendedFunction<T>& source, const Box & a_region, const Real&  a_scale);
  int& nComp() const{ return m_val->nComp()};
  int& direction() const{ return m_val->direction()};
  Box& getRegion() const{ return m_val->getRegion()};
  

  T* m_val;
  bool m_valAllocated;
  differential* m_diff;
  bool m_diffAllocated;
  vector<T*> m_FAB;
  bool m_FABAllocated;
  //for listing purposes
  int m_numVofs; // this shoul be equal to m_diff->rowSize
  int m_blockSize;
};
template<class T> T& ExtendedFunction<T>::getVal() {return *m_val};
template<class T> BaseEBCellFAB<int> & ExtendedFunction<T>::getCELLmap() {return *(m_diff->fromCELLtoIb);};
template<class T> BaseEBFaceFAB<int> & ExtendedFunction<T>::getFACEmap() {return *(m_diff->fromFACEtoIb);};
template<class T> Mat& ExtendedFunction<T>::getDeriv(const int& a_i) {return m_diff->grad[a_i];};
template<class T> Mat& ExtendedFunction<T>::getDeriv(const string& a_str) 
{
  int num = size();
  string PnameF = a_str;
  std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
  PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
  for  (int j=0;j<num;j++)
    {
      string Pname = m_diff->dependent[j];
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      
      if(Pname == PnameF)
	return m_diff->grad[j];
    }
  //failure
   MayDay::Error("ExtendedFunction::getDeriv: !found name");
};
template<class T> T& ExtendedFunction<T>::getFAB(const string& a_str) 
{
  int num = size();
  string PnameF = a_str;
  std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
  PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
  for  (int j=0;j<num;j++)
    {
      string Pname = m_diff->dependent[j];
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      
      if(Pname == PnameF)
	return *m_FAB[j];
    }
  //failure
   MayDay::Error("ExtendedFunction::getDeriv: !found name");
};

template<class T> BaseFab<Real>& ExtendedFunction<T>::getSingleValuedFAB(const int& a_n) 
{

  if(!m_FABAllocated)
    MayDay::Error("ExtendedFunction::getSingleValuedFAB: FAB were not allocated");

  return (*m_FAB)[a_n].getSingleValuedFAB();
};
template<class T> int ExtendedFunction<T>::getCol(const string& a_str) 
{
  int num = size();
  string PnameF = a_str;
  std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
  PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
  for  (int j=0;j<num;j++)
    {
      string Pname = m_diff->dependent[j];
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      
      if(Pname == PnameF)
	return m_diff->colSize[j];
    }
  //failure
   MayDay::Error("ExtendedFunction::getDeriv: !found name");
};
template<class T> IVect& ExtendedFunction<T>::getNNZ(const string& a_str) 
{
  int num = size();
  string PnameF = a_str;
  std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
  PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
  for  (int j=0;j<num;j++)
    {
      string Pname = m_diff->dependent[j];
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      
      if(Pname == PnameF)
	return m_diff->nnzs[j];
    }
  //failure
   MayDay::Error("ExtendedFunction::getNNZ: !found name");
};
template<class T> BaseEBFaceFAB<int>& ExtendedFunction<T>::getFACEmap(const string& a_str) 
{
  int num = size();
  string PnameF = a_str;
  std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
  PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
  for  (int j=0;j<num;j++)
    {
      string Pname = m_diff->dependent[j];
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      
      if(Pname == PnameF)
	return *(m_diff->ColsfromFACEtoIb[j]);
    }
  //failure
  MayDay::Error("ExtendedFunction::getFACEmap: !found name");
};

template<class T> differential& ExtendedFunction<T>::getDifferential() {return *m_diff;};
template<class T> differential* ExtendedFunction<T>::getDiffPtr() {return m_diff;};
template<class T> string& ExtendedFunction<T>::getDependentName(const int& a_i) {return m_diff->dependent[a_i];};

//operators
template<class T> const Real& ExtendedFunction<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc) const
{
  return(*m_val(a_ndin,a_nVarLoc));
}
template<class T> Real& ExtendedFunction<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc)
{
  return(*m_val(a_ndin,a_nVarLoc));
}
template<class T> const Real& ExtendedFunction<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc,int  a_nDeriv) const
{
  CH_assert(m_FABAllocated);
  CH_assert(a_nDeriv<m_FAB.size());
  return((*m_FAB[a_nDeriv])(a_ndin,a_nVarLoc));
}
template<class T> Real& ExtendedFunction<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc,int  a_nDeriv)
{
  CH_assert(m_FABAllocated);
  CH_assert(a_nDeriv<m_FAB.size());
  return((*m_FAB[a_nDeriv])(a_ndin,a_nVarLoc));
}


template<class T> void ExtendedFunction<T>::toggleCreated(const string& a_str) 
{
  int num = size();
  string PnameF = a_str;
  std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
  PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
  for  (int j=0;j<num;j++)
    {
      string Pname = name();
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      
      if(Pname == PnameF)
	{
	 m_diff->isMatCreated[j] = true;
	 return;
	}
    }

};

template<class T> void ExtendedFunction<T>::convertFaceFABtoMAT() 
{
  if(!m_FABAllocated) MayDay::Error("ExtendedFunction::convertFaceFABtoMAT: !m_FABAllocated");
    
  int blkSize = m_blockSize;
  int nvars = blkSize;
  PetscScalar *data = new PetscScalar[nvars*nvars];
  //the MAts shouldnt have been alocated at this point
  allocateBlockMatrices(1);
  for (int k=0;k<m_FAB.size();k++)
    {
      if(m_FAB[k] != NULL)
	{
	  Mat& PTSC_mat = m_diff->grad[k];
	  EBFaceFAB& coeff = *m_FAB[k];
	  
	  int dir = coeff.direction();
	  const EBGraph& graph = coeff.getEBISBox().getEBGraph();
	  Box cellBox = coeff.getCellRegion();
	  //fill the matrix
	  IntVectSet ivs(cellBox);
	  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
	  for (FaceIterator faceit(ivs, graph, dir, stopCrit); faceit.ok(); ++faceit)
	    {
	      const FaceIndex& face = faceit(); //would be nice to add a check it this is in the faces, but I wasnt able to think of a good one
	      if(! (m_diff->rowType == FACE) && !(m_diff->colType[k] == FACE))
		MayDay::Error("ExtendedFunction::convertFaceFABtoMAT: restricting it to face-face interactions");
	      int ibRow = (*m_diff->fromFACEtoIb)(face,0);
	      int ibCol = (*m_diff->ColsfromFACEtoIb[k])(face,0);
	      Real maxdata = 0e0;
	      for (int ivar = 0; ivar < nvars*nvars; ivar++)
		{
		  data[ivar] = coeff(face, ivar);
		  maxdata = Max(maxdata, abs(data[ivar]));
		}
	      if(maxdata > 1e-16)
		MatSetValuesBlocked(PTSC_mat, 1, &ibRow, 1, &ibCol, data, INSERT_VALUES );
	    }
	  int ierr = MatAssemblyBegin(PTSC_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(PTSC_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
	  delete m_FAB[k];
	  m_FAB[k]=NULL;
	}
    }

  m_FABAllocated = false;
  delete [] data;
 
  
};



template<class T> void ExtendedFunction<T>::convertIVFABtoMAT() 
{
  if(!m_FABAllocated) MayDay::Error("ExtendedFunction::convertIVFABtoMAT: !m_FABAllocated");
    
  int blkSize = m_blockSize;
  int nvars = blkSize;
  PetscScalar *data = new PetscScalar[nvars*nvars];
  //the MAts shouldnt have been alocated at this point
  allocateBlockMatrices(1);
  for (int k=0;k<m_FAB.size();k++)
    {
      if(m_FAB[k] != NULL)
	{
	  Mat& PTSC_mat = m_diff->grad[k];
	  T& coeff = *m_FAB[k];
	  
	  //fill the matrix
	  const IntVectSet& ivs = coeff.getIVS();

	  for (VoFIterator vofit(ivs, coeff.getEBGraph());vofit.ok(); ++vofit)
	    {
	      const VolIndex& cell = vofit();
	      if(! (m_diff->rowType == CELL) && !(m_diff->colType[k] == CELL))
		MayDay::Error("ExtendedFunction::convertCellFABtoMAT: restricting it to cell-cell interactions");
	      int ibRow = *(m_diff->fromCELLtoIb)(cell,0);
	      int ibCol = *(m_diff->ColsfromCELLtoIb[k])(cell,0);
	      Real maxdata = 0e0;
	      for (int ivar = 0; ivar < nvars*nvars; ivar++)
		{
		  data[ivar] = coeff(cell, ivar);
		  maxdata = Max(maxdata, abs(data[ivar]));
		}
	      if(maxdata > 1e-16)
		MatSetValuesBlocked(PTSC_mat, 1, &ibRow, 1, &ibCol, data, INSERT_VALUES );
	    }
	  int ierr = MatAssemblyBegin(PTSC_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(PTSC_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
	  delete m_FAB[k];
	  m_FAB[k]=NULL;
	}
    }

  m_FABAllocated = false;
  delete [] data;
 
  
};


template<class T> void ExtendedFunction<T>::convertCellFABtoMAT() 
{
  if(!m_FABAllocated) MayDay::Error("ExtendedFunction::convertCellFABtoMAT: !m_FABAllocated");
    
  int blkSize = m_blockSize;
  int nvars = blkSize;
  PetscScalar *data = new PetscScalar[nvars*nvars];
  //the MAts shouldnt have been alocated at this point
  allocateBlockMatrices(1);
  for (int k=0;k<m_FAB.size();k++)
    {
      if(m_FAB[k] != NULL)
	{
	  Mat& PTSC_mat = m_diff->grad[k];
	  T& coeff = *m_FAB[k];
	  
	  //fill the matrix
	  Box wholeBox = coeff.getRegion();
	  IntVectSet ivs(wholeBox);

	  for (VoFIterator vofit(ivs, coeff.getEBGraph());vofit.ok(); ++vofit)
	    {
	      const VolIndex& cell = vofit();
	      if(! (m_diff->rowType == CELL) && !(m_diff->colType[k] == CELL))
		MayDay::Error("ExtendedFunction::convertCellFABtoMAT: restricting it to cell-cell interactions");
	      int ibRow = *(m_diff->fromCELLtoIb)(cell,0);
	      int ibCol = *(m_diff->ColsfromCELLtoIb[k])(cell,0);
	      Real maxdata = 0e0;
	      for (int ivar = 0; ivar < nvars*nvars; ivar++)
		{
		  data[ivar] = coeff(cell, ivar);
		  maxdata = Max(maxdata, abs(data[ivar]));
		}
	      if(maxdata > 1e-16)
		MatSetValuesBlocked(PTSC_mat, 1, &ibRow, 1, &ibCol, data, INSERT_VALUES );
	    }
	  int ierr = MatAssemblyBegin(PTSC_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(PTSC_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
	  delete m_FAB[k];
	  m_FAB[k]=NULL;
	}
    }

  m_FABAllocated = false;
  delete [] data;
 
  
};


template<class T> void ExtendedFunction<T>::registerEF(const int num, ...)
{
  va_list args;
  
  if(m_diff == NULL)
    { m_diff = new differential;
      m_diffAllocated = true;
    }
  else
    {
      for(int k=0;k<size();k++)
	{
	  if(isMatCreated(k))
	    {
	      MatDestroy(&getDeriv(k));
	    }
	}
    }
  
  deque<Mat>& myGrad = m_diff->grad;
  vector<string>& myDeps = m_diff->dependent;
  vector<int>& myCols = m_diff->colSize;
  vector<EFdata_t>& myColTypes = m_diff->colType;
  vector<BaseEBCellFAB<int>* >& myColsfromCELLtoIb = m_diff->ColsfromCELLtoIb;
  vector<BaseEBCellFAB<int>* >& myColsfromFACEtoIb = m_diff->ColsfromFACEtoIb;  
  vector<IVect>& myNNZ = m_diff->nnzs;
  
  int myrowSize = m_numVofs;
  
  myDeps.resize(num);
  myGrad.resize(num);
  myCols.resize(num);
  myColsfromCELLtoIb.resize(num,NULL);
  myColsfromFACEtoIb.resize(num,NULL);
  myColTypes.resize(num);
  myNNZ.resize(num);
  m_diff->isMatCreated.resize(num,false);
  m_diff->rowSize = myrowSize;
  

  
  va_start(args, num);
  for  (int j=0;j<num;j++)
    {
      differential* var = va_arg(args, differential*);
      
      string Pname = var->name;
      std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
      Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );
      getDependentName(j) = Pname;
      myCols[j] = var->rowSize;
      myColsfromCELLtoIb[j] = var->fromCELLtoIb;
      myColsfromFACEtoIb[j] = var->fromFACEtoIb;
      myColTypes[j] = var->colType;
      myNNZ[j] = IVect(m_numVofs,0,0);
    }
  va_end(args);
}


template<class T> void ExtendedFunction<T>::chainEF(const int num, ...)
{
 va_list args;  
   
    int origSize =size();

    IVect eraseList(0);

    deque<Mat>& myGrad = m_diff->grad;
    vector<string>& myDeps = m_diff->dependent;
    
    for (int i=0;i<origSize;i++)
      {
	string PnameF = getDependentName(i);
	std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
	PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );

	Mat & dF_dP = getDeriv(i);
	MatConvert(dF_dP,  MATSEQAIJ, MAT_REUSE_MATRIX, &dF_dP);

	va_start(args, num);
	for  (int j=0;j<num;j++)
	  {
	    differential* var = va_arg(args, differential*);
	    string Pname = var->name;
	    std::transform(Pname.begin(), Pname.end(), Pname.begin(), ::tolower);
	    Pname.erase( std::remove_if( Pname.begin(), Pname.end(), ::isspace ), Pname.end() );

	    if(Pname == PnameF)
	      {
		for  (int k=0;k<var->size();k++)
		  {
		    string Xname = var->dependent[k];
		    std::transform(Xname.begin(), Xname.end(), Xname.begin(), ::tolower);
		    Xname.erase( std::remove_if( Xname.begin(), Xname.end(), ::isspace ), Xname.end() );		    
		    
		    Mat & dP_dX = var->grad[k];
		    MatConvert(dP_dX,  MATSEQAIJ, MAT_REUSE_MATRIX, &dP_dX);
                    
		    bool found = false;
		    int mfound = -1;
		    for  (int m=0;m< size();m++)
		      {
			if(m==i) continue;
			string XnameF = getDependentName(m);
			std::transform(XnameF.begin(), XnameF.end(), XnameF.begin(), ::tolower);
			XnameF.erase( std::remove_if( XnameF.begin(), XnameF.end(), ::isspace ), XnameF.end() );
			found = XnameF == Xname;
			if(found) {mfound=m; break;}
		      }
		    if(!found)
		      {
			mfound = size();
			myGrad.resize(mfound+1);
			myDeps.push_back(var->dependent[k]);
			Mat& dF_dX =  getDeriv(mfound);
			MatMatMult(dF_dP, dP_dX, MAT_INITIAL_MATRIX, 1., &dF_dX);
		      }
		    else
		      {
			Mat dF_dX2;
			MatMatMult(dF_dP, dP_dX, MAT_INITIAL_MATRIX, 1., &dF_dX2);
			Mat& dF_dX =  getDeriv(mfound);
			MatConvert(dF_dX,  MATSEQAIJ, MAT_REUSE_MATRIX, &dF_dX);
			MatAXPY(dF_dX,1.0, dF_dX2, DIFFERENT_NONZERO_PATTERN);
			MatDestroy(&dF_dX2);
		      }
		  }
		MatDestroy(&dF_dP);
		eraseList.push_back(i);
		break;
	      }
	  }
	va_end(args);
      }

    //fix-up the My(this) variables
    for (int i=0;i<eraseList.size();i++)
      {
	myGrad.erase(myGrad.begin()+ eraseList[i]);
	myDeps.erase(myDeps.begin()+ eraseList[i]);
	eraseList-=1;
	//uniquify dependent names
	std::vector<string>::iterator it;
	it = std::unique (myDeps.begin(), myDeps.end()); 
	myDeps.resize( std::distance(myDeps.begin(),it) );
      }
}


template<class T> void ExtendedFunction<T>::setVal(const Real& a_val)
{
  //reset an Extended Function
  m_val->setval(a_val);
  if(m_diff != NULL)
    {
      for(int k=0;k<size();k++)
	{
	  if(isMatCreated(k))
	    {
	      MatDestroy(&getDeriv(k));
	    }
	}
      m_diff->dependent.resize(0);
      m_diff->isMatCreated.resize(0);
      if(m_diffAllocated) 
	{
	  delete m_diff;
	  m_diff=NULL;
	  m_diffAllocated = false;
	}
    }
  if(m_FAB.size() > 0)
    {
      if(m_FABAllocated) 
	{
	  for (int k=0;k<m_FAB.size();k++)
	    {
	      if(m_FAB[k] != NULL)
		{
		  delete m_FAB[k];
		  m_FAB[k]=NULL;
		}
	    }
	}
      m_FAB.resize(0);
      m_FABAllocated = false;
    }
}

template <class T>
ExtendedFunction<T>&
ExtendedFunction<T>::copy(const ExtendedFunction<T>& a_src)
{
  m_val.copy(a_src.getVal());
  m_valAllocated = false;

  if(m_diff != NULL)
    {
      if(m_diffAllocated)
	delete m_diff;
      m_diff = NULL;
    }  
  m_diff = new differential;
  m_diffAllocated = true;

  // copy vectors
  int nele = a_src.size();
  m_diff->grad.resize(nele);
  for (int k=0;k<nele;k++)
    {
      if(a_src.isMatCreated(k))
	MatConvert(a_src.getDeriv(k), MATSAME, MAT_INITIAL_MATRIX , &m_diff->grad[k])
    }
  m_diff->dependent = a_src.m_diff->dependent;
  m_diff->colSize = a_src.m_diff->colSize;
  m_diff->colType = a_src.m_diff->colType;
  m_diff->isMatCreated = a_src.m_diff->isMatCreated;
  m_diff->ColsfromCELLtoIb = a_src.m_diff->ColsfromCELLtoIb;
  m_diff->ColsfromFACEtoIb = a_src.m_diff->ColsfromFACEtoIb;

  m_diff->name = a_src.name() + string("copy");
  m_diff->rowSize = a_src.m_diff->rowSize;
  m_diff->rowType = a_src.m_diff->rowType;
  m_diff->fromCELLtoIb = a_src.m_diff->fromCELLtoIb;
  m_diff->fromFACEtoIb = a_src.m_diff->fromFACEtoIb;

  m_FAB.resize(0,NULL);
  m_FABAllocated = false;

  m_numVofs = a_src.m_numVofs;
  m_blockSize = a_src.m_blockSize;
  
}



template <class T>
ExtendedFunction<T>&
ExtendedFunction<T>::incr( const ExtendedFunction<T>& a_source, const Box & a_region, const Real&  a_scale)
{
  m_val->plus(a_source.getVal(),a_region,0,0,m_val->nComp(), a_scale);
  
  deque<Mat>& myGrad = m_diff->grad;
  const deque<Mat>& soGrad = a_source.m_diff->grad;

  vector<string> myDeps = m_diff->dependent;
  const vector<string> soDeps = a_source.m_diff->dependent;


  for(int k =0; k< a_source.size(); k++)
    {

	string PnameF = a_source.getDependentName(i);
	std::transform(PnameF.begin(), PnameF.end(), PnameF.begin(), ::tolower);
	PnameF.erase( std::remove_if( PnameF.begin(), PnameF.end(), ::isspace ), PnameF.end() );
	const Mat& dS_dX=soGrad[k];
	for  (int m=0;m< size();m++)
	  {
	    string XnameF = getDependentName(m);
	    std::transform(XnameF.begin(), XnameF.end(), XnameF.begin(), ::tolower);
	    XnameF.erase( std::remove_if( XnameF.begin(), XnameF.end(), ::isspace ), XnameF.end() );
	    found = XnameF == PnameF;
	    if(found) {mfound=m; break;}
	  }
	if(!found)
	  {
	    mfound = size();
	    myDeps.push_back(soDeps[k]);	    
	    myGrad.push_back(dS_dX);
	    Mat& dF_dX =  getDeriv(mfound);
	    MatConvert(dF_dX,  MATSEQAIJ, MAT_REUSE_MATRIX, &dF_dX);
	    MatScale(dF_dX, a_scale);
	    // fill other vectors
	    m_diff->colSize.push_back(a_source.m_diff->colSize[k]);
	    m_diff->colType.push_back(a_source.m_diff->colType[k]);
	    m_diff->isMatCreated.push_back(a_source.m_diff->isMatCreated[k]);
	    m_diff->ColsfromCELLtoIb.push_back(a_source.m_diff->ColsfromCELLtoIb[k]);
	    m_diff->ColsfromFACEtoIb.push_back(a_source.m_diff->ColsfromFACEtoIb[k]);
	  }
	else
	  {
	    Mat& dF_dX = getDeriv(mfound);
	    MatConvert(dF_dX,  MATSEQAIJ, MAT_REUSE_MATRIX, &dF_dX);
	    // to maintain constancy of a_src
	    Mat dS_dX2;
	    MatConvert(dS_dX, MATSEQAIJ, MAT_INITIAL_MATRIX , &dS_dX2)
	    MatAXPY(dF_dX,a_scale, dS_dX2, DIFFERENT_NONZERO_PATTERN);
	    MatDestroy(&dS_dX2);
	  }
	
    }
  
};
#include "NamespaceFooter.H"
#endif
