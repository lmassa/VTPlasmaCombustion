#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBAMR.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
void EBAMR::setupForGivenDBLRun(const Vector<DisjointBoxLayout>& a_grids)
{

  CH_assert(isDefined());

  m_isSetUp = true;

  int numLevels = a_grids.size();
  m_finest_level = numLevels -1;
  m_finest_level_old = m_finest_level;

  // set to zero to start
  m_finest_level=0;

  m_use_meshrefine = false;

  // set this to -1 initially
  m_dt_base = -1;

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  m_amr_grids.resize(numLevels);
  for (int ilev = 0; ilev < numLevels; ilev++)  m_amr_grids[ilev] = a_grids[ilev].boxArray();

  m_regrid_intervals.resize(m_max_level, -1);
  m_cur_step = 0;
  s_step = m_cur_step;
  m_finest_level = m_amr_grids.size() - 1;

  if (m_finest_level > m_max_level)
    {
      cerr << "EBAMR::setupForFixedHierarchy: too many levels in input grids" << endl;
      cerr << "max_level from define function = " << m_max_level << endl;;
      cerr << "finest level(top level of grids in setup) = = " << m_finest_level << endl;
      MayDay::Error("EBAMR::setupForFixedHierarchy: define call and setup call inconsistent in class EBAMR");
    }

  // lucacancel::This will have to be revisited in order to do adaptation on the perturbation
  /*{
    CH_TIME("FixGridCancel");
    for (int level = m_finest_level; level <= m_max_level; ++level) 
      {
	m_amrlevels[level]->finerLevelPtr(NULL);
      }
    for (int level = m_finest_level + 1; level <= m_max_level; ++level)
      {
	m_amrlevels[level]->initialGrid(Vector<Box>());
	m_amrlevels[level]->regrid(Vector<Box>());
	//m_amrlevels[level] = NULL;
      }
    m_max_level=m_finest_level;
    }*/

  {
    CH_TIME("Initialize");
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_amrlevels[level]->initialGrid(a_grids[level]);
      }
    for (int level = m_finest_level; level >= 0; --level)
      {
        m_amrlevels[level]->postInitialGrid(false);
      }
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_amrlevels[level]->initialData();
      }
  }
      
  {
    CH_TIME("PostInitialize");
    // call post-initialize once all the levels have been defined
    for (int level = m_finest_level; level >= 0; --level)
      {
        m_amrlevels[level]->postInitialize();
      }
  }

  {
    CH_TIME("InitialDt");
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_dt_new[level] = m_amrlevels[level]->computeInitialDt();
        m_dt_cur[level] = m_dt_new[level];
      }

    assignDt();


    //for (int level = m_finest_level; level <= m_max_level; ++level) 
    //{
      //m_amrlevels[level]->finerLevelPtr(NULL);
    //}
    for (int level = m_finest_level + 1; level <= m_max_level; ++level)
      {
	m_amrlevels[level]->initialGrid(Vector<Box>());
	m_amrlevels[level]->regrid(Vector<Box>());
	//m_amrlevels[level] = NULL;
      }
  }
}

//-----------------------------------------------------------------------
void EBAMR::setupForGivenICRun(const Vector<DisjointBoxLayout>& a_grids)
{

  CH_assert(isDefined());

  m_isSetUp = true;

  int numLevels = a_grids.size();
  m_finest_level = numLevels -1;
  m_finest_level_old = m_finest_level;


  m_use_meshrefine = false;

  // set this to -1 initially
  m_dt_base = -1;
  m_old_dt_base = 1e-5; // luca just a number

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  m_amr_grids.resize(numLevels);
  for (int ilev = 0; ilev < numLevels; ilev++)  m_amr_grids[ilev] = a_grids[ilev].boxArray();

  m_regrid_intervals.resize(m_max_level, -1);
  m_cur_step = 0;
  m_cur_time = 0;
  s_step = m_cur_step;

  if (m_finest_level > m_max_level)
    {
      cerr << "EBAMR::setupForFixedHierarchy: too many levels in input grids" << endl;
      cerr << "max_level from define function = " << m_max_level << endl;;
      cerr << "finest level(top level of grids in setup) = = " << m_finest_level << endl;
      MayDay::Error("EBAMR::setupForFixedHierarchy: define call and setup call inconsistent in class EBAMR");
    }

  // lucacancel::This will have to be revisited in order to run adaptation on the perturbation
  {
    CH_TIME("FixGridCancel");
    for (int level = m_finest_level; level <= m_max_level; ++level) 
      {
	m_amrlevels[level]->finerLevelPtr(NULL);
      }
    for (int level = m_finest_level + 1; level <= m_max_level; ++level)
      {
	/*m_amrlevels[level]->initialGrid(Vector<Box>());
	  m_amrlevels[level]->regrid(Vector<Box>());*/
	m_amrlevels[level] = NULL;
      }
    m_max_level=m_finest_level;
  }

  {
    CH_TIME("Initialize");
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_amrlevels[level]->initialGrid(a_grids[level]);
      }
    for (int level = m_finest_level; level >= 0; --level)
      {
        m_amrlevels[level]->postInitialGrid(false);
      }
  }

  {
    CH_TIME("PostInitialize");
    // call post-initialize once all the levels have been defined
    for (int level = m_finest_level; level >= 0; --level)
      {
        m_amrlevels[level]->postInitialize();
      }
  }

  {
    CH_TIME("InitialDt");
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_dt_new[level] = m_amrlevels[level]->computeInitialDt();
        m_dt_cur[level] = m_dt_new[level];
      }

    assignDt();
  }
}

//-----------------------------------------------------------------------
void EBAMR::writePlotFile(const string& a_prefix, const int a_ival) const
{
  CH_TIME("EBAMR::writePlotFile");

  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "EBAMR::writePlotFile" << endl;
    }

#ifdef CH_USE_HDF5
  string iter_str = a_prefix;

  char suffix[100];
  sprintf(suffix,"%06d.%dd.hdf5",a_ival,SpaceDim);

  iter_str += suffix;

  if (m_verbosity >= 2)
    {
      pout() << "plot file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str.c_str(), HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = a_ival;
  header.m_real["time"]       = m_cur_time;

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if (m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // write physics class header data
  m_amrlevels[0]->writePlotHeader(handle);

  // write physics class per-level data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writePlotLevel(handle);
    }

  handle.close();
#endif

  // Here's an extra hook for doing custom plots. :-P -JNJ
  m_amrlevels[0]->writeCustomPlotFile(m_plotfile_prefix, a_ival);
}
void EBAMR::setFinestLevel(int a_finest_level)
{
  m_finest_level = a_finest_level;   
}
void EBAMR::writeCheckpointData(HDF5Handle& a_handle) const
{
  // write amr data
  HDF5HeaderData header;
  header.m_int ["Smax_level"]  = m_max_level;
  header.m_int ["Snum_levels"] = m_finest_level + 1;
  header.m_int ["Siteration"]  = m_cur_step;
  header.m_real["Stime"]       = m_cur_time;
  header.writeToFile(a_handle);

  // write physics class data
  m_amrlevels[0]->writeCheckpointHeader(a_handle);

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writeCheckpointLevel(a_handle);
    }

}
#ifdef CH_USE_HDF5
// same as the original in AMR.cpp, but it reads the S data, to allow for a different number of levels between fluids and chemistry
void EBAMR::setupForRestartS(HDF5Handle& a_handle)
{
  CH_TIME("EBAMR::setupForRestartS");

  CH_assert(m_isDefined);

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::restart" << endl;
    }

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (m_verbosity >= 3)
    {
      pout() << "hdf5 header data: " << endl;
      pout() << header << endl;
    }

  // read max level
  if (header.m_int.find("Smax_level") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain max_level");
    }

  // note that this check should result in a warning rather than an error,
  // because you should be able to restart with a different number of levels
  // (DFM 2/27/02)
  int max_level_check = header.m_int ["Smax_level"];
  if (max_level_check != m_max_level)
    {
      pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "max level input to define = " << m_max_level << endl;
      pout() << "max level in checkpoint = " << max_level_check << endl;
      ///MayDay::Warning("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }

  if (m_verbosity >= 2)
    {
      pout() << "read max_level = " << m_max_level << endl;
    }

  // read finest level
  if (header.m_int.find("Snum_levels") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain num_levels");
    }
  int num_levels = header.m_int ["Snum_levels"];

  if (m_verbosity >= 2)
    {
      pout() << "read num_levels = " << num_levels << endl;
    }

  m_finest_level = num_levels - 1;

  if (m_finest_level > m_max_level)
    {
      pout() << "EBAMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "numlevels input to define = " << m_max_level + 1<< endl;
      pout() << "numlevels in checkpoint = " << num_levels << endl;
      m_finest_level = min(m_finest_level,m_max_level);
      //MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }
  m_finest_level_old = m_finest_level;

  if (m_verbosity >= 2)
    {
      pout() << "set finest_level = " << m_finest_level << endl;
    }

  if (m_finest_level > m_max_level)
    {
      MayDay::Error("AMR::restart: finest_level > max_level");
    }

  if (header.m_int.find("Siteration") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
    }

  m_cur_step = header.m_int ["Siteration"];
  s_step = m_cur_step;

  if (m_verbosity >= 2)
    {
      pout() << "read cur_step = " << m_cur_step << endl;
    }

  m_restart_step = m_cur_step;

  if (m_verbosity >= 2)
    {
      pout() << "set restart_step = " << m_restart_step << endl;
    }

  if (header.m_real.find("Stime") == header.m_real.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain time");
    }

  m_cur_time = header.m_real["Stime"];
  m_old_dt_base = 1e-5;

  if (m_verbosity >= 2)
    {
      pout() << "read cur_time = " << m_cur_time << endl;
    }

  // read physics class header data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      // reset to root group to read header info
      a_handle.setGroup("/");
      m_amrlevels[level]->readCheckpointHeader(a_handle);
      m_amrlevels[level]->readCheckpointLevel(a_handle);
    }

  for (int level = 0; level < m_finest_level; ++level)
    {
      int refratio_test = m_amrlevels[level]->refRatio();
      if (refratio_test != m_ref_ratios[level])
        {
          pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
          pout() << "for level " << level << endl;
          pout() << "refratio input to define = " << m_ref_ratios[level] << endl;
          pout() << "refratio in checkpoint = "  << refratio_test << endl;
          MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
        }
    }

  // fine to coarse transversal of levels for grid setup
  for (int level = m_finest_level; level >= 0; --level)
    {
      m_amrlevels[level]->postInitialGrid(true);
    }

  // maintain time steps
  m_dt_new.resize(m_max_level+1);
  m_dt_cur.resize(m_max_level+1);

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_dt_new[level] = m_amrlevels[level]->dt();
      m_dt_cur[level] = m_dt_new[level];
      pout() << "  maintain time steps m_dt_new of " << level << " = " << m_dt_new[level] << endl;
    }

  assignDt();

  // maintain steps since regrid
  m_steps_since_regrid.resize(m_max_level+1,0);

  //restart cell updates(we could also output them to the chk file)
  m_cell_updates.resize(m_max_level+1,0);

  // final thing to do -- call initialGrid and initialData on undefined levels
  // (just in case there are setup things which need to be done there
  // (DFM 2/27/02)
  for (int level = m_finest_level+1; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->initialGrid(Vector<Box>());
      m_amrlevels[level]->initialData();
    }
}
#endif
//-----------------------------------------------------------------------
void EBAMR::setTime(const Real& a_time)
{
  
  for (int level = 0; level <= m_max_level; ++level)
    m_amrlevels[level]->time(a_time);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
#include "NamespaceFooter.H"
