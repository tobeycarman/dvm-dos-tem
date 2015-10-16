#ifndef COHORT_H_
#define COHORT_H_

#include "../../include/Climate.h"

#include "../ecodomain/Ground.h"
#include "../ecodomain/Vegetation.h"

#include "../vegetation/Vegetation_Env.h"
#include "../vegetation/Vegetation_Bgc.h"

#include "../snowsoil/Snow_Env.h"
#include "../snowsoil/Soil_Env.h"
#include "../snowsoil/SoilParent_Env.h"
#include "../snowsoil/Soil_Bgc.h"

#include "../disturb/WildFire.h"

#include "../data/CohortData.h"

#include "../data/EnvData.h"
#include "../data/BgcData.h"
#include "../data/FirData.h"

#include "../data/RestartData.h"

#include "../lookup/CohortLookup.h"

#include "Integrator.h"

// headers for run
#include "ModelData.h"
#include "OutRetrive.h"

class Cohort {
public :
  Cohort();
  Cohort(int y, int x, ModelData* modeldatapointer);
  ~Cohort();

  ///< A lookup object that serves as a bridge between parameters stored in files
  ///< and setting the parameter values inside the model's run time.
  CohortLookup chtlu;

  /*
    Note: FRI is a member of CohortData because it is checked in
    Soil_Bgc::prepareintegration(...), and at that point there is no access to
    the members/fields of a Cohort...
  */

  
  /** @name Location / Spatial Reference */
  ///@{
  int y;       ///< pixel coordinate, (row)
  int x;       ///< pixel coordinate, (col)
  float lon;   ///< degrees W
  float lat;   ///< degres N
  ///@}


  /** @name Deprecated? */
  ///@{
  int errorid;    ///< model running status
  bool failed;    ///< when an exception is caught, set failed to be true
  OutRetrive outbuffer; ///< output...

  //double pfsize[NUM_FSIZE];
  //double pfseason[NUM_FSEASON];

  ///@}

  /** @name Domain 
   * Not sure why these objects/structs are grouped separately from the
   * obects in the Data section below??
  */
  ///@{
  Vegetation veg;
  Ground ground;
  Climate climate;
  ///@}

  /** @name Processes
   * Need a better description here...
  */
  ///@{
  Vegetation_Env vegenv[NUM_PFT];
  Vegetation_Bgc vegbgc[NUM_PFT];
  Snow_Env snowenv;
  Soil_Env soilenv;
  SoilParent_Env solprntenv;
  Soil_Bgc soilbgc;
  WildFire fire;
  ///@}


  /** @name Data
   *  Each of these objects contains structs with state, diagnostic, and flux
   *  values. Each object also has methods assorted methods for clearing data
   *  at the begining and end of each year and some mechanisms for summing and
   *  or averaging the state, diagnostic and flux values.
  */
  ///@{
  EnvData ed[NUM_PFT];   ///< explicitly for PFTs
  BgcData bd[NUM_PFT];   ///< explicitly for PFTs
  EnvData * edall;       ///< summed across PFTs?
  BgcData * bdall;       ///< summed across PFTs?
  FirData * fd;

  // may need a new grouping for these?
  CohortData cd;           ///< ? different..??
  ModelData * md;          ///< ? should be called Config
  RestartData restartdata; ///< ? the model state for saving/restarting
  ///@}

  /** @name Setup functions. */
  ///@{
  void initialize_internal_pointers();

  void setModelData(ModelData* md);
  void setProcessData(EnvData * alledp, BgcData * allbdp, FirData *fdp);

  void initialize_state_parameters();
  ///@}

  /** @name Driving functions. */
  ///@{
  void updateMonthly(const int yrcnt, const int currmind, const std::string& stage);
  
  void set_state_from_restartdata();
  void set_restartdata_from_state();
  ///@}

private:

  Integrator vegintegrator[NUM_PFT];
  Integrator solintegrator;


  void updateMonthly_DIMveg(const int & currmind, const bool & dvmmodule);
  void updateMonthly_DIMgrd(const int & currmind, const bool & dslmodule);

  void updateMonthly_Env(const int & currmind);
  void updateMonthly_Bgc(const int & currmind);
  void updateMonthly_Dsb(const int & yrcnt, const int & currmind, const std::string& stage);

  // Fire is a type of disturbance?
  void updateMonthly_Fir(const int & year, const int & midx, const std::string& stage);

  // update root distribution
  void getSoilFineRootFrac_Monthly();
  double assignSoilLayerRootFrac(const double & topz, const double & botz,
                                 const double csumrootfrac[MAX_ROT_LAY],
                                 const double dzrotlay[MAX_ROT_LAY]);

  // distribute values amongst PFTs...???
  void assignAtmEd2pfts_daily();
  void assignGroundEd2pfts_daily();
  void assignSoilBd2pfts_monthly();

  // aggregate values from all PFTs...???
  void getSoilTransfactor4all_daily();
  void getEd4allveg_daily();
  void getEd4land_daily();
  void getBd4allveg_monthly();

};
#endif /*COHORT_H_*/
