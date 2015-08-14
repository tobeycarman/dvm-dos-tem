#ifndef COHORT_H_
#define COHORT_H_

#include "../Climate.h"

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
  
  int y;
  int x;

  float lon;
  float lat;

  // model running status
  int errorid;
  bool failed;    // when an exception is caught, set failed to be true


  // 1) for eq and sp we need a fire recurrance interval
  // varies with space and *NOT* time
  // generated by some statictical munging of the SNAP Fire History tif files
  // will result in a single band tiff
  //int fri; /* PUT THIS IN COHORT DATA?? */

  // 2) for tr and sc FRI is based on ALF outputs or historical records
  // (for AK from 1952) year_of_burn - explicit replacement for FRI in tr and
  // sc (maybe from ALF)
  std::vector<int> fire_years;
  std::vector<int> fire_sizes;

  // old? can I deprecate these??
  double pfsize[NUM_FSIZE];
  double pfseason[NUM_FSEASON];
  
  //inputs
  CohortLookup chtlu;

  // domain
  Vegetation veg;
  Ground ground;
  
  // new domain
  Climate climate;

  // processes
  Vegetation_Env vegenv[NUM_PFT];
  Snow_Env snowenv;
  Soil_Env soilenv;
  SoilParent_Env solprntenv;

  Vegetation_Bgc vegbgc[NUM_PFT];
  Soil_Bgc soilbgc;

  WildFire fire;

  // output
  OutRetrive outbuffer;

  // data
  EnvData ed[NUM_PFT];
  BgcData bd[NUM_PFT];
  EnvData * edall;
  BgcData * bdall;

  FirData * fd;   // this for all PFTs and their soil

  ModelData * md;

  CohortData cd;
  RestartData resid;    //for input
  RestartData restartdata;
  

//  void NEW_load_climate_from_file(int y, int x);
//  void NEW_load_veg_class_from_file(int y, int x);
//  void NEW_load_fire_from_file(int y, int x);

  void initialize_internal_pointers();

  void setModelData(ModelData* md);
  void setProcessData(EnvData * alledp, BgcData * allbdp, FirData *fdp);

  void initialize_state_parameters();
//  void prepareAllDrivingData();
  void prepareDayDrivingData(const int & yrcnt, const int &usedatmyr);
  void updateMonthly(const int & yrcnt, const int & currmind,
                     const int & dinmcurr);
  
  void sync_state_to_restartdata();


private:

  Integrator vegintegrator[NUM_PFT];
  Integrator solintegrator;


  void updateMonthly_DIMveg(const int & currmind, const bool & dvmmodule);
  void updateMonthly_DIMgrd(const int & currmind, const bool & dslmodule);

  void updateMonthly_Env(const int & currmind, const int & dinmcurr);
  void updateMonthly_Bgc(const int & currmind);
  void updateMonthly_Fir(const int & yrcnt, const int & currmind);

  // update root distribution
  void getSoilFineRootFrac_Monthly();
  double assignSoilLayerRootFrac(const double & topz, const double & botz,
                                 const double csumrootfrac[MAX_ROT_LAY],
                                 const double dzrotlay[MAX_ROT_LAY]);

  //
  void assignAtmEd2pfts_daily();
  void assignGroundEd2pfts_daily();
  void getSoilTransfactor4all_daily();
  void getEd4allveg_daily();
  void getEd4land_daily();

  void assignSoilBd2pfts_monthly();
  void getBd4allveg_monthly();

};
#endif /*COHORT_H_*/
