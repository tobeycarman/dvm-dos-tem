#include <string>
#include <algorithm>
#include <json/writer.h>


#ifdef WITHNCPAR
#include <netcdf_par.h>
#else
#include <netcdf.h>
#endif

#ifdef WITHMPI
#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <stdio.h>
#endif

#include "../include/Runner.h"
#include "../include/Cohort.h"
#include "../include/TEMUtilityFunctions.h"
#include "../include/TEMLogger.h"
#include "../include/tbc-debug-util.h"

extern src::severity_logger< severity_level > glg;

Runner::Runner(ModelData mdldata, bool cal_mode, int y, int x):
    calibrationMode(false), y(y), x(x) {

  BOOST_LOG_SEV(glg, note) << "RUNNER Constructing a Runner, new style, with ctor-"
                           << "injected ModelData, and for explicit (y,x) "
                           << "position w/in the input data region.";
  this->md = mdldata;
  this->cohort = Cohort(y, x, &mdldata); // explicitly constructed cohort...

  BOOST_LOG_SEV(glg, info) << "Calibration mode?: " << cal_mode;
  if ( cal_mode ) {
    this->calcontroller_ptr.reset( new CalController(&this->cohort) );
    this->calcontroller_ptr->clear_and_create_json_storage();
  } // else null??


  // within-grid cohort-level aggregated 'ed' (i.e. 'edall in 'cht')
  BOOST_LOG_SEV(glg, debug) << "Create some empty containers for 'cohort-level "
                            << "aggregations of 'ed', (i.e. 'edall in 'cohort')";
  this->chted = EnvData();
  this->chtbd = BgcData();
  this->chtfd = FirData();
  
  // Now give the cohort pointers to these containers.
  this->cohort.setProcessData(&this->chted, &this->chtbd, &this->chtfd);

}


Runner::~Runner() {
};

void Runner::run_years(int start_year, int end_year, const std::string& stage) {

  /** YEAR TIMESTEP LOOP */
  BOOST_LOG_NAMED_SCOPE("Y") {
  for (int iy = start_year; iy < end_year; ++iy) {
    BOOST_LOG_SEV(glg, debug) << "(Beginning of year loop) " << cohort.ground.layer_report_string("depth thermal CN");
    BOOST_LOG_SEV(glg, err) << "y: "<<this->y<<" x: "<<this->x<<" Year: "<<iy;

    /* Interpolate all the monthly values...? */
    if( (stage.find("eq") != std::string::npos
           || stage.find("pre") != std::string::npos) ){
      this->cohort.climate.prepare_daily_driving_data(iy, stage);
    }

    else if(stage.find("sp") != std::string::npos){
      //FIX - 30 should not be hardcoded
      this->cohort.climate.prepare_daily_driving_data(iy%30, stage);
    }

    else if(stage.find("tr") != std::string::npos
              || stage.find("sc") != std::string::npos){
      this->cohort.climate.prepare_daily_driving_data(iy, stage);
    }


    if (this->calcontroller_ptr) { // should be null unless we are in "calibration mode"

      this->output_debug_daily_drivers(iy, this->calcontroller_ptr->daily_json);

      // Run any pre-configured directives
      this->calcontroller_ptr->run_config(iy, stage);

      // See if a signal has arrived (possibly from user
      // hitting Ctrl-C) and if so, stop the simulation
      // and drop into the calibration "shell".
      this->calcontroller_ptr->check_for_signals();

    }

    /** MONTH TIMESTEP LOOP */
    BOOST_LOG_NAMED_SCOPE("M") {
      for (int im = 0; im < 12; ++im) {
        BOOST_LOG_SEV(glg, note) << "(Beginning of month loop, iy:"<<iy<<", im:"<<im<<") " << cohort.ground.layer_report_string("depth thermal CN desc");

        this->cohort.updateMonthly(iy, im, DINM[im], stage);

        this->monthly_output(iy, im, stage);

      } // end month loop
    } // end named scope

    this->yearly_output(iy, stage, start_year, end_year);

    BOOST_LOG_SEV(glg, note) << "(END OF YEAR) " << cohort.ground.layer_report_string("depth thermal CN ptr");

    BOOST_LOG_SEV(glg, note) << "Completed year " << iy << " for cohort/cell (row,col): (" << this->y << "," << this->x << ")";

  }} // end year loop (and named scope
}

void Runner::monthly_output(const int year, const int month, const std::string& runstage) {

  if (md.output_monthly) {

    // Calibration json files....
    if(this->calcontroller_ptr) {
      BOOST_LOG_SEV(glg, debug) << "Write monthly data to json files...";
      this->output_caljson_monthly(year, month, runstage, this->calcontroller_ptr->monthly_json);
    }

  } else {
    BOOST_LOG_SEV(glg, debug) << "Monthly output turned off in config settings.";
  }

  // NetCDF output is not controlled by the monthly output flag in
  // the config file. TODO Semi-kludgy
  if(   (runstage.find("eq")!=std::string::npos && md.nc_eq)
     || (runstage.find("sp")!=std::string::npos && md.nc_sp)
     || (runstage.find("tr")!=std::string::npos && md.nc_tr)
     || (runstage.find("sc")!=std::string::npos && md.nc_sc) ){
    BOOST_LOG_SEV(glg, debug) << "Monthly NetCDF output function call, runstage: "<<runstage<<" year: "<<year<<" month: "<<month;
    output_netCDF_monthly(year, month, runstage);
  }

}

void Runner::yearly_output(const int year, const std::string& stage,
    const int startyr, const int endyr) {

  if(this->calcontroller_ptr) {
    if ( -1 == md.last_n_json_files ) {
      this->output_caljson_yearly(year, stage, this->calcontroller_ptr->yearly_json);
    }

    if ( year >= (endyr - md.last_n_json_files) ) {
      this->output_caljson_yearly(year, stage, this->calcontroller_ptr->yearly_json);
    }
  }

  // NetCDF output. Semi-kludgy
  if(   (stage.find("eq")!=std::string::npos && md.nc_eq)
     || (stage.find("sp")!=std::string::npos && md.nc_sp)
     || (stage.find("tr")!=std::string::npos && md.nc_tr)
     || (stage.find("sc")!=std::string::npos && md.nc_sc) ){
    BOOST_LOG_SEV(glg, debug) << "Yearly NetCDF output function call, runstage: "<<stage<<" year: "<<year;
    output_netCDF_yearly(year, stage);
  }


}

std::string Runner::report_not_equal(const std::string& a_desc,
                                     const std::string& b_desc,
                                     int PFT,
                                     double A, double B) {
  std::stringstream ss;
  if ( !temutil::AlmostEqualRelative(A, B) ) {
    ss << "PFT:" << PFT << " " << a_desc << " and " << b_desc
       << " not summing correctly!" << " A: "<< A <<" B: "<< B
       << " (A-B: "<< A - B <<")";
  }
  return ss.str(); // empty string if no error
}

std::string Runner::report_not_equal(double A, double B, const std::string& msg) {
  std::stringstream ss;
  if ( !temutil::AlmostEqualRelative(A, B) ) {
    ss << msg <<" A: "<< A <<" B: "<< B <<" (A-B: "<< A - B <<")";
  }
  return ss.str(); // empty string if no error
}

/** Used to check that sums across PFT compartments match the 
    corresponding 'all' container.

*/
std::list<std::string> Runner::check_sum_over_compartments() {

  std::list<std::string> errlist;

  for (int ip = 0; ip < NUM_PFT; ++ip) {

    errlist.push_back(report_not_equal(
                      "whole plant C", "plant C PART", ip,
                      cohort.bd[ip].m_vegs.call,
                      cohort.bd[ip].m_vegs.c[I_leaf] +
                      cohort.bd[ip].m_vegs.c[I_stem] +
                      cohort.bd[ip].m_vegs.c[I_root]));

    errlist.push_back(report_not_equal("whole plant C", "plant C PART", ip,
                  cohort.bd[ip].m_vegs.call,
                  cohort.bd[ip].m_vegs.c[I_leaf] +
                  cohort.bd[ip].m_vegs.c[I_stem] +
                  cohort.bd[ip].m_vegs.c[I_root]));

    errlist.push_back(report_not_equal("whole plant strn", "plant strn PART", ip,
                  cohort.bd[ip].m_vegs.strnall,
                  cohort.bd[ip].m_vegs.strn[I_leaf] +
                  cohort.bd[ip].m_vegs.strn[I_stem] +
                  cohort.bd[ip].m_vegs.strn[I_root]));

    errlist.push_back(report_not_equal("whole plant ingpp", "plant ingpp PART", ip,
                  cohort.bd[ip].m_a2v.ingppall,
                  cohort.bd[ip].m_a2v.ingpp[I_leaf] +
                  cohort.bd[ip].m_a2v.ingpp[I_stem] +
                  cohort.bd[ip].m_a2v.ingpp[I_root]));


    errlist.push_back(report_not_equal("whole plant gpp", "plant gpp PART", ip,
                  cohort.bd[ip].m_a2v.gppall,
                  cohort.bd[ip].m_a2v.gpp[I_leaf] +
                  cohort.bd[ip].m_a2v.gpp[I_stem] +
                  cohort.bd[ip].m_a2v.gpp[I_root]));

    errlist.push_back(report_not_equal("whole plant npp", "plant npp PART", ip,
                  cohort.bd[ip].m_a2v.nppall,
                  cohort.bd[ip].m_a2v.npp[I_leaf] +
                  cohort.bd[ip].m_a2v.npp[I_stem] +
                  cohort.bd[ip].m_a2v.npp[I_root]));

    errlist.push_back(report_not_equal("whole plant innpp", "plant innpp PART", ip,
                  cohort.bd[ip].m_a2v.innppall,
                  cohort.bd[ip].m_a2v.innpp[I_leaf] +
                  cohort.bd[ip].m_a2v.innpp[I_stem] +
                  cohort.bd[ip].m_a2v.innpp[I_root]));

    errlist.push_back(report_not_equal("whole plant rm", "plant rm PART", ip,
                  cohort.bd[ip].m_v2a.rmall,
                  cohort.bd[ip].m_v2a.rm[I_leaf] +
                  cohort.bd[ip].m_v2a.rm[I_stem] +
                  cohort.bd[ip].m_v2a.rm[I_root]));

    errlist.push_back(report_not_equal("whole plant rg", "plant rg PART", ip,
                  cohort.bd[ip].m_v2a.rgall,
                  cohort.bd[ip].m_v2a.rg[I_leaf] +
                  cohort.bd[ip].m_v2a.rg[I_stem] +
                  cohort.bd[ip].m_v2a.rg[I_root]));

    errlist.push_back(report_not_equal("whole plant N litterfall", "plant N litterfall PART", ip,
                  cohort.bd[ip].m_v2soi.ltrfalnall + cohort.bd[ip].m_v2soi.mossdeathn,
                  cohort.bd[ip].m_v2soi.ltrfaln[I_leaf] +
                  cohort.bd[ip].m_v2soi.ltrfaln[I_stem] +
                  cohort.bd[ip].m_v2soi.ltrfaln[I_root]));

    errlist.push_back(report_not_equal("whole plant C litterfall", "plant C litterfall PART", ip,
                  cohort.bd[ip].m_v2soi.ltrfalcall + cohort.bd[ip].m_v2soi.mossdeathc,
                  cohort.bd[ip].m_v2soi.ltrfalc[I_leaf] +
                  cohort.bd[ip].m_v2soi.ltrfalc[I_stem] +
                  cohort.bd[ip].m_v2soi.ltrfalc[I_root]));

    errlist.push_back(report_not_equal("whole plant snuptake", "plant snuptake PART", ip,
                  cohort.bd[ip].m_soi2v.snuptakeall,
                  cohort.bd[ip].m_soi2v.snuptake[I_leaf] +
                  cohort.bd[ip].m_soi2v.snuptake[I_stem] +
                  cohort.bd[ip].m_soi2v.snuptake[I_root]));

    errlist.push_back(report_not_equal("whole plant nmobil", "plant nmobil PART", ip,
                  cohort.bd[ip].m_v2v.nmobilall,
                  cohort.bd[ip].m_v2v.nmobil[I_leaf] +
                  cohort.bd[ip].m_v2v.nmobil[I_stem] +
                  cohort.bd[ip].m_v2v.nmobil[I_root]));

    errlist.push_back(report_not_equal("whole plant nresorb", "plant nresorb PART", ip,
                  cohort.bd[ip].m_v2v.nresorball,
                  cohort.bd[ip].m_v2v.nresorb[I_leaf] +
                  cohort.bd[ip].m_v2v.nresorb[I_stem] +
                  cohort.bd[ip].m_v2v.nresorb[I_root]));

  } // end loop over PFTS

  // The way this works, we end up with a bunch of empty items in the list,
  // so here we remove them.
  errlist.erase(std::remove_if(errlist.begin(), errlist.end(), temutil::emptyContainer<std::string>), errlist.end());

  return errlist;
}

/** Sum across PFTs, compare with ecosystem totals (eg data from 'bdall').

Used to add up across all pfts (data held in cohort's bd array of BgcData objects)
and compare with the data held in cohort's bdall BgcData object.
*/
std::list<std::string> Runner::check_sum_over_PFTs(){

  double ecosystem_C = 0;
  double ecosystem_C_by_compartment = 0;

  double ecosystem_N = 0;
  double ecosystem_strn = 0;
  double ecosystem_strn_by_compartment = 0;
  double ecosystem_labn = 0;

  double ecosystem_ingpp = 0;
  double ecosystem_gpp = 0;
  double ecosystem_innpp = 0;
  double ecosystem_npp = 0;

  double ecosystem_rm = 0;
  double ecosystem_rg = 0;

  double ecosystem_ltrfalc = 0;
  double ecosystem_ltrfaln = 0;
  double ecosystem_snuptake = 0;
  double ecosystem_nmobil = 0;
  double ecosystem_nresorb = 0;

  // sum various quantities over all PFTs
  for (int ip = 0; ip < NUM_PFT; ++ip) {
    ecosystem_C += this->cohort.bd[ip].m_vegs.call;
    ecosystem_C_by_compartment += (this->cohort.bd[ip].m_vegs.c[I_leaf] +
                                   this->cohort.bd[ip].m_vegs.c[I_stem] +
                                   this->cohort.bd[ip].m_vegs.c[I_root]);

    ecosystem_strn += this->cohort.bd[ip].m_vegs.strnall;
    ecosystem_strn_by_compartment += (this->cohort.bd[ip].m_vegs.strn[I_leaf] +
                                      this->cohort.bd[ip].m_vegs.strn[I_stem] +
                                      this->cohort.bd[ip].m_vegs.strn[I_root]);
    ecosystem_labn += this->cohort.bd[ip].m_vegs.labn;

    ecosystem_N += (this->cohort.bd[ip].m_vegs.strnall + this->cohort.bd[ip].m_vegs.labn);

    ecosystem_ingpp += this->cohort.bd[ip].m_a2v.ingppall;
    ecosystem_gpp += this->cohort.bd[ip].m_a2v.gppall;
    ecosystem_innpp += this->cohort.bd[ip].m_a2v.innppall;
    ecosystem_npp += this->cohort.bd[ip].m_a2v.nppall;

    ecosystem_rm += this->cohort.bd[ip].m_v2a.rmall;
    ecosystem_rg += this->cohort.bd[ip].m_v2a.rgall;

    ecosystem_ltrfalc += this->cohort.bd[ip].m_v2soi.ltrfalcall;
    ecosystem_ltrfaln += this->cohort.bd[ip].m_v2soi.ltrfalnall;

    ecosystem_snuptake += this->cohort.bd[ip].m_soi2v.snuptakeall;

    ecosystem_nmobil += this->cohort.bd[ip].m_v2v.nmobilall;
    ecosystem_nresorb += this->cohort.bd[ip].m_v2v.nresorball;

  }

  std::list<std::string> errlist;

  // Check that the sums are equal to the Runner level containers (ecosystem totals)
  errlist.push_back(report_not_equal(this->cohort.bdall->m_vegs.call, ecosystem_C, "Runner:: ecosystem veg C not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_vegs.call, ecosystem_C_by_compartment, "Runner:: ecosystem veg C not matching sum over compartments"));

  errlist.push_back(report_not_equal(this->cohort.bdall->m_vegs.nall, ecosystem_N, "Runner:: ecosystem nall not matching sum over PFTs (of strn and nall)!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_vegs.labn, ecosystem_labn, "Runner:: ecosystem labn not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_vegs.strnall, ecosystem_strn, "Runner:: ecosystem strn not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_vegs.strnall, ecosystem_strn_by_compartment, "Runner:: ecosystem strn not matching sum over compartments!"));

  errlist.push_back(report_not_equal(this->cohort.bdall->m_a2v.ingppall, ecosystem_ingpp, "Runner:: ecosystem npp not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_a2v.gppall, ecosystem_gpp, "Runner:: ecosystem npp not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_a2v.innppall, ecosystem_innpp, "Runner:: ecosystem innpp not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_a2v.nppall, ecosystem_npp, "Runner:: ecosystem npp not matching sum over PFTs!"));

  errlist.push_back(report_not_equal(this->cohort.bdall->m_v2a.rmall, ecosystem_rm, "Runner:: ecosystem rm not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_v2a.rgall, ecosystem_rg, "Runner:: ecosystem rg not matching sum over PFTs!"));

  errlist.push_back(report_not_equal(this->cohort.bdall->m_v2soi.ltrfalcall, ecosystem_ltrfalc, "Runner:: ecosystem ltrfalc not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_v2soi.ltrfalnall, ecosystem_ltrfaln, "Runner:: ecosystem ltrfaln not matching sum over PFTs!"));

  errlist.push_back(report_not_equal(this->cohort.bdall->m_soi2v.snuptakeall, ecosystem_snuptake, "Runner:: ecosystem snuptake not matching sum over PFTs!"));

  errlist.push_back(report_not_equal(this->cohort.bdall->m_v2v.nmobilall, ecosystem_nmobil, "Runner:: ecosystem nmobil not matching sum over PFTs!"));
  errlist.push_back(report_not_equal(this->cohort.bdall->m_v2v.nresorball, ecosystem_nresorb, "Runner:: ecosystem nresorb not matching sum over PFTs!"));

  // The way this works, we end up with a bunch of empty items in the list,
  // so here we remove them.
  errlist.erase(std::remove_if(errlist.begin(), errlist.end(), temutil::emptyContainer<std::string>), errlist.end());

  return errlist;

}

void Runner::output_caljson_monthly(int year, int month, std::string stage, boost::filesystem::path p){

  std::list<std::string> compartment_err_report = check_sum_over_compartments();
  std::list<std::string> pft_err_report = check_sum_over_PFTs();


  if ( !compartment_err_report.empty() || !pft_err_report.empty() ) {
    BOOST_LOG_SEV(glg, err) << "========== MONTHLY CHECKSUM ERRORS ============";
    while (!compartment_err_report.empty()) {
      BOOST_LOG_SEV(glg, err) << compartment_err_report.front();
      compartment_err_report.pop_front();
    }
    while ( !(pft_err_report.empty()) ){
      BOOST_LOG_SEV(glg, err) << pft_err_report.front();
      pft_err_report.pop_front();
    }
    BOOST_LOG_SEV(glg, err) << "========== END MONTHLY CHECKSUMMING month: " << month << " year: " << year << " ============";
  }


  // CAUTION: this->md and this->cohort.md are different instances!

  Json::Value data;
  std::ofstream out_stream;

  /* Not PFT dependent */
  data["Runstage"] = stage;
  data["Year"] = year;
  data["Month"] = month;
  data["CMT"] = this->cohort.chtlu.cmtcode;
  data["Lat"] = this->cohort.lat;
  data["Lon"] = this->cohort.lon;

  data["Nfeed"] = this->cohort.md->get_nfeed();
  data["AvlNFlag"] = this->cohort.md->get_avlnflg();
  data["Baseline"] = this->cohort.md->get_baseline();
  data["EnvModule"] = this->cohort.md->get_envmodule();
  data["BgcModule"] = this->cohort.md->get_bgcmodule();
  data["DvmModule"] = this->cohort.md->get_dvmmodule();
  data["DslModule"] = this->cohort.md->get_dslmodule();
  data["DsbModule"] = this->cohort.md->get_dsbmodule();

  data["TAir"] = cohort.edall->m_atms.ta;
  data["Snowfall"] = cohort.edall->m_a2l.snfl;
  data["Rainfall"] = cohort.edall->m_a2l.rnfl;
  data["WaterTable"] = cohort.edall->m_sois.watertab;
  data["ActiveLayerDepth"] = cohort.edall->m_soid.ald;
  data["CO2"] = cohort.edall->m_atms.co2;
  data["VPD"] = cohort.edall->m_atmd.vpd;
  data["EET"] = cohort.edall->m_l2a.eet;
  data["PET"] = cohort.edall->m_l2a.pet;
  data["PAR"] = cohort.edall->m_a2v.pardown;            // <-- from edall
  data["PARAbsorb"] = cohort.edall->m_a2v.parabsorb;    // <-- from edall

  data["VWCShlw"] = cohort.edall->m_soid.vwcshlw;
  data["VWCDeep"] = cohort.edall->m_soid.vwcdeep;
  data["VWCMineA"] = cohort.edall->m_soid.vwcminea;
  data["VWCMineB"] = cohort.edall->m_soid.vwcmineb;
  data["VWCMineC"] = cohort.edall->m_soid.vwcminec;
  data["TShlw"] = cohort.edall->m_soid.tshlw;
  data["TDeep"] = cohort.edall->m_soid.tdeep;
  data["TMineA"] = cohort.edall->m_soid.tminea;
  data["TMineB"] = cohort.edall->m_soid.tmineb;
  data["TMineC"] = cohort.edall->m_soid.tminec;


  data["NMobilAll"] = cohort.bdall->m_v2v.nmobilall;
  data["NResorbAll"] = cohort.bdall->m_v2v.nresorball;

  data["StNitrogenUptakeAll"] = cohort.bdall->m_soi2v.snuptakeall;
  data["InNitrogenUptakeAll"] = cohort.bdall->m_soi2v.innuptake;
  data["AvailableNitrogenSum"] = cohort.bdall->m_soid.avlnsum;
  data["OrganicNitrogenSum"] = cohort.bdall->m_soid.orgnsum;
  data["CarbonShallow"] = cohort.bdall->m_soid.shlwc;
  data["CarbonDeep"] = cohort.bdall->m_soid.deepc;
  data["CarbonMineralSum"] = cohort.bdall->m_soid.mineac
                             + cohort.bdall->m_soid.minebc
                             + cohort.bdall->m_soid.minecc;
  // pools
  data["StandingDeadC"] = cohort.bdall->m_vegs.deadc;
  data["StandingDeadN"] = cohort.bdall->m_vegs.deadn;
  data["WoodyDebrisC"] = cohort.bdall->m_sois.wdebrisc;
  data["WoodyDebrisN"] = cohort.bdall->m_sois.wdebrisn;
  // fluxes
  data["MossDeathC"] = cohort.bdall->m_v2soi.mossdeathc;
  data["MossdeathNitrogen"] = cohort.bdall->m_v2soi.mossdeathn;
  data["D2WoodyDebrisC"] = cohort.bdall->m_v2soi.d2wdebrisc;
  data["D2WoodyDebrisN"] = cohort.bdall->m_v2soi.d2wdebrisn;

  double t = 0.0;
  for (int il=0; il < MAX_SOI_LAY; il++) {
    t += cohort.bdall->m_soi2v.nextract[il];
  }
  data["NExtract"] = t;
  data["NetNMin"] = cohort.bdall->m_soi2soi.netnminsum;
  data["NetNImmob"] = cohort.bdall->m_soi2soi.nimmobsum;
  data["OrgNInput"] = cohort.bdall->m_a2soi.orgninput;
  data["AvlNInput"] = cohort.bdall->m_a2soi.avlninput;
  data["AvlNLost"] = cohort.bdall->m_soi2l.avlnlost;
  data["RHraw"] = cohort.bdall->m_soi2a.rhrawcsum;
  data["RHsoma"] = cohort.bdall->m_soi2a.rhsomasum;
  data["RHsompr"] = cohort.bdall->m_soi2a.rhsomprsum;
  data["RHsomcr"] = cohort.bdall->m_soi2a.rhsomcrsum;
  data["RHwdeb"] = cohort.bdall->m_soi2a.rhwdeb;
  data["RH"] = cohort.bdall->m_soi2a.rhtot;

  data["YearsSinceDisturb"] = cohort.cd.yrsdist;
  data["Burnthick"] = cohort.year_fd[month].fire_soid.burnthick;
  data["BurnVeg2AirC"] = cohort.year_fd[month].fire_v2a.orgc;
  data["BurnVeg2AirN"] = cohort.year_fd[month].fire_v2a.orgn;
  data["BurnVeg2SoiAbvVegC"] = cohort.year_fd[month].fire_v2soi.abvc;
  data["BurnVeg2SoiBlwVegC"] = cohort.year_fd[month].fire_v2soi.blwc;
  data["BurnVeg2SoiAbvVegN"] = cohort.year_fd[month].fire_v2soi.abvn;
  data["BurnVeg2SoiBlwVegN"] = cohort.year_fd[month].fire_v2soi.blwn;
  data["BurnSoi2AirC"] = cohort.year_fd[month].fire_soi2a.orgc;
  data["BurnSoi2AirN"] = cohort.year_fd[month].fire_soi2a.orgn;
  data["BurnAir2SoiN"] = cohort.year_fd[month].fire_a2soi.orgn;
  data["BurnAbvVeg2DeadC"] = cohort.year_fd[month].fire_v2dead.vegC;
  data["BurnAbvVeg2DeadN"] = cohort.year_fd[month].fire_v2dead.strN;
  data["RawCSum"] = cohort.bdall->m_soid.rawcsum;
  data["SomaSum"] = cohort.bdall->m_soid.somasum;
  data["SomcrSum"] = cohort.bdall->m_soid.somcrsum;
  data["SomprSum"] = cohort.bdall->m_soid.somprsum;


  /* PFT dependent variables */

  // calculated ecosystem summary values
  double litterfallCsum = 0;
  double litterfallNsum = 0;
  double parDownSum = 0;
  double parAbsorbSum = 0;

  for(int pft=0; pft<NUM_PFT; pft++) {
    char pft_chars[5];
    sprintf(pft_chars, "%d", pft);
    std::string pft_str = std::string(pft_chars);
    // c++0x equivalent: std::string pftvalue = std::to_string(pft);
    data["PFT" + pft_str]["VegCarbon"]["Leaf"] = cohort.bd[pft].m_vegs.c[I_leaf];
    data["PFT" + pft_str]["VegCarbon"]["Stem"] = cohort.bd[pft].m_vegs.c[I_stem];
    data["PFT" + pft_str]["VegCarbon"]["Root"] = cohort.bd[pft].m_vegs.c[I_root];
    data["PFT" + pft_str]["VegStructuralNitrogen"]["Leaf"] = cohort.bd[pft].m_vegs.strn[I_leaf];
    data["PFT" + pft_str]["VegStructuralNitrogen"]["Stem"] = cohort.bd[pft].m_vegs.strn[I_stem];
    data["PFT" + pft_str]["VegStructuralNitrogen"]["Root"] = cohort.bd[pft].m_vegs.strn[I_root];
    data["PFT" + pft_str]["VegLabileNitrogen"] = cohort.bd[pft].m_vegs.labn;

    data["PFT" + pft_str]["NAll"] = cohort.bd[pft].m_vegs.nall; // <-- Sum of labn and strn
    data["PFT" + pft_str]["StandingDeadC"] = cohort.bd[pft].m_vegs.deadc;
    data["PFT" + pft_str]["StandingDeadN"] = cohort.bd[pft].m_vegs.deadn;

    data["PFT" + pft_str]["NMobil"] = cohort.bd[pft].m_v2v.nmobilall; // <- the all denotes multi-compartment
    data["PFT" + pft_str]["NResorb"] = cohort.bd[pft].m_v2v.nresorball;

    data["PFT" + pft_str]["GPPAll"] = cohort.bd[pft].m_a2v.gppall;
    data["PFT" + pft_str]["GPP"]["Leaf"] = cohort.bd[pft].m_a2v.gpp[I_leaf];
    data["PFT" + pft_str]["GPP"]["Stem"] = cohort.bd[pft].m_a2v.gpp[I_stem];
    data["PFT" + pft_str]["GPP"]["Root"] = cohort.bd[pft].m_a2v.gpp[I_root];

    data["PFT" + pft_str]["NPPAll"] = cohort.bd[pft].m_a2v.nppall;
    data["PFT" + pft_str]["NPP"]["Leaf"] = cohort.bd[pft].m_a2v.npp[I_leaf];
    data["PFT" + pft_str]["NPP"]["Stem"] = cohort.bd[pft].m_a2v.npp[I_stem];
    data["PFT" + pft_str]["NPP"]["Root"] = cohort.bd[pft].m_a2v.npp[I_root];

    data["PFT" + pft_str]["GPPAllIgnoringNitrogen"] = cohort.bd[pft].m_a2v.ingppall;
    data["PFT" + pft_str]["NPPAllIgnoringNitrogen"] = cohort.bd[pft].m_a2v.innppall;

    data["PFT" + pft_str]["LitterfallCarbonAll"] = cohort.bd[pft].m_v2soi.ltrfalcall;
    data["PFT" + pft_str]["LitterfallCarbon"]["Leaf"] = cohort.bd[pft].m_v2soi.ltrfalc[I_leaf];
    data["PFT" + pft_str]["LitterfallCarbon"]["Stem"] = cohort.bd[pft].m_v2soi.ltrfalc[I_stem];
    data["PFT" + pft_str]["LitterfallCarbon"]["Root"] = cohort.bd[pft].m_v2soi.ltrfalc[I_root];

    data["PFT" + pft_str]["RespGrowth"]["Leaf"] = cohort.bd[pft].m_v2a.rg[I_leaf];
    data["PFT" + pft_str]["RespGrowth"]["Stem"] = cohort.bd[pft].m_v2a.rg[I_stem];
    data["PFT" + pft_str]["RespGrowth"]["Root"] = cohort.bd[pft].m_v2a.rg[I_root];

    data["PFT" + pft_str]["RespMaint"]["Leaf"] = cohort.bd[pft].m_v2a.rm[I_leaf];
    data["PFT" + pft_str]["RespMaint"]["Stem"] = cohort.bd[pft].m_v2a.rm[I_stem];
    data["PFT" + pft_str]["RespMaint"]["Root"] = cohort.bd[pft].m_v2a.rm[I_root];

    data["PFT" + pft_str]["LitterfallNitrogenPFT"] = cohort.bd[pft].m_v2soi.ltrfalnall;
    data["PFT" + pft_str]["LitterfallNitrogen"]["Leaf"] = cohort.bd[pft].m_v2soi.ltrfaln[I_leaf];
    data["PFT" + pft_str]["LitterfallNitrogen"]["Stem"] = cohort.bd[pft].m_v2soi.ltrfaln[I_stem];
    data["PFT" + pft_str]["LitterfallNitrogen"]["Root"] = cohort.bd[pft].m_v2soi.ltrfaln[I_root];

    data["PFT" + pft_str]["StNitrogenUptake"] = cohort.bd[pft].m_soi2v.snuptakeall;
    data["PFT" + pft_str]["InNitrogenUptake"] = cohort.bd[pft].m_soi2v.innuptake;
    data["PFT" + pft_str]["LabNitrogenUptake"] = cohort.bd[pft].m_soi2v.lnuptake;
    data["PFT" + pft_str]["TotNitrogenUptake"] = cohort.bd[pft].m_soi2v.snuptakeall + cohort.bd[pft].m_soi2v.lnuptake;
    data["PFT" + pft_str]["MossDeathC"] = cohort.bd[pft].m_v2soi.mossdeathc;

    litterfallCsum += cohort.bd[pft].m_v2soi.ltrfalcall;
    litterfallNsum += cohort.bd[pft].m_v2soi.ltrfalnall;

    data["PFT" + pft_str]["PARDown"] = cohort.ed[pft].m_a2v.pardown;
    data["PFT" + pft_str]["PARAbsorb"] = cohort.ed[pft].m_a2v.parabsorb;

    parDownSum += cohort.ed[pft].m_a2v.pardown;
    parAbsorbSum += cohort.ed[pft].m_a2v.parabsorb;

  }

  data["LitterfallCsum"] = litterfallCsum;
  data["LitterfallNsum"] = litterfallNsum;

  data["PARAbsorbSum"] = parAbsorbSum;
  data["PARDownSum"] = parDownSum;
  data["GPPSum"] = cohort.bdall->m_a2v.gppall;
  data["NPPSum"] = cohort.bdall->m_a2v.nppall;

  // Writes files like this:
  //  0000000.json, 0000001.json, 0000002.json, ...

  std::stringstream filename;
  filename.fill('0');
  filename << std::setw(7) << (12 * year) + month << ".json";

  // Add the file name to the path
  p /= filename.str();

  // write out the data
  out_stream.open(p.string().c_str(), std::ofstream::out);
  out_stream << data << std::endl;
  out_stream.close();

}


void Runner::output_caljson_yearly(int year, std::string stage, boost::filesystem::path p) {

  std::list<std::string> compartment_err_report = check_sum_over_compartments();
  std::list<std::string> pft_err_report = check_sum_over_PFTs();

  if ( !compartment_err_report.empty() || !pft_err_report.empty() ) {
    BOOST_LOG_SEV(glg, err) << "========== YEARLY CHECKSUM ERRORS ============";
    while (!compartment_err_report.empty()) {
      BOOST_LOG_SEV(glg, err) << compartment_err_report.front();
      compartment_err_report.pop_front();
    }
    while ( !(pft_err_report.empty()) ){
      BOOST_LOG_SEV(glg, err) << pft_err_report.front();
      pft_err_report.pop_front();
    }
    BOOST_LOG_SEV(glg, err) << "========== END YEARLY CHECKSUMMING year: " << year << " ============";
  }

  // CAUTION: this->md and this->cohort.md are different instances!

  Json::Value data;
  std::ofstream out_stream;

  /* Not PFT dependent */
  data["Runstage"] = stage;
  data["Year"] = year;
  data["CMT"] = this->cohort.chtlu.cmtcode;
  data["Lat"] = this->cohort.lat;
  data["Lon"] = this->cohort.lon;

  data["Nfeed"] = this->cohort.md->get_nfeed();
  data["AvlNFlag"] = this->cohort.md->get_avlnflg();
  data["Baseline"] = this->cohort.md->get_baseline();
  data["EnvModule"] = this->cohort.md->get_envmodule();
  data["BgcModule"] = this->cohort.md->get_bgcmodule();
  data["DvmModule"] = this->cohort.md->get_dvmmodule();
  data["DslModule"] = this->cohort.md->get_dslmodule();
  data["DsbModule"] = this->cohort.md->get_dsbmodule();

  data["TAir"] = cohort.edall->y_atms.ta;
  data["Snowfall"] = cohort.edall->y_a2l.snfl;
  data["Rainfall"] = cohort.edall->y_a2l.rnfl;
  data["WaterTable"] = cohort.edall->y_sois.watertab;
  data["ActiveLayerDepth"]= cohort.edall-> y_soid.ald;
  data["CO2"] = cohort.edall->y_atms.co2;
  data["VPD"] = cohort.edall->y_atmd.vpd;
  data["EET"] = cohort.edall->y_l2a.eet;
  data["PET"] = cohort.edall->y_l2a.pet;
  data["PAR"] = cohort.edall->y_a2l.par;
  data["PARAbsorb"] = cohort.edall->y_a2v.parabsorb;

  data["VWCShlw"] = cohort.edall->y_soid.vwcshlw;
  data["VWCDeep"] = cohort.edall->y_soid.vwcdeep;
  data["VWCMineA"] = cohort.edall->y_soid.vwcminea;
  data["VWCMineB"] = cohort.edall->y_soid.vwcmineb;
  data["VWCMineC"] = cohort.edall->y_soid.vwcminec;
  data["TShlw"] = cohort.edall->y_soid.tshlw;
  data["TDeep"] = cohort.edall->y_soid.tdeep;
  data["TMineA"] = cohort.edall->y_soid.tminea;
  data["TMineB"] = cohort.edall->y_soid.tmineb;
  data["TMineC"] = cohort.edall->y_soid.tminec;

  data["NMobilAll"] = cohort.bdall->y_v2v.nmobilall;
  data["NResorbAll"] = cohort.bdall->y_v2v.nresorball;

  data["StNitrogenUptakeAll"] = cohort.bdall->y_soi2v.snuptakeall;
  data["InNitrogenUptakeAll"] = cohort.bdall->y_soi2v.innuptake;
  data["AvailableNitrogenSum"] = cohort.bdall->y_soid.avlnsum;
  data["OrganicNitrogenSum"] = cohort.bdall->y_soid.orgnsum;
  data["CarbonShallow"] = cohort.bdall->y_soid.shlwc;
  data["CarbonDeep"] = cohort.bdall->y_soid.deepc;
  data["CarbonMineralSum"] = cohort.bdall->y_soid.mineac
                             + cohort.bdall->y_soid.minebc
                             + cohort.bdall->y_soid.minecc;
  // pools
  data["StandingDeadC"] = cohort.bdall->y_vegs.deadc;
  data["StandingDeadN"] = cohort.bdall->y_vegs.deadn;
  data["WoodyDebrisC"] = cohort.bdall->y_sois.wdebrisc;
  data["WoodyDebrisN"] = cohort.bdall->y_sois.wdebrisn;
  // fluxes
  data["MossDeathC"] = cohort.bdall->y_v2soi.mossdeathc;
  data["MossdeathNitrogen"] = cohort.bdall->y_v2soi.mossdeathn;
  data["D2WoodyDebrisC"] = cohort.bdall->y_v2soi.d2wdebrisc;
  data["D2WoodyDebrisN"] = cohort.bdall->y_v2soi.d2wdebrisn;

  double t = 0.0;
  for (int il=0; il < MAX_SOI_LAY; il++) {
    t += cohort.bdall->m_soi2v.nextract[il];
  }
  data["NExtract"] = t;
  data["NetNMin"] = cohort.bdall->y_soi2soi.netnminsum;
  data["NetNImmob"] = cohort.bdall->y_soi2soi.nimmobsum;
  data["OrgNInput"] = cohort.bdall->y_a2soi.orgninput;
  data["AvlNInput"] = cohort.bdall->y_a2soi.avlninput;
  data["AvlNLost"] = cohort.bdall->y_soi2l.avlnlost;
  data["RHraw"] = cohort.bdall->y_soi2a.rhrawcsum;
  data["RHsoma"] = cohort.bdall->y_soi2a.rhsomasum;
  data["RHsompr"] = cohort.bdall->y_soi2a.rhsomprsum;
  data["RHsomcr"] = cohort.bdall->y_soi2a.rhsomcrsum;
  data["RH"] = cohort.bdall->y_soi2a.rhtot;
 
  //Placeholders for summing fire variables for the entire year
  double burnthick = 0.0, veg2airc = 0.0, veg2airn = 0.0, veg2soiabvvegc=0.0, veg2soiabvvegn = 0.0, veg2soiblwvegc = 0.0, veg2soiblwvegn = 0.0, veg2deadc = 0.0, veg2deadn = 0.0, soi2airc = 0.0, soi2airn = 0.0, air2soin = 0.0;
 
  for(int im=0; im<12; im++){
    char mth_chars[2];
    sprintf(mth_chars, "%02d", im);
    std::string mth_str = std::string(mth_chars);
    data["Fire"][mth_str]["Burnthick"] = cohort.year_fd[im].fire_soid.burnthick;
    data["Fire"][mth_str]["Veg2AirC"] = cohort.year_fd[im].fire_v2a.orgc;
    data["Fire"][mth_str]["Veg2AirN"] = cohort.year_fd[im].fire_v2a.orgn;
    data["Fire"][mth_str]["Veg2SoiAbvVegC"] = cohort.year_fd[im].fire_v2soi.abvc;
    data["Fire"][mth_str]["Veg2SoiBlwVegC"] = cohort.year_fd[im].fire_v2soi.blwc;
    data["Fire"][mth_str]["Veg2SoiAbvVegN"] = cohort.year_fd[im].fire_v2soi.abvn;
    data["Fire"][mth_str]["Veg2SoiBlwVegN"] = cohort.year_fd[im].fire_v2soi.blwn;
    data["Fire"][mth_str]["Veg2DeadC"] = cohort.year_fd[im].fire_v2dead.vegC;
    data["Fire"][mth_str]["Veg2DeadN"] = cohort.year_fd[im].fire_v2dead.strN;
    data["Fire"][mth_str]["Soi2AirC"] = cohort.year_fd[im].fire_soi2a.orgc;
    data["Fire"][mth_str]["Soi2AirN"] = cohort.year_fd[im].fire_soi2a.orgn;
    data["Fire"][mth_str]["Air2SoiN"] = cohort.year_fd[im].fire_a2soi.orgn/12;

    //Summed data for the entire year
    burnthick += cohort.year_fd[im].fire_soid.burnthick;
    veg2airc += cohort.year_fd[im].fire_v2a.orgc;
    veg2airn += cohort.year_fd[im].fire_v2a.orgn;
    veg2soiabvvegc += cohort.year_fd[im].fire_v2soi.abvc;
    veg2soiblwvegc += cohort.year_fd[im].fire_v2soi.blwc;
    veg2soiabvvegn += cohort.year_fd[im].fire_v2soi.abvn;
    veg2soiblwvegn += cohort.year_fd[im].fire_v2soi.blwn;
    veg2deadc += cohort.year_fd[im].fire_v2dead.vegC;
    veg2deadn += cohort.year_fd[im].fire_v2dead.strN;
    soi2airc += cohort.year_fd[im].fire_soi2a.orgc;
    soi2airn += cohort.year_fd[im].fire_soi2a.orgn;
    air2soin += cohort.year_fd[im].fire_a2soi.orgn/12;

  }
  data["Burnthick"] = burnthick;
  data["BurnVeg2AirC"] = veg2airc;
  data["BurnVeg2AirN"] = veg2airn;
  data["BurnVeg2SoiAbvVegC"] = veg2soiabvvegc;
  data["BurnVeg2SoiAbvVegN"] = veg2soiabvvegn;
  data["BurnVeg2SoiBlwVegN"] = veg2soiblwvegn;
  data["BurnVeg2SoiBlwVegC"] = veg2soiblwvegc;
  data["BurnSoi2AirC"] = soi2airc;
  data["BurnSoi2AirN"] = soi2airn;
  data["BurnAir2SoiN"] = air2soin;
  data["BurnAbvVeg2DeadC"] = veg2deadc;
  data["BurnAbvVeg2DeadN"] = veg2deadn;

  //Calculated sums
  double litterfallCsum = 0;
  double litterfallNsum = 0;

  for(int pft=0; pft<NUM_PFT; pft++) { //NUM_PFT
    char pft_chars[5];
    sprintf(pft_chars, "%d", pft);
    std::string pft_str = std::string(pft_chars);
    //c++0x equivalent: std::string pftvalue = std::to_string(pft);
    data["PFT" + pft_str]["VegCarbon"]["Leaf"] = cohort.bd[pft].y_vegs.c[I_leaf];
    data["PFT" + pft_str]["VegCarbon"]["Stem"] = cohort.bd[pft].y_vegs.c[I_stem];
    data["PFT" + pft_str]["VegCarbon"]["Root"] = cohort.bd[pft].y_vegs.c[I_root];
    data["PFT" + pft_str]["VegStructuralNitrogen"]["Leaf"] = cohort.bd[pft].y_vegs.strn[I_leaf];
    data["PFT" + pft_str]["VegStructuralNitrogen"]["Stem"] = cohort.bd[pft].y_vegs.strn[I_stem];
    data["PFT" + pft_str]["VegStructuralNitrogen"]["Root"] = cohort.bd[pft].y_vegs.strn[I_root];
    data["PFT" + pft_str]["VegLabileNitrogen"] = cohort.bd[pft].y_vegs.labn;

    data["PFT" + pft_str]["NAll"] = cohort.bd[pft].y_vegs.nall; // <-- Sum of labn and strn
    data["PFT" + pft_str]["StandingDeadC"] = cohort.bd[pft].y_vegs.deadc;
    data["PFT" + pft_str]["StandingDeadN"] = cohort.bd[pft].y_vegs.deadn;

    data["PFT" + pft_str]["NMobil"] = cohort.bd[pft].y_v2v.nmobilall; // <- the all denotes multi-compartment
    data["PFT" + pft_str]["NResorb"] = cohort.bd[pft].y_v2v.nresorball;

    data["PFT" + pft_str]["GPPAll"] = cohort.bd[pft].y_a2v.gppall;
    data["PFT" + pft_str]["GPP"]["Leaf"] = cohort.bd[pft].y_a2v.gpp[I_leaf];
    data["PFT" + pft_str]["GPP"]["Stem"] = cohort.bd[pft].y_a2v.gpp[I_stem];
    data["PFT" + pft_str]["GPP"]["Root"] = cohort.bd[pft].y_a2v.gpp[I_root];

    data["PFT" + pft_str]["NPPAll"] = cohort.bd[pft].y_a2v.nppall;
    data["PFT" + pft_str]["NPP"]["Leaf"] = cohort.bd[pft].y_a2v.npp[I_leaf];
    data["PFT" + pft_str]["NPP"]["Stem"] = cohort.bd[pft].y_a2v.npp[I_stem];
    data["PFT" + pft_str]["NPP"]["Root"] = cohort.bd[pft].y_a2v.npp[I_root];

    data["PFT" + pft_str]["GPPAllIgnoringNitrogen"] = cohort.bd[pft].y_a2v.ingppall;
    data["PFT" + pft_str]["NPPAllIgnoringNitrogen"] = cohort.bd[pft].y_a2v.innppall;

    data["PFT" + pft_str]["LitterfallCarbonAll"] = cohort.bd[pft].y_v2soi.ltrfalcall;
    data["PFT" + pft_str]["LitterfallCarbon"]["Leaf"] = cohort.bd[pft].y_v2soi.ltrfalc[I_leaf];
    data["PFT" + pft_str]["LitterfallCarbon"]["Stem"] = cohort.bd[pft].y_v2soi.ltrfalc[I_stem];
    data["PFT" + pft_str]["LitterfallCarbon"]["Root"] = cohort.bd[pft].y_v2soi.ltrfalc[I_root];

    data["PFT" + pft_str]["RespGrowth"]["Leaf"] = cohort.bd[pft].y_v2a.rg[I_leaf];
    data["PFT" + pft_str]["RespGrowth"]["Stem"] = cohort.bd[pft].y_v2a.rg[I_stem];
    data["PFT" + pft_str]["RespGrowth"]["Root"] = cohort.bd[pft].y_v2a.rg[I_root];

    data["PFT" + pft_str]["RespMaint"]["Leaf"] = cohort.bd[pft].y_v2a.rm[I_leaf];
    data["PFT" + pft_str]["RespMaint"]["Stem"] = cohort.bd[pft].y_v2a.rm[I_stem];
    data["PFT" + pft_str]["RespMaint"]["Root"] = cohort.bd[pft].y_v2a.rm[I_root];

    data["PFT" + pft_str]["LitterfallNitrogenPFT"] = cohort.bd[pft].y_v2soi.ltrfalnall;
    data["PFT" + pft_str]["LitterfallNitrogen"]["Leaf"] = cohort.bd[pft].y_v2soi.ltrfaln[I_leaf];
    data["PFT" + pft_str]["LitterfallNitrogen"]["Stem"] = cohort.bd[pft].y_v2soi.ltrfaln[I_stem];
    data["PFT" + pft_str]["LitterfallNitrogen"]["Root"] = cohort.bd[pft].y_v2soi.ltrfaln[I_root];

    data["PFT" + pft_str]["StNitrogenUptake"] = cohort.bd[pft].y_soi2v.snuptakeall;
    data["PFT" + pft_str]["InNitrogenUptake"] = cohort.bd[pft].y_soi2v.innuptake;
    data["PFT" + pft_str]["LabNitrogenUptake"] = cohort.bd[pft].y_soi2v.lnuptake;
    data["PFT" + pft_str]["TotNitrogenUptake"] = cohort.bd[pft].y_soi2v.snuptakeall + cohort.bd[pft].y_soi2v.lnuptake;

    data["PFT" + pft_str]["PARDown"] = cohort.ed[pft].y_a2v.pardown;
    data["PFT" + pft_str]["PARAbsorb"] = cohort.ed[pft].y_a2v.parabsorb;

    litterfallCsum += cohort.bd[pft].y_v2soi.ltrfalcall; 
    litterfallNsum += cohort.bd[pft].y_v2soi.ltrfalnall;

  }

  data["LitterfallCsum"] = litterfallCsum;
  data["LitterfallNsum"] = litterfallNsum;

  // Writes files like this:
  //  00000.json, 00001.json, 00002.json

  std::stringstream filename;
  filename.fill('0');
  filename << std::setw(5) << year << ".json";

  // Add the filename to the path
  p /= filename.str();

  out_stream.open(p.string().c_str(), std::ofstream::out);
  out_stream << data << std::endl;
  out_stream.close();

}

void Runner::output_debug_daily_drivers(int iy, boost::filesystem::path p) {

  // Writes files like this:
  //  year_00000_daily_drivers.text
  //  year_00001_daily_drivers.text
  //  year_00002_daily_drivers.text

  std::stringstream filename;
  filename.fill('0');
  filename << "year_" << std::setw(5) << iy << "_daily_drivers.text";

  // NOTE: (FIX?) This may not be the best place for these files as they are
  // not exactly the same format/layout as the normal "calibration" json files

  // Add the file name to the path
  p /= filename.str();

  std::ofstream out_stream;
  out_stream.open(p.string().c_str(), std::ofstream::out);

  out_stream << "tair_d = [" << temutil::vec2csv(cohort.climate.tair_d) << "]" << std::endl;
  out_stream << "nirr_d = [" << temutil::vec2csv(cohort.climate.nirr_d) << "]" << std::endl;
  out_stream << "vapo_d = [" << temutil::vec2csv(cohort.climate.vapo_d) << "]" << std::endl;
  out_stream << "prec_d = [" << temutil::vec2csv(cohort.climate.prec_d) << "]" << std::endl;
  out_stream << "rain_d = [" << temutil::vec2csv(cohort.climate.rain_d) << "]" << std::endl;
  out_stream << "snow_d = [" << temutil::vec2csv(cohort.climate.snow_d) << "]" << std::endl;
  out_stream << "svp_d = [" << temutil::vec2csv(cohort.climate.svp_d) << "]" << std::endl;
  out_stream << "vpd_d = [" << temutil::vec2csv(cohort.climate.vpd_d) << "]" << std::endl;
  out_stream << "girr_d = [" << temutil::vec2csv(cohort.climate.girr_d) << "]" << std::endl;
  out_stream << "cld_d = [" << temutil::vec2csv(cohort.climate.cld_d) << "]" << std::endl;
  out_stream << "par_d = [" << temutil::vec2csv(cohort.climate.par_d) << "]" << std::endl;

  out_stream.close();
}


void Runner::output_netCDF_monthly(int year, int month, std::string stage){
    BOOST_LOG_SEV(glg, debug)<<"NetCDF monthly output, year: "<<year<<" month: "<<month;
    output_netCDF(md.monthly_netcdf_outputs, year, month, stage);

    BOOST_LOG_SEV(glg, debug)<<"Outputting accumulated daily data on the monthly timestep";
    output_netCDF(md.daily_netcdf_outputs, year, month, stage);
}

void Runner::output_netCDF_yearly(int year, std::string stage){
    BOOST_LOG_SEV(glg, debug)<<"NetCDF yearly output, year: "<<year;
    output_netCDF(md.yearly_netcdf_outputs, year, 0, stage);
}

//The following two functions are the beginning of an attempt to
// generalize the output of different variables. The end goal is
// for the output_netCDF() function to be shorter, more readable,
// and have fewer redundant pieces of code.
//20181006 Currently there are only a few variables using these
// general functions, but since they include the front outputs
// we are merging this despite it being incomplete.
void Runner::output_nc_soil_layer(int ncid, int cv, int *data, int max_var_count, int start_timestep, int timesteps){

  //timestep, layer, row, col
  size_t soilstart[4];
  soilstart[0] = start_timestep;
  soilstart[1] = 0;
  soilstart[2] = this->y;
  soilstart[3] = this->x;

  size_t soilcount[4];
  //soilcount[0] = 1;
  soilcount[0] = timesteps;
  soilcount[1] = max_var_count;
  soilcount[2] = 1;
  soilcount[3] = 1;

  temutil::nc( nc_put_vara_int(ncid, cv, soilstart, soilcount, data) );
}

void Runner::output_nc_soil_layer(int ncid, int cv, double *data, int max_var_count, int start_timestep, int timesteps){

  //timestep, layer, row, col
  size_t soilstart[4];
  soilstart[0] = start_timestep;
  soilstart[1] = 0;
  soilstart[2] = this->y;
  soilstart[3] = this->x;

  size_t soilcount[4];
  //soilcount[0] = 1;
  soilcount[0] = timesteps;
  soilcount[1] = max_var_count;
  soilcount[2] = 1;
  soilcount[3] = 1;

  temutil::nc( nc_put_vara_double(ncid, cv, soilstart, soilcount, data) );
}

/** The wrapper makes decisions about whether to send the data to an IO slave
    or whether to write the data out here and now. */
template<typename T>
void Runner::io_wrapper(const std::string& vname,
                        const std::string& curr_filename,
                        const std::vector<size_t>& starts,
                        const std::vector<size_t>& counts,
                        const T& values) {

#ifdef WITHMPI

  int rank = MPI::COMM_WORLD.Get_rank();
  int ntasks = MPI::COMM_WORLD.Get_size();

  if (rank == 0) {
    //std::cout << "Normal old write...\n";
    temutil::write_var_to_netcdf(vname, curr_filename, starts, counts, values);
  } else {
    add_to_package_for_IO_slave(vname, curr_filename, starts, counts, values);
    //std::cout << "rank " << rank << " sending data to IO slave for " << curr_filename << "\n";
    // try this for starters...maybe the processes will overwrite eachother??
    //temutil::write_var_to_netcdf(vname, curr_filename, starts, counts, values);
  }
#else
  temutil::write_var_to_netcdf(vname, curr_filename, starts, counts, values);
#endif

}

template<typename T>
void Runner::add_to_package_for_IO_slave(const std::string & vname, 
                                         const std::string & curr_filename,
                                         const std::vector<size_t> & starts, 
                                         const std::vector<size_t> & counts, 
                                         const T & values) {
#ifdef WITHMPI

  // int id = MPI::COMM_WORLD.Get_rank();
  // int ntasks = MPI::COMM_WORLD.Get_size();

  OutputDataNugget odn = OutputDataNugget(curr_filename, vname, starts, counts, values);

  // Get a unique tag for this message (so boost serializing can break it into
  // multiple messages and be sure to correctly re-build it based on tag).
  int tag = temutil::get_uid(this->md.io_data_comm_ptr->rank());

  // Could set this based on this processes rank??
  // OR could set this based on the variable name??
  int designated_IO_process = 0;

  // SEND IT!
  this->md.io_data_comm_ptr->send(designated_IO_process, tag, odn);

  //BOOST_LOG_SEV(glg, debug) << "id: " << id << " (of " << ntasks << ") is sending an MPI message --to--> " << designated_io_slave << "\n";
  //std::cout << "id: " << id << " (of " << ntasks << ") is sending an MPI message --to--> " << designated_io_slave << "\n";
 
#else
  // pass - empty function body, this should never get called w/o MPI
#endif

}

void Runner::output_netCDF(std::map<std::string, OutputSpec> &netcdf_outputs, int year, int month, std::string stage){
  int month_timestep = year*12 + month;

  int day_timestep = year*365;
  for(int im=0; im<month; im++){
    day_timestep += DINM[im];
  }

  //For outputting subsets of driving data arrays
  int doy = temutil::day_of_year(month, 0); 

  std::string file_stage_suffix;
  if(stage.find("eq")!=std::string::npos){
    file_stage_suffix = "_eq.nc";
  }
  else if(stage.find("sp")!=std::string::npos){
    file_stage_suffix = "_sp.nc";
  }
  else if(stage.find("tr")!=std::string::npos){
    file_stage_suffix = "_tr.nc";
  }
  else if(stage.find("sc")!=std::string::npos){
    file_stage_suffix = "_sc.nc";
  }

  std::string curr_filename;

  int dinm = DINM[month];

  int rowidx = this->y;
  int colidx = this->x;

  OutputSpec curr_spec;
  int ncid;
  int timeD; //unlimited dimension
  int xD;
  int yD;
  int pftD;
  int pftpartD;
  int layerD;
  int cv; //reusable variable handle

  std::map<std::string, OutputSpec>::iterator map_itr;

  /*** 3D variables ***/
  std::vector<size_t> start3(3, 0);
  // Index 0 is set later
  start3[1] = rowidx;
  start3[2] = colidx;

  // For single value variables, an empty count vector
  std::vector<size_t> count0;

  // For daily-level variables
  std::vector<size_t> count3(3, 0);
  count3[0] = dinm;
  count3[1] = 1;
  count3[2] = 1;

  /*** Soil Variables ***/
  std::vector<size_t> soilstart4(4, 0);
  // Index 0 is set later
  soilstart4[1] = 0;
  soilstart4[2] = rowidx;
  soilstart4[3] = colidx;

  std::vector<size_t> soilcount4(4, 0);
  soilcount4[0] = 1;
  soilcount4[1] = MAX_SOI_LAY;
  soilcount4[2] = 1;
  soilcount4[3] = 1;

  /*** Fronts ***/
  std::vector<size_t> frontstart4(4, 0);
  // Index 0 is set later 
  frontstart4[1] = 0;
  frontstart4[2] = rowidx;
  frontstart4[3] = colidx;

  std::vector<size_t> frontcount4(4, 0);
  frontcount4[0] = 1; 
  frontcount4[1] = MAX_NUM_FNT;
  frontcount4[2] = 1;
  frontcount4[3] = 1;

  std::vector<size_t> frontstart5(5, 0);
  //Index 0 is set later   // time
  frontstart5[1] = 0;      // day of month
  frontstart5[2] = 0;      // front layer (index?)
  frontstart5[3] = rowidx;
  frontstart5[4] = colidx;

  std::vector<size_t> frontcount5(5, 0);
  frontcount5[0] = 1;
  frontcount5[1] = dinm;
  frontcount5[2] = MAX_NUM_FNT;
  frontcount5[3] = 1;
  frontcount5[4] = 1;

  /*** PFT variables ***/
  std::vector<size_t> PFTstart4(4, 0);
  // Index 0 is set later
  PFTstart4[1] = 0; // PFT
  PFTstart4[2] = rowidx;
  PFTstart4[3] = colidx;

  std::vector<size_t> PFTcount4(4, 0);
  PFTcount4[0] = 1;
  PFTcount4[1] = NUM_PFT;
  PFTcount4[2] = 1;
  PFTcount4[3] = 1;

  /*** Compartment variables ***/
  std::vector<size_t> CompStart4(4, 0);
  // Index 0 is set later
  CompStart4[1] = 0; // PFT compartment
  CompStart4[2] = rowidx;
  CompStart4[3] = colidx;

  std::vector<size_t> CompCount4(4, 0);
  CompCount4[0] = 1;
  CompCount4[1] = NUM_PFT_PART;
  CompCount4[2] = 1;
  CompCount4[3] = 1;

  /*** PFT and PFT compartment variables ***/
  std::vector<size_t> start5(5, 0);
  // Index 0 is set later
  start5[1] = 0; // PFT Compartment
  start5[2] = 0; // PFT
  start5[3] = rowidx;
  start5[4] = colidx;

  std::vector<size_t> count5(5, 0);
  count5[0] = 1;
  count5[1] = NUM_PFT_PART;
  count5[2] = NUM_PFT;
  count5[3] = 1;
  count5[4] = 1;

  /*** Single option vars: (year) ***/
  map_itr = netcdf_outputs.find("ALD");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputALD)
    {

    start3[0] = year;

    std::vector<double> values(1, cohort.edall->y_soid.ald); 
    io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputALD)
  }//end ALD
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("DEEPDZ");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDEEPDZ)
    {

      start3[0] = year;

      double deepdz = 0;
      Layer* currL = cohort.ground.toplayer;
      while(currL!=NULL){
        if(currL->isHumic){
          deepdz += currL->dz;
        }
        currL = currL->nextl;
      }

      std::vector<double> values(1, deepdz);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputDEEPDZ)
  }//end DEEPDZ
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("GROWEND");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputGROWEND)
    {

      start3[0] = year;

      std::vector<double> values(1, cohort.edall->y_soid.rtdpGEoutput);

      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputGROWEND)
  }//end GROWEND
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("GROWSTART");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputGROWSTART)
    {

      start3[0] = year;

      std::vector<double> growstart(1, cohort.edall->y_soid.rtdpGSoutput);
      io_wrapper(svname, curr_filename, start3, count0, growstart);

    }//end critical(outputGROWSTART)
  }//end GROWSTART
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("MOSSDZ");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputMOSSDZ)
    {
      start3[0] = year;

      std::vector<double> mossdz(1, 0);
      Layer* currL = cohort.ground.toplayer;
      while(currL!=NULL){
        if(currL->isMoss){
          mossdz[0] += currL->dz;
        }
        currL = currL->nextl;
      }
      //The following may never get set to anything useful?
      //y_soil.mossthick;

      io_wrapper(svname, curr_filename, start3, count0, mossdz);

    }//end critical(outputMOSSDZ)
  }//end MOSSDZ
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("ROLB");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputROLB)
    {
      double rolb;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        rolb = cohort.year_fd[month].fire_soid.rolb;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        //FIX. This will not work if there is more than one fire per year
        // What does yearly ROLB even mean with multiple fires?
        rolb = cohort.fd->fire_soid.rolb;
      }

      std::vector<double> values(1, rolb);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputROLB)
  }//end ROLB
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("PERMAFROST");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputPERMAFROST)
    {
      start3[0] = year;
      std::vector<double> permafrost(1, cohort.edall->y_soid.permafrost);
      io_wrapper(svname, curr_filename, start3, count0, permafrost);
    }//end critical(outputPERMAFROST)
  }//end PERMAFROST
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("SHLWDZ");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSHLWDZ)
    {
      start3[0] = year;

      double shlwdz = 0;
      Layer* currL = cohort.ground.toplayer;
      while(currL!=NULL){
        if(currL->isFibric){
          shlwdz += currL->dz;
        }
        currL = currL->nextl;
      }
      std::vector<double> values(1, shlwdz);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputSHLWDZ)
  }//end SHLWDZ
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("SNOWEND");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSNOWEND)
    {
      start3[0] = year;
      std::vector<double> snowend(1, cohort.edall->y_snws.snowend);
      io_wrapper(svname, curr_filename, start3, count0, snowend);
    }//end critical(outputSNOWEND)
  }//end SNOWEND
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("SNOWSTART");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSNOWSTART)
    {
      start3[0] = year;
      std::vector<double> snowstart(1, cohort.edall->y_snws.snowstart);
      io_wrapper(svname, curr_filename, start3, count0, snowstart);

    }//end critical(outputSNOWSTART)
  }//end SNOWSTART
  map_itr = netcdf_outputs.end();


  //Years since disturbance
  map_itr = netcdf_outputs.find("YSD");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputYSD)
    {
      start3[0] = year;
      std::vector<double> ysd(1, cohort.cd.yrsdist);
      io_wrapper(svname, curr_filename, start3, count0, ysd);
    }//end critical(outputYSD)
  }//end YSD
  map_itr = netcdf_outputs.end();


  /*** Two combination vars: (month, year) ***/
  map_itr = netcdf_outputs.find("BURNAIR2SOIN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNAIR2SOIN)
    {

      double burnair2soin;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnair2soin = cohort.year_fd[month].fire_a2soi.orgn;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnair2soin = 0;
        for(int im=0; im<12; im++){
          burnair2soin += cohort.year_fd[im].fire_a2soi.orgn;
        }
      }

      std::vector<double> values(1, burnair2soin);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNAIR2SOIN)
  }//end BURNAIR2SOIN
  map_itr = netcdf_outputs.end();


  //Burn thickness
  map_itr = netcdf_outputs.find("BURNTHICK");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNTHICK)
    {

      double burnthick;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnthick = cohort.year_fd[month].fire_soid.burnthick;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnthick = 0;
        for(int im=0; im<12; im++){
          burnthick += cohort.year_fd[im].fire_soid.burnthick;
        }
      }
      std::vector<double> values(1, burnthick);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputBURNTHICK)
  }//end BURNTHICK
  map_itr = netcdf_outputs.end();


  //Standing dead C
  map_itr = netcdf_outputs.find("DEADC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDEADC)
    {
      double deadc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        deadc = cohort.bdall->m_vegs.deadc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        deadc = cohort.bdall->y_vegs.deadc;
      }
      std::vector<double> values(1, deadc);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputDEADC)
  }//end DEADC
  map_itr = netcdf_outputs.end();


  //Standing dead N
  map_itr = netcdf_outputs.find("DEADN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDEADN)
    {
      double deadn;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        deadn = cohort.bdall->m_vegs.deadn;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        deadn = cohort.bdall->y_vegs.deadn;
      }
      std::vector<double> values(1, deadn);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputDEADN)
  }//end DEADN
  map_itr = netcdf_outputs.end();


  //Deep C
  map_itr = netcdf_outputs.find("DEEPC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDEEPC)
    {
      double deepc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        deepc = cohort.bdall->m_soid.deepc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        deepc = cohort.bdall->y_soid.deepc;
      }
      std::vector<double> values(1, deepc);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputDEEPC)
  }//end DEEPC
  map_itr = netcdf_outputs.end();


  //Woody debris C
  map_itr = netcdf_outputs.find("DWDC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDWDC)
    {

      double woodyc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        woodyc = cohort.bdall->m_sois.wdebrisc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        woodyc = cohort.bdall->y_sois.wdebrisc;
      }
      std::vector<double> values(2, woodyc);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputDWDC)
  }//end DWDC
  map_itr = netcdf_outputs.end();


  //Woody debris N
  map_itr = netcdf_outputs.find("DWDN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDWDN)
    {

      double woodyn;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        woodyn = cohort.bdall->m_sois.wdebrisn;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        woodyn = cohort.bdall->y_sois.wdebrisn;
      }
      std::vector<double> values(1, woodyn);
      io_wrapper(svname, curr_filename, start3, count0, values);
      
    }//end critical(outputDWDN)
  }//end DWDN
  map_itr = netcdf_outputs.end();


  //Mineral C
  map_itr = netcdf_outputs.find("MINEC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputMINEC)
    {
      double minec;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        minec = cohort.bdall->m_soid.mineac
                + cohort.bdall->m_soid.minebc
                + cohort.bdall->m_soid.minecc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        minec = cohort.bdall->y_soid.mineac
                + cohort.bdall->y_soid.minebc
                + cohort.bdall->y_soid.minecc;
      }
      std::vector<double> values(1, minec);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputMINEC)
  }//end MINEC
  map_itr = netcdf_outputs.end();


  //MOSSDEATHC
  map_itr = netcdf_outputs.find("MOSSDEATHC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputMOSSDEATHC)
    {
      double mossdeathc = 0;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        mossdeathc = cohort.bdall->m_v2soi.mossdeathc; 
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        mossdeathc = cohort.bdall->y_v2soi.mossdeathc; 
      }
      std::vector<double> values(1, mossdeathc);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputMOSSDEATHC)
  }//end MOSSDEATHC
  map_itr = netcdf_outputs.end();


  //MOSSDEATHN
  map_itr = netcdf_outputs.find("MOSSDEATHN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputMOSSDEATHN)
    {
      double mossdeathn = 0;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        mossdeathn = cohort.bdall->m_v2soi.mossdeathn; 
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        mossdeathn = cohort.bdall->y_v2soi.mossdeathn; 
      }
      std::vector<double> values(1, mossdeathn);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputMOSSDEATHN)
  }//end MOSSDEATHN
  map_itr = netcdf_outputs.end();


  //QDRAINAGE
  map_itr = netcdf_outputs.find("QDRAINAGE");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputQDRAINAGE)
    {
      double qdrainage = 0;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        qdrainage = cohort.edall->m_soi2l.qdrain; 
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        qdrainage = cohort.edall->y_soi2l.qdrain; 
      }
      std::vector<double> values(1, qdrainage);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputQDRAINAGE)
  }//end QDRAINAGE 
  map_itr = netcdf_outputs.end();


  //QINFILTRATION
  map_itr = netcdf_outputs.find("QINFILTRATION");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputQINFILTRATION)
    {
      double qinfil = 0;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        qinfil = cohort.edall->m_soi2l.qinfl; 
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        qinfil = cohort.edall->y_soi2l.qinfl; 
      }
      std::vector<double> values(1, qinfil);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputQINFILTRATION)
  }//end QINFILTRATION
  map_itr = netcdf_outputs.end();


  //QRUNOFF
  map_itr = netcdf_outputs.find("QRUNOFF");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputQRUNOFF)
    {
      double qrunoff = 0;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        qrunoff = cohort.edall->m_soi2l.qover; 
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        qrunoff = cohort.edall->y_soi2l.qover;
      }
      std::vector<double> values(1, qrunoff);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputQRUNOFF)
  }//end QRUNOFF 
  map_itr = netcdf_outputs.end();


  //Shallow C
  map_itr = netcdf_outputs.find("SHLWC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSHLWC)
    {
      double shlwc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        shlwc = cohort.bdall->m_soid.shlwc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        shlwc = cohort.bdall->y_soid.shlwc;
      }
      std::vector<double> values(1, shlwc);
      io_wrapper(svname, curr_filename, start3, count0, values);
       
    }//end critical(outputSHLWC)
  }//end SHLWC 
  map_itr = netcdf_outputs.end();


  //Woody debris RH
  map_itr = netcdf_outputs.find("WDRH");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputWDRH)
    {
      double woodyrh;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        woodyrh = cohort.bdall->m_soi2a.rhwdeb;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        woodyrh = cohort.bdall->y_soi2a.rhwdeb;
      }
      std::vector<double> values(1, woodyrh);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputWDRH)
  }//end WDRH
  map_itr = netcdf_outputs.end();


  /*** Three combination vars: (year, month, day) ***/
  //HKDEEP
  map_itr = netcdf_outputs.find("HKDEEP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputHKDEEP)
    {
      if(curr_spec.daily){      
        start3[0] = day_timestep; 
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_hkdeep[id];
        } 
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.hkdeep);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.hkdeep);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }


    }//end critical(outputHKDEEP)
  }//end HKDEEP 
  map_itr = netcdf_outputs.end();


  //HKLAYER
  map_itr = netcdf_outputs.find("HKLAYER");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputHKLAYER)
    {
      std::vector<double> values(MAX_SOI_LAY, -99999);

      for(int il=0; il<MAX_SOI_LAY; il++){
        if(curr_spec.monthly){
          soilstart4[0] = month_timestep;
          values[il] = cohort.edall->m_soid.hcond[il];
        }
        else if(curr_spec.yearly){
          soilstart4[0] = year;
            values[il] = cohort.edall->y_soid.hcond[il];
        }
      }

      io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);

    }//end critical(outputHKLAYER)
  }//end HKLAYER
  map_itr = netcdf_outputs.end();


  //HKMINEA
  map_itr = netcdf_outputs.find("HKMINEA");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputHKMINEA)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_hkminea[id];
        } 
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.hkminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.hkminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputHKMINEA)
  }//end HKMINEA
  map_itr = netcdf_outputs.end();


  //HKMINEB
  map_itr = netcdf_outputs.find("HKMINEB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputHKMINEB)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_hkmineb[id];
        } 
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.hkmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.hkmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputHKMINEB)
  }//end HKMINEB
  map_itr = netcdf_outputs.end();


  //HKMINEC
  map_itr = netcdf_outputs.find("HKMINEC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputHKMINEC)
    {
      double hkminec;
      if(curr_spec.daily){
        start3[0] = day_timestep; 
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_hkminec[id];
        } 
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.hkminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.hkminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputHKMINEC)
  }//end HKMINEC
  map_itr = netcdf_outputs.end();


  //HKSHLW
  map_itr = netcdf_outputs.find("HKSHLW");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputHKSHLW)
    {
      double hkshlw;
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for(int id=0; id<dinm; id++){
          values[id] = cohort.edall->daily_hkshlw[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.hkshlw);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.hkshlw);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }
  }//end HKSHLW 
  map_itr = netcdf_outputs.end();


  //Snowthick - a snapshot of the time when output is called
  map_itr = netcdf_outputs.find("SNOWTHICK");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSNOWTHICK)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_snowthick[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        double snowthick=0.0;
        Layer* currL = cohort.ground.toplayer;
        while(currL->isSnow){
          snowthick += currL->dz;
          currL = currL->nextl;
        }
        std::vector<double> values(1, snowthick);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        double snowthick=0.0;
        Layer* currL = cohort.ground.toplayer;
        while(currL->isSnow){
          snowthick += currL->dz;
          currL = currL->nextl;
        }
        std::vector<double> values(1, snowthick);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputSNOWTHICK)
  }//end SNOWTHICK
  map_itr = netcdf_outputs.end();


  //SNOWFALL
  map_itr = netcdf_outputs.find("SNOWFALL");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSNOWFALL)
    {
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_a2l.snfl);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_a2l.snfl);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

    }//end critical(outputSNOWFALL)
  }//end SNOWFALL
  map_itr = netcdf_outputs.end();


  //DRIVINGSNOWFALL (this is not really driving data as it is calculated??)
  map_itr = netcdf_outputs.find("DRIVINGSNOWFALL");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDRIVINGSNOWFALL)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        int z = 0;
        for (int id=doy; id<(doy+dinm); id++) {
          values[z] = cohort.climate.snow_d[id];
          z += 1;
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
    }//end critical(outputDRIVINGSNOWFALL)
  }//end DRIVINGSNOWFALL
  map_itr = netcdf_outputs.end();


  //DRIVINGRAINFALL
  map_itr = netcdf_outputs.find("DRIVINGRAINFALL");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDRIVINGRAINFALL)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        int z = 0;
        for(int id=doy; id<(doy+dinm); id++) {
          values[z] = cohort.climate.rain_d[id];
          z += 1;
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
    }//end critical(outputDRIVINGRAINFALL)
  }//end DRIVINGRAINFALL
  map_itr = netcdf_outputs.end();


  //DRIVINGTAIR
  map_itr = netcdf_outputs.find("DRIVINGTAIR");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDRIVINGTAIR)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        int z = 0;
        for(int id=doy; id<(doy+dinm); id++) {
          values[z] = cohort.climate.tair_d[id];
          z += 1;
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
    }//end critical(outputDRIVINGTAIR)
  }//end DRIVINGTAIR
  map_itr = netcdf_outputs.end();


  //DRIVINGVAPO
  map_itr = netcdf_outputs.find("DRIVINGVAPO");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDRIVINGVAPO)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        int z = 0;
        for(int id=doy; id<(doy+dinm); id++) {
          values[z] = cohort.climate.vapo_d[id];
          z += 1;
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
    }//end critical(outputDRIVINGVAPO)
  }//end DRIVINGVAPO
  map_itr = netcdf_outputs.end();


  //DRIVINGNIRR
  map_itr = netcdf_outputs.find("DRIVINGNIRR");

  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputDRIVINGNIRR)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        int z = 0;
        for(int id=doy; id<(doy+dinm); id++) {
          values[z] = cohort.climate.nirr_d[id];
          z += 1;
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
    }//end critical(outputDRIVINGNIRR)
  }//end DRIVINGNIRR
  map_itr = netcdf_outputs.end();


  //RAINFALL
  map_itr = netcdf_outputs.find("RAINFALL");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputRAINFALL)
    {
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_a2l.rnfl);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_a2l.rnfl);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputRAINFALL)
  }//end RAINFALL
  map_itr = netcdf_outputs.end();


  //Snow water equivalent - a snapshot of the time when output is called
  map_itr = netcdf_outputs.find("SWE");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSWE)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++){
          values[id] = cohort.edall->daily_swesum[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }

      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_snws.swesum);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_snws.swesum);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

    }//end critical(outputSWE)
  }//end SWE
  map_itr = netcdf_outputs.end();


  //TCDEEP
  map_itr = netcdf_outputs.find("TCDEEP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTCDEEP)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tcdeep[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tcdeep);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tcdeep);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

    }//end critical(outputTCDEEP)
  }//end TCDEEP
  map_itr = netcdf_outputs.end();


  //TCLAYER
  map_itr = netcdf_outputs.find("TCLAYER");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTCLAYER)
    {

      std::vector<double> values(MAX_SOI_LAY, -99999);

      if(curr_spec.monthly){
        soilstart4[0] = month_timestep;
        for (int il=0; il<MAX_SOI_LAY; il++){
          values[il] = cohort.edall->m_soid.tcond[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
      else if(curr_spec.yearly){
        soilstart4[0] = year;
        for (int il=0; il<MAX_SOI_LAY; il++){
          values[il] = cohort.edall->y_soid.tcond[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }

    }//end critical(outputTCLAYER)
  }//end TCLAYER
  map_itr = netcdf_outputs.end();


  //TCMINEA
  map_itr = netcdf_outputs.find("TCMINEA");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTCMINEA)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tcminea[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tcminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tcminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTCMINEA)
  }//end TCMINEA
  map_itr = netcdf_outputs.end();


  //TCMINEB
  map_itr = netcdf_outputs.find("TCMINEB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTCMINEB)
    {
      double tcmineb;
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tcmineb[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tcmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tcmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTCMINEB)
  }//end TCMINEB
  map_itr = netcdf_outputs.end();


  //TCMINEC
  map_itr = netcdf_outputs.find("TCMINEC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTCMINEC)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tcminec[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tcminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tcminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTCMINEC)
  }//end TCMINEC
  map_itr = netcdf_outputs.end();


  //TCSHLW
  map_itr = netcdf_outputs.find("TCSHLW");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTCSHLW)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tcshlw[id];
        } 
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> tcshlw(1, cohort.edall->m_soid.tcshlw);
        io_wrapper(svname, curr_filename, start3, count0, tcshlw);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> tcshlw(1, cohort.edall->y_soid.tcshlw);
        io_wrapper(svname, curr_filename, start3, count0, tcshlw);
      }
      
    }//end critical(outputTCSHLW)
  }//end TCSHLW
  map_itr = netcdf_outputs.end();


  //TDEEP
  map_itr = netcdf_outputs.find("TDEEP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTDEEP)
    {
      double tdeep;
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tdeep[id];
        } 
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> tdeep(1, cohort.edall->m_soid.tdeep);
        io_wrapper(svname, curr_filename, start3, count0, tdeep);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> tdeep(1, cohort.edall->y_soid.tdeep);
        io_wrapper(svname, curr_filename, start3, count0, tdeep);
      }
      
    }//end critical(outputTDEEP)
  }//end TDEEP
  map_itr = netcdf_outputs.end();


  //TLAYER
  map_itr = netcdf_outputs.find("TLAYER");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTLAYER)
    {
      std::vector<double> values(MAX_SOI_LAY, -99999);
      if(curr_spec.monthly){
        soilstart4[0] = month_timestep;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.edall->m_sois.ts[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
      else if(curr_spec.yearly){
        soilstart4[0] = year;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.edall->y_sois.ts[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
    }//end critical(outputTLAYER)
  }//end TLAYER
  map_itr = netcdf_outputs.end();


  // FRONTSTYPE
  map_itr = netcdf_outputs.find("FRONTSTYPE");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputFRONTSTYPE)
    {

    // This uses the summary structs, but might be more accurate
    // if the deque of fronts was checked directly.
    if(curr_spec.daily){

      // BROKEN!!!! Gives netcdf "start+count exceeds dimension bounds" 
      frontstart5[0] = day_timestep;
      std::vector<double> values(dinm*MAX_NUM_FNT, -99999);

      for (int day_of_month=0; day_of_month<dinm; day_of_month++) {
        for (int fnt_idx=0; fnt_idx<MAX_NUM_FNT; fnt_idx++){
          int idx = day_of_month*MAX_NUM_FNT + fnt_idx;
          values[idx] = cohort.edall->daily_frontstype[day_of_month][fnt_idx];
        }
      }
      //output_nc_soil_layer(ncid, cv, &cohort.edall->daily_frontstype[0][0], MAX_NUM_FNT, day_timestep, dinm);
      io_wrapper(svname, curr_filename, frontstart5, frontcount5, values);

    }
    else if(curr_spec.monthly){
      frontstart4[0] = month_timestep;
      std::vector<double> values(MAX_NUM_FNT, -99999);
      for (int i=0; i < MAX_NUM_FNT; i++){
        values[i] = cohort.ground.frnttype[i];
      }
      io_wrapper(svname, curr_filename, frontstart4, frontcount4, values);
      //temutil::nc( nc_put_vara_int(ncid, cv, soilstart4, frontcount4, &cohort.ground.frnttype[0]) );
      //nc_output_soil(ncid, cv, &cohort.ground.frnttype[0], MAX_NUM_FNT, month_timestep)
    }
    else if(curr_spec.yearly){
      frontstart4[0] = year;
      std::vector<double> values(MAX_NUM_FNT, -99999);
      for (int i=0; i < MAX_NUM_FNT; i++){
        values[i] = cohort.ground.frnttype[i];
      }
      io_wrapper(svname, curr_filename, frontstart4, frontcount4, values);
      //temutil::nc( nc_put_vara_int(ncid, cv, soilstart4, frontcount4, &cohort.ground.frnttype[0]) );
      //nc_output_soil(ncid, cv, &cohort.ground.frnttype[0], MAX_NUM_FNT, year)
    }

    }//end critical(outputFRONTSTYPE)
  }//end FRONTSTYPE
  map_itr = netcdf_outputs.end();


  // FRONTSDEPTH
  map_itr = netcdf_outputs.find("FRONTSDEPTH");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputFRONTSDEPTH)
    {

    // This uses the summary structs, but might be more accurate
    // if the deque of fronts was checked directly.
    if(curr_spec.daily){
      output_nc_soil_layer(ncid, cv, &cohort.edall->daily_frontsdepth[0][0], MAX_NUM_FNT, day_timestep, dinm);

      // Likely BROKEN!!!! Gives netcdf "start+count exceeds dimension bounds" 
      frontstart5[0] = day_timestep;
      std::vector<double> values(dinm*MAX_NUM_FNT, -99999);
      for (int day_of_month=0; day_of_month<dinm; day_of_month++) {
        for (int fnt_idx=0; fnt_idx<MAX_NUM_FNT; fnt_idx++){
          int idx = day_of_month*MAX_NUM_FNT + fnt_idx;
          values[idx] = cohort.edall->daily_frontsdepth[day_of_month][fnt_idx];
        }
      }
      io_wrapper(svname, curr_filename, frontstart5, frontcount5, values);
    }
    else if(curr_spec.monthly){
      frontstart4[0] = month_timestep;
      std::vector<double> values(MAX_NUM_FNT, -99999);
      for (int i=0; i<MAX_NUM_FNT; i++){
        values[i] = cohort.ground.frntz[i];
      }
      io_wrapper(svname, curr_filename, frontstart4, frontcount4, values);
    }
    else if(curr_spec.yearly){
      frontstart4[0] = year;
      std::vector<double> values(MAX_NUM_FNT, -99999);
      for (int i=0; i<MAX_NUM_FNT; i++){
        values[i] = cohort.ground.frntz[i];
      }
      io_wrapper(svname, curr_filename, frontstart4, frontcount4, values);
    }
   }//end critical(outputFRONTSDEPTH)
  }//end FRONTSDEPTH
  map_itr = netcdf_outputs.end();


  //TMINEA
  map_itr = netcdf_outputs.find("TMINEA");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTMINEA)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tminea[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTMINEA)
  }//end TMINEA
  map_itr = netcdf_outputs.end();


  //TMINEB
  map_itr = netcdf_outputs.find("TMINEB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTMINEB)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tmineb[id];
        }

        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTMINEB)
  }//end TMINEB
  map_itr = netcdf_outputs.end();


  //TMINEC
  map_itr = netcdf_outputs.find("TMINEC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTMINEC)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tminec[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTMINEC)
  }//end TMINEC
  map_itr = netcdf_outputs.end();


  //TSHLW
  map_itr = netcdf_outputs.find("TSHLW");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputTSHLW)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_tshlw[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.tshlw);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.tshlw);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputTSHLW)
  }//end TSHLW
  map_itr = netcdf_outputs.end();


  //VWCDEEP
  map_itr = netcdf_outputs.find("VWCDEEP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVWCDEEP)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_vwcdeep[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.vwcdeep);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.vwcdeep);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputVWCDEEP)
  }//end VWCDEEP
  map_itr = netcdf_outputs.end();


  //VWCLAYER
  map_itr = netcdf_outputs.find("VWCLAYER");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVWCLAYER)
    {
      std::vector<double> values(MAX_SOI_LAY, -99999);
      if(curr_spec.monthly){
        soilstart4[0] = month_timestep;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.edall->m_soid.vwc[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
      else if(curr_spec.yearly){
        soilstart4[0] = year;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.edall->y_soid.vwc[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
    }//end critical(outputVWCLAYER)
  }//end VWCLAYER
  map_itr = netcdf_outputs.end();


  //VWCMINEA
  map_itr = netcdf_outputs.find("VWCMINEA");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVWCMINEA)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_vwcminea[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.vwcminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.vwcminea);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputVWCMINEA)
  }//end VWCMINEA
  map_itr = netcdf_outputs.end();


  //VWCMINEB
  map_itr = netcdf_outputs.find("VWCMINEB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVWCMINEB)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_vwcmineb[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.vwcmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.vwcmineb);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputVWCMINEB)
  }//end VWCMINEB
  map_itr = netcdf_outputs.end();


  //VWCMINEC
  map_itr = netcdf_outputs.find("VWCMINEC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVWCMINEC)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_vwcminec[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.vwcminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.vwcminec);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputVWCMINEC)
  }//end VWCMINEC
  map_itr = netcdf_outputs.end();


  //VWCSHLW
  map_itr = netcdf_outputs.find("VWCSHLW");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVWCSHLW)
    {
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_vwcshlw[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_soid.vwcshlw);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_soid.vwcshlw);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputVWCSHLW)
  }//end VWCSHLW
  map_itr = netcdf_outputs.end();


  //WATERTAB
  map_itr = netcdf_outputs.find("WATERTAB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputWATERTAB)
    {
      double watertab;
      if(curr_spec.daily){
        start3[0] = day_timestep;
        std::vector<double> values(dinm, -99999);
        for (int id=0; id<dinm; id++) {
          values[id] = cohort.edall->daily_watertab[id];
        }
        io_wrapper(svname, curr_filename, start3, count3, values);
      }
      else if(curr_spec.monthly){
        start3[0] = month_timestep;
        std::vector<double> values(1, cohort.edall->m_sois.watertab);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        std::vector<double> values(1, cohort.edall->y_sois.watertab);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputWATERTAB)
  }//end WATERTAB
  map_itr = netcdf_outputs.end();


  /*** Three combination vars. (year, month)x(layer) ***/
  //LAYERDEPTH
  map_itr = netcdf_outputs.find("LAYERDEPTH");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputLAYERDEPTH)
    {
      std::vector<double> values(MAX_SOI_LAY, -99999);
      if(curr_spec.monthly){
        soilstart4[0] = month_timestep;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.cd.m_soil.z[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
      else if(curr_spec.yearly){
        soilstart4[0] = year;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.cd.y_soil.z[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
    }//end critical(outputLAYERDEPTH)
  }//end LAYERDEPTH 
  map_itr = netcdf_outputs.end();


  //LAYERDZ
  map_itr = netcdf_outputs.find("LAYERDZ");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputLAYERDZ)
    {
      std::vector<double> values(MAX_SOI_LAY, -99999);
      if(curr_spec.monthly){
        soilstart4[0] = month_timestep;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.cd.m_soil.dz[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
      else if(curr_spec.yearly){
        soilstart4[0] = year;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.cd.y_soil.dz[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
    }//end critical(outputLAYERDZ)
  }//end LAYERDZ 
  map_itr = netcdf_outputs.end();


  //LAYERTYPE
  map_itr = netcdf_outputs.find("LAYERTYPE");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputLAYERTYPE)
    {
      std::vector<double> values(MAX_SOI_LAY, -99999);
      if(curr_spec.monthly){
        soilstart4[0] = month_timestep;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.cd.m_soil.type[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
      else if(curr_spec.yearly){
        soilstart4[0] = year;
        for (int il=0; il<MAX_SOI_LAY; il++) {
          values[il] = cohort.cd.y_soil.type[il];
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
      }
    }//end critical(outputLAYERTYPE)
  }//end LAYERTYPE 
  map_itr = netcdf_outputs.end();


  /*** Four combination vars. (year, month)x(layer, tot)  ***/
  //AVLN
  map_itr = netcdf_outputs.find("AVLN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputAVLN)
    {
      if(curr_spec.layer){
        std::vector<double> avln(MAX_SOI_LAY, -99999);
        if(curr_spec.monthly){
          soilstart4[0] = month_timestep;
          for(int il=0; il<MAX_SOI_LAY; il++){
            avln[il] = cohort.bdall->m_sois.avln[il];
          }
        }
        else if(curr_spec.yearly){
          soilstart4[0] = year;
          for(int il=0; il<MAX_SOI_LAY; il++){
            avln[il] = cohort.bdall->y_sois.avln[il];
          }
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, avln);
      }
      else if(!curr_spec.layer){

        double avln;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          avln = cohort.bdall->m_soid.avlnsum;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          avln = cohort.bdall->y_soid.avlnsum;
        }
        std::vector<double> values(1, avln);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputAVLN)
  }//end AVLN
  map_itr = netcdf_outputs.end();


  //Burned soil carbon
  map_itr = netcdf_outputs.find("BURNSOIC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNSOIC)
    {
      if(curr_spec.layer){
        /*** STUB ***/
        //By-layer may not be feasible yet.
      }
      else if(!curr_spec.layer){

        double burnsoilC = 0.0;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          burnsoilC = cohort.year_fd[month].fire_soi2a.orgc;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          burnsoilC = 0;
          for(int im=0; im<12; im++){
            burnsoilC += cohort.year_fd[im].fire_soi2a.orgc;
          }
        }
        std::vector<double> values(1, burnsoilC);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputBURNSOIC)
  }//end BURNSOIC
  map_itr = netcdf_outputs.end();


  //Burned soil nitrogen 
  map_itr = netcdf_outputs.find("BURNSOILN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNSOILN)
    {
      if(curr_spec.layer){
        /*** STUB ***/
      }
      else if(!curr_spec.layer){

        double burnSoilN = 0.0;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          burnSoilN = cohort.year_fd[month].fire_soi2a.orgn;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          burnSoilN = 0;
          for(int im=0; im<12; im++){
            burnSoilN += cohort.year_fd[im].fire_soi2a.orgn;
          }
        }
        std::vector<double> values(1, burnSoilN);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputBURNSOILN)
  }//end BURNSOILN
  map_itr = netcdf_outputs.end();


  //NDRAIN
  map_itr = netcdf_outputs.find("NDRAIN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNDRAIN)
    {
      if(curr_spec.layer){

        if(curr_spec.monthly){
          soilstart4[0] = month_timestep;
          std::vector<double> values(MAX_SOI_LAY, -99999);
          for (int il=0; il<MAX_SOI_LAY; il++) {
            /*** STUB ??? ***/
            //values[il] = cohort.soilbgc.bd->m_soi2l.ndrain[il];
          }
          io_wrapper(svname, curr_filename, soilstart4, soilcount4, values);
        }
        else if(curr_spec.yearly){
          /*** STUB ***/
        }
      }
      else if(!curr_spec.layer){

        double ndrain = 0;
        if(curr_spec.monthly){
          start3[0] = month_timestep;

          for(int il=0; il<MAX_SOI_LAY; il++){
            /*** STUB ***/
            //ndrain += bd->m_soi2l.ndrain[il];
          }

        }
        if(curr_spec.yearly){
          start3[0] = year;
          /*** STUB ***/
        }
        std::vector<double> values(1, ndrain);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputNDRAIN)
  }//end NDRAIN
  map_itr = netcdf_outputs.end();


  //NETNMIN
  map_itr = netcdf_outputs.find("NETNMIN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNETNMIN)
    {
      if(curr_spec.layer){

        std::vector<double> netnmin(MAX_SOI_LAY, -99999);
        for(int il=0; il<MAX_SOI_LAY; il++){
          if(curr_spec.monthly){
            soilstart4[0] = month_timestep;
            netnmin[il] = cohort.bdall->m_soi2soi.netnmin[il];
          }
          else if(curr_spec.yearly){
            soilstart4[0] = year;
            netnmin[il] = cohort.bdall->y_soi2soi.netnmin[il];
          }
        }

        io_wrapper(svname, curr_filename, soilstart4, soilcount4, netnmin);
      }
      //Total, instead of by layer
      else if(!curr_spec.layer){

        double netnmin;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          netnmin = cohort.bdall->m_soi2soi.netnminsum;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          netnmin = cohort.bdall->y_soi2soi.netnminsum;
        }
        std::vector<double> values(1, netnmin);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputNETNMIN)
  }//end NETNMIN
  map_itr = netcdf_outputs.end();


  //NIMMOB
  map_itr = netcdf_outputs.find("NIMMOB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNIMMOB)
    {
      if(curr_spec.layer){

        std::vector<double> nimmob(MAX_SOI_LAY, -99999);
        for(int il=0; il<MAX_SOI_LAY; il++){
          if(curr_spec.monthly){
            soilstart4[0] = month_timestep;
            nimmob[il] = cohort.bdall->m_soi2soi.nimmob[il];
          }
          else if(curr_spec.yearly){
            soilstart4[0] = year;
            nimmob[il] = cohort.bdall->y_soi2soi.nimmob[il];
          }
        }

        io_wrapper(svname, curr_filename, soilstart4, soilcount4, nimmob);
      }
      else if(!curr_spec.layer){

        double nimmob;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          nimmob = cohort.bdall->m_soi2soi.nimmobsum;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          nimmob = cohort.bdall->y_soi2soi.nimmobsum;
        }
        std::vector<double> values(1, nimmob);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

    }//end critical(outputNIMMOB)
  }//end NIMMOB
  map_itr = netcdf_outputs.end();


  //NINPUT
  map_itr = netcdf_outputs.find("NINPUT");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNINPUT)
    {
      if(curr_spec.layer){
        /*** STUB ***/
      }
      else if(!curr_spec.layer){

        double ninput;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          ninput = cohort.bdall->m_a2soi.avlninput;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          ninput = cohort.bdall->y_a2soi.avlninput;
        }
        std::vector<double> values(1, ninput);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputNINPUT)
  }//end NINPUT
  map_itr = netcdf_outputs.end();


  //NLOST
  map_itr = netcdf_outputs.find("NLOST");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNLOST)
    {
      double nlost;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        nlost = cohort.bdall->m_soi2l.avlnlost
              + cohort.bdall->m_soi2l.orgnlost;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        nlost = cohort.bdall->y_soi2l.avlnlost
              + cohort.bdall->y_soi2l.orgnlost;
      }
      std::vector<double> values(1, nlost);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputNLOST)
  }//end NLOST
  map_itr = netcdf_outputs.end();


  //ORGN
  map_itr = netcdf_outputs.find("ORGN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputORGN)
    {
      if(curr_spec.layer){

        if(curr_spec.monthly){
          soilstart4[0] = month_timestep;
        }
        else if(curr_spec.yearly){
          soilstart4[0] = year;
        }

        std::vector<double> orgn(MAX_SOI_LAY, -99999);
        int il = 0;
        Layer* currL = this->cohort.ground.fstsoill;
        while(currL != NULL) {

          if ( (currL == this->cohort.ground.lstsoill) && ( !(currL->prevl->isSoil) || (currL->nextl->isSoil) ) ) {
            BOOST_LOG_SEV(glg, fatal) << "ERROR - Something is wrong with the soil layer pointers!!";
          }
          if (il >= MAX_SOI_LAY) {
            BOOST_LOG_SEV(glg, err) << "PROBLEM! More than " << MAX_SOI_LAY
                                    << " layers.! il[" << il << "] "
                                    << " y:" << rowidx << " x:" << colidx
                                    << " year: " << year
                                    << " currL->orgn is:" << currL->orgn;
          }

          orgn.at(il) = currL->orgn; // .at() has bounds checking and will 
                                     // throw exception when there is an 
                                     // error

          if ( (currL == this->cohort.ground.lstsoill) && (currL->prevl->isSoil) && (!(currL->nextl->isSoil)) ) {
            // Can't rely completely on NULL to break the loop. There can
            // be instances where there are rock layers (non NULL) below the 
            // soil layer. So we check that the previous layer is soil, the
            // current layer is the last soil layer, and the next layer is not
            // soil.
            break;
          }

          il++;
          currL = currL->nextl;
        }

        io_wrapper(svname, curr_filename, soilstart4, soilcount4, orgn);
      }
      // Total, instead of by layer
      else if(!curr_spec.layer){

        double orgn;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          orgn = cohort.bdall->m_soid.orgnsum;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          orgn = cohort.bdall->y_soid.orgnsum;
        }
        std::vector<double> values(1, orgn);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputORGN)
  }//end ORGN
  map_itr = netcdf_outputs.end();


  map_itr = netcdf_outputs.find("RH");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputRH)
    {
      if(curr_spec.layer){

        std::vector<double> rh(MAX_SOI_LAY, -99999);
        if(curr_spec.monthly){
          soilstart4[0] = month_timestep;
          for(int il=0; il<MAX_SOI_LAY; il++){
            rh[il] = cohort.bdall->m_soi2a.rhrawc[il]
                   + cohort.bdall->m_soi2a.rhsoma[il]
                   + cohort.bdall->m_soi2a.rhsompr[il]
                   + cohort.bdall->m_soi2a.rhsomcr[il];
          }
        }
        else if(curr_spec.yearly){
          soilstart4[0] = year;
          for(int il=0; il<MAX_SOI_LAY; il++){
            rh[il] = cohort.bdall->y_soi2a.rhrawc[il]
                   + cohort.bdall->y_soi2a.rhsoma[il]
                   + cohort.bdall->y_soi2a.rhsompr[il]
                   + cohort.bdall->y_soi2a.rhsomcr[il];
          }
        }

        io_wrapper(svname, curr_filename, soilstart4, soilcount4, rh);
      }

      else if(!curr_spec.layer){

        double rh;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          rh = cohort.bdall->m_soi2a.rhtot;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          rh = cohort.bdall->y_soi2a.rhtot;
        }
        std::vector<double> values(1, rh);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

      
    }//end critical(outputRH)
  }//end RH 
  map_itr = netcdf_outputs.end();


  //SOC
  map_itr = netcdf_outputs.find("SOC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputSOC)
    {
      if(curr_spec.layer){

        if(curr_spec.monthly){
          soilstart4[0] = month_timestep;
        }
        else if(curr_spec.yearly){
          soilstart4[0] = year;
        }

        std::vector<double> soilc(MAX_SOI_LAY, -99999);
        int il = 0;

        Layer* currL = this->cohort.ground.fstsoill;
        while(currL != NULL) {

          if ( (currL == this->cohort.ground.lstsoill) && ( !(currL->prevl->isSoil) || (currL->nextl->isSoil) ) ) {
            BOOST_LOG_SEV(glg, fatal) << "ERROR - Something is wrong with the soil layer pointers!!";
          }
          if (il >= MAX_SOI_LAY) {
            BOOST_LOG_SEV(glg, err) << "PROBLEM! More than " << MAX_SOI_LAY
                                    << " layers.! il[" << il << "] "
                                    << " y:" << rowidx << " x:" << colidx
                                    << " year: " << year
                                    << " currL->rawc is:" << currL->rawc;
          }

          soilc.at(il) = currL->rawc; // .at() has bounds checking and will 
                                      // throw exception when there is an 
                                      // error

          if ( (currL == this->cohort.ground.lstsoill) && (currL->prevl->isSoil) && (!(currL->nextl->isSoil)) ) {
            // Can't rely completely on NULL to break the loop. There can
            // be instances where there are rock layers (non NULL) below the 
            // soil layer. So we check that the previous layer is soil, the
            // current layer is the last soil layer, and the next layer is not
            // soil.
            break;
          }

          il++;
          currL = currL->nextl;
        }
        io_wrapper(svname, curr_filename, soilstart4, soilcount4, soilc);
      }
      //Total, instead of by layer
      else if(!curr_spec.layer){

        double soilc;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          soilc = cohort.bdall->m_soid.rawcsum;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          soilc = cohort.bdall->y_soid.rawcsum;
        }
        std::vector<double> values(1, soilc);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      
    }//end critical(outputSOC)
  }//end SOC
  map_itr = netcdf_outputs.end();


  /*** Six combination vars: (year,month)x(PFT,Comp,Both)***/
  //BURNVEG2AIRC
  map_itr = netcdf_outputs.find("BURNVEG2AIRC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2AIRC)
    {
      double burnveg2airc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2airc = cohort.year_fd[month].fire_v2a.orgc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2airc = cohort.fd->fire_v2a.orgc;
      }
      std::vector<double> values(1, burnveg2airc);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNVEG2AIRC)
  }//end BURNVEG2AIRC
  map_itr = netcdf_outputs.end();


  //BURNVEG2AIRN
  map_itr = netcdf_outputs.find("BURNVEG2AIRN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2AIRN)
    {
      double burnveg2airn;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2airn = cohort.year_fd[month].fire_v2a.orgn;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2airn = cohort.fd->fire_v2a.orgn;
      }

      std::vector<double> values(1, burnveg2airn);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNVEG2AIRN)
  }//end BURNVEG2AIRN
  map_itr = netcdf_outputs.end();


  //BURNVEG2DEADC
  map_itr = netcdf_outputs.find("BURNVEG2DEADC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2DEADC)
    {
      double burnveg2deadc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2deadc = cohort.year_fd[month].fire_v2dead.vegC;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2deadc = cohort.fd->fire_v2dead.vegC;
      }
      std::vector<double> values(1, burnveg2deadc);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNVEG2DEADC)
  }//end BURNVEG2DEADC
  map_itr = netcdf_outputs.end();


  //BURNVEG2DEADN
  map_itr = netcdf_outputs.find("BURNVEG2DEADN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2DEADN)
    {
      double burnveg2deadn;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2deadn = cohort.year_fd[month].fire_v2dead.strN;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2deadn = cohort.fd->fire_v2dead.strN;
      }
      std::vector<double> values(1, burnveg2deadn);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputBURNVEG2DEADN)
  }//end BURNVEG2DEADN
  map_itr = netcdf_outputs.end();


  //BURNVEG2SOIABVC
  map_itr = netcdf_outputs.find("BURNVEG2SOIABVC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2SOIABVC)
    {
      double burnveg2soiabvc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2soiabvc = cohort.year_fd[month].fire_v2soi.abvc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2soiabvc = cohort.fd->fire_v2soi.abvc;
      }
      std::vector<double> values(1, burnveg2soiabvc);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNVEG2SOIABVC)
  }//end BURNVEG2SOIABVC
  map_itr = netcdf_outputs.end();


  //BURNVEG2SOIABVN
  map_itr = netcdf_outputs.find("BURNVEG2SOIABVN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2SOIABVN)
    {
      double burnveg2soiabvn;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2soiabvn = cohort.year_fd[month].fire_v2soi.abvn;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2soiabvn = cohort.fd->fire_v2soi.abvn;
      }
      std::vector<double> values(1, burnveg2soiabvn);
      io_wrapper(svname, curr_filename, start3, count0, values);
    }//end critical(outputBURNVEG2SOIABVN)
  }//end BURNVEG2SOIABVN
  map_itr = netcdf_outputs.end();


  //BURNVEG2SOIBLWC
  map_itr = netcdf_outputs.find("BURNVEG2SOIBLWC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2SOIBLWC)
    {
      double burnveg2soiblwc;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2soiblwc = cohort.year_fd[month].fire_v2soi.blwc;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2soiblwc = cohort.fd->fire_v2soi.blwc;
      }
      std::vector<double> values(1, burnveg2soiblwc);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNVEG2SOIBLWC)
  }//end BURNVEG2SOIBLWC
  map_itr = netcdf_outputs.end();


  //BURNVEG2SOIBLWN
  map_itr = netcdf_outputs.find("BURNVEG2SOIBLWN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputBURNVEG2SOIBLWN)
    {
      double burnveg2soiblwn;
      if(curr_spec.monthly){
        start3[0] = month_timestep;
        burnveg2soiblwn = cohort.year_fd[month].fire_v2soi.blwn;
      }
      else if(curr_spec.yearly){
        start3[0] = year;
        burnveg2soiblwn = cohort.fd->fire_v2soi.blwn;
      }
      std::vector<double> values(1, burnveg2soiblwn);
      io_wrapper(svname, curr_filename, start3, count0, values);

    }//end critical(outputBURNVEG2SOIBLWN)
  }//end BURNVEG2SOIBLWN
  map_itr = netcdf_outputs.end();


  //GPP
  map_itr = netcdf_outputs.find("GPP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputGPP)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd gpp(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double gpp[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              start5[0] = month_timestep;
              gpp[ipp][ip] = cohort.bd[ip].m_a2v.gpp[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              gpp[ipp][ip] = cohort.bd[ip].y_a2v.gpp[ipp];
            }
          }
        }
        std::vector<double> gppF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            gppF[(i*NUM_PFT)+j] = gpp[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, gppF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> gpp(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            gpp[ip] = cohort.bd[ip].m_a2v.gppall; 
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            gpp[ip] = cohort.bd[ip].y_a2v.gppall; 
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, gpp);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> gpp(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){

            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              gpp[ipp] += cohort.bd[ip].m_a2v.gpp[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              gpp[ipp] += cohort.bd[ip].m_a2v.gpp[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, gpp);
      }
      //Neither PFT nor Compartment - total instead
      else if(!curr_spec.pft && !curr_spec.compartment){

        double gpp;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          gpp = cohort.bdall->m_a2v.gppall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          gpp = cohort.bdall->y_a2v.gppall;
        }
        std::vector<double> values(1, gpp);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      
    }//end critical(outputGPP)
  }//end GPP
  map_itr = netcdf_outputs.end();


  //INGPP
  map_itr = netcdf_outputs.find("INGPP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputINGPP)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd ingpp(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double ingpp[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              start5[0] = month_timestep;
              ingpp[ipp][ip] = cohort.bd[ip].m_a2v.ingpp[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              ingpp[ipp][ip] = cohort.bd[ip].y_a2v.ingpp[ipp];
            }
          }
        }
        std::vector<double> ingppF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            ingppF[(i*NUM_PFT)+j] = ingpp[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, ingppF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> ingpp(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            ingpp[ip] = cohort.bd[ip].m_a2v.ingppall; 
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            ingpp[ip] = cohort.bd[ip].y_a2v.ingppall; 
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, ingpp);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> ingpp(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){

            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              ingpp[ipp] += cohort.bd[ip].m_a2v.ingpp[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              ingpp[ipp] += cohort.bd[ip].m_a2v.ingpp[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, ingpp);
      }
      //Neither PFT nor Compartment - total instead
      else if(!curr_spec.pft && !curr_spec.compartment){

        double ingpp;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          ingpp = cohort.bdall->m_a2v.ingppall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          ingpp = cohort.bdall->y_a2v.ingppall;
        }
        std::vector<double> values(1, ingpp);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputINGPP)
  }//end INGPP
  map_itr = netcdf_outputs.end();


  //INNPP
  map_itr = netcdf_outputs.find("INNPP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputINNPP)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd innpp(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double innpp[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              start5[0] = month_timestep;
              innpp[ipp][ip] = cohort.bd[ip].m_a2v.innpp[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              innpp[ipp][ip] = cohort.bd[ip].y_a2v.innpp[ipp];
            }
          }
        }
        std::vector<double> innppF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            innppF[(i*NUM_PFT)+j] = innpp[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, innppF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> innpp(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            innpp[ip] = cohort.bd[ip].m_a2v.innppall; 
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            innpp[ip] = cohort.bd[ip].y_a2v.innppall; 
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, innpp);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> innpp(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){

            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              innpp[ipp] += cohort.bd[ip].m_a2v.innpp[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              innpp[ipp] += cohort.bd[ip].m_a2v.innpp[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, innpp);
      }
      //Neither PFT nor Compartment - total instead
      else if(!curr_spec.pft && !curr_spec.compartment){

        double innpp;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          innpp = cohort.bdall->m_a2v.innppall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          innpp = cohort.bdall->y_a2v.innppall;
        }
        std::vector<double> values(1, innpp);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputINNPP)
  }//end INNPP
  map_itr = netcdf_outputs.end();


  //LAI
  map_itr = netcdf_outputs.find("LAI");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputLAI)
    {

      //PFT
      if(curr_spec.pft){
        std::vector<double> values(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){
          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            values[ip] = cohort.cd.m_veg.lai[ip];
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            values[ip] = cohort.cd.y_veg.lai[ip];
          }
        }
        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, values);

      }
      //Total
      else if(!curr_spec.pft){

        double lai = 0;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          for(int ip=0; ip<NUM_PFT; ip++){
            if(cohort.cd.m_veg.vegcov[ip]>0.){
              lai += cohort.cd.m_veg.lai[ip];
            }
          }
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          for(int ip=0; ip<NUM_PFT; ip++){
            if(cohort.cd.y_veg.vegcov[ip]>0.){
              lai += cohort.cd.y_veg.lai[ip];
            }
          }
        }
        std::vector<double> values(1, lai);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputLAI)
  }//end LAI
  map_itr = netcdf_outputs.end();


  //LTRFALC
  map_itr = netcdf_outputs.find("LTRFALC");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputLTRFALC)
    {

      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd ltrfalc(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double ltrfalc[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              start5[0] = month_timestep;
              ltrfalc[ipp][ip] = cohort.bd[ip].m_v2soi.ltrfalc[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              ltrfalc[ipp][ip] = cohort.bd[ip].y_v2soi.ltrfalc[ipp];
            }
          }
        }
        std::vector<double> ltrfalcF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            ltrfalcF[(i*NUM_PFT)+j] = ltrfalc[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, ltrfalcF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> ltrfalc(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            ltrfalc[ip] = cohort.bd[ip].m_v2soi.ltrfalcall;
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            ltrfalc[ip] = cohort.bd[ip].y_v2soi.ltrfalcall;
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, ltrfalc);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> ltrfalc(NUM_PFT_PART, 0);
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              ltrfalc[ipp] += cohort.bd[ip].m_v2soi.ltrfalc[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              ltrfalc[ipp] += cohort.bd[ip].y_v2soi.ltrfalc[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, ltrfalc);
      }
      //Neither PFT nor compartment - totals
      else if(!curr_spec.pft && !curr_spec.compartment){

        double ltrfalc = 0;
        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            start3[0] = month_timestep;
            ltrfalc += cohort.bd[ip].m_v2soi.ltrfalcall;
          }
          else if(curr_spec.yearly){
            start3[0] = year;
            ltrfalc += cohort.bd[ip].y_v2soi.ltrfalcall;
          }
        }
        std::vector<double> values(1, ltrfalc);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputLTRFALC)
  }//end LTRFALC
  map_itr = netcdf_outputs.end();


  //LTRFALN
  map_itr = netcdf_outputs.find("LTRFALN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputLTRFALN)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd ltrfaln(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double ltrfaln[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              start5[0] = month_timestep;
              ltrfaln[ipp][ip] = cohort.bd[ip].m_v2soi.ltrfaln[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              ltrfaln[ipp][ip] = cohort.bd[ip].y_v2soi.ltrfaln[ipp];
            }
          }
        }
        std::vector<double> ltrfalnF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            ltrfalnF[(i*NUM_PFT)+j] = ltrfaln[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, ltrfalnF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> ltrfaln(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){
          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            ltrfaln[ip] = cohort.bd[ip].m_v2soi.ltrfalnall;
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            ltrfaln[ip] = cohort.bd[ip].y_v2soi.ltrfalnall;
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, ltrfaln);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> ltrfaln(NUM_PFT_PART, 0);
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){

            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              ltrfaln[ipp] += cohort.bd[ip].m_v2soi.ltrfaln[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              ltrfaln[ipp] += cohort.bd[ip].y_v2soi.ltrfaln[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, ltrfaln);
      }
      //Neither PFT nor compartment - totals
      else if(!curr_spec.pft && !curr_spec.compartment){

        double ltrfaln = 0;
        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            start3[0] = month_timestep;
            ltrfaln += cohort.bd[ip].m_v2soi.ltrfalnall;
          }
          else if(curr_spec.yearly){
            start3[0] = year;
            ltrfaln += cohort.bd[ip].y_v2soi.ltrfalnall;
          }
        }
        std::vector<double> values(1, ltrfaln);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputLTRFALN)
  }//end LTRFALN
  map_itr = netcdf_outputs.end();


  //NPP
  map_itr = netcdf_outputs.find("NPP");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNPP)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd npp(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double npp[NUM_PFT_PART][NUM_PFT];
        if(curr_spec.monthly){
          start5[0] = month_timestep;
          for(int ip=0; ip<NUM_PFT; ip++){
            for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
              npp[ipp][ip] = cohort.bd[ip].m_a2v.npp[ipp];
            }
          }
        }
        else if(curr_spec.yearly){
          start5[0] = year;
          for(int ip=0; ip<NUM_PFT; ip++){
            for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
              npp[ipp][ip] = cohort.bd[ip].y_a2v.npp[ipp];
            }
          }
        }
        std::vector<double> nppF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            nppF[(i*NUM_PFT)+j] = npp[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, nppF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> npp(NUM_PFT, -99999);
        if(curr_spec.monthly){
          PFTstart4[0] = month_timestep;
          for(int ip=0; ip<NUM_PFT; ip++){
            npp[ip] = cohort.bd[ip].m_a2v.nppall; 
          }
        }
        else if(curr_spec.yearly){
          PFTstart4[0] = year;
          for(int ip=0; ip<NUM_PFT; ip++){
            npp[ip] = cohort.bd[ip].y_a2v.nppall; 
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, npp);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> npp(NUM_PFT_PART, 0);
        if(curr_spec.monthly){
          CompStart4[0] = month_timestep;
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
            for(int ip=0; ip<NUM_PFT; ip++){
              npp[ipp] += cohort.bd[ip].m_a2v.npp[ipp];
            }
          }
        }
        else if(curr_spec.yearly){
          CompStart4[0] = year;
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
            for(int ip=0; ip<NUM_PFT; ip++){
              npp[ipp] += cohort.bd[ip].y_a2v.npp[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, npp);
      }
      //Neither PFT nor Compartment - total instead
      else if(!curr_spec.pft && !curr_spec.compartment){

        double npp;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          npp = cohort.bdall->m_a2v.nppall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          npp = cohort.bdall->y_a2v.nppall;
        }
        std::vector<double> values(1, npp);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
      
    }//end critical(outputNPP)
  }//end NPP
  map_itr = netcdf_outputs.end();


  //NUPTAKEIN
  map_itr = netcdf_outputs.find("NUPTAKEIN");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNUPTAKEIN)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){
        /*** STUB ***/
        //Currently unavailable. N uptake will need to be made accessible
        // by PFT compartment.
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){
        
        std::vector<double> innuptake(NUM_PFT, -99999);

        if(curr_spec.monthly){
          PFTstart4[0] = month_timestep;

          for(int ip=0; ip<NUM_PFT; ip++){
            innuptake[ip] = cohort.bd[ip].m_soi2v.innuptake;
          }
        }
        else if(curr_spec.yearly){
          PFTstart4[0] = year;

          for(int ip=0; ip<NUM_PFT; ip++){
            innuptake[ip] = cohort.bd[ip].y_soi2v.innuptake;
          }
        }
        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, innuptake);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){
        /*** STUB ***/
      }
      //Neither PFT nor compartment
      else if(!curr_spec.pft && !curr_spec.compartment){
        double innuptake = 0;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          innuptake = cohort.bdall->m_soi2v.innuptake;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          innuptake = cohort.bdall->y_soi2v.innuptake;
        }
        std::vector<double> values(1, innuptake);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

      
    }//end critical(outputNUPTAKEIN)
  }//end NUPTAKEIN
  map_itr = netcdf_outputs.end();


  //NUPTAKELAB
  map_itr = netcdf_outputs.find("NUPTAKELAB");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNUPTAKELAB)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){
        /*** STUB ***/
        //Currently unavailable. Labile N uptake will need to be made
        // accessible by PFT compartment.
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){
        std::vector<double> labnuptake(NUM_PFT, -99999);
        if(curr_spec.monthly){
          PFTstart4[0] = month_timestep;

          for(int ip=0; ip<NUM_PFT; ip++){
            labnuptake[ip] = cohort.bd[ip].m_soi2v.lnuptake;
          }
        }
        else if(curr_spec.yearly){
          PFTstart4[0] = year;

          for(int ip=0; ip<NUM_PFT; ip++){
            labnuptake[ip] = cohort.bd[ip].y_soi2v.lnuptake;
          }
        }
        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, labnuptake);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){
        /*** STUB ***/
      }
      //Neither PFT nor compartment
      else if(!curr_spec.pft && !curr_spec.compartment){
        double labnuptake = 0;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          labnuptake = cohort.bdall->m_soi2v.lnuptake;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          labnuptake = cohort.bdall->y_soi2v.lnuptake;
        }
        std::vector<double> values(1, labnuptake);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

      
    }//end critical(outputNUPTAKELAB)
  }//end NUPTAKELAB
  map_itr = netcdf_outputs.end();


  //NUPTAKEST
  map_itr = netcdf_outputs.find("NUPTAKEST");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputNUPTAKEST)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){
        
        //ma2dd snuptake(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double snuptake[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
            if(curr_spec.monthly){
              start5[0] = month_timestep;
              snuptake[ipp][ip] = cohort.bd[ip].m_soi2v.snuptake[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              snuptake[ipp][ip] = cohort.bd[ip].y_soi2v.snuptake[ipp];
            }
          }
        }
        std::vector<double> snuptakeF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            snuptakeF[(i*NUM_PFT)+j] = snuptake[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, snuptakeF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){
        std::vector<double> snuptake(NUM_PFT, 0);
        if(curr_spec.monthly){
          PFTstart4[0] = month_timestep;
          for(int ip=0; ip<NUM_PFT; ip++){
            snuptake[ip] = cohort.bd[ip].m_soi2v.snuptakeall;
          }
        }
        else if(curr_spec.yearly){
          PFTstart4[0] = year;
          for(int ip=0; ip<NUM_PFT; ip++){
            snuptake[ip] = cohort.bd[ip].y_soi2v.snuptakeall;
          }
        }
        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, snuptake);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){
        std::vector<double> snuptake(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){
            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              snuptake[ipp] += cohort.bd[ip].m_soi2v.snuptake[ipp]; 
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              snuptake[ipp] += cohort.bd[ip].y_soi2v.snuptake[ipp]; 
            }
          }
        }
        io_wrapper(svname, curr_filename, CompStart4, CompCount4, snuptake);
      }
      //Neither PFT nor compartment
      else if(!curr_spec.pft && !curr_spec.compartment){
        double snuptakeall;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          snuptakeall = cohort.bdall->m_soi2v.snuptakeall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          snuptakeall = cohort.bdall->y_soi2v.snuptakeall;
        }
        std::vector<double> values(1, snuptakeall);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

    }//end critical(outputNUPTAKEST)
  }//end NUPTAKEST
  map_itr = netcdf_outputs.end();


  //RG
  map_itr = netcdf_outputs.find("RG");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputRG)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd rg(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double rg[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
            if(curr_spec.monthly){
              start5[0] = month_timestep;
              rg[ipp][ip] = cohort.bd[ip].m_v2a.rg[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              rg[ipp][ip] = cohort.bd[ip].y_v2a.rg[ipp];
            }
          }
        }
        std::vector<double> rgF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            rgF[(i*NUM_PFT)+j] = rg[i][j];
          }
        }

        io_wrapper(svname, curr_filename, start5, count5, rgF); 
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> rg(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){
          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            rg[ip] = cohort.bd[ip].m_v2a.rgall;
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            rg[ip] = cohort.bd[ip].y_v2a.rgall;
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, rg); 
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> rg(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){
            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              rg[ipp] += cohort.bd[ip].m_v2a.rg[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              rg[ipp] += cohort.bd[ip].y_v2a.rg[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, rg);
   
      }
      //Neither PFT nor compartment - Total
      else if(!curr_spec.pft && !curr_spec.compartment){

        double rg;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          rg = cohort.bdall->m_v2a.rgall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          rg = cohort.bdall->y_v2a.rgall;
        }
        std::vector<double> values(1, rg);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }

    }//end critical(outputRG)
  }//end RG
  map_itr = netcdf_outputs.end();


  //RM
  map_itr = netcdf_outputs.find("RM");
  if(map_itr != netcdf_outputs.end()){
    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputRM)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd rm(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double rm[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
            if(curr_spec.monthly){
              start5[0] = month_timestep;
              rm[ipp][ip] = cohort.bd[ip].m_v2a.rm[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              rm[ipp][ip] = cohort.bd[ip].y_v2a.rm[ipp];
            }
          }
        }
        std::vector<double> rmF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            rmF[(i*NUM_PFT)+j] = rm[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, rmF); 
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> rm(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){
          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            rm[ip] = cohort.bd[ip].m_v2a.rmall;
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            rm[ip] = cohort.bd[ip].y_v2a.rmall;
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, rm); 
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> rm(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){
            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              rm[ipp] += cohort.bd[ip].m_v2a.rm[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              rm[ipp] += cohort.bd[ip].y_v2a.rm[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, rm);
   
      }
      //Neither PFT nor compartment - Total
      else if(!curr_spec.pft && !curr_spec.compartment){

        double rm;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          rm = cohort.bdall->m_v2a.rmall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          rm = cohort.bdall->y_v2a.rmall;
        }
        std::vector<double> values(1, rm);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    }//end critical(outputRM)
  }//end RM
  map_itr = netcdf_outputs.end();


  //VEGC
  map_itr = netcdf_outputs.find("VEGC");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVEGC)
    {

      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd vegc(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double vegc[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){

          for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
            if(curr_spec.monthly){
              start5[0] = month_timestep;
              vegc[ipp][ip] = cohort.bd[ip].m_vegs.c[ipp];
            }
            else if(curr_spec.yearly){
              start5[0] = year;
              vegc[ipp][ip] = cohort.bd[ip].y_vegs.c[ipp];
            }

          }
        }
        std::vector<double> vegcF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            vegcF[(i*NUM_PFT)+j] = vegc[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, vegcF);

      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> vegc(NUM_PFT, -99999);

        for(int ip=0; ip<NUM_PFT; ip++){

          if(curr_spec.monthly){
            PFTstart4[0] = month_timestep;
            vegc[ip] = cohort.bd[ip].m_vegs.call;
          }
          else if(curr_spec.yearly){
            PFTstart4[0] = year;
            vegc[ip] = cohort.bd[ip].y_vegs.call;
          }
        }

        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, vegc);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        // does this result in something like
        // "ecosystem wide leaf", "ecosuystem wide wood", "ecosystem wide root"s?
        // summed across PFTs?

        std::vector<double> vegc(NUM_PFT_PART, 0);

        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){

            if(curr_spec.monthly){
              CompStart4[0] = month_timestep;
              vegc[ipp] += cohort.bd[ip].m_vegs.c[ipp];
            }
            else if(curr_spec.yearly){
              CompStart4[0] = year;
              vegc[ipp] += cohort.bd[ip].y_vegs.c[ipp];
            }
          }
        }

        io_wrapper(svname, curr_filename, CompStart4, CompCount4, vegc);
      }
      //Neither PFT nor compartment
      else if(!curr_spec.pft && !curr_spec.compartment){

        std::vector<double> vegc(1, -99999);
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          vegc[0] = cohort.bdall->m_vegs.call;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          vegc[0] = cohort.bdall->y_vegs.call;
        }
        io_wrapper(svname, curr_filename, start3, count0, vegc);
      }
    }//end critical(outputVEGC)
  }//end VEGC
  map_itr = netcdf_outputs.end();


  //VEGN
  map_itr = netcdf_outputs.find("VEGN");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputVEGN)
    {
      //PFT and compartment
      if(curr_spec.pft && curr_spec.compartment){

        //ma2dd vegn(boost::extents[NUM_PFT_PART][NUM_PFT]);
        double vegn[NUM_PFT_PART][NUM_PFT];
        for(int ip=0; ip<NUM_PFT; ip++){
          if(cohort.cd.m_veg.vegcov[ip]>0.){//only check PFTs that exist

            for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
              if(curr_spec.monthly){
                start5[0] = month_timestep;
                vegn[ipp][ip] = cohort.bd[ip].m_vegs.strn[ipp];
              }
              else if(curr_spec.yearly){
                start5[0] = year;
                vegn[ipp][ip] = cohort.bd[ip].y_vegs.strn[ipp];
              }
            }

          }
        }
        std::vector<double> vegnF(NUM_PFT_PART*NUM_PFT, -99999);
        for(int i = 0; i < NUM_PFT_PART; i++) {
          for(int j = 0; j < NUM_PFT; j++) {
            vegnF[(i*NUM_PFT)+j] = vegn[i][j];
          }
        }
        io_wrapper(svname, curr_filename, start5, count5, vegnF);
      }
      //PFT only
      else if(curr_spec.pft && !curr_spec.compartment){

        std::vector<double> vegn(NUM_PFT, -99999);
        for(int ip=0; ip<NUM_PFT; ip++){
          if(cohort.cd.m_veg.vegcov[ip]>0.){//only check PFTs that exist

            if(curr_spec.monthly){
              PFTstart4[0] = month_timestep;
              vegn[ip] = cohort.bd[ip].m_vegs.strnall;
            }
            else if(curr_spec.yearly){
              PFTstart4[0] = year;
              vegn[ip] = cohort.bd[ip].y_vegs.strnall;
            }

          }
        }
        io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, vegn);
      }
      //Compartment only
      else if(!curr_spec.pft && curr_spec.compartment){

        std::vector<double> vegn(NUM_PFT_PART, 0);
        for(int ipp=0; ipp<NUM_PFT_PART; ipp++){
          for(int ip=0; ip<NUM_PFT; ip++){
            if(cohort.cd.m_veg.vegcov[ip]>0.){//only check PFTs that exist

              if(curr_spec.monthly){
                CompStart4[0] = month_timestep;
                vegn[ipp] += cohort.bd[ip].m_vegs.strn[ipp];
              }
              else if(curr_spec.yearly){
                CompStart4[0] = year;
                vegn[ipp] += cohort.bd[ip].y_vegs.strn[ipp];
              }

            }
          }
        }
        io_wrapper(svname, curr_filename, CompStart4, CompCount4, vegn);
      }
      //Neither PFT nor compartment
      else if(!curr_spec.pft && !curr_spec.compartment){

        double vegn;
        if(curr_spec.monthly){
          start3[0] = month_timestep;
          vegn = cohort.bdall->m_vegs.nall;
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          vegn = cohort.bdall->y_vegs.nall;
        }
        std::vector<double> values(1, vegn);
        io_wrapper(svname, curr_filename, start3, count0, values);
      }
    } //end critical(outputVEGN)
  }//end VEGN
  map_itr = netcdf_outputs.end();


  /*** Six combination vars: (year,month,day)x(PFT,total) ***/

  //EET
  map_itr = netcdf_outputs.find("EET");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputEET)
    {

      //PFT
      if(curr_spec.pft){

        if(curr_spec.daily){
          PFTstart4[0] = day_timestep;
          PFTcount4[0] = dinm;

          //ma2dd EET(boost::extents[dinm][NUM_PFT]);
          double EET[dinm][NUM_PFT];
          for(int ip=0; ip<NUM_PFT; ip++){
            for(int id=0; id<dinm; id++){
              EET[id][ip] = cohort.ed[ip].daily_eet[id];
            }
          }
          std::vector<double> EETf(NUM_PFT_PART*NUM_PFT, -99999);
          for(int i = 0; i < NUM_PFT_PART; i++) {
            for(int j = 0; j < NUM_PFT; j++) {
              EETf[(i*NUM_PFT)+j] = EET[i][j];
            }
          }
          io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, EETf);
        }
        else if(curr_spec.monthly){
          PFTstart4[0] = month_timestep;
          std::vector<double> EET(NUM_PFT, -99999);
          for(int ip=0; ip<NUM_PFT; ip++){
            EET[ip] = cohort.ed[ip].m_l2a.eet;
          }
          io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, EET);
        }
        else if(curr_spec.yearly){
          PFTstart4[0] = year;
          std::vector<double> EET(NUM_PFT, -99999);
          for(int ip=0; ip<NUM_PFT; ip++){
            EET[ip] = cohort.ed[ip].y_l2a.eet;
          }
          io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, EET);
        }
      }
      //Not PFT. Total
      else if(!curr_spec.pft){

        if(curr_spec.daily){
          start3[0] = day_timestep;
          std::vector<double> eet(31, 0);
          for(int ii=0; ii<31; ii++){
            for(int ip=0; ip<NUM_PFT; ip++){
              eet[ii] += cohort.ed[ip].daily_eet[ii];
            }
          }
          io_wrapper(svname, curr_filename, start3, count3, eet);
        }
        else if(curr_spec.monthly){
          start3[0] = month_timestep;
          std::vector<double> eet(1, cohort.edall->m_l2a.eet);
          io_wrapper(svname, curr_filename, start3, count0, eet);
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          std::vector<double> eet(1, cohort.edall->y_l2a.eet);
          io_wrapper(svname, curr_filename, start3, count0, eet);
        }
      }
    }//end critical(outputEET)
  }//end EET
  map_itr = netcdf_outputs.end();


  //PET
  map_itr = netcdf_outputs.find("PET");
  if(map_itr != netcdf_outputs.end()){

    std::string svname = map_itr->first;
    curr_spec = map_itr->second;
    curr_filename = curr_spec.file_path + curr_spec.filename_prefix + file_stage_suffix;

    #pragma omp critical(outputPET)
    {
      //PFT
      if(curr_spec.pft){

        if(curr_spec.daily){
          PFTstart4[0] = day_timestep;
          PFTcount4[0] = dinm;
          //ma2dd PET(boost::extents[dinm][NUM_PFT]);
          double PET[dinm][NUM_PFT];
          for(int ip=0; ip<NUM_PFT; ip++){
            for(int id=0; id<dinm; id++){
              PET[id][ip] = cohort.ed[ip].daily_pet[id];
            }
          }
          std::vector<double> PETf(NUM_PFT_PART*NUM_PFT, -99999);
          for(int i = 0; i < NUM_PFT_PART; i++) {
            for(int j = 0; j < NUM_PFT; j++) {
              PETf[(i*NUM_PFT)+j] = PET[i][j];
            }
          }
          io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, PETf);
        }
        else if(curr_spec.monthly){
          PFTstart4[0] = month_timestep;
          std::vector<double> PET(NUM_PFT, -99999);
          for(int ip=0; ip<NUM_PFT; ip++){
            PET[ip] = cohort.ed[ip].m_l2a.pet;
          }
          io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, PET);
        }
        else if(curr_spec.yearly){
          PFTstart4[0] = year;
          std::vector<double> PET(NUM_PFT, -99999);
          for(int ip=0; ip<NUM_PFT; ip++){
            PET[ip] = cohort.ed[ip].y_l2a.pet;
          }
          io_wrapper(svname, curr_filename, PFTstart4, PFTcount4, PET);
        }
      }
      //Not PFT. Total
      else if(!curr_spec.pft){

        if(curr_spec.daily){
          start3[0] = day_timestep;
          std::vector<double> pet(31, 0);
          for(int ii=0; ii<31; ii++){
            for(int ip=0; ip<NUM_PFT; ip++){
              pet[ii] += cohort.ed[ip].daily_pet[ii];
            }
          }
          io_wrapper(svname, curr_filename, start3, count3, pet);
        }
        else if(curr_spec.monthly){
          start3[0] = month_timestep;
          std::vector<double> pet(1, cohort.edall->m_l2a.pet);
          io_wrapper(svname, curr_filename, start3, count0, pet);
        }
        else if(curr_spec.yearly){
          start3[0] = year;
          std::vector<double> pet(1, cohort.edall->y_l2a.pet);
          io_wrapper(svname, curr_filename, start3, count0, pet);
        }
      }
    }//end critical(outputPET)
  }//end PET
  map_itr = netcdf_outputs.end();

}


