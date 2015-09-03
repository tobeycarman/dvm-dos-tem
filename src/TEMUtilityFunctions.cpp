//
//  TEMUtilityFunctions.cpp
//  dvm-dos-tem
//
//  Created by Tobey Carman on 4/10/14.
//  Copyright (c) 2014 Spatial Ecology Lab. All rights reserved.
//

#include <string>
#include <stdexcept>
#include <fstream>
#include <cerrno>
#include <sstream>

#include <netcdfcpp.h>

#include <json/reader.h>
#include <json/value.h>

#include "inc/physicalconst.h" // for PI
#include "inc/timeconst.h" // for mapping from first day of month -> day of year

#include "TEMLogger.h"
#include "TEMUtilityFunctions.h"

extern src::severity_logger< severity_level > glg;

namespace temutil {
  
  /** Sets a negative number to 0.0.
  Can be used with std::for_each to make sure the contents of a
  container does not contain positive numbers.
  */
  void force_negative2zero(float& i) {
    if (i < 0) {
      i = 0.0;
    } else {
      // do nothing...
    }
  }


  /** Maybe useful for preventing divide by zero errors?
      - probably not very effecient for production runs, but maybe helpful for 
        debugging
      - sign parameter lets you force the result to a given sign
  */
  double NON_ZERO(const double val, const int sign) {
    assert((sign == 1 || sign == -1) && "Invalid parameter for sign. Must be 1, or -1");

    if (val != 0) {
      return val; // nothing to do, value is already non-zero.
    }

    // otherwise, return some arbitrary very small number
    double vsm = 0.0000000000000001 * sign;
    BOOST_LOG_SEV(glg, debug) << "NON_ZERO: setting to a very small number: " << vsm;
    return vsm;

  }

  /** Return the day of year based on month and day, everything is zero based. */
  int day_of_year(int month, int day) {
    return DOYINDFST[month] + day;
  }

  /** Length of day as a function of latitude (degrees) and day of year.
  */
  float length_of_day(float lat, int doy) {
    // OLD COMMENTS - should update.
    // the following are the original algorithm, and
    //  modified as below by Yi (2013 Feb):
    //  double ampl;
    //  ampl = exp(7.42 +0.045 *gd.lat)/3600.;
    //  gd.alldaylengths[id] = ampl * (sin ((id -79) *0.01721)) +12.0;
    // make sure all arguments in sin, cos and tan are
    //  in unit of arc (not degree)
    //http://www.jgiesen.de/astro/solarday.htm
    //http://www.gandraxa.com/length_of_day.xml

    double m = 1 - tan(lat*PI/180.0) * tan(23.45 * cos(doy*PI/182.625) * PI/180.0);
    m = fmax(m, 0.0);
    m = fmin(m, 2.0);
    double b = acos(1-m)/PI;
    double daylength = b * 24;

    return daylength;
  }

  /** Takes an integer number and returns a string like "CMT01".
  * Inserts leading zeros if needed. Works if 0 <= cmtnumber <= 99.
  */
  std::string cmtnum2str(int cmtnumber) {

    // get string representation of number
    std::stringstream cmtnumber_ss;
    cmtnumber_ss << cmtnumber;

    // take care of leading zero...
    std::string prefix = "";
    if (cmtnumber < 10) {
      prefix =  "CMT0";
    } else {
      prefix = "CMT";
    }

    return prefix + cmtnumber_ss.str();
  }


  /** Returns true for 'on' and false for 'off'.
   * Throws exception if s is not "on" or "off".
   * might want to inherit from std exception or do something else?
   */
  bool onoffstr2bool(const std::string &s) {
    if (s.compare("on") == 0) {
      return true;
    } else if (s.compare("off") == 0) {
      return false;
    } else {
      throw std::runtime_error("Invalid string! Must be 'on' or 'off'.");
    }
  }


  /** Read a file into a string. Return the string. 
  * Throws exceptions if there is an error reading the file. Poached from:
  * http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
  */
  std::string file2string(const char *filename) {
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in) {
      std::string contents;
      in.seekg(0, std::ios::end);
      contents.resize(in.tellg());
      in.seekg(0, std::ios::beg);
      in.read(&contents[0], contents.size());
      in.close();
      return(contents);
    }
    throw(errno);
  }

  /** Reads a json foramtted "control file", returning a json data object.
  */
  Json::Value parse_control_file(const std::string &filepath) {

    BOOST_LOG_SEV(glg, debug) << "Read the control file ('" << filepath << "') into a string...";
    std::string datastring = file2string(filepath.c_str());

    BOOST_LOG_SEV(glg, debug) << "Creating Json Value and Reader objects...";
    Json::Value root;   // will contain the root value after parsing
    Json::Reader reader;

    BOOST_LOG_SEV(glg, debug) << "Trying to parse the json data string...";

    bool parsingSuccessful = reader.parse( datastring, root );

    BOOST_LOG_SEV(glg, debug) << "Parsing successful?: " << parsingSuccessful;

    if ( !parsingSuccessful ) {
        BOOST_LOG_SEV(glg, fatal) << "Failed to parse configuration file: " << filepath;
        BOOST_LOG_SEV(glg, fatal) << reader.getFormatedErrorMessages();
        exit(-1);
    }

    return root;
  }



  /** Opens a netcdf file and collects a bunch of info into a string for printing.
  */
  std::string report_on_netcdf_file(const std::string& fname, const std::string& varname) {

    std::stringstream ss;

    int ncid;
    temutil::nc( nc_open(fname.c_str(), NC_NOWRITE, &ncid) );

    // lookup variable by name
    int vid;
    temutil::nc( nc_inq_varid(ncid, varname.c_str(), &vid) );

    // stuff to report
    nc_type v_type;                      /* variable type */
    int v_ndims;                         /* number of dims */
    int v_dimids[NC_MAX_VAR_DIMS];       /* dimension IDs */
    int v_natts;                         /* number of attributes */

    temutil::nc( nc_inq_var (ncid, vid, 0 /*NC_MAX_NAME*/, &v_type,
                             &v_ndims, v_dimids, &v_natts ) );

    temutil::nc( nc_close(ncid) );

    ss << "varname: " << varname.c_str()
       << " id: " << vid
       << " type: " << v_type      // crude - prints number
       << " ndims: " << v_ndims
       << " dimids: " << v_dimids  // crude - prints address
       << " natts: " << v_natts;

    return ss.str();

  }

  /** Opens a netcdf file, assumed to be setup for dvmdostem, and reports y,x
  * dimension lengths.
  */
  std::string report_yx_pixel_dims2str(const std::string& fname) {

    std::stringstream ss;

    int ncid;
    temutil::nc( nc_open(fname.c_str(), NC_NOWRITE, &ncid) );

    int xD, yD;
    size_t xD_len, yD_len;

    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );
    temutil::nc( nc_inq_dimlen(ncid, xD, &xD_len) );

    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    temutil::nc( nc_inq_dimlen(ncid, yD, &yD_len) );

    temutil::nc( nc_close(ncid) );

    ss << "Y len: " << yD_len << " X len: " << xD_len << " (" << fname << ")";

    return ss.str();
  }


  /** Opens a netcdf file for reading, returns NcFile object.
  * 
  * NetCDF library error mode is set to silent (no printing to std::out), and
  * non-fatal. If the file open fails, it logs a message and exits the program
  * with a non-zero exit code.
  */
  NcFile open_ncfile(std::string filename) {
    
    BOOST_LOG_SEV(glg, info) << "Opening NetCDF file: " << filename;

    NcError err(NcError::silent_nonfatal);
    NcFile file(filename.c_str(), NcFile::ReadOnly);
    
    if( !file.is_valid() ) {
      BOOST_LOG_SEV(glg, fatal) << "Problem opening/reading " << filename;
      exit(-1);
    }

    return file;

  }

  /** Given an NcFile object and dimension name, reutrns a pointer to the NcDim.
  * 
  * If the dimension-read is not valid, then an error message is logged and 
  * the program exits with a non-zero status.
  */
  NcDim* get_ncdim(const NcFile& file, std::string dimname) {
  
    BOOST_LOG_SEV(glg, debug) << "Looking for dimension '" << dimname << "' in NetCDF file...";
    NcDim* dim = file.get_dim(dimname.c_str());
    
    BOOST_LOG_SEV(glg, debug) << "'" << dimname <<"' is valid?: " << dim->is_valid();
    if ( !dim->is_valid() ) {
      BOOST_LOG_SEV(glg, fatal) << "Problem with '" << dimname << "' in NetCDF file.";
      exit(-1);
    }

    return dim;

  }

  /** Given an NcFile object and a variable name, returns a pointer to the NcVar.
  *
  * If getting the variable fails, then an error message is logged, and the
  * the program exits with a non-zero status.
  */
  NcVar* get_ncvar(const NcFile& file, std::string varname) {
    BOOST_LOG_SEV(glg, debug) << "Looking for variable '" << varname << "' in NetCDF file...";
    NcVar* var = file.get_var(varname.c_str());
    if (var == NULL) {
      BOOST_LOG_SEV(glg, fatal) << "Problem with '" << varname << "' variable in NetCDF file!";
      exit(-1);
    }
    return var;
  }


  /** Look up a lat-lon pair in a NetCDF file given a rec_id.
  *
  * Note: rec_id - the order (from ZERO) in the .nc file,
  *       gridid - the grid id user-defined in the dataset
  */
  std::pair<float, float> get_location(std::string gridfilename, int rec_id) {

    float lat;
    float lon;

    NcFile grid_file = temutil::open_ncfile(gridfilename);

    NcVar* latV = temutil::get_ncvar(grid_file, "LAT");
    latV->set_cur(rec_id);
    latV->get(&lat, 1);

    NcVar* lonV = temutil::get_ncvar(grid_file, "LON");
    lonV->set_cur(rec_id);
    lonV->get(&lon, 1);

    return std::pair<float, float> (lat, lon);

  }
  
  /** Handles NetCDF errors by printing message and exiting. */
  void handle_error(int status) {
    if (status != NC_NOERR) {
      fprintf(stderr, "%s\n", nc_strerror(status));
      BOOST_LOG_SEV(glg, fatal) << nc_strerror(status);
      exit(-1);
    }
  }
  
  /** Two letter alias for handle_error() */
  void nc(int status) {
    handle_error(status);
  }

  /** rough draft for reading a timeseries for a single location from a
  *   new-style input file
  */
  template <typename DTYPE>
  std::vector<DTYPE> get_timeseries(const std::string &filename,
                                    const std::string &var,
                                    const int y, const int x) {

    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;
    BOOST_LOG_SEV(glg, debug) << "Getting variable: " << var;

    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

    int timeD;
    size_t timeD_len;

    temutil::nc( nc_inq_dimid(ncid, "time", &timeD) );
    temutil::nc( nc_inq_dimlen(ncid, timeD, &timeD_len) );

    int timeseries_var;
    temutil::nc( nc_inq_varid(ncid, var.c_str(), &timeseries_var) );

    BOOST_LOG_SEV(glg, note) << "Getting value for pixel(y,x): ("<< y <<","<< x <<").";
    int yD, xD;
    size_t yD_len, xD_len;

    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    temutil::nc( nc_inq_dimlen(ncid, yD, &yD_len) );

    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );
    temutil::nc( nc_inq_dimlen(ncid, xD, &xD_len) );

    // specify start indices for each dimension (y, x)
    size_t start[3];
    start[0] = 0;     // from begining of time
    start[1] = y;     // specified location
    start[2] = x;     // specified location

    // specify counts for each dimension
    size_t count[3];
    count[0] = timeD_len;     // all time
    count[1] = 1;             // one location
    count[2] = 1;             // one location

    // might need to add a call to nc_inq_var so we can find the type and call
    // the right nc_get_var_type(...) function..
    char vname[NC_MAX_NAME+1];
    nc_type the_type;
    int num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];
    int num_atts;
    temutil::nc( nc_inq_var(ncid, timeseries_var, vname, &the_type, &num_dims, dim_ids, &num_atts) );

    std::vector<DTYPE> data2;

    BOOST_LOG_SEV(glg, debug) << "Grab the data from the netCDF file...";
    if (the_type == NC_INT64) {
      int dataI[timeD_len];
      temutil::nc( nc_get_vara_int(ncid, timeseries_var, start, count, &dataI[0]) );
      unsigned dataArraySize = sizeof(dataI) / sizeof(DTYPE);
      data2.insert(data2.end(), &dataI[0], &dataI[dataArraySize]);
    }
    if (the_type == NC_FLOAT) {
      float dataF[timeD_len];
      temutil::nc( nc_get_vara_float(ncid, timeseries_var, start, count, &dataF[0]) );
      unsigned dataArraySize = sizeof(dataF) / sizeof(DTYPE);
      data2.insert(data2.end(), &dataF[0], &dataF[dataArraySize]);

    } else {
      BOOST_LOG_SEV(glg, err) << "Unknown datatype: '" << the_type << "'. Returning empty vector.";
    }
    return data2;
  }

  /** rough draft for reading a timeseries of co2 data from a new-style co2 file.
  */
  std::vector<float> get_timeseries(const std::string &filename,
                                    const std::string& var) {

    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;
    BOOST_LOG_SEV(glg, debug) << "Getting variable: " << var;

    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

    int timeseries_var;
    temutil::nc( nc_inq_varid(ncid, var.c_str(), &timeseries_var) );

    int yearD;
    size_t yearD_len;

    temutil::nc( nc_inq_dimid(ncid, "year", &yearD) );
    temutil::nc( nc_inq_dimlen(ncid, yearD, &yearD_len) );

    BOOST_LOG_SEV(glg, debug) << "Allocate a vector with enough space for the whole timeseries (" << yearD_len << " timesteps)";
    std::vector<float> data(yearD_len);

    size_t start[1];
    start[0] = 0;         // from beginning of time

    size_t count[1];
    count[0] = yearD_len; // all time

    BOOST_LOG_SEV(glg, debug) << "Grab the data from the netCDF file...";
    temutil::nc( nc_get_vara_float(ncid, timeseries_var, start, count, &data[0]) );

    return data;
  }

  /** rough draft - look up lon/lat in nc file from y,x coordinates. 
      Assumes that the file has some coordinate dimensions...
  */
  std::pair<float, float> get_latlon(const std::string& filename, int y, int x) {
    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;
    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid ) );

    int yD, xD;
    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );

    int latV;
    int lonV;
    temutil::nc( nc_inq_varid(ncid, "lat", &latV) );
    temutil::nc( nc_inq_varid(ncid, "lon", &lonV) );

    size_t start[2];
    start[0] = y;
    start[1] = x;

    float lon_value;
    float lat_value;
    temutil::nc( nc_get_var1_float(ncid, latV, start, &lat_value));
    temutil::nc( nc_get_var1_float(ncid, lonV, start, &lon_value));

    return std::pair<float, float>(lat_value, lon_value);
  }

  /** rough draft for reading a fri for a single location */
  int get_fri(const std::string &filename, int y, int x) {
    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;
    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid ) );

    int yD, xD;
    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );

    int friV;
    temutil::nc( nc_inq_varid(ncid, "fri", &friV));

    size_t start[2];
    start[0] = y;
    start[1] = x;

    int fri_value;
    temutil::nc( nc_get_var1_int(ncid, friV, start, &fri_value)  );

    return fri_value;

  }

  /** rough draft for reading a single location's, list of 'fire years' 
    (explicit replacement for FRI) */
  std::vector<int> get_fire_years(const std::string &filename, int y, int x) {
  // FIX: implement this!
/*
    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;
    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid ) );

    int yD, xD;
    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );
  
  
    int fyV;
    temutil::nc( nc_inq_varid, "fire_years", fyV);

       int nc_get_vara       (int ncid, int varid, const size_t start[],
                            const size_t count[], void *ip);
    size_t start[2];
    start[0] = y;
    start[1] = x;
    nc_get_vara(ncid, fyV, start, &)
    
    nc_vlen_t fyrs[?];
    nc_vlen_t fszs[?];

    
//ncid
//The ncid of the file that contains the VLEN type.
//xtype
//The type of the VLEN to inquire about. 
//name
//A pointer for storage for the types name. The name will be NC_MAX_NAME characters or less. 
//datum_sizep
//A pointer to a size_t, this will get the size of one element of this vlen. 
//base_nc_typep
//A pointer to an nc_type, this will get the type of the VLEN base type. (In other words, what type is this a VLEN of?)

    temutil::nc( nc_inq_vlen(ncid, NC_INT, "fire_years", ??,??))
  
  
    if (nc_inq_vlen(ncid, typeid, name_in, &size_in, &base_nc_type_in)) ERR;
    
    std::vector<int> fy = ??
*/
  }

  /** rough draft for reading a single location's list of 'fire sizes'. 
      This should parallel fire_years. In otherwords for every year mentioned
      in fire-years, there must be a corresponding fire size?
      Does this imply that the fire_years must be sorted?
  */
  std::vector<int> get_fire_sizes(const std::string &filename, int y, int x){
    // FIX: implement this!
  }

  /** rough draft for reading a single location, veg classification
  */
  int get_veg_class(const std::string &filename, int y, int x) {

    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;
    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid ) );

    int xD, yD;

    //size_t yD_len, xD_len;

    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    //temutil::nc( nc_inq_dimlen(ncid, yD, &yD_len) );

    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );
    //temutil::nc( nc_inq_dimlen(ncid, xD, &xD_len) );

    int veg_classificationV;
    temutil::nc( nc_inq_varid(ncid, "veg_class", &veg_classificationV) );

    size_t start[2];
    start[0] = y;
    start[1] = x;

    int veg_class_value;
    temutil::nc( nc_get_var1_int(ncid, veg_classificationV, start, &veg_class_value)  );

    return veg_class_value;
  }

  /** Looks up a pixel's drainage class from a netcdf file */
  int get_drainage_class(const std::string& filename, int y, int x) {

    BOOST_LOG_SEV(glg, debug) << "Opening dataset: " << filename;

    int ncid;
    temutil::nc( nc_open(filename.c_str(), NC_NOWRITE, &ncid ) );

    int xD, yD;

    temutil::nc( nc_inq_dimid(ncid, "Y", &yD) );
    temutil::nc( nc_inq_dimid(ncid, "X", &xD) );

    int drainage_classV;
    temutil::nc( nc_inq_varid(ncid, "drainage_class", &drainage_classV) );

    size_t start[2];
    start[0] = y;
    start[1] = x;

    int drainage_class_value;
    temutil::nc( nc_get_var1_int(ncid, drainage_classV, start, &drainage_class_value)  );

    return drainage_class_value;
  }


  /** Parses a string, looking for a community code.
   Reads the string, finds the first occurrence of the characters "CMT", and
   returns a string consisting of CMT and the following two characters.

   Returns something like "CMT01".
  */
  std::string read_cmt_code(std::string s) {
    int pos = s.find("CMT");
    return s.substr(pos, 5);
  }

  /** Parses a string, looking for a community code, returns an integer.
  */
  int cmtcode2num(std::string s) {
    int pos = s.find("CMT");
    
    return std::atoi( s.substr(pos+3, 2).c_str() );
  }

  /** Reads a file, returning a contiguous section of lines surrounded by "CMT".
  * Each line from the file is an element in the vector. 
  */  
  std::vector<std::string> get_cmt_data_block(std::string filename, int cmtnum) {

    BOOST_LOG_SEV(glg, note) << "Opening file: " << filename;
    std::ifstream par_file(filename.c_str(), std::ifstream::in);

    if ( !par_file.is_open() ) {
      BOOST_LOG_SEV(glg, fatal) << "Problem opening " << filename << "!";
      exit(-1);
    }

    std::string cmtstr = cmtnum2str(cmtnum);

    // create a place to store lines making up the community data "block"
    std::vector<std::string> cmt_block_vector; 
    BOOST_LOG_SEV(glg, note) << "Searching file for community: " << cmtstr;
    for (std::string line; std::getline(par_file, line); ) {
      int pos = line.find(cmtstr);
      if ( pos != std::string::npos ) {

        // add the 'header line' to the data block
        cmt_block_vector.push_back(line);			

      for (std::string block_line; std::getline(par_file, block_line); ) {

        int block_line_pos = block_line.find("CMT");
        if ( block_line_pos != std::string::npos ) {
          //std::cout << "Whoops - line contains 'CMT'. Must be first line of next community data block; breaking loop.\n";
          break;
        } else {
          //std::cout << "Add line to cmt_block_vector: " << block_line << std::endl;
          cmt_block_vector.push_back(block_line);
        }
      }

      }
    }
    return cmt_block_vector;
  }

  /** Takes a cmt data "block" and strips any comments. */
  std::list<std::string> strip_comments(std::vector<std::string> idb) {
    
    std::list<std::string> l;

    for (std::vector<std::string>::iterator it = idb.begin(); it != idb.end(); ++it ) {

      std::string line = *it;

      // // strip comment and everthing after
      // size_t pos = line.find("//");
      // line = line.substr(0, pos);

      // Split into data and comment (everything after '//')
      size_t pos = line.find("//");
      std::string data = line.substr(0, pos);
      std::string comment = "";
      if (pos != std::string::npos) {
        comment = line.substr(pos+2, std::string::npos);
      }
      //std::cout << "Data: " << data << " Comment: " << comment << std::endl;

      if (data.size() == 0) {
        // pass
      } else {
        l.push_back(line);
      }

    }

    return l;
  }


  /** Given a file name, a community number and a number for expected lines of
   * data, returns a list of strings with just that data, after stripping 
   * comments.
   */
  std::list<std::string> parse_parameter_file(
      const std::string& fname, int cmtnumber, int linesofdata) {
    
    BOOST_LOG_SEV(glg, note) << "Parsing '"<< fname << "', "
                             << "for community number: " << cmtnumber;

    // get a vector of strings for that cmt "block". includes comments.
    std::vector<std::string> v(get_cmt_data_block(fname, cmtnumber));
    
    // strip the comments and turn it into a list
    std::list<std::string> datalist(strip_comments(v));

    // handy for debugging...
    //std::list<std::string>::iterator it = datalist.begin();
    //for (it; it != datalist.end(); ++it) {
    //  std::cout << "list item: " << *it << std::endl;
    //}

    // check the size
    if (datalist.size() != linesofdata) {
      BOOST_LOG_SEV(glg, err) << "Expected " << linesofdata << ". "
                              << "Found " << datalist.size() << ". "
                              << "(" << fname << ", community " << cmtnumber << ")";
      exit(-1);
    }
    
    return datalist;
  }


  // Explicit instatiation of different template types...??
  // Not sure if this is necessary ??
  template void pfll2data(std::list<std::string> &l, double &data);
  template void pfll2data(std::list<std::string> &l, float &data);
  
  template void pfll2data_pft(std::list<std::string> &l, double *data);
  template void pfll2data_pft(std::list<std::string> &l, float *data);

  // inorder to keep the template function definition out of the header
  // we have to explicitly instantiate it here...
  template std::vector<int> get_timeseries<int>(const std::string &filename,
      const std::string &var, const int y, const int x);
  template std::vector<float> get_timeseries<float>(const std::string &filename,
      const std::string &var, const int y, const int x);

}
