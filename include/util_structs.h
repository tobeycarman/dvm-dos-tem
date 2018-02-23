/*
 * Provides simple utility structs.
 */
#ifndef UTIL_STRUCTS_H_
#define UTIL_STRUCTS_H_

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

//template <typename T>
struct OutputDataNugget {
  std::string file_path;
  std::string vname;
  std::vector<size_t> starts;
  std::vector<size_t> counts;
  std::vector<double> data;

  OutputDataNugget(){}

  OutputDataNugget(
      const std::string file_path, 
      const std::string vname, 
      const std::vector<size_t> starts, 
      const std::vector<size_t> counts, 
      const std::vector<double> data):
    file_path(file_path),
    vname(vname),
    starts(starts),
    counts(counts),
    data(data) {}
};


struct OutputSpec {
  std::string file_path; // Subjective file path, not including filename
  std::string filename_prefix; //example: ALD_monthly
  int dim_count;

  // Which dimensions to define
  bool pft;
  bool compartment;
  bool layer;
  bool yearly;
  bool monthly;
  bool daily;
};


namespace boost {
  namespace serialization {

  template<class Archive>
  void serialize(Archive & ar, OutputDataNugget & odn, const unsigned int version) {
    ar & odn.file_path;
    ar & odn.vname;
    ar & odn.starts;
    ar & odn.counts;
    ar & odn.data;
  }

  template<class Archive>
  void serialize(Archive & ar, OutputSpec & os, const unsigned int version){

   ar & os.file_path;
   ar & os.filename_prefix;
   ar & os.dim_count;

   ar & os.pft;
   ar & os.compartment;
   ar & os.layer;
   ar & os.yearly;
   ar & os.monthly;
   ar & os.daily;
  }

  } // namespace serialization
} // namespace boost
 

#endif /* UTIL_STRUCTS_H_ */
