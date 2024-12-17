#ifndef __OPTIONS
#define __OPTIONS

#include <string>
#include <vector>
#include <unordered_map>
#include <deque>

/////////////////////////////////////////////////////////////////////////////////////////////////////
#define		BIORAD_VERSION			"0.6, September 6, 2024"

// Version 0.6
//	- Added linear regression
//	- Added k-nearest neighbor regression
// Version 0.5
//	- Added a scoring function as a surrogate for experimental affinity data
// Version 0.4
//	- Numerous improvements
//	- Atom adjacency based features
// Version 0.3
//	- Added the ability to specify an arbitry feature vector from a CSV file.
// Version 0.2
//	- Added the Department of Defense software distribution statement: 
//	  "Distribution Statement A: Approved for public release: distribution is unlimited."
// Version 0.1
//	- First C++ version

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

// Forward declaration of mpi helper functions to keep the compiler happy
template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

struct Options
{
	// Clustering methods
	enum{
		CLUSTER_BY_FILE,
		CLUSTER_BY_COMPOSITION,
		CLUSTER_BY_SEQ_ALIGN	// Not implemented yet
	};

	// Feature selection method
	enum{
		FEATURE_DEFAULT,
		FEATURE_GREEDY,
		FEATURE_NOT_SET
	};

	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define BIORAD_OPTIONS_MEMBERS \
        VARIABLE(std::string, output_file) \
		VARIABLE(std::string, affinity_file) \
		VARIABLE(std::string, cluster_file) \
		VARIABLE(std::string, rosetta_score_file) \
		VARIABLE(std::string, score_file) \
		VARIABLE(std::string, pdb_dir) \
		VARIABLE(std::string, hist_dir) \
		VARIABLE(unsigned int, cluster_by) \
		VARIABLE(unsigned int, num_fold) \
		VARIABLE(unsigned int, seed) \
		VARIABLE(unsigned int, num_trial) \
		VARIABLE(unsigned int, forest_size) \
		VARIABLE(float, forest_data_bag) \
		VARIABLE(float, forest_feature_bag) \
		VARIABLE(unsigned int, forest_leaf) \
		VARIABLE(unsigned int, k_nearest_neighbors) \
		VARIABLE(unsigned int, num_permute) \
		VARIABLE(unsigned int, feature_type) \
		VARIABLE(float, feature_hist_min) \
		VARIABLE(float, feature_hist_max) \
		VARIABLE(unsigned int, feature_hist_num_bin) \
		VARIABLE(unsigned int, feature_self_max_count) \
		VARIABLE(unsigned int, feature_hetero_max_count) \
		VARIABLE(unsigned int, affinity_min_year) \
		VARIABLE(unsigned int, affinity_max_year) \
		VARIABLE(float, affinity_min) \
		VARIABLE(float, affinity_max) \
		VARIABLE(float, feature_self_distance) \
		VARIABLE(float, feature_hetero_distance) \
		VARIABLE(bool, ignore_hetatom) \
		VARIABLE(bool, ignore_charge) \
		VARIABLE(bool, force_dimer) \
		VARIABLE(bool, random_affinity) \
		VARIABLE(bool, surrogate_affinity) \
		VARIABLE(bool, append_rosetta_scores) \
		VARIABLE(bool, append_scores) \
		VARIABLE(bool, ignore_kd) \
		VARIABLE(bool, ignore_ki) \
		VARIABLE(bool, ignore_ic50) \
		VARIABLE(bool, use_linear_regression) \
		VARIABLE(bool, use_forest_weight) \
		VARIABLE(bool, use_weighted_stats) \
		VARIABLE(bool, quit)
	
	#define VARIABLE(A, B) A B;
		BIORAD_OPTIONS_MEMBERS
	#undef VARIABLE
	
	Options() {};

	void load(int argc, char* argv[]);
};

template<> size_t mpi_size(const Options &m_opt);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_opt);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_opt);

#endif // __OPTIONS