#ifndef __REGRESSION
#define __REGRESSION

#include "pdb.h"
#include "options.h"
#include <deque>
#include <vector>
#include <math.h>
#include <stdlib.h>

#include <iostream> //debug

// Forward declaration of mpi helper functions to keep the compiler happy
//template <class T> size_t mpi_size(const T &m_obj);
//template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
//template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

struct WeightedValue
{
    float value;
    float weight;
    unsigned int cluster_index;

    WeightedValue() :
        value(0.0), weight(1.0), cluster_index(0)
    {

    };

    WeightedValue(const float &m_value) :
        value(m_value), weight(1.0), cluster_index(0)
    {

    };

    WeightedValue(const float &m_value, const float &m_weight) :
        value(m_value), weight(m_weight), cluster_index(0)
    {

    };

    WeightedValue(const float &m_value, const float &m_weight, const unsigned int &m_cluster_index) :
        value(m_value), weight(m_weight), cluster_index(m_cluster_index)
    {

    };
};

typedef float FeatureType;
typedef std::vector<FeatureType> Features;

class RandomForest
{	
	public:

		struct TreeNode
		{

			// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
			// ensure that structure variable are correctly serialized
			#define TREE_NODE_MEMBERS \
        		VARIABLE(SINGLE_ARG(std::pair<unsigned int /*index*/, FeatureType /*threshold*/>), boundary) \
				VARIABLE(unsigned int, left) \
				VARIABLE(unsigned int, right) \
				VARIABLE(float, prediction) \
                VARIABLE(float, weighted_mse)

			#define VARIABLE(A, B) A B;
				TREE_NODE_MEMBERS
			#undef VARIABLE

			TreeNode() :
                left(0), right(0), prediction(0.0), weighted_mse(0.0)
			{
			};

			inline bool is_leaf() const
			{
				return left == right;
			};
		};

	private:
    
		typedef std::deque<TreeNode> Tree; 

		// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
		// ensure that structure variables are correctly serialized
		#define TREE_MEMBERS \
			VARIABLE(std::vector<Tree>, forest) \
			VARIABLE(size_t, forest_size) \
			VARIABLE(size_t, forest_leaf) \
			VARIABLE(float, forest_data_bag) \
            VARIABLE(float, forest_feature_bag) \
			VARIABLE(struct drand48_data *, rand_ptr)

		#define VARIABLE(A, B) A B;
			TREE_MEMBERS
		#undef VARIABLE

        // A lightweight structure to enable random sampling of
        // features and depended variables without modifying the
        // input structures.
        struct XY
        {
            const Features *x;
            const WeightedValue *y;
            bool is_left;

            XY() : x(NULL), y(NULL), is_left(false)
            {
            };

            XY(const Features *m_x, const WeightedValue *m_y) : x(m_x), y(m_y), is_left(false)
            {
                if(m_x == NULL){
                    throw __FILE__ ":XY::XY(): m_x == NULL";
                }

                if(m_y == NULL){
                    throw __FILE__ ":XY::XY(): m_y == NULL";
                }
            };
        };

        struct IsLeft
        {
                inline bool operator()(const XY &m_xy) const
                {
                    return (m_xy.is_left);
                };
        };

        std::pair<float, float> average_and_weighted_mse(
            const std::vector<XY>::const_iterator m_begin,
            const std::vector<XY>::const_iterator m_end);	

        bool best_split(std::pair<unsigned int, FeatureType> &m_boundary, std::vector<unsigned int> &m_left, 
            const size_t &m_leaf, 
            const std::vector<XY>::const_iterator &m_begin, 
            const std::vector<XY>::const_iterator &m_end,
            const std::vector<unsigned int>::const_iterator &m_feature_begin,
            const std::vector<unsigned int>::const_iterator &m_feature_end);

		void build_tree(Tree &m_tree, 
			std::vector<XY>::iterator m_data_begin /*copy*/,
			std::vector<XY>::iterator m_data_end /*copy*/,
            const std::vector<unsigned int>::const_iterator &m_feature_begin,
            const std::vector<unsigned int>::const_iterator &m_feature_end);
			
		float predict_tree(const Tree &m_tree, const Features &m_features) const;
		
	public:
	
		RandomForest(const size_t &m_forest_size, 
			const size_t &m_forest_leaf, 
			const float &m_forest_data_bag,
            const float &m_forest_feature_bag,
			struct drand48_data *m_rand_ptr) :
			forest_size(m_forest_size),
			forest_leaf(m_forest_leaf),
			forest_data_bag(m_forest_data_bag),
            forest_feature_bag(m_forest_feature_bag),
			rand_ptr(m_rand_ptr)
		{
		};
		
		void build(const std::vector<WeightedValue> &m_data,
			const std::vector<const Features*> &m_features, const bool m_verbose = false);
		
		float predict(const Features &m_features) const;

        float residual() const
        {
            float ret = 0.0;

            for(std::vector<Tree>::const_iterator i = forest.begin();i != forest.end();++i){

                for(Tree::const_iterator j = i->begin();j != i->end();++j){
                    if( j->is_leaf() ){
                        ret += j->weighted_mse;
                    }
                }
            }

            return ret;
        };

        // The sum of all splits in all trees
        size_t total_num_split() const
        {
            size_t ret = 0;

            for(std::vector<Tree>::const_iterator i = forest.begin();i != forest.end();++i){

                for(Tree::const_iterator j = i->begin();j != i->end();++j){
                    ret += (j->is_leaf() == false);
                }
            }

            return ret;
        };
};

struct RowColumn
{
    unsigned int row;
    unsigned int col;

    RowColumn() : row(0), col(0)
    {
    };

    RowColumn(const unsigned int &m_row, const unsigned int &m_col) :
        row(m_row), col(m_col)
    {
    };

    //inline bool operator<(const RowColumn &m_rhs) const
    //{
    //    if(row == m_rhs.row){
    //        return col < m_rhs.col;
    //    }
    //
    //    return row < m_rhs.row;
    //};
};

typedef std::vector<RowColumn> PairwiseDistanceFeature;

// In features.cpp
Features composition_features(const PDBComplex &m_pdb);

#ifdef USE_GSL
std::vector<double> linear_features(const PDBComplex &m_pdb, const std::unordered_map<std::string, unsigned int> &m_atom_table);
std::vector<double> compute_linear_regression(const std::vector<WeightedValue> &m_data, 
    const std::vector<const std::vector<double>*> &m_features, const bool m_use_bias);
double linear_prediction(const std::vector<double> &m_model, const std::vector<double> &m_features);
#endif // USE_GSL

Features extract_features(const PDBComplex &m_pdb, const std::unordered_map<std::string, unsigned int> &m_atom_table, 
    const std::string &m_id, const Options &m_opt);

void write_histogram(std::ofstream &m_fout, const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const Options &m_opt);

float nn_prediction(const std::string &m_label, const Features &m_x, const std::vector<std::string> &m_db_labels, 
    std::unordered_map<std::string, Features> &m_db_x, std::vector<WeightedValue> &m_db_y, 
    const unsigned int &m_k_nearest_neighbors, const float &m_ground_truth = 0.0);

#endif // __REGRESSION
