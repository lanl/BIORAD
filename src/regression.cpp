#include "regression.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <list>
#include <unordered_map>
#include <math.h>

#include "shuffle.h"

using namespace std;

// A way to potentially speed up the sorting and partitioning routines (which can consume a non-trivial
// amount of CPU time) is to use the g++-specific __gnu_parallel:: funcions instead
// of the std:: functions. Since __gnu_parallel:: is a drop in replacement for
// std::, include:
//#include <parallel/algorithm>
// which defines _GLIBCXX_PARALLEL_PARALLEL_H and enables parallel sorting
// Please note that approach will not work on non-g++ compilers (i.e. clang or intel).

#ifdef _GLIBCXX_PARALLEL_PARALLEL_H

	#define	SORT		__gnu_parallel::sort
	#define	PARTITION	__gnu_parallel::partition
#else

	#define	SORT		std::sort
	#define	PARTITION	std::partition
#endif // _GLIBCXX_PARALLEL_PARALLEL_H

struct sort_by_first
{
	inline bool operator()(const pair<FeatureType, unsigned int> &m_a, const pair<FeatureType, unsigned int> &m_b) const
	{
		return m_a.first < m_b.first;
	};
};
			
void RandomForest::build(const vector<WeightedValue> &m_data, 
	const vector< const vector<FeatureType>* > &m_features,
    const bool m_verbose)
{
	if(forest_size == 0){
		
		// when number of trees is less than the number if workers, an
		// individual worker may have nothing to do.
		return;
	}
	
	if( (forest_data_bag <= 0.0) || (forest_data_bag > 1.0) ){
		throw __FILE__ ":RandomForest::build: Please specify 0 < forest_data_bag <= 1.0";
	}
	
	if( (forest_feature_bag <= 0.0) || (forest_feature_bag > 1.0) ){
		throw __FILE__ ":RandomForest::build: Please specify 0 < forest_feature_bag <= 1.0";
	}

	const unsigned int num_data = m_data.size();
	
    if( num_data != m_features.size() ){
        throw __FILE__ ":RandomForest::build: Number of data points != number of feature vectors";
    }

    if(num_data == 0){
        throw __FILE__ ":RandomForest::build: No training data!";
    }

	// Make sure that all of the data have the same number of features
	// and the same number of response values
	const unsigned int num_features = m_features.front()->size();

	if(num_features == 0){
		throw __FILE__ ":RandomForest::build: No features!";
	}

	const unsigned int num_bagged_data = max(1U, (unsigned int)(forest_data_bag*num_data) );
	const unsigned int num_bagged_features = max(1U, (unsigned int)(forest_feature_bag*num_features) );

	vector<XY> data_ptr;

	// The STRATIFIED_DATA_SAMPLING code randomly selects data in a two step process:
	//	- Randomly shuffle the cluster bins
	//	- Randomly shuffle the order of protein complexes within each bin
	//	- Select data by visiting randomized bins ("round-robin") and selecting a
	//    single protein complex until num_bagged_data protein complexes have been selected
	// The non-STRATIFIED_DATA_SAMPLING code randomly selects data by:
	//	- Randomly select num_bagged_data protein complexes without replacement
	//
	// ** Update 02/15/2023 **: Stratified sampling (as implemented below) does not appear to improve
	// prediction accuracy for predicting PDBBind binding affinity for the pairwise distance features.
	// The performance is *worse* for the baseline model.
	
	//#define STRATIFIED_DATA_SAMPLING

	#ifdef STRATIFIED_DATA_SAMPLING
	unordered_map< unsigned int, deque<XY> > clustered_data;
	data_ptr.resize(num_bagged_data);
	#else
	data_ptr.resize(num_data);
	#endif // STRATIFIED_DATA_SAMPLING

	for(size_t i = 0;i < num_data;++i){
		
		if( num_features != m_features[i]->size() ){
			throw __FILE__ ":RandomForest::build: Variable number of features not allowed";
		}

        // Pointers to the features and data that can be sorted without modifying the
		// input order
		#ifdef STRATIFIED_DATA_SAMPLING
		clustered_data[m_data[i].cluster_index].push_back( XY(m_features[i], &(m_data[i])) );
		#else
		data_ptr[i].x = m_features[i];
		data_ptr[i].y = &(m_data[i]);
		#endif // STRATIFIED_DATA_SAMPLING
	}

	forest.resize(forest_size);

	vector<unsigned int> feature_index(num_features);

	for(unsigned int i = 0;i < num_features;++i){
		feature_index[i] = i;
	}

	string info;

    if(m_verbose){
	    cout << "\tBuilding trees: ";
    }

	time_t profile = time(NULL);

	for(size_t i = 0;i < forest_size;++i){

		if(m_verbose){
			
			profile = time(NULL) - profile;

			for(string::const_iterator j = info.begin();j != info.end();++j){
				cout << '\b';
			}

			for(string::const_iterator j = info.begin();j != info.end();++j){
				cout << ' ';
			}

			for(string::const_iterator j = info.begin();j != info.end();++j){
				cout << '\b';
			}

			stringstream ssout;

			ssout << (100.0*i)/forest_size << "% in " << profile << " sec";

			info = ssout.str();

			cout << info;

			cout.flush();

			profile = time(NULL);
		}


		//randomize(data_ptr.begin(), data_ptr.end(), rand_ptr);

		//vector<XY> data_ptr(num_bagged_data);

		#ifdef STRATIFIED_DATA_SAMPLING
		// Make a copy of the clustered data pointers
		vector< deque<XY> > local_data;
		
		local_data.reserve(clustered_data.size());

		// Shuffle the order of the data within each cluster bin
		for(unordered_map< unsigned int, deque<XY> >::const_iterator j = clustered_data.begin();j != clustered_data.end();++j){
			
			local_data.push_back(j->second);

			// Randomize the order of the data within each cluster
			randomize(local_data.back().begin(), local_data.back().end(), rand_ptr);
		}

		// Randomize the order of the clusters
		randomize(local_data.begin(), local_data.end(), rand_ptr);

		vector< deque<XY> >::iterator iter = local_data.begin();

		// Stratify the random data sampling by cluster
		for(unsigned int j = 0;j < num_bagged_data;++j){

			while(true){

				// "Round robin" sampling of the clusters
				if(iter == local_data.end()){
					iter = local_data.begin();
				}

				// Can't sample from empty clusters
				if(iter->empty()){
					++iter;
				}
				else{
					// This is a valid cluster to sample from
					break;
				}
			}

			// Since the data within each cluster bin has already been randomized,
			// we can select the first element
			data_ptr[j] = iter->back();

			iter->pop_back();
		}
		#else
		randomize(data_ptr.begin(), data_ptr.end(), rand_ptr);
		#endif // STRATIFIED_DATA_SAMPLING

		randomize(feature_index.begin(), feature_index.end(), rand_ptr);

		// Remove any existing tree data (for the case when we are reusing a RandomForest object
		// for multiple training runs)
		forest[i].clear();

		RandomForest::build_tree(forest[i], 
			data_ptr.begin(), data_ptr.begin() + num_bagged_data,
			feature_index.begin(), feature_index.begin() + num_bagged_features);
	}

	if(m_verbose){
		cout << endl;
	}
}

void RandomForest::build_tree(Tree &m_tree, 
	vector<XY>::iterator m_data_begin /*copy*/,
	vector<XY>::iterator m_data_end /*copy*/,
	const vector<unsigned int>::const_iterator &m_feature_begin,
	const vector<unsigned int>::const_iterator &m_feature_end)
{
	if(m_data_end <= m_data_begin){
		throw __FILE__ ":RandomForest::build_tree: No data!";
	}
	
	if(m_feature_end <= m_feature_begin){
		throw __FILE__ ":RandomForest::build_tree: No features!";
	}

	const unsigned int num_data = m_data_end - m_data_begin;
	
	const pair<float, float> ave_mse = average_and_weighted_mse(m_data_begin, m_data_end);

	TreeNode local;
	
	m_tree.push_back(local);
	
	// The weighted mean squared error of this node does depend on
	// whether the node is a split or a leaf
	m_tree.back().weighted_mse = ave_mse.second;

	// Do we have enough data to make a split?
	if(num_data <= forest_leaf){
		
		// Return the average over all leaf members
		m_tree.back().prediction = ave_mse.first;
		return;
	}
	
	// Search for the partition that obtains the smallest *weighted* mean square error
	pair<unsigned int, float> best_boundary;
	vector<unsigned int> best_left;

	if( best_split(best_boundary, best_left, forest_leaf, m_data_begin, m_data_end,
		m_feature_begin, m_feature_end) == false){
	
		// We could not find a valid split, return the average over all leaf members
		m_tree.back().prediction = ave_mse.first;
		return;
	}

	// Partition the data into left and right branches. Since a non-trivial amount of time
	// is spent paritioning the data, this code has some admittedly kludgy hacks to make it run as
	// fast as posible. A single bit (borrowed from the high bit of the index member 
	// variable) is used to indicate membership in the left hand set.	
	for(vector<unsigned int>::iterator i = best_left.begin();i != best_left.end();++i){
        (m_data_begin + *i)->is_left = true;
	}
	
	PARTITION( m_data_begin, m_data_end, IsLeft() );
	
	vector<XY>::iterator boundary_iter = m_data_begin + best_left.size();
	
	// Unset the sign bit so we can use the weight variable normally
	for(vector<XY>::iterator i = m_data_begin;i != boundary_iter;++i){
		i->is_left = false;
	}
	
	const unsigned int node = m_tree.size() - 1;
	
	// DEBUG
	//if(best_weight < 0.0){
	//	cerr << "best_weight = " << best_weight << endl;
	//	throw __FILE__ ": weight out of bounds!";
	//}

	m_tree[node].boundary = best_boundary;
	m_tree[node].left = m_tree.size();

	build_tree(m_tree, m_data_begin, boundary_iter,
		m_feature_begin, m_feature_end);

	m_tree[node].right = m_tree.size();

	build_tree(m_tree, boundary_iter, m_data_end,
		m_feature_begin, m_feature_end);
}

// Each rank computes the total predicted value for all of the trees that
// belong the rank. The ranks then share this total will each other (via
// AllReduce) and *all* ranks compute the average predicted value over all
// trees in the entire, distributed forest.
float RandomForest::predict(
	const vector<FeatureType> &m_features) const
{
	if( forest_size != forest.size() ){
		throw __FILE__ ":RandomForest::predict: forest_size != forest.size()";
	}
	
    if(forest_size == 0){
        throw __FILE__ ":RandomForest::predict: No trees in the forest!";
    }

	double sum = 0.0;
	
	// A quick bench mark suggests that parallelizing this for loop
	// is slower (by a factor of two) than the serial version.
	//#pragma omp parallel for reduction(+:sum)
	for(size_t i = 0;i < forest_size;++i){
		sum += predict_tree(forest[i], m_features);
	}
	
	return sum/forest_size;
}

float RandomForest::predict_tree( const Tree &m_tree,
	const vector<FeatureType> &m_features) const
{
	if( m_tree.empty() ){
		throw __FILE__ ":RandomForest::predict_tree: Empty tree!";
	}
	
	const unsigned int num_features = m_features.size();
	
	unsigned int index = 0;

	while(true){

		const TreeNode &node = m_tree[index];

		if( node.is_leaf() ){
			return node.prediction;
		}
		
		if(num_features <= node.boundary.first){
			throw __FILE__ ":RandomForest::predict_tree: Feature index out of bounds!";
		}

		index = (m_features[node.boundary.first] < node.boundary.second) ?
			node.left :
			node.right;
	}
	
	throw __FILE__ ":RandomRegressionForest::predict_tree: Should never get here!";
	
	return 0.0f;
}

pair<float /*ave*/, float /*mse*/> RandomForest::average_and_weighted_mse(
	const vector<XY>::const_iterator m_begin,
	const vector<XY>::const_iterator m_end)
{	
	if(m_end <= m_begin){
		throw __FILE__ ":average: No data!";
	}
	
	// Accumulate as double to avoid loss of precision for
	// large datasets
	double ave = 0.0;
	double mse = 0.0;
	double norm = 0.0;

	for(vector<XY>::const_iterator i = m_begin;i != m_end;++i){

		const WeightedValue &ref = *(i->y);
		
		ave += ref.value*ref.weight;
		mse += ref.value*ref.value*ref.weight;
		norm += ref.weight;
	}

	ave /= norm;

	// Prevent round off error from making mse < 0
	mse = fabs(mse/norm - ave*ave);

	return make_pair(ave, mse);
}

// Return true if we found a valid split, false otherwise
bool RandomForest::best_split(pair<unsigned int, FeatureType> &m_boundary, vector<unsigned int> &m_left,
	const size_t &m_leaf, 
	const vector<XY>::const_iterator &m_data_begin, 
	const vector<XY>::const_iterator &m_data_end,
	const std::vector<unsigned int>::const_iterator &m_feature_begin,
    const std::vector<unsigned int>::const_iterator &m_feature_end)
{
	bool ret = false;

	float best_score = std::numeric_limits<float>::max(); // <-- A really large value!
	
	if(m_data_end <= m_data_begin){
		throw __FILE__ ":best_split: No data!";
	}
	
	if(m_feature_end <= m_feature_begin){
		throw __FILE__ ":best_split: No features!";
	}

	const unsigned int num_data = m_data_end - m_data_begin;
	const unsigned int num_features_to_test = m_feature_end - m_feature_begin;
	
	if(m_leaf < 1){
		throw __FILE__ ":best_split: m_leaf < 1";
	}

	#pragma omp parallel
	{
		// Store the values for the i^th independent feature. We can allocate
		// this memory outside the for loop, since the size does not change
		vector< pair<FeatureType, unsigned int> > feature_slice(num_data);

		float local_best_score = std::numeric_limits<float>::max(); // <-- A really large value!
		pair<unsigned int, float> local_boundary;
		vector<unsigned int> local_left;
		
		// Test each feature and every possible boundary value within a feature 
		#pragma omp for
		for(unsigned int f = 0;f < num_features_to_test;++f){

			for(unsigned int i = 0;i < num_data;++i){
				feature_slice[i] = 
					make_pair( (*((m_data_begin + i)->x))[*(m_feature_begin + f)], i);
			}

			// Sort the feature values in ascending order. The sort_by_first() function is an optimization
			// to avoid the default comparison behavior of std::pair (which compares both elements of pair, even
			// though we only need to compare by the first element)
			sort( feature_slice.begin(), feature_slice.end(), sort_by_first() );

			// To make the calculation efficient, track the running sum of y values in the left and 
			// right branches. Use double to accumulate the floating point moments
			double sum_left = 0.0;
			double sum_right = 0.0;

			double sum_left2 = 0.0;
			double sum_right2 = 0.0;

			float num_left = 0;
			float num_right = 0;
			
			for(unsigned int i = 0;i < num_data;++i){

				const WeightedValue &ref = *((m_data_begin + i)->y);
				
				sum_right += ref.value*ref.weight;
				sum_right2 += ref.value*ref.value*ref.weight;

				// By default, all of the data starts in the *right* branch
				num_right += ref.weight;
			}
			
			//const float total_weight = num_right;

			for(unsigned int i = 0;i < num_data;++i){

				// Move data point i from the right branch to the left branch
				const WeightedValue &ref = *((m_data_begin + feature_slice[i].second)->y);

				const float y = ref.value*ref.weight;

				// Since the sums are unnormalized, we can simply remove a point from the right and
				// add it to the left
				sum_right -= y;
				sum_left += y;

				const float y2 = ref.value*ref.value*ref.weight;
				
				sum_right2 -= y2;
				sum_left2 += y2;

				num_right -= ref.weight;
				num_left += ref.weight;

				if( (i < m_leaf) || (i >= (num_data - m_leaf) ) ){
					continue;
				}

				// Don't split on equal values! We can access element i + 1 of the
				// feature slice vector since m_leaf must be greater than 0 (and the
				// above test on i will trigger a "continue" at the end of the vector range)
				if(feature_slice[i].first == feature_slice[i + 1].first){
					continue;
				}
				
				if( (num_left <= 0.0) || (num_right <= 0.0) ){

					// Don't throw an error, since round-off error in computing num_left and num_right
					// can make the values <= 0.0
					//throw __FILE__ ":best_split: Unable to normalize split variance";
					continue;
				}

				const float left_mean_square_error = (sum_left2 - sum_left*sum_left/num_left);
				const float right_mean_square_error = (sum_right2 - sum_right*sum_right/num_right);

				// Here is the original, more readable code:
				//left_mean_square_error /= num_left;
				//right_mean_square_error /= num_right;
				//
				//const float trial_mean_square_error = (num_left*left_mean_square_error + 
				//	num_right*right_mean_square_error)/num_data;

				// And here is the code when we cancel the factors of L and R
				//const float trial_mean_square_error = 
				//	(left_mean_square_error + right_mean_square_error)/total_weight;

				// Do not normalize the MSE. All of the split boundaries that are being
				// compared have the same total_weight (so the normalization does not
				// effect the choice of best split) and we want to accumulate the weighted
				// sum of MSE values for each split in a tree.
				const float trial_mean_square_error = (left_mean_square_error + right_mean_square_error);

				if(trial_mean_square_error < local_best_score){

					local_best_score = trial_mean_square_error;

					// DEBUG
					//if(local_best_weight < 0.0){

					//	cerr << "** Warning **: trial_mean_square_error = " << trial_mean_square_error << endl;
					//	cerr << "\tleft_mean_square_error = " << left_mean_square_error << endl;
					//	cerr << "\tright_mean_square_error = " << right_mean_square_error << endl;
					//	cerr << "\ttotal_weight = " << total_weight << endl;
					//}
					
					// Place the boundary at the midpoint between the current and the
					// next feature value.
					local_boundary = make_pair(*(m_feature_begin + f), 
						0.5*(feature_slice[i].first + feature_slice[i + 1].first) );

					local_left.resize(i + 1);

					for(unsigned int j = 0;j <= i;++j){
						local_left[j] = feature_slice[j].second;
					}
					
					// If a data point is not in the left branch, it must be in the
					// right branch (so we don't need to explicitly store the data points
					// that belong to the right hand branch).
				}
			}
		}
		
		#pragma omp critical
		if(local_best_score < best_score){
			
			best_score = local_best_score;
			m_boundary = local_boundary;
			m_left = local_left;
			
			ret = true;
		}
	}
	
	return ret;
}

float nn_prediction(const string &m_label, const Features &m_x, const vector<string> &m_db_labels, 
    unordered_map<string, Features> &m_db_x, vector<WeightedValue> &m_db_y, 
    const unsigned int &m_k_nearest_neighbors, const float &m_ground_truth /*= 0.0*/)
{
	const unsigned int num_data = m_db_labels.size();

	if( num_data != m_db_y.size() ){
		throw __FILE__ ":nn_prediction: |x| != |y|";
	}

	const size_t num_features = m_x.size();
	vector< pair<float /*distance*/, unsigned int /*index*/> > distance_index(num_data);

	for(unsigned int i = 0;i < num_data;++i){

		unordered_map<string, Features>::const_iterator iter = m_db_x.find(m_db_labels[i]);

		if( iter == m_db_x.end() ){
			throw __FILE__ ":nn_prediction: Unable to look up label";
		}

		float d = 0.0;
		
		for(size_t j = 0;j < num_features;++j){

			const float delta = m_x[j] - iter->second[j];

			d += delta*delta;
		}

		distance_index[i] = make_pair(d, i);
	}

	// Sort in ascending order
	sort(distance_index.begin(), distance_index.end());

	const unsigned int k = min(m_k_nearest_neighbors, num_data);

	float ret = 0.0;

	for(unsigned int i = 0;i < k;++i){
		ret += m_db_y[distance_index[i].second].value;
	}

	#ifdef DEBUG
	// DEBUG -- skip self matches (from training set residual calculations)
	if(m_label != m_db_labels[distance_index[0].second]){

		//cout << distance_index[0].first << '\t' << ret << '\t' 
		//	<< m_label << '\t' << m_db_labels[distance_index[0].second] << endl;
		
		cout << distance_index[0].first << '\t' << fabs(m_ground_truth - m_db_y[distance_index[0].second].value) << '\t' 
			<< m_label << '\t' << m_db_labels[distance_index[0].second] << endl;
	}
	#endif // DEBUG

	// Return the average over the k nearest neighbors
	return ret/k;
}