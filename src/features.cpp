#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm> // DEBUG

#include "regression.h"
#include "symmetric_matrix.h"
#include "shuffle.h"
#include "tensor.h"

// Use the GSL for feature parameter fitting
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_multifit_nlinear.h>

using namespace std;

// Global MPI variables
extern int mpi_rank;
extern int mpi_numtasks;

#define     INVALID_ATOM_PAIR    0xFFFFFFFF
#define     IS_ROOT             (mpi_rank == 0)
//#define     DISTANCE_INDEX(AA_PROPERTY, ATOM_TYPE)  (ATOM_TYPE*AminoAcidData::NUM_AMINO_ACID_PROPERTY + AA_PROPERTY)
#define     DISTANCE_INDEX(AA_PROPERTY, ATOM_TYPE)  (ATOM_TYPE)

#define     MIN_INTERFACE_HBOND_DISTANCE    1.0
#define     MAX_INTERFACE_HBOND_DISTANCE    5.0

#define     MIN_INTERFACE_HBOND_ANGLE    0.5
#define     MAX_INTERFACE_HBOND_ANGLE    M_PI

#define     NUM_ANGLE_BINS              10
#define     NUM_DIST_BINS               10

// 3D coordinate indicies
enum{
    X, Y, Z
};

template<class T>
T operator*(const vector<T> &m_a, const vector<T> &m_b)
{
    if(m_a.size() != m_b.size()){
        throw __FILE__ ":operator*: vector dot product require equal sized vectors";
    }

    const size_t N = m_a.size();

    T ret = 0.0;

    for(size_t i = 0;i < N;++i){
        ret += m_a[i]*m_b[i];
    }

    return ret;
};

struct descending_abudance
{
    inline bool operator()(const pair<float, RowColumn> &m_a, const pair<float, RowColumn> &m_b) const
    {
        return m_a.first > m_b.first;
    };
};

Features distance_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt);
Features distance_sum_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt);
Features adjacency_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt);
Features hydrogen_bond_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt);
Features interacting_residue_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt);
Features interacting_atom_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt);

template <class T>
vector<T> unroll_symmetric_matrix(const SymmetricMatrix< vector<T> > &m_matrix)
{
    // Pack the feature vector by "un-rolling" the symmetrix matrix of vectors (i.e., tensor of 
    // pairwise distances) into a 1D vector
    vector<T> ret;

    if( m_matrix.empty() ){
        return ret;
    }

    const size_t num_elem = m_matrix.front().size();

    ret.reserve(m_matrix.size()*num_elem);

    for(typename SymmetricMatrix< vector<T> >::const_iterator i = m_matrix.begin();i != m_matrix.end();++i){

        if(i->size() != num_elem){
            throw __FILE__ ":unroll_symmetric_matrix: Variable length sub-vectors are not allowed";
        }

        for(typename vector<T>::const_iterator j = i->begin();j != i->end();++j){
            ret.push_back(*j);
        }
    }

    return ret;
}

template <class T>
vector<T> concatinate(const vector<T> &m_a, const vector<T> &m_b)
{
    vector<T> ret = m_a;

    ret.insert( ret.end(), m_b.begin(), m_b.end() );

    return ret;
}

Features extract_features(const PDBComplex &m_pdb, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt)
{
    //return hydrogen_bond_features(m_pdb, m_atom_table, m_id, m_opt);

    //return adjacency_features(m_pdb, m_atom_table, m_id, m_opt);
    return distance_features(m_pdb, m_atom_table, m_id, m_opt);
    //return interacting_atom_features(m_pdb, m_atom_table, m_id, m_opt);
    //return interacting_residue_features(m_pdb, m_atom_table, m_id, m_opt);
    //Features adj = adjacency_features(m_pdb, m_atom_table, m_id, m_opt);

    //Features dist_sum = distance_sum_features(m_pdb, m_atom_table, m_id, m_opt);

    //adj.insert(adj.end(), dist_sum.begin(), dist_sum.end()); //  Concatinate the feature vectors

    //return adj;

    //const Features dst = distance_features(m_pdb, m_atom_table, m_id, m_opt);

    //adj.insert(adj.end(), dst.begin(), dst.end()); //  Concatinate the feature vectors

    //return adj;
}

// The composition-based features are simply the raw counts of each amino acid
// type in a protein complex
Features composition_features(const PDBComplex &m_pdb)
{
    #define SINGLE_AMINO_ACID

    #ifdef SINGLE_AMINO_ACID // *Single* amino acid counts are features
    Features ret(AminoAcidData::NUM_REAL_AMINO_ACID_TYPES);

    #define AMINO_ACID_COUNT(PROTEIN) \
        for(vector<AminoAcidData>::const_iterator i = PROTEIN.begin();i != PROTEIN.end();++i){ \
            if(i->type < AminoAcidData::NUM_REAL_AMINO_ACID_TYPES){ \
                ++ret[i->type]; \
            } \
        }

    #else // Amino acid *pair* count are features

    unsigned int last_aa = AminoAcidData::NUM_AMINO_ACID_TYPES;

    #define AMINO_ACID_COUNT(PROTEIN) \
        for(vector<AminoAcidData>::const_iterator i = PROTEIN.begin();i != PROTEIN.end();++i){ \
            if(i->type < AminoAcidData::NUM_REAL_AMINO_ACID_TYPES){ \
                if(last_aa < AminoAcidData::NUM_REAL_AMINO_ACID_TYPES){ \
                    ++ret[last_aa*AminoAcidData::NUM_REAL_AMINO_ACID_TYPES + i->type]; \
                } \
                last_aa = i->type; \
            } \
        }

    Features ret(AminoAcidData::NUM_REAL_AMINO_ACID_TYPES*AminoAcidData::NUM_REAL_AMINO_ACID_TYPES);


    #endif // SINGLE_AMINO_ACID

    AMINO_ACID_COUNT(m_pdb.first);
    AMINO_ACID_COUNT(m_pdb.second);

    return ret;
}

Features distance_features(const PDBComplex &m_pdb, const unordered_map<string, unsigned int> &m_atom_table, const string &m_id, 
    const Options &m_opt)
{
    unsigned int num_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        num_atom_type = max(num_atom_type, i->second);
    }

    ++num_atom_type;

    //const size_t num_aa = AminoAcidData::NUM_AMINO_ACID_PROPERTY;
    
    // Use an explicit histogram so that we can allow for non-uniform bin sizes in the future
    vector<float> distance_bins(m_opt.feature_hist_num_bin);

    #define LINEAR_HISTOGRAM
    //#define NONLINEAR_HISTOGRAM
    //#define LOG_HISTOGRAM

    #ifdef LINEAR_HISTOGRAM
    const float delta = (m_opt.feature_hist_max - m_opt.feature_hist_min)/m_opt.feature_hist_num_bin;
    
    for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
        distance_bins[i] = m_opt.feature_hist_min + i*delta;
    }
    #endif // LINEAR_HISTOGRAM

    #ifdef NONLINEAR_HISTOGRAM

    if(m_opt.feature_hist_max <= 0.0){
        throw __FILE__ ":distance_features: m_opt.feature_hist_max <= 0.0";
    }

    if(m_opt.feature_hist_min <= 0.0){
        throw __FILE__ ":distance_features: m_opt.feature_hist_min <= 0.0";
    }

    const float delta = 1.5; // <-- what value should be used?

    if(delta <= 1.0){
        throw __FILE__ ":distance_features: delta must be > 1.0";
    }

    // The "start" and "stop" values are are choosen to make the histogram boundaries [m_opt.feature_hist_min, m_opt.feature_hist_max]
    const float start = (m_opt.feature_hist_min*powf(delta, m_opt.feature_hist_num_bin - 1.0) - m_opt.feature_hist_max)/
        (powf(delta, m_opt.feature_hist_num_bin - 1.0) - 1.0);
    const float stop = m_opt.feature_hist_max;

    for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
        // distance_bins must be sorted in ascending order
        //distance_bins[i] = expf(log(m_opt.feature_hist_min) + i*delta);
        distance_bins[ (m_opt.feature_hist_num_bin - 1) - i] = start + (stop - start)/powf(delta, i);
    }

    // DEBUG
    //for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
    //    cerr << i << '\t' << distance_bins[i] << endl;
    //}

    #endif // NONLINEAR_HISTOGRAM

    #ifdef LOG_HISTOGRAM

    if(m_opt.feature_hist_max <= 0.0){
        throw __FILE__ ":distance_features: m_opt.feature_hist_max <= 0.0";
    }

    if(m_opt.feature_hist_min <= 0.0){
        throw __FILE__ ":distance_features: m_opt.feature_hist_min <= 0.0";
    }

    // The "start" and "stop" values are are choosen to make the histogram boundaries [m_opt.feature_hist_min, m_opt.feature_hist_max]
    const float start = log(m_opt.feature_hist_min);
    const float stop = log(m_opt.feature_hist_max);
    const float delta = (stop - start)/m_opt.feature_hist_num_bin;

    for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
        // distance_bins must be sorted in ascending order
        //distance_bins[i] = expf(log(m_opt.feature_hist_min) + i*delta);
        distance_bins[i] = expf(start + i*delta);
    }

    // DEBUG
    //for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
    //    cerr << i << '\t' << distance_bins[i] << endl;
    //}

    #endif // LOG_HISTOGRAM


    //SymmetricMatrix< vector<float> > hist(num_aa*num_atom_type, vector<float>(m_opt.feature_hist_num_bin));
    SymmetricMatrix< vector<float> > hist(num_atom_type, vector<float>(m_opt.feature_hist_num_bin));

    //#define MINIMUM_DISTANCE
    #define ALL_DISTANCES

    #ifdef MINIMUM_DISTANCE

    // Find the atoms in protein 2 that are the nearest neighbor of each atom in protein 1
    for(vector<AminoAcidData>::const_iterator aa_1 = m_pdb.first.begin();aa_1 != m_pdb.first.end();++aa_1){
        for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){
            
            vector<AminoAcidData>::const_iterator min_aa = m_pdb.second.begin();
            vector<AtomData>::const_iterator min_atom = min_aa->atoms.begin();
            float min_d = atom_1->distance(*min_atom);

            for(vector<AminoAcidData>::const_iterator aa_2 = m_pdb.second.begin();aa_2 != m_pdb.second.end();++aa_2){
                for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){

                    const float d = atom_1->distance(*atom_2);
                    
                    if(d < min_d){

                        min_d = d;
                        min_aa = aa_2;
                        min_atom = atom_2;
                    }
                }
            }

            // Clamp the distance
            min_d = min(max(min_d, m_opt.feature_hist_min), m_opt.feature_hist_max);

            vector<float>::const_iterator iter = lower_bound(distance_bins.begin(), distance_bins.end(), min_d);

            size_t index = ( iter - distance_bins.begin() );

            if(index == 0){
                //cerr << "min_distance = " << min_distance << endl;
                //cerr << "min_dist2 = " << min_dist2 << endl;
                //cerr << "max_dist2 = " << max_dist2 << endl;
                //cerr << "index = " << index << endl;

                // Allow index=0 distances (corresponding to min_distance <= min_dist2)
                //throw __FILE__ ":histogram_distance_features: Index == 0";
            }
            else{
                --index;
            }

            hist(DISTANCE_INDEX(aa_1->property(), atom_1->type), 
                 DISTANCE_INDEX(min_aa->property(), min_atom->type))[index] += 1.0;
        }
    }

    // Find the atoms in protein 1 that are the nearest neighbor of each atom in protein 2
    for(vector<AminoAcidData>::const_iterator aa_2 = m_pdb.second.begin();aa_2 != m_pdb.second.end();++aa_2){
        for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){
            
            vector<AminoAcidData>::const_iterator min_aa = m_pdb.first.begin();
            vector<AtomData>::const_iterator min_atom = min_aa->atoms.begin();
            float min_d = atom_2->distance(*min_atom);

            for(vector<AminoAcidData>::const_iterator aa_1 = m_pdb.first.begin();aa_1 != m_pdb.first.end();++aa_1){
                for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){

                    const float d = atom_2->distance(*atom_1);
                    
                    if(d < min_d){

                        min_d = d2;
                        min_aa = aa_1;
                        min_atom = atom_1;
                    }
                }
            }

            // Clamp the distance
            min_d = min(max(min_d, m_opt.feature_hist_min), m_opt.feature_hist_max);

            vector<float>::const_iterator iter = lower_bound(distance_bins.begin(), distance_bins.end(), min_d);

            size_t index = ( iter - distance_bins.begin() );

            if(index == 0){
                // Allow index=0 distances (corresponding to min_distance <= min_dist2)
                //throw __FILE__ ":histogram_distance_features: Index == 0";
            }
            else{
                --index;
            }

            hist(DISTANCE_INDEX(aa_2->property(), atom_2->type), 
                 DISTANCE_INDEX(min_aa->property(), min_atom->type))[index] += 1.0;
        }
    }

    #endif // MINIMUM_DISTANCE

    #ifdef ALL_DISTANCES

    #define ACCUMULATE_DISTANCES(FIRST, SECOND) \
        for(vector<AminoAcidData>::const_iterator aa_1 = FIRST.begin();aa_1 != FIRST.end();++aa_1){ \
            for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){ \
                \
                for(vector<AminoAcidData>::const_iterator aa_2 = SECOND.begin();aa_2 != SECOND.end();++aa_2){ \
                    for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){ \
                        \
                        const float d = atom_1->distance(*atom_2); \
                        \
                        if( (d > m_opt.feature_hist_max) || (d < m_opt.feature_hist_min) ){ \
                            /* Skip the expensive lower_bound operation if this distance is out of bounds */ \
                            continue; \
                        } \
                        \
                        vector<float>::const_iterator iter = lower_bound(distance_bins.begin(), distance_bins.end(), d); \
                        \
                        size_t index = ( iter - distance_bins.begin() ); \
                        \
                        if(index == 0){ \
                            /* Skip distances that are too short */ \
                            continue; \
                        } \
                        \
                        --index; \
                        \
                        hist(DISTANCE_INDEX(aa_1->property(), atom_1->type), \
                            DISTANCE_INDEX(aa_2->property(), atom_2->type))[index] += 1.0; \
                    } \
                } \
            } \
        }
    
    ACCUMULATE_DISTANCES(m_pdb.first, m_pdb.second);
    //ACCUMULATE_DISTANCES(m_pdb.first, m_pdb.first);
    //ACCUMULATE_DISTANCES(m_pdb.second, m_pdb.second);

    #undef ACCUMULATE_DISTANCES

    #endif // ALL_DISTANCES

    return unroll_symmetric_matrix(hist);
}

void write_histogram(ofstream &m_fout, const PDBComplex &m_pdb, 
    const unordered_map<string, unsigned int> &m_atom_table, const Options &m_opt)
{
    unsigned int num_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        num_atom_type = max(num_atom_type, i->second);
    }

    ++num_atom_type;

    SymmetricMatrix< vector<float> > hist(num_atom_type, vector<float>(m_opt.feature_hist_num_bin));

    // Distance-base histrogram bins
    //const vector<float> distance_bins = {POW2(2.0), POW2(4.0), POW2(6.0), POW2(8.0)};
    //const vector<float> distance2_bins = {POW2(1.5), POW2(2.0), POW2(2.5), POW2(3.0), POW2(4.0), POW2(5.0), POW2(6.0), POW2(7.0), POW2(8.0)};
    
    // Use an explicit histogram so that we can allow for non-uniform bin sizes in the future
    vector<float> distance2_bins(m_opt.feature_hist_num_bin);
    const float delta = (m_opt.feature_hist_max - m_opt.feature_hist_min)/m_opt.feature_hist_num_bin;
    const float max_dist2 = m_opt.feature_hist_max*m_opt.feature_hist_max;
    const float min_dist2 = m_opt.feature_hist_min*m_opt.feature_hist_min;

    for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
        distance2_bins[i] = powf(m_opt.feature_hist_min + i*delta, 2.0);
    }

    for(vector<AminoAcidData>::const_iterator aa_1 = m_pdb.first.begin();aa_1 != m_pdb.first.end();++aa_1){
        for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){
            
            for(vector<AminoAcidData>::const_iterator aa_2 = m_pdb.second.begin();aa_2 != m_pdb.second.end();++aa_2){
                for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){
                    
                    const float d2 = atom_1->distance2(*atom_2);
                    
                    if( (d2 > max_dist2) || (d2 < min_dist2) ){

                        // Skip the expensive lower_bound operation if this distance is out of bounds
                        continue;
                    }

                    vector<float>::const_iterator iter = lower_bound(distance2_bins.begin(), distance2_bins.end(), d2);

                    size_t index = ( iter - distance2_bins.begin() );

                    if(index == 0){
                        // Skip distances that are too close
                        //throw __FILE__ ":histogram_distance_features: Index == 0";
                        continue;
                    }

                    --index;

                    hist(atom_1->type, atom_2->type)[index] += 1.0;

                    // Weighting the histogram count by the atomic temperature factors does not seem to
                    // help. Is the functional form for the dependence on the individual temperature
                    // factors incorrect? This likely needs additional investigation!
                    //hist(atom_1->type, atom_2->type)[index] += 1.0/sqrt(atom_1->temperature + atom_2->temperature);
                }
            }
        }
    }

    m_fout << "atom_1,atom_2";

    for(unsigned int i = 0;i < m_opt.feature_hist_num_bin;++i){
        m_fout << ",bin" << i;
    }

    m_fout << endl;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        for(unordered_map<string, unsigned int>::const_iterator j = i;j != m_atom_table.end();++j){

            const vector<float> &pair_hist = hist(i->second, j->second);

            // Only write a histogram for a given pairwise distance if it has at least one bin with a count > 0
            vector<float>::const_iterator iter = max_element(pair_hist.begin(), pair_hist.end());

            if( iter == pair_hist.end() ){
                throw __FILE__ ":write_histogram: Unexpected empty histogram";
            }

            if(*iter == 0.0){
                continue;
            }

            m_fout << i->first << ',' << j->first;

            for(vector<float>::const_iterator h = pair_hist.begin();h != pair_hist.end();++h){
                m_fout << ',' << *h;
            }

            m_fout << endl;
            
        }
    }
}

// Return the unrolled adjacency matrix, which counts the number of atoms that have a given number of neighbors
// within specified distance bounds: d_low <= di <= d_high
//                 ___________________
//  #heterodimer 0 |__|__|__|__|__|__|
//    neighbors  1 |__|__|__|__|__|__|
//               2 |__|__|__|__|__|__|
//               3 |__|__|__|__|__|__|
//               4 |__|__|__|__|__|__|
//                  1   2  3  4  5 6
//                   #self neighbors
//
// Each feature vector can be comprised of multiple neighbor count histograms, each conditioned on
//  - Distance range
//  - Source amino acid type
//  - Source atom type
Features adjacency_features(const PDBComplex &m_pdb, 
    const unordered_map<string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt)
{
    const unsigned int max_self_count = m_opt.feature_self_max_count;
    const unsigned int max_hetero_count = m_opt.feature_hetero_max_count;

    // The min and max *squared* distances
    const float min_self_dist2 = 0.0;
    const float max_self_dist2 = powf(m_opt.feature_self_distance, 2.0);
    const float min_hetero_dist2 = 0.0;
    const float max_hetero_dist2 = powf(m_opt.feature_hetero_distance, 2.0);

    unsigned int num_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        num_atom_type = max(num_atom_type, i->second);
    }

    ++num_atom_type;

    //const size_t num_aa = AminoAcidData::NUM_AMINO_ACID_PROPERTY;
    
    vector<size_t> hist_dim;
    
    //hist_dim.reserve( (num_atom_type) * (num_aa) * (max_self_count) * (max_hetero_count) );
    hist_dim.reserve( (num_atom_type) * (max_self_count) * (max_hetero_count) );

    hist_dim.push_back(num_atom_type);
    //hist_dim.push_back(num_aa);
    hist_dim.push_back(max_self_count);
    hist_dim.push_back(max_hetero_count);
    
    Tensor<FeatureType> hist(hist_dim);

    #define ACCUMULATE_HIST(SRC, NEIGHBOR) \
        for(vector<AminoAcidData>::const_iterator aa_1 = SRC.begin();aa_1 != SRC.end();++aa_1){ \
            for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){ \
                \
                Tensor<uint8_t> self_count(1); \
                \
                for(vector<AminoAcidData>::const_iterator aa_self = SRC.begin();aa_self != SRC.end();++aa_self){ \
                    for(vector<AtomData>::const_iterator atom_self = aa_self->atoms.begin();atom_self != aa_self->atoms.end();++atom_self){ \
                        \
                        if(atom_1 != atom_self){ \
                            \
                            const float d2 = atom_1->distance2(*atom_self); \
                            \
                            self_count(0) += ( (d2 <= max_self_dist2) && (d2 >= min_self_dist2) ); \
                        } \
                    } \
                } \
                \
                Tensor<uint8_t> hetero_count(1); \
                \
                for(vector<AminoAcidData>::const_iterator aa_2 = NEIGHBOR.begin();aa_2 != NEIGHBOR.end();++aa_2){ \
                    for(vector<AtomData>::const_iterator atom_hetero = aa_2->atoms.begin();atom_hetero != aa_2->atoms.end();++atom_hetero){ \
                        \
                        const float d2 = atom_1->distance2(*atom_hetero); \
                        \
                        hetero_count(0) += ( (d2 <= max_hetero_dist2) && (d2 >= min_hetero_dist2) ); \
                    } \
                } \
                \
                if(self_count.sum() == 0){ \
                    /* Skip "orphaned" atoms (most likely a result of incomplete data) */ \
                    continue; \
                } \
                \
                --self_count(0); \
                \
                /* Clamp the neighbor counts to the maximum allowed values*/ \
                self_count.clamp(max_self_count - 1); \
                hetero_count.clamp(max_hetero_count - 1); \
                vector<size_t> index; \
                index.reserve( hist_dim.size() ); \
                index.insert( index.end(), atom_1->type ); \
                /*index.insert( index.end(), aa_1->property() );*/ \
                index.insert( index.end(), self_count.begin(), self_count.end() ); \
                index.insert( index.end(), hetero_count.begin(), hetero_count.end() ); \
                ++hist(index); \
            } \
        }

    // self_count = min(self_count, max_self_count - 1);
    // hetero_count = min(hetero_count, max_hetero_count - 1);

    ACCUMULATE_HIST(m_pdb.first, m_pdb.second);
    ACCUMULATE_HIST(m_pdb.second, m_pdb.first);

    // DEBUG
    //cout << "count histogram" << endl;

    //for(unsigned int s = 0;s < max_self_count;++s){
    //    cout << "s(" << s + 1 << ")";
    
    //    for(unsigned int h = 0;h < max_hetero_count;++h){
    //        cout << '\t' << hist[INDEX(s, h)];
    //    }
    //    cout << endl;
    //}

    //#define CHECK_PROXIMITY
    #ifdef CHECK_PROXIMITY
    // Check to make sure that the two molecules are within the proximity cutoff.
    // This can catch potential problems with PDB coordinate transformations.
    size_t total_hetero_count = 0;

    for(size_t i = 0;i < num_atom_type;++i){
        for(size_t j = 0;j < num_aa;++j){
            for(size_t k = 0;k < max_self_count;++k){

                // Total number of atoms with one or more hetero-neighbors
                for(size_t m = 1;m < max_hetero_count;++m){
                    total_hetero_count += hist(i, j, k, m);
                }
            }
        }
    }

    if(total_hetero_count == 0){
        cerr << "Warning: heterodimers are outside interaction cutoff for " << m_id << endl;
    }
    #endif // CHECK_PROXIMITY

    // Normalizing the adjacency features *decreases* the prediction accuracy!
    // hist.normalize();

    // Since the tensor class is derived from vector<>, C++ will implicity cast
    // to vector<> when we return a tensor!
    return hist;
}

Features distance_sum_features(const PDBComplex &m_pdb, const unordered_map<string, unsigned int> &m_atom_table, const string &m_id, 
    const Options &m_opt)
{
    unsigned int num_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        num_atom_type = max(num_atom_type, i->second);
    }

    ++num_atom_type;
    
    // Compute the sum of dist^-n for each type of atom pair
    SymmetricMatrix<float> inv_dist(num_atom_type, 0.0f);
    SymmetricMatrix<float> inv_dist2(num_atom_type, 0.0f);

    for(vector<AminoAcidData>::const_iterator aa_1 = m_pdb.first.begin();aa_1 != m_pdb.first.end();++aa_1){
        for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){
            
            for(vector<AminoAcidData>::const_iterator aa_2 = m_pdb.second.begin();aa_2 != m_pdb.second.end();++aa_2){
                for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){

                    const float d2 = atom_1->distance2(*atom_2);
                    
                    inv_dist2(atom_1->type, atom_2->type) += 1.0/d2;
                    inv_dist(atom_1->type, atom_2->type) += 1.0/sqrt(d2);
                }
            }
        }
    }

    //return vector<float>(inv_dist.begin(), inv_dist.end());
    Features ret(inv_dist.begin(), inv_dist.end());

    ret.insert(ret.end(), inv_dist2.begin(), inv_dist2.end());

    return ret;
}

// Return the unrolled hydrogen bond paramter distribution matrix, 
// which counts the number of interface hydrogen bonds that have a given angle and distance
//               ___________________
//               |__|__|__|__|__|__|
//    bond       |__|__|__|__|__|__|
//    angle      |__|__|__|__|__|__|
//               |__|__|__|__|__|__|
//               |__|__|__|__|__|__|
//               
//                 bond distance
//
Features hydrogen_bond_features(const PDBComplex &m_pdb, 
    const unordered_map<string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt)
{   
    unsigned int num_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        num_atom_type = max(num_atom_type, i->second);
    }

    ++num_atom_type;

    //Tensor<FeatureType> hist(num_atom_type /*donor*/, num_atom_type /*acceptor*/, NUM_ANGLE_BINS, NUM_DIST_BINS);
    Tensor<FeatureType> hist(AminoAcidData::NUM_AMINO_ACID_TYPES /*donor*/, AminoAcidData::NUM_AMINO_ACID_TYPES /*acceptor*/, NUM_ANGLE_BINS, NUM_DIST_BINS);

    // Compute the hydrogen bond distances and angle
    #define     MAKE_HIST(MOL_A, MOL_B) \
        for(vector<AminoAcidData>::const_iterator i = MOL_A.begin();i != MOL_A.end();++i){ \
            for(vector<AtomData>::const_iterator j = i->atoms.begin();j != i->atoms.end();++j){ \
                \
                if(j->has_interface_hydrogen_bond() == false){ \
                    continue; \
                } \
                \
                const AtomData& hydrogen_atom = *j; \
                const AtomData& donor_atom = MOL_A[j->h_donor_index.first].atoms[j->h_donor_index.second]; \
                const AtomData& acceptor_atom = MOL_B[j->interface_h_acceptor_index.first].atoms[j->interface_h_acceptor_index.second]; \
                \
                /*Distance between donor and acceptor atoms*/ \
                const float d = donor_atom.distance(acceptor_atom); \
                \
                vector<float> DH(3); /* Vector from donor to hydrogen */ \
                \
                DH[X] = hydrogen_atom.x - donor_atom.x; \
                DH[Y] = hydrogen_atom.y - donor_atom.y; \
                DH[Z] = hydrogen_atom.z - donor_atom.z; \
                \
                vector<float> HA(3); /* Vector from hydrogen to acceptor */ \
                \
                HA[X] = acceptor_atom.x - hydrogen_atom.x; \
                HA[Y] = acceptor_atom.y - hydrogen_atom.y; \
                HA[Z] = acceptor_atom.z - hydrogen_atom.z; \
                \
                const float angle = M_PI - acos(DH*HA/sqrt( (DH*DH) * (HA*HA))); \
                \
                if( (d >= MIN_INTERFACE_HBOND_DISTANCE) && (d <= MAX_INTERFACE_HBOND_DISTANCE) && \
                    (angle >= MIN_INTERFACE_HBOND_ANGLE) && (angle <= MAX_INTERFACE_HBOND_ANGLE) ){ \
                    \
                    int angle_index = NUM_ANGLE_BINS*(angle - MIN_INTERFACE_HBOND_ANGLE)/(MAX_INTERFACE_HBOND_ANGLE - MIN_INTERFACE_HBOND_ANGLE); \
                    int d_index = NUM_DIST_BINS*(d - MIN_INTERFACE_HBOND_DISTANCE)/(MAX_INTERFACE_HBOND_DISTANCE - MIN_INTERFACE_HBOND_DISTANCE); \
                    \
                    angle_index = min(angle_index, NUM_ANGLE_BINS - 1); \
                    d_index = min(d_index, NUM_DIST_BINS - 1); \
                    \
                    /*++hist(donor_atom.type, acceptor_atom.type, angle_index, d_index); */ \
                    ++hist(MOL_A[j->h_donor_index.first].type, MOL_B[j->interface_h_acceptor_index.first].type, angle_index, d_index); \
                } \
            } \
        }

    MAKE_HIST(m_pdb.first, m_pdb.second);
    MAKE_HIST(m_pdb.second, m_pdb.first);

    return hist;
}

Features interacting_atom_features(const PDBComplex &m_pdb, 
    const unordered_map<string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt)
{
    // Stratify by amino acid type
    const size_t num_aa = AminoAcidData::NUM_AMINO_ACID_TYPES;

    SymmetricMatrix<float> count(num_aa);

    #define COUNT_INTERACTING_ATOMS(FIRST, SECOND) \
        for(vector<AminoAcidData>::const_iterator aa_1 = FIRST.begin();aa_1 != FIRST.end();++aa_1){ \
            \
            for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){ \
                \
                for(vector<AminoAcidData>::const_iterator aa_2 = SECOND.begin();aa_2 != SECOND.end();++aa_2){ \
                    for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){ \
                        \
                        const float d = atom_1->distance(*atom_2); \
                        \
                        count(aa_1->type, aa_2->type) += expf(-d/m_opt.feature_hetero_distance); \
                    } \
                } \
            } \
        }

    COUNT_INTERACTING_ATOMS(m_pdb.first, m_pdb.second);

    #undef COUNT_INTERACTING_ATOMS

    return vector<float>( count.begin(), count.end() );
}

Features interacting_residue_features(const PDBComplex &m_pdb, 
    const unordered_map<string, unsigned int> &m_atom_table, const string &m_id, const Options &m_opt)
{
    // Stratify by amino acid type
    const size_t num_aa = AminoAcidData::NUM_AMINO_ACID_TYPES;

    SymmetricMatrix<float> count(num_aa);

    #define COUNT_INTERACTING_AA(FIRST, SECOND) \
        for(vector<AminoAcidData>::const_iterator aa_1 = FIRST.begin();aa_1 != FIRST.end();++aa_1){ \
            \
            for(vector<AminoAcidData>::const_iterator aa_2 = SECOND.begin();aa_2 != SECOND.end();++aa_2){ \
                \
                bool interacting = false; \
                \
                for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){ \
                    \
                    for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){ \
                        \
                        const float d = atom_1->distance(*atom_2); \
                        \
                        if(d <= m_opt.feature_hetero_distance){ \
                            interacting = true; \
                        } \
                    } \
                } \
                \
                if(interacting){ \
                    count(aa_1->type, aa_2->type) += 1.0; \
                } \
            } \
        }

    COUNT_INTERACTING_AA(m_pdb.first, m_pdb.second);

    #undef COUNT_INTERACTING_AA
    
    return vector<float>( count.begin(), count.end() );
}

vector<double> linear_features(const PDBComplex &m_pdb, const unordered_map<string, unsigned int> &m_atom_table)
{
    unsigned int num_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        num_atom_type = max(num_atom_type, i->second);
    }

    ++num_atom_type;

    // For each pair of atom types (e.g., H-O, C-C, P-S, H-H, ...), compute the sum of pairwise distances raised
    // to different powers
    //const vector<float> powers = {-12.0, -6.0, -2.0, -1.0};
    const vector<float> powers = {-2.0, -1.0, 0.0, 1.0};
    
    const size_t num_powers = powers.size();

    if(num_powers == 0){
        throw __FILE__ ":linear_features: Please use a least one power of the pairwise distance";
    }

    SymmetricMatrix< vector<double> > between( num_atom_type, vector<double>(num_powers) );
    SymmetricMatrix< vector<double> > within( num_atom_type, vector<double>(num_powers) );

    #define ACCUMULATE_DISTANCES(FIRST, SECOND, DST) \
        for(vector<AminoAcidData>::const_iterator aa_1 = FIRST.begin();aa_1 != FIRST.end();++aa_1){ \
            for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){ \
                \
                for(vector<AminoAcidData>::const_iterator aa_2 = SECOND.begin();aa_2 != SECOND.end();++aa_2){ \
                    for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){ \
                        \
                        const float d = atom_1->distance(*atom_2); \
                        \
                        if(d <= 0.0){ \
                            continue; /*don't compare the same atoms from the same molecule*/ \
                        } \
                        \
                        for(size_t i = 0;i < num_powers;++i){ \
                            DST(atom_1->type, atom_2->type)[i] += pow(d, powers[i]); \
                        } \
                    } \
                } \
            } \
        }
    
    ACCUMULATE_DISTANCES(m_pdb.first, m_pdb.second, between);
    //ACCUMULATE_DISTANCES(m_pdb.first, m_pdb.first, within);
    //ACCUMULATE_DISTANCES(m_pdb.second, m_pdb.second, within);

    #undef ACCUMULATE_DISTANCES

    //return concatinate( unroll_symmetric_matrix(between), unroll_symmetric_matrix(within) );

    return unroll_symmetric_matrix(between);
}