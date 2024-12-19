#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>

#include "pdb.h"
#include "options.h"
#include "mpi_util.h"
#include "parse_util.h"
#include "cluster.h"
#include "shuffle.h"
#include "regression.h"
#include "rosetta.h"
#include "affinity.h"
#include "symmetric_matrix.h"

// For creating histogram directories if needed
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

// Global variables for MPI
int mpi_numtasks = 1; // Default for non-MPI operation
int mpi_rank = 0; // Default for non-MPI operation
double start_time;

#define     IS_ROOT     (mpi_rank == 0)

struct find_by_accession
{
    inline bool operator()(const pair<string, Affinity> &m_a, const string &m_b) const
    {
        return m_a.first < m_b;
    };
};

typedef vector< pair<float /*observed*/, float /*predicted*/> > ClusterPredictions;

void terminate_program(int m_sig);

vector< pair<string, Affinity> > parse_affinity(const string &m_affinity_file);
vector< pair<string, Affinity> >::const_iterator affinity_lookup(const vector< pair<string, Affinity> > &m_data, 
    const string &m_pdb, const char *m_file_name, const size_t &m_line_number);
Affinity surrogate_affinity_score(const PDBComplex &m_mol, const SymmetricMatrix<float> &m_param);
SymmetricMatrix<float> surrogate_affinity_param(const unordered_map<string, unsigned int> &m_atom_table, struct drand48_data *m_rand_ptr);

#define AFFINITY_LOOKUP(DATA,  PDB) affinity_lookup(DATA, PDB, __FILE__, __LINE__)

float pearson_r(const deque<ClusterPredictions> &m_clusters, const bool m_use_weights);
float weighted_pearson_r(const deque<ClusterPredictions> &m_clusters);
float pearson_pvalue(float m_r /*copy*/, const deque<ClusterPredictions> &m_clusters, 
    const size_t &m_num_permute, struct drand48_data *m_rand_ptr, const bool m_use_weights);

float compute_rmse(const deque<ClusterPredictions> &m_clusters, const bool m_use_weights);
float compute_R2(const deque<ClusterPredictions> &m_predict, const bool m_use_weights);
float compute_R2(const deque<ClusterPredictions> &m_predict, const deque<ClusterPredictions> &m_reference, 
    const bool m_use_weights);

vector< pair<string, Affinity> > random_affinity(const vector< pair<float, vector<Cluster> > > &m_clusters, 
    struct drand48_data *m_rand_ptr);
string cluster_id(Cluster m_cluster /*copy*/);

void remove_atoms(vector<AminoAcidData> &m_chain, const unsigned int &m_atom_type);

// DEBUG
float feature_similarity(const Features &m_a, const Features &m_b);

int main(int argc, char* argv[])
{
    try{

        #ifdef USE_MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        signal( SIGINT, terminate_program );
        signal( SIGTERM, terminate_program );
        signal( SIGSEGV, terminate_program );

        start_time = MPI_Wtime();
        #endif // USE_MPI

        Options opt;
        
        if(IS_ROOT){

            // Per Department of Defense software release rules, this statement must appear every time the
	        // program is run
            cerr << "Distribution Statement A: Approved for public release: distribution is unlimited." << endl;

            opt.load(argc, argv);
        }

        #ifdef USE_MPI
        broadcast(opt, 0);
        #endif // USE_MPI

        if(opt.quit){

            #ifdef USE_MPI
            MPI_Finalize();
            #endif // USE_MPI

            return EXIT_SUCCESS;
        }

        vector< pair<string, Affinity> > affinity_data;
        vector< pair<float /*min distance between clusters*/, vector<Cluster> > > clusters;
        unordered_map<string, PDBComplex> pdb_data;
        unordered_map<string, vector<float> > rosetta_scores;
        unordered_map<string, vector<float> > scores;
        unordered_map<string, unsigned int> atom_table;

        ofstream fout( (IS_ROOT && !opt.output_file.empty()) ? opt.output_file.c_str() : "/dev/null");

        if(!fout){
            throw __FILE__ ":main: Unable to open output file for writing";
        }

        ostream &out = (IS_ROOT && opt.output_file.empty()) ? cout : fout; 
        
        // Every MPI rank gets a different random number seed
        opt.seed += mpi_rank;

        // Convert the seed to a format suitable for the drand48_r random number generator
        // (this is Linux-specific reentrant version of drand48).
        struct drand48_data rand_buffer;

        srand48_r(opt.seed, &rand_buffer);

        // Parse the affinity data and list of PDB files to analyze
        if(IS_ROOT){

            out << "Command line:";

            for(int i = 0;i < argc;++i){
                out << " " << argv[i];
            }

            out << endl;

            if(opt.use_weighted_stats){
                out << "Reporting the ** weighted **  statistical metrics" << endl;
            }
            else{
                out << "Reporting the ** unweighted ** statistical metrics" << endl;
            }

            if(opt.random_affinity){
                out << "** Random affinity data! **" << endl;
            }
            else if(opt.surrogate_affinity){
                out << "** Score-based surrogate affinity data! **" << endl;
            }
            else{
                out << "Affinity input data = " << opt.affinity_file << endl;
            }
            
            out << "PDB data directory = " << opt.pdb_dir << endl;

            if( !opt.rosetta_score_file.empty() ){

                out << "Rosetta score file = " << opt.rosetta_score_file << endl;

                out << "Rosetta scores will " << (opt.append_rosetta_scores ? 
                    "be appended to structure-based features" :
                    "be the only features") << endl;
            }

            if( !opt.score_file.empty() ){

                out << "Score file = " << opt.score_file << endl;

                out << "Scores will " << (opt.append_scores ? 
                    "be appended to structure-based features" :
                    "be the only features") << endl;
            }

            out << "Clustering file = " << opt.cluster_file << endl;
            out << "Force dimer? = " << (opt.force_dimer ? "true" : "false") << endl;
            out << "Ignore hetatom? = " << (opt.ignore_hetatom ? "true" : "false") << endl;
            out << "seed = " << opt.seed << endl;
            out << "Num permutations (p-value) = " << opt.num_permute << endl;
            out << "Num fold = " << opt.num_fold << endl;
            out << "Num trial = " << opt.num_trial << endl;
            out << "Forest size = " << opt.forest_size << endl;
            out << "Forest data bag = " << opt.forest_data_bag << endl;
            out << "Forest feature bag = " << opt.forest_feature_bag << endl;
            out << "Forest leaf = " << opt.forest_leaf << endl;
            
            if(opt.use_forest_weight){
                out << "Applying inverse cluster size training data weights" << endl;
            }
            else{
                out << "All training points are equally weighted" << endl;
            }

            if( opt.affinity_min_year != std::numeric_limits<unsigned int>::lowest() ){
                out << "Removing affinity data earlier than " << opt.affinity_min_year << endl;
            }

            if( opt.affinity_max_year != std::numeric_limits<unsigned int>::max() ){
                out << "Removing affinity data later than " << opt.affinity_max_year << endl;
            }

            if( opt.affinity_min != 0.0 ){
                out << "Removing affinity data < " << opt.affinity_min << endl;
            }

            if( opt.affinity_max != 1.0 ){
                out << "Removing affinity data > " << opt.affinity_max << endl;
            }

            out << "Using histogram features" << endl;
            out << "Min histogram distance = " << opt.feature_hist_min << endl;
            out << "Max histogram distance = " << opt.feature_hist_max << endl;
            out << "Num histogram bins = " << opt.feature_hist_num_bin << endl;
            
            out << "There are " << AminoAcidData::NUM_AMINO_ACID_TYPES << " amino acid types:";
            
            for(unsigned int i = 0;i < AminoAcidData::NUM_AMINO_ACID_TYPES;++i){
                out << " '" << amino_acid_names[i] << '\'';
            }

            out << endl;

            // Parse the PDB data
            parse_pdb_data(pdb_data, opt.pdb_dir, atom_table, opt.ignore_hetatom, opt.ignore_charge, opt.force_dimer, out);

            out << "Found " << pdb_data.size() << " protein complexes in " << opt.pdb_dir << endl;

            //#define REMOVE_HYDROGENS
            #ifdef REMOVE_HYDROGENS
            // Unless hydrogens have been added to X-ray structures, remove all hydrogens to put
            // NMR and X-ray structures on an equal footing.
            unordered_map<string, unsigned int>::const_iterator remove_hydrogen = atom_table.find("H");

            if( remove_hydrogen != atom_table.end() ){

                for(unordered_map<string, PDBComplex>::iterator i = pdb_data.begin();i != pdb_data.end();++i){
                    remove_atoms(i->second.first, remove_hydrogen->second);
                    remove_atoms(i->second.second, remove_hydrogen->second);
                }

                atom_table.erase(remove_hydrogen);

                out << "Removed all hydrogen atoms" << endl;
            }
            #endif // REMOVE_HYDROGENS

            out << "Atom table has " << atom_table.size() << " entries" << endl;

            #define PRINT_ATOM_INVENTORY
            #ifdef PRINT_ATOM_INVENTORY
            // Count the number of atom types
            unordered_map<unsigned int /*type*/, size_t /*count*/> atom_count;

            for(unordered_map<string, PDBComplex>::const_iterator i = pdb_data.begin();i != pdb_data.end();++i){

                for(vector<AminoAcidData>::const_iterator j = i->second.first.begin();j != i->second.first.end();++j){
                    for(vector<AtomData>::const_iterator k = j->atoms.begin();k != j->atoms.end();++k){
                        ++atom_count[k->type];
                    }
                }

                for(vector<AminoAcidData>::const_iterator j = i->second.second.begin();j != i->second.second.end();++j){
                    for(vector<AtomData>::const_iterator k = j->atoms.begin();k != j->atoms.end();++k){
                        ++atom_count[k->type];
                    }
                }
            }

            out << "Atom type inventory:" << endl;

            deque< pair<unsigned int, string> > atom_inventory;

            for(unordered_map<string, unsigned int>::const_iterator i = atom_table.begin();i != atom_table.end();++i){

                if( atom_count.find(i->second) != atom_count.end() ){
                    atom_inventory.push_back( make_pair(atom_count[i->second], i->first) );
                }
            }

            sort(atom_inventory.begin(), atom_inventory.end());

            for(deque< pair<unsigned int, string> >::const_reverse_iterator i = atom_inventory.rbegin();i != atom_inventory.rend();++i){
                out << '\t' << i->second << '\t' << i->first << endl;
            }
            #endif // PRINT_ATOM_INVENTORY

            out << "Computing hydrogen bond donor atoms ...";

            for(unordered_map<string, PDBComplex>::iterator i = pdb_data.begin();i != pdb_data.end();++i){

                inventory_hydrogen_bond_donors(i->second.first, atom_table, i->first + ":A");
                inventory_hydrogen_bond_donors(i->second.second, atom_table, i->first + ":B");

                inventory_interface_hydrogen_bond_acceptors(i->second, atom_table, i->first);
            }
            
            out << "done." << endl;

            out << "Computing interface hydrogen bond acceptor atoms ...";

            out << "done." << endl;

            if(opt.hist_dir.empty() == false){

                out << "Writing histograms to " << opt.hist_dir << endl;

                // Create the directory if needed
                struct stat dir_info;

                if (stat(opt.hist_dir.c_str(), &dir_info) == -1) {
                    mkdir(opt.hist_dir.c_str(), 0700);
                }

                for(unordered_map<string, PDBComplex>::const_iterator i = pdb_data.begin();i != pdb_data.end();++i){

                    ofstream fout(opt.hist_dir + "/" + i->first + ".hist");

                    if(!fout){
                        throw __FILE__ ":main: Unable to open histogram file for writing";
                    }

                    write_histogram(fout, i->second, atom_table, opt);
                }
            }

            // If the user has provided affinity values in a file, read them now so we can remove unneeded PDB
            // files prior to clustering
            if(!opt.random_affinity && !opt.surrogate_affinity){
                
                affinity_data = parse_affinity(opt.affinity_file);

                // Convert the Kd values into ln(Kd). Leave score value unchanged
                for(vector< pair<string, Affinity> >::iterator i = affinity_data.begin();i != affinity_data.end();++i){

                    if(i->second.affinity_type != Affinity::AFFINITY_SCORE){

                        if(i->second.value <= 0.0){
                            throw __FILE__ ":main: Unable to compute ln(Kd <= 0.0)";
                        }

                        i->second.value = log(i->second.value);
                    }
                }

                // Remove PDB files that do not have an associate affinity value prior to clustering
                deque<string> reaper;

                for(unordered_map<string, PDBComplex>::const_iterator i = pdb_data.begin();i != pdb_data.end();++i){
                    
                    // Note that affinity_data (as returned by the parse_affinity function) is sorted for fast lookup
                    vector< pair<string, Affinity> >::const_iterator iter = 
                        lower_bound( affinity_data.begin(), affinity_data.end(), i->first, find_by_accession() );

                    if( ( iter == affinity_data.end() ) || (iter->first != i->first) ){
                        reaper.push_back(i->first);
                    }
                }

                out << "Removing " << reaper.size() << " PDB records that do not have corresponding affinity values" << endl;

                while( !reaper.empty() ){

                    pdb_data.erase( reaper.back() );
                    reaper.pop_back();
                }
            }

            switch(opt.cluster_by){
                case Options::CLUSTER_BY_COMPOSITION:
                    cluster_by_composition(opt.cluster_file, pdb_data, opt);
                    break;
                case Options::CLUSTER_BY_SEQ_ALIGN:
                    cluster_by_seq_align(opt.cluster_file, pdb_data, opt);
                    break;
            };

            // Load the clusters from disk
            clusters = read_clusters(opt.cluster_file);

            out << "Read a total of " << clusters.size() << " cluster configurations" << endl;

            if(opt.random_affinity){

                // Generate random "affinity" values that reflect the sequence-based clustering of the protein complexes
                affinity_data = random_affinity(clusters, &rand_buffer);
            }
            else if(opt.surrogate_affinity){

                const SymmetricMatrix<float> param = surrogate_affinity_param(atom_table, &rand_buffer);

                affinity_data.clear();

                for(unordered_map<string, PDBComplex>::const_iterator i = pdb_data.begin();i != pdb_data.end();++i){
                    affinity_data.push_back( make_pair( i->first, surrogate_affinity_score(i->second, param) ) );
                }

                // The affinity data must be sorted by PDB accession for fast lookup
                sort(affinity_data.begin(), affinity_data.end());
            }

            out << "Found " << affinity_data.size() << " PDB files with associated affinity values" << endl;

            if( !opt.rosetta_score_file.empty() ){

                parse_rosetta_scores(opt.rosetta_score_file, rosetta_scores);

                // Append the default PDB id, ".0", to the accessions listed in the rosetta score file. This 
                // is required to allow processing of the SKEMPI dataset, which has multiple sequence variants for
                // each "wildtype" PDB file.
                if(false){
                    
                    cerr << "Appending the default PDB id \".0\" to Rosetta score PDB accessions" << endl;

                    unordered_map<string, vector<float> > new_rosetta_scores;

                    for(unordered_map<string, vector<float> >::const_iterator i = rosetta_scores.begin();i != rosetta_scores.end();++i){
                        new_rosetta_scores[i->first + ".0" ] = i->second;
                    }

                    new_rosetta_scores.swap(rosetta_scores);
                }
            }

            if( !opt.score_file.empty() ){

                parse_csv_scores(opt.score_file, scores);

                // Append the default PDB id, ".0", to the accessions listed in the score file. This 
                // is reuqired to allow processing of the SKEMPI dataset, which has multiple sequence variants for
                // each "wildtype" PDB file.
                if(false){
                    
                    cerr << "Appending the default PDB id \".0\" to score PDB accessions" << endl;

                    unordered_map<string, vector<float> > new_scores;

                    for(unordered_map<string, vector<float> >::const_iterator i = scores.begin();i != scores.end();++i){
                        new_scores[i->first + ".0" ] = i->second;
                    }

                    new_scores.swap(scores);
                }
            }

            // Only include the affinity data that has a valid associated molecular structure
            // and is not explicitly ignored by the user
            vector< pair<string, Affinity> > new_affinity_data;
            
            new_affinity_data.reserve( affinity_data.size() );

            for(vector< pair<string, Affinity> >::const_iterator i = affinity_data.begin();i != affinity_data.end();++i){

                bool valid_pdb = true;

                if( pdb_data.find(i->first) == pdb_data.end() ){
                    valid_pdb = false;
                }

                if(!rosetta_scores.empty()  && (rosetta_scores.find(i->first) == rosetta_scores.end() ) ){
                    valid_pdb = false;
                }

                if(!scores.empty()  && (scores.find(i->first) == scores.end() ) ){
                    valid_pdb = false;
                }

                if( opt.ignore_kd && (i->second.affinity_type == Affinity::AFFINITY_KD) ){
                    valid_pdb = false;
                }

                if( opt.ignore_ki && (i->second.affinity_type == Affinity::AFFINITY_KI) ){
                    valid_pdb = false;
                }

                if( opt.ignore_ic50 && (i->second.affinity_type == Affinity::AFFINITY_IC50) ){
                    valid_pdb = false;
                }

                if( (i->second.year < opt.affinity_min_year) || (i->second.year > opt.affinity_max_year) ){
                    valid_pdb = false;
                }

                if( (opt.affinity_min > 0.0) && (i->second.value < log(opt.affinity_min)) ){
                    valid_pdb = false;
                }

                if( (opt.affinity_max < 1.0) && (i->second.value > log(opt.affinity_max)) ){
                    valid_pdb = false;
                }

                if(valid_pdb){
                    new_affinity_data.push_back(*i);
                }
                else{
                    out << "Skipping PDB structure: " << i->first << endl;
                }
            }

            if(new_affinity_data.size() < affinity_data.size()){

                out << "Removed " << ( affinity_data.size() - new_affinity_data.size() ) 
                    << " affinity measurements that do not have associated structural or score information" << endl;
                
                out << "There are " << new_affinity_data.size() 
                    << " valid protein complexes with associated affinity values and feature values" << endl;
            }

            new_affinity_data.swap(affinity_data);

            // Summarize the affinity metadata
            vector<size_t> structure_summary(Affinity::NUM_STRUCTURE);
            vector<size_t> affinity_summary(Affinity::NUM_AFFINITY);

            for(vector< pair<string, Affinity> >::const_iterator i = affinity_data.begin();i != affinity_data.end();++i){

                ++structure_summary[i->second.structure_type];
                ++affinity_summary[i->second.affinity_type];
            }

            out << "Affinity measurement summary:" << endl;
            out << "\t" << affinity_summary[Affinity::AFFINITY_KD] << " Kd affinities" << endl;
            out << "\t" << affinity_summary[Affinity::AFFINITY_KI] << " Ki affinities" << endl;
            out << "\t" << affinity_summary[Affinity::AFFINITY_IC50] << " IC50 affinities" << endl;
            out << "\t" << affinity_summary[Affinity::AFFINITY_SCORE] << " computational scores" << endl;

            out << "Structure summary:" << endl;
            out << "\t" << structure_summary[Affinity::STRUCTURE_X_RAY] << " X-ray structures" << endl;
            out << "\t" << structure_summary[Affinity::STRUCTURE_NMR] << " NMR structures" << endl;
            out << "\t" << structure_summary[Affinity::STRUCTURE_PREDICTED] << " predicted structures" << endl;

            // Remove any proteins from the cluster list that do not have associated affinity data
            // If the clustering was computed for all protein structures, this can include structures
            // that don't have valid affinity measurements.
            vector< pair<float, vector<Cluster> > > new_clusters;

            new_clusters.reserve(clusters.size());

            for(vector< pair<float, vector<Cluster> > >::const_iterator i = clusters.begin();i != clusters.end();++i){

                new_clusters.push_back( make_pair(i->first, vector<Cluster>() ) );

                for(vector<Cluster>::const_iterator j = i->second.begin();j != i->second.end();++j){

                    new_clusters.back().second.push_back( Cluster() );

                    for(Cluster::const_iterator k = j->begin();k != j->end();++k){

                        vector< pair<string, Affinity> >::const_iterator iter = 
                            lower_bound(affinity_data.begin(), affinity_data.end(), *k, find_by_accession());

                        if( ( iter != affinity_data.end() ) && (*k == iter->first) ){
                            // This protein accession has a corresponding valid affinity value
                            new_clusters.back().second.back().push_back(*k);
                        }
                    }

                    // If the current cluster is empty, remove it
                    if( new_clusters.back().second.back().empty() ){
                        new_clusters.back().second.pop_back();
                    }
                }

                if( new_clusters.back().second.empty() ){
                    new_clusters.pop_back();
                }
            }

            // Swap out the clusters so that "cluster" will now only contain protein accessions with valid
            // affinity values
            new_clusters.swap(clusters);
        }

        #ifdef USE_MPI
        // Share the affinity, (optional) rosetta scores, structural and cluster data with all nodes
        broadcast(affinity_data, 0);
        broadcast(rosetta_scores, 0);
        broadcast(scores, 0);
        broadcast(pdb_data, 0);
        broadcast(clusters, 0);
        broadcast(atom_table, 0);
        #endif // USE_MPI

        // Predict protein-complex affinity for each clustering
        const size_t num_clusterings = clusters.size();
        const size_t num_affinity = affinity_data.size();
        
        /////////////////////////////////////////////////////////////////
        #define DEFINE_AND_ALLOCATE(VAR, LEN) \
            float *VAR = new float[LEN]; \
            if(VAR == NULL){ \
                throw __FILE__ ":main: Unable to allocate " #VAR; \
            } \
            memset(VAR, 0x0, sizeof(float)*LEN);

        /////////////////////////////////////////////////////////////////
        DEFINE_AND_ALLOCATE(baseline_ave_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(baseline_stdev_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(baseline_ave_p_value, num_clusterings);

        DEFINE_AND_ALLOCATE(model_ave_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(model_stdev_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(model_ave_p_value, num_clusterings);

        DEFINE_AND_ALLOCATE(ave_prediction, num_affinity);
        DEFINE_AND_ALLOCATE(stdev_prediction, num_affinity);
        //////////////////////////////////////////////////////////////////
        // baseline RMSE
        DEFINE_AND_ALLOCATE(baseline_ave_rmse, num_clusterings);
        DEFINE_AND_ALLOCATE(baseline_stdev_rmse, num_clusterings);

        // R2 of baseline relative to average
        DEFINE_AND_ALLOCATE(baseline_ave_R2, num_clusterings);
        DEFINE_AND_ALLOCATE(baseline_stdev_R2, num_clusterings);

        // model RMSE
        DEFINE_AND_ALLOCATE(model_ave_rmse, num_clusterings);
        DEFINE_AND_ALLOCATE(model_stdev_rmse, num_clusterings);

        // R2 of model relative to average
        DEFINE_AND_ALLOCATE(model_ave_R2, num_clusterings);
        DEFINE_AND_ALLOCATE(model_stdev_R2, num_clusterings);

        // R2 of model relative to baseline
        DEFINE_AND_ALLOCATE(model_ave_R2_baseline, num_clusterings);
        DEFINE_AND_ALLOCATE(model_stdev_R2_baseline, num_clusterings);
        //////////////////////////////////////////////////////////////////

        #ifdef USE_GSL
        DEFINE_AND_ALLOCATE(linear_ave_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(linear_stdev_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(linear_ave_p_value, num_clusterings);
        #endif // USE_GSL

        //////////////////////////////////////////////////////////////////

        DEFINE_AND_ALLOCATE(nn_ave_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(nn_stdev_correlation, num_clusterings);
        DEFINE_AND_ALLOCATE(nn_ave_p_value, num_clusterings);

        //////////////////////////////////////////////////////////////////
        // Precompute all of the features for the structures that have associated affinity data
        unordered_map<string, Features> baseline_features;
        unordered_map<string, Features> structural_features;

        #ifdef USE_GSL
        unordered_map< string, vector<double> > linear_regression_features;
        #endif // USE_GSL

        time_t profile = time(NULL);

        if(IS_ROOT){
            cerr << "Computing features for " << affinity_data.size() << " structures ... ";
        }

        #pragma omp parallel for
        for(size_t i = 0;i < num_affinity;++i){
            
            unordered_map<string, PDBComplex>::const_iterator j = pdb_data.find(affinity_data[i].first);

            if( j == pdb_data.end() ){
                throw __FILE__ ":main: Unable to lookup protein complex by id";
            }

            /*const*/ Features baseline = composition_features(j->second);

            // DEBUG
            //baseline.push_back(affinity_data[i].second.structure_ph);
            //baseline.push_back(affinity_data[i].second.structure_temperature);

            Features structural;

            if( ( rosetta_scores.empty() || (opt.append_rosetta_scores == true) ) &&
                ( scores.empty() || (opt.append_scores == true) ) ){
                
                structural = extract_features(j->second, atom_table, affinity_data[i].first, opt);

                // DEBUG
                //structural.push_back(affinity_data[i].second.structure_ph);
                //structural.push_back(affinity_data[i].second.structure_temperature);
            }

            unordered_map< string, vector<float> >::const_iterator score_iter = rosetta_scores.find(affinity_data[i].first);

            if( score_iter != rosetta_scores.end() ){
                
                if(opt.append_rosetta_scores){
                    
                    // Append the rosetta scores to the distance features
                    structural.insert( structural.end(), score_iter->second.begin(), score_iter->second.end() );
                }
                else{
                    structural = score_iter->second;
                }
            }

            score_iter = scores.find(affinity_data[i].first);

            if( score_iter != scores.end() ){
                
                if(opt.append_scores){
                    
                    // Append the scores to the distance features
                    structural.insert( structural.end(), score_iter->second.begin(), score_iter->second.end() );
                }
                else{
                    structural = score_iter->second;
                }
            }

            #ifdef USE_GSL
            const vector<double> linear = opt.use_linear_regression ? 
                linear_features(j->second, atom_table) : 
                vector<double>();
            #endif // USE_GSL

            #pragma omp critical
            {
                baseline_features[affinity_data[i].first] = baseline;
                structural_features[affinity_data[i].first] = structural;

                #ifdef USE_GSL
                linear_regression_features[affinity_data[i].first] = linear;
                #endif // USE_GSL
            }
        }

        profile = time(NULL) - profile;

        if(IS_ROOT){
            cerr << "done in " << profile << " sec." << endl;
        }

        for(size_t c = mpi_rank;c < num_clusterings;c += mpi_numtasks){

            vector<Cluster> &curr_clustering = clusters[c].second;

            const size_t num_cluster = curr_clustering.size();

            // For a given clustering, we need at least num_fold clusters to perform cross validation
            if(num_cluster < opt.num_fold){
                break;
            }

            // Repeat the cross validation proceedure num_trial times and retain the cluster structure of
            // the input data
            for(size_t trial = 0;trial < opt.num_trial;++trial){

                // Predict the test set
                deque<ClusterPredictions> baseline_predictions;

                #ifdef USE_GSL
                deque<ClusterPredictions> linear_predictions;
                #endif // USE_GSL

                deque<ClusterPredictions> nn_predictions;
                deque<ClusterPredictions> model_predictions;

                // Predict the *training* set to compute prediction residual
                deque<ClusterPredictions> residual_baseline_predictions;

                #ifdef USE_GSL
                deque<ClusterPredictions> residual_linear_predictions;
                #endif // USE_GSL

                deque<ClusterPredictions> residual_nn_predictions;
                deque<ClusterPredictions> residual_model_predictions;

                // Randomly shuffle the order of the protein-complex clusters
                randomize(curr_clustering.begin(), curr_clustering.end(), &rand_buffer);

                // Perform cross validation
                for(size_t fold = 0;fold < opt.num_fold;++fold){

                    // Partition the clusters into test and training
                    deque<Cluster> test_cluster;
                    deque<Cluster> train_cluster;

                    for(size_t i = 0;i < num_cluster;++i){

                        if(i%opt.num_fold == fold){
                            test_cluster.push_back( curr_clustering[i] );
                        }
                        else{
                            train_cluster.push_back( curr_clustering[i] );
                        }
                    }

                    if(IS_ROOT){
                        
                        cerr << num_cluster << " clusters; iteration " << trial << "; fold " << fold << endl;
                        cerr << "\t|test| = " << test_cluster.size() << endl;
                        cerr << "\t|train| = " << train_cluster.size() << endl;
                    }

                    vector<WeightedValue> train_y;
                    float ave_train_y = 0.0;
                    vector<string> train_x_label;

                    for(deque<Cluster>::const_iterator i = train_cluster.begin();i != train_cluster.end();++i){

                        if( i->empty() ){
                            throw __FILE__ ":main: Empty cluster";
                        }

                        const float cluster_weight = opt.use_forest_weight ? 1.0/i->size() : 1.0;

                        const unsigned int cluster_index = i - train_cluster.begin();

                        for(Cluster::const_iterator j = i->begin();j != i->end();++j){

                            vector< pair<string, Affinity> >::const_iterator affinity_iter = AFFINITY_LOOKUP(affinity_data, *j);

                            train_y.push_back( WeightedValue(affinity_iter->second.value, cluster_weight, cluster_index) );
                            ave_train_y += affinity_iter->second.value*cluster_weight;
                            train_x_label.push_back(affinity_iter->first);
                        }
                    }

                    ave_train_y /= train_cluster.size();

                    // Select atomic-level features in the form of a histogram weighting function.
                    const size_t num_x = train_x_label.size();

                    vector<const Features*> baseline_ptr_x(num_x);

                    #ifdef USE_GSL
                    vector<const vector<double>*> linear_ptr_x(num_x);
                    #endif // USE_GSL

                    vector<const Features*> structural_ptr_x(num_x);

                    for(size_t i = 0;i < num_x;++i){

                        unordered_map<string, Features>::const_iterator j = baseline_features.find(train_x_label[i]);

                        if( j == baseline_features.end() ){
                            throw __FILE__ ":main: Unable to lookup baseline features by id";
                        }

                        baseline_ptr_x[i] = &(j->second);

                        j = structural_features.find(train_x_label[i]);

                        if( j == structural_features.end() ){
                            throw __FILE__ ":main: Unable to lookup structural features by id";
                        }

                        structural_ptr_x[i] = &(j->second);

                        #ifdef USE_GSL
                        unordered_map<string, vector<double> >::const_iterator linear_iter = 
                            linear_regression_features.find(train_x_label[i]);

                        if( linear_iter == linear_regression_features.end() ){
                            throw __FILE__ ":main: Unable to lookup linear features by id";
                        }

                        linear_ptr_x[i] = &(linear_iter->second);
                        #endif // USE_GSL
                    }

                    // Train the models
                    RandomForest baseline_model(opt.forest_size, opt.forest_leaf, opt.forest_data_bag, 
                        opt.forest_feature_bag, &rand_buffer);

                    baseline_model.build(train_y, baseline_ptr_x, false /*verbose*/);

                    #ifdef USE_GSL
                    // The linear_model is a vector of feature weights. The *first* element is the bias
                    const vector<double> linear_model = opt.use_linear_regression ? 
                        compute_linear_regression(train_y, linear_ptr_x, false /*no bias*/) :
                        vector<double>(); 
                    #endif // USE_GSL

                    RandomForest model(opt.forest_size, opt.forest_leaf, opt.forest_data_bag, 
                        opt.forest_feature_bag, &rand_buffer);

                    model.build(train_y, structural_ptr_x, false /*verbose*/);

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Test set evaluation
                    /////////////////////////////////////////////////////////////////////////////////////////
                    float ave_test_y = 0.0;

                    for(deque<Cluster>::const_iterator i = test_cluster.begin();i != test_cluster.end();++i){

                        if( i->empty() ){
                            throw __FILE__ ":main: Empty test cluster";
                        }

                        ClusterPredictions local_baseline;
                        
                        #ifdef USE_GSL
                        ClusterPredictions local_linear;
                        #endif // USE_GSL

                        ClusterPredictions local_nn;
                        ClusterPredictions local_model;

                        local_baseline.reserve( i->size() );

                        #ifdef USE_GSL
                        local_linear.reserve( i->size() );
                        #endif // USE_GSL

                        local_nn.reserve( i->size() );
                        local_model.reserve( i->size() );

                        for(Cluster::const_iterator j = i->begin();j != i->end();++j){

                            vector< pair<string, Affinity> >::const_iterator affinity_iter = 
                                AFFINITY_LOOKUP(affinity_data, *j);

                            unordered_map<string, Features>::const_iterator baseline_iter = 
                                baseline_features.find(*j);

                            if( baseline_iter == baseline_features.end() ){
                                throw __FILE__ ":main: Unable to lookup baseline test features";
                            }

                            unordered_map<string, Features>::const_iterator structural_iter = 
                                structural_features.find(*j);

                            if( structural_iter == structural_features.end() ){
                                throw __FILE__ ":main: Unable to lookup structural test features";
                            }

                            #ifdef USE_GSL
                            unordered_map< string, vector<double> >::const_iterator linear_iter = 
                                linear_regression_features.find(*j);

                            if( linear_iter == linear_regression_features.end() ){
                                throw __FILE__ ":main: Unable to lookup linear test features";
                            }
                            #endif // USE_GSL

                            ave_test_y += affinity_iter->second.value/i->size();

                            local_baseline.push_back( make_pair(affinity_iter->second.value, 
                                baseline_model.predict(baseline_iter->second) ) );
                            
                            #ifdef USE_GSL
                            if(opt.use_linear_regression){
                                local_linear.push_back( make_pair(affinity_iter->second.value, 
                                    linear_prediction(linear_model, linear_iter->second) ) );
                            }
                            #endif // USE_GSL

                            if(opt.k_nearest_neighbors > 0){

                                // Use structural features for k-nearest neighbor classification
                                local_nn.push_back( make_pair(affinity_iter->second.value, 
                                    nn_prediction(*j, structural_iter->second, train_x_label, 
                                        structural_features, train_y, opt.k_nearest_neighbors, affinity_iter->second.value) ));
                            }

                            const float pred = model.predict(structural_iter->second);

                            local_model.push_back( make_pair(affinity_iter->second.value, pred ) );

                            #ifdef DEBUG_RANDOM_SAMPLE
                            // DEBUG
                            out << "DEBUG: " << affinity_iter->first << '\t' << affinity_iter->second.value << '\t' << pred << endl;

                            // Add random noise to the atomic coordinates -- how much to the predictions changes?
                            const float noise = 0.1;

                            unordered_map<string, PDBComplex>::const_iterator debug_iter = pdb_data.find(affinity_iter->first);

                            if( debug_iter == pdb_data.end() ){
                                throw __FILE__ ":main: Unable to lookup protein complex by id";
                            }
                            
                            for(size_t k = 0;k < 10;++k){

                                PDBComplex tmp = debug_iter->second;

                                // Add random noise to all atoms
                                #define ADD_NOISE(X) \
                                    for(vector<AminoAcidData>::iterator aa_iter = X.begin();aa_iter != X.end();++aa_iter){ \
                                        for(vector<AtomData>::iterator atom_iter = aa_iter->atoms.begin();atom_iter != aa_iter->atoms.end();++atom_iter){ \
                                            double r = 0.0; \
                                            drand48_r(&rand_buffer, &r); \
                                            atom_iter->x += 2.0*(r - 0.5)*noise; \
                                            drand48_r(&rand_buffer, &r); \
                                            atom_iter->y += 2.0*(r - 0.5)*noise; \
                                            drand48_r(&rand_buffer, &r); \
                                            atom_iter->x += 2.0*(r - 0.5)*noise; \
                                        } \
                                    }

                                ADD_NOISE(tmp.first);
                                ADD_NOISE(tmp.second);

                                const Features debug_structural = extract_features(tmp, atom_table, affinity_iter->first, opt);

                                out <<  "DEBUG: " << affinity_iter->first << '\t' << affinity_iter->second.value << '\t' 
                                    << model.predict(debug_structural) << endl;
                            }
                            #endif // DEBUG_RANDOM_SAMPLE

                            const size_t index = affinity_iter - affinity_data.begin();
                            
                            ave_prediction[index] += pred;
                            stdev_prediction[index] += pred*pred;
                        }

                        baseline_predictions.push_back(local_baseline);

                        #ifdef USE_GSL
                        linear_predictions.push_back(local_linear);
                        #endif // USE_GSL

                        nn_predictions.push_back(local_nn);
                        model_predictions.push_back(local_model);
                    }

                    ave_test_y /= test_cluster.size();

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Training set evaluation for residual calculation
                    /////////////////////////////////////////////////////////////////////////////////////////
                    for(deque<Cluster>::const_iterator i = train_cluster.begin();i != train_cluster.end();++i){

                        if( i->empty() ){
                            throw __FILE__ ":main: Empty train cluster";
                        }

                        ClusterPredictions local_baseline;

                        #ifdef USE_GSL
                        ClusterPredictions local_linear;
                        #endif // USE_GSL

                        ClusterPredictions local_nn;
                        ClusterPredictions local_model;

                        local_baseline.reserve( i->size() );

                        #ifdef USE_GSL
                        local_linear.reserve( i->size() );
                        #endif // USE_GSL

                        local_nn.reserve( i->size() );
                        local_model.reserve( i->size() );

                        for(Cluster::const_iterator j = i->begin();j != i->end();++j){

                            vector< pair<string, Affinity> >::const_iterator affinity_iter = 
                                AFFINITY_LOOKUP(affinity_data, *j);

                            unordered_map<string, Features>::const_iterator baseline_iter = 
                                baseline_features.find(*j);

                            if( baseline_iter == baseline_features.end() ){
                                throw __FILE__ ":main: Unable to lookup baseline features (residual)";
                            }

                            unordered_map<string, Features>::const_iterator structural_iter = 
                                structural_features.find(*j);

                            if( structural_iter == structural_features.end() ){
                                throw __FILE__ ":main: Unable to lookup structural features (residual)";
                            }

                            #ifdef USE_GSL
                            unordered_map< string, vector<double> >::const_iterator linear_iter = 
                                linear_regression_features.find(*j);

                            if( linear_iter == linear_regression_features.end() ){
                                throw __FILE__ ":main: Unable to lookup linear regression features (residual)";
                            }
                            #endif // USE_GSL

                            local_baseline.push_back( make_pair(affinity_iter->second.value, 
                                baseline_model.predict(baseline_iter->second) ) );
                            
                            local_model.push_back( make_pair(affinity_iter->second.value, 
                                model.predict(structural_iter->second) ) );

                            #ifdef USE_GSL
                            if(opt.use_linear_regression){
                                local_linear.push_back( make_pair(affinity_iter->second.value, 
                                    linear_prediction(linear_model, linear_iter->second) ) );
                            }
                            #endif // USE_GSL

                            if(opt.k_nearest_neighbors > 0){

                                // Use structural features for k-nearest neighbor classification
                                local_nn.push_back( make_pair(affinity_iter->second.value, 
                                    nn_prediction(*j, structural_iter->second, train_x_label,
                                        structural_features, train_y, opt.k_nearest_neighbors) ));
                            }
                        }

                        residual_baseline_predictions.push_back(local_baseline);

                        #ifdef USE_GSL
                        residual_linear_predictions.push_back(local_linear);
                        #endif // USE_GSL

                        residual_nn_predictions.push_back(local_nn);
                        residual_model_predictions.push_back(local_model);
                    }

                    if(IS_ROOT){
                    
                        cerr << "\t<training affinity> = " << ave_train_y << endl;
                        cerr << "\t<testing affinity> = " << ave_test_y << endl;

                        // Update the Pearson correlation coefficient based on the partial
                        // sampling performed by rank 0
                        cerr << "\tincremental baseline r = " << pearson_r(baseline_predictions, opt.use_weighted_stats) 
                            << "; training set residual r = " << pearson_r(residual_baseline_predictions, opt.use_weighted_stats) << endl;
                        
                        #ifdef USE_GSL
                        if(opt.use_linear_regression){
                            cerr << "\tincremental linear r = " << pearson_r(linear_predictions, opt.use_weighted_stats) 
                                << "; training set residual r = " << pearson_r(residual_linear_predictions, opt.use_weighted_stats) << endl;
                        }
                        #endif // USE_GSL

                        if(opt.k_nearest_neighbors > 0){

                            cerr << "\tincremental k-nn r = " << pearson_r(nn_predictions, opt.use_weighted_stats) 
                                << "; training set residual r = " << pearson_r(residual_nn_predictions, opt.use_weighted_stats) << endl;
                        }

                        cerr << "\tincremental model r = " << pearson_r(model_predictions, opt.use_weighted_stats) 
                            << "; training set residual r = " << pearson_r(residual_model_predictions, opt.use_weighted_stats) << endl;
                    }
                }

                if(IS_ROOT){
                
                    // Update the Pearson correlation coefficient based on the partial
                    // sampling performed by rank 0
                    cerr << "number of clusters = " << curr_clustering.size() << "; trial = " << trial << endl;
                    cerr << "\tbaseline r = " << pearson_r(baseline_predictions, opt.use_weighted_stats) << endl;

                    #ifdef USE_GSL
                    if(opt.use_linear_regression){
                        cerr << "\tlinear r = " << pearson_r(linear_predictions, opt.use_weighted_stats) << endl;
                    }
                    #endif // USE_GSL

                    if(opt.k_nearest_neighbors > 0){
                        cerr << "\tk-nn r = " << pearson_r(nn_predictions, opt.use_weighted_stats) << endl;
                    }

                    cerr << "\tmodel = " << pearson_r(model_predictions, opt.use_weighted_stats) << endl;
                }

                // Compute the correlation p-value
                const float r_baseline = pearson_r(baseline_predictions, opt.use_weighted_stats);

                const float p_baseline = (opt.num_permute == 0) ? -1.0 : pearson_pvalue(r_baseline, 
                        baseline_predictions, opt.num_permute, &rand_buffer, opt.use_weighted_stats);

                const float r_model = pearson_r(model_predictions, opt.use_weighted_stats);

                const float p_model = (opt.num_permute == 0) ? -1.0 : pearson_pvalue(r_model, 
                        model_predictions, opt.num_permute, &rand_buffer, opt.use_weighted_stats);

                baseline_ave_correlation[c] += r_baseline;
                baseline_stdev_correlation[c] += r_baseline*r_baseline;
                baseline_ave_p_value[c] += p_baseline;

                #ifdef USE_GSL
                float r_linear = 0.0;
                float p_linear = 0.0;

                if(opt.use_linear_regression){

                    r_linear = pearson_r(linear_predictions, opt.use_weighted_stats);

                    p_linear = (opt.num_permute == 0) ? -1.0 : pearson_pvalue(r_linear, 
                            linear_predictions, opt.num_permute, &rand_buffer, opt.use_weighted_stats);

                    linear_ave_correlation[c] += r_linear;
                    linear_stdev_correlation[c] += r_linear*r_linear;
                    linear_ave_p_value[c] += p_linear;
                }
                #endif // USE_GSL

                float r_nn = 0.0;
                float p_nn = 0.0;

                if(opt.k_nearest_neighbors > 0){

                    r_nn = pearson_r(nn_predictions, opt.use_weighted_stats);

                    p_nn = (opt.num_permute == 0) ? -1.0 : pearson_pvalue(r_nn, 
                            nn_predictions, opt.num_permute, &rand_buffer, opt.use_weighted_stats);

                    nn_ave_correlation[c] += r_nn;
                    nn_stdev_correlation[c] += r_nn*r_nn;
                    nn_ave_p_value[c] += p_nn;
                }

                model_ave_correlation[c] += r_model;
                model_stdev_correlation[c] += r_model*r_model;
                model_ave_p_value[c] += p_model;

                if(IS_ROOT){

                    cerr << curr_clustering.size() << '\t';
                    cerr << r_baseline << '\t' << p_baseline << '\t';

                    #ifdef USE_GSL
                    if(opt.use_linear_regression){
                        cerr << r_linear << '\t' << p_linear << '\t';
                    }
                    #endif // USE_GSL

                    cerr << r_model << '\t' << p_model << endl;                    
                }

                // Compute the RMSD and R2 metrics
                const float rmse_baseline = compute_rmse(baseline_predictions, opt.use_weighted_stats);
                
                baseline_ave_rmse[c] += rmse_baseline;
                baseline_stdev_rmse[c] += rmse_baseline*rmse_baseline;

                const float R2_baseline = compute_R2(baseline_predictions, opt.use_weighted_stats);
                
                baseline_ave_R2[c] += R2_baseline;
                baseline_stdev_R2[c] += R2_baseline*R2_baseline;

                const float rmse_model = compute_rmse(model_predictions, opt.use_weighted_stats);
                
                model_ave_rmse[c] += rmse_model;
                model_stdev_rmse[c] += rmse_model*rmse_model;

                const float R2_model = compute_R2(model_predictions, opt.use_weighted_stats);
                
                model_ave_R2[c] += R2_model;
                model_stdev_R2[c] += R2_model*R2_model;

                const float R2_model_baseline = compute_R2(model_predictions, baseline_predictions, opt.use_weighted_stats);
                
                model_ave_R2_baseline[c] += R2_model_baseline;
                model_stdev_R2_baseline[c] += R2_model_baseline*R2_model_baseline;
            }
        }

        #ifdef USE_MPI

        // Collect all of the results at rank 0 and write the results to the output file
        #define REDUCE(VAR, LEN) \
            MPI_Reduce( (mpi_rank == 0) ? MPI_IN_PLACE : VAR, VAR, LEN, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD)
        
        REDUCE(baseline_ave_correlation, num_clusterings);
        REDUCE(model_ave_correlation, num_clusterings);

        #ifdef USE_GSL
        REDUCE(linear_ave_correlation, num_clusterings);
        #endif // USE_GSL

        REDUCE(nn_ave_correlation, num_clusterings);

        REDUCE(baseline_stdev_correlation, num_clusterings);
        REDUCE(model_stdev_correlation, num_clusterings);

        #ifdef USE_GSL
        REDUCE(linear_stdev_correlation, num_clusterings);
        #endif // USE_GSL

        REDUCE(nn_stdev_correlation, num_clusterings);

        REDUCE(baseline_ave_p_value, num_clusterings);
        REDUCE(model_ave_p_value, num_clusterings);

        #ifdef USE_GSL
        REDUCE(linear_ave_p_value, num_clusterings);
        #endif // USE_GSL

        REDUCE(nn_ave_p_value, num_clusterings);

        REDUCE(baseline_ave_rmse, num_clusterings);
        REDUCE(baseline_stdev_rmse, num_clusterings);
        REDUCE(baseline_ave_R2, num_clusterings);
        REDUCE(baseline_stdev_R2, num_clusterings);
        REDUCE(model_ave_rmse, num_clusterings);
        REDUCE(model_stdev_rmse, num_clusterings);
        REDUCE(model_ave_R2, num_clusterings);
        REDUCE(model_stdev_R2, num_clusterings);
        REDUCE(model_ave_R2_baseline, num_clusterings);
        REDUCE(model_stdev_R2_baseline, num_clusterings);

        REDUCE(ave_prediction, num_affinity);
        REDUCE(stdev_prediction, num_affinity);
        
        #endif // USE_MPI

        if(mpi_rank == 0){
            
            if(opt.num_trial > 1){
                out << "min_cluster_distance\tnumber_of_clusters\tave_baseline_r\tstdev_baseline_r\tbaseline_p_value"
                    << "\tave_baseline_RMSE\tstdev_baseline_RMSE\tave_baseline_R2\tstdev_baseline_R2";

                #ifdef USE_GSL
                if(opt.use_linear_regression){
                    out << "\tave_linear_r\tstdev_linear_r\tlinear_p_value";
                }
                #endif // USE_GSL

                if(opt.k_nearest_neighbors > 0){
                    out << "\tave_nn_r\tstdev_nn_r\tnn_p_value";
                }

                out << "\tave_model_r\tstdev_model_r\tmodel_p_value" 
                    << "\tave_model_RMSE\tstdev_model_RMSE\tave_model_R2\tstdev_model_R2"
                    << "\tave_model_R2(baseline)\tstdev_model_R2(baseline)" << endl;
            }
            else{
                out << "min_cluster_distance\tnumber_of_clusters\tave_baseline_r"
                    << "\tave_baseline_RMSE\tave_baseline_R2";

                #ifdef USE_GSL
                if(opt.use_linear_regression){
                    out << "\tave_linear_r";
                }
                #endif // USE_GSL

                if(opt.k_nearest_neighbors){
                    out << "\tave_nn_r";
                }

                out << "\tave_model_r\tave_model_RMSE\tave_model_R2\tave_model_R2(baseline)" << endl;
            }

            for(size_t c = 0;c < num_clusterings;++c){

                #define COMPUTE_STDEV(AVE, STDEV) \
                    sqrt(fabs(STDEV[c]/opt.num_trial - AVE[c]*AVE[c]))

                baseline_ave_correlation[c] /= opt.num_trial;
                model_ave_correlation[c] /= opt.num_trial;
                
                baseline_stdev_correlation[c] = COMPUTE_STDEV(baseline_ave_correlation, baseline_stdev_correlation);
                model_stdev_correlation[c] = COMPUTE_STDEV(model_ave_correlation, model_stdev_correlation);

                baseline_ave_p_value[c] /= opt.num_trial;
                model_ave_p_value[c] /= opt.num_trial;

                baseline_ave_rmse[c] /= opt.num_trial;
                baseline_ave_R2[c] /= opt.num_trial;

                baseline_stdev_rmse[c] = COMPUTE_STDEV(baseline_ave_rmse, baseline_stdev_rmse);
                baseline_stdev_R2[c] = COMPUTE_STDEV(baseline_ave_R2, baseline_stdev_R2);

                model_ave_rmse[c] /= opt.num_trial;
                model_ave_R2[c] /= opt.num_trial;
                model_ave_R2_baseline[c] /= opt.num_trial;

                model_stdev_rmse[c] = COMPUTE_STDEV(model_ave_rmse, model_stdev_rmse);
                model_stdev_R2[c] = COMPUTE_STDEV(model_ave_R2, model_stdev_R2);
                model_stdev_R2_baseline[c] = COMPUTE_STDEV(model_ave_R2_baseline, model_stdev_R2_baseline);

                #ifdef USE_GSL
                if(opt.use_linear_regression){

                    linear_ave_correlation[c] /= opt.num_trial;

                    linear_stdev_correlation[c] = COMPUTE_STDEV(linear_ave_correlation, linear_stdev_correlation);
                    linear_ave_p_value[c] /= opt.num_trial;
                }
                #endif // USE_GSL

                if(opt.k_nearest_neighbors > 0){
                    
                    nn_ave_correlation[c] /= opt.num_trial;
                    nn_stdev_correlation[c] = COMPUTE_STDEV(nn_ave_correlation, nn_stdev_correlation);

                    nn_ave_p_value[c] /= opt.num_trial;
                }

                if(opt.num_trial > 1){

                    out << clusters[c].first << '\t' << clusters[c].second.size() << '\t';
                    out << baseline_ave_correlation[c] << '\t' << baseline_stdev_correlation[c] << '\t' << baseline_ave_p_value[c] << '\t';
                    out << baseline_ave_rmse[c] << '\t' << baseline_stdev_rmse[c] << '\t' << baseline_ave_R2[c] << '\t' << baseline_stdev_R2[c] << '\t';

                    #ifdef USE_GSL
                    if(opt.use_linear_regression){
                        out << linear_ave_correlation[c] << '\t' << linear_stdev_correlation[c] << '\t' << linear_ave_p_value[c] << '\t';
                    }
                    #endif // USE_GSL

                    if(opt.k_nearest_neighbors > 0){
                        out << nn_ave_correlation[c] << '\t' << nn_stdev_correlation[c] << '\t' << nn_ave_p_value[c] << '\t';
                    }

                    out << model_ave_correlation[c] << '\t' << model_stdev_correlation[c] << '\t' << model_ave_p_value[c] << '\t';
                    out << model_ave_rmse[c] << '\t' << model_stdev_rmse[c] << '\t' << model_ave_R2[c] << '\t' << model_stdev_R2[c] << '\t';
                    out << model_ave_R2_baseline[c] << '\t' << model_stdev_R2_baseline[c] << endl;
                }
                else{

                    out << clusters[c].first << '\t' << clusters[c].second.size() << '\t';
                    out << baseline_ave_correlation[c] << '\t' << baseline_ave_rmse[c] << '\t' << baseline_ave_R2[c];

                    #ifdef USE_GSL
                    if(opt.use_linear_regression){
                        out << linear_ave_correlation[c] << '\t';
                    }
                    #endif // USE_GSL

                    if(opt.k_nearest_neighbors > 0){
                        out << nn_ave_correlation[c] << '\t';
                    }

                    out << model_ave_correlation[c] << '\t' << model_ave_rmse[c] << '\t' << model_ave_R2[c] << '\t' << model_ave_R2_baseline[c] << endl;
                }
            }

            // Print the per-structure predicition results
            if(true){

                if(opt.num_trial*num_clusterings > 1){
                    out << "\naccession\tlog(Kd)\tave predict log(Kd)\tstdev predict log(Kd)" << endl;
                }
                else{
                    out << "\naccession\tlog(Kd)\tave predict log(Kd)" << endl;
                }

                for(size_t i = 0;i < num_affinity;++i){

                    ave_prediction[i] /= opt.num_trial*num_clusterings;

                    if(opt.num_trial*num_clusterings > 1){
                        stdev_prediction[i] = sqrt(stdev_prediction[i]/(opt.num_trial*num_clusterings) - ave_prediction[i]*ave_prediction[i]);

                        out << affinity_data[i].first << '\t' << affinity_data[i].second.value << '\t' 
                            << ave_prediction[i] << '\t' << stdev_prediction[i];
                    }
                    else{
                        out << affinity_data[i].first << '\t' << affinity_data[i].second.value << '\t' << ave_prediction[i];
                    }

                    #ifdef OUTPUT_CHAIN_LENGTHS
                    unordered_map<string, PDBComplex>::const_iterator structure_iter = pdb_data.find(affinity_data[i].first);

                    if( structure_iter == pdb_data.end() ){
                        throw __FILE__ ":main: Unable to lookup protein complex by id";
                    }

                    out << '\t' << (structure_iter->second.first.size() + structure_iter->second.second.size()) 
                        << '\t' << min(structure_iter->second.first.size(), structure_iter->second.second.size())
                        << '\t' << max(structure_iter->second.first.size(), structure_iter->second.second.size()) << endl;
                    #else
                    out << endl;
                    #endif // OUTPUT_CHAIN_LENGTHS
                }
            }
        }
        
        #define DELETE(VAR) \
            if(VAR != NULL){ \
                \
                delete [] VAR; \
                VAR = NULL; \
            }

        DELETE(baseline_ave_correlation);
        DELETE(model_ave_correlation);

        #ifdef USE_GSL
        DELETE(linear_ave_correlation);
        #endif // USE_GSL

        DELETE(nn_ave_correlation);
        DELETE(baseline_stdev_correlation);
        DELETE(model_stdev_correlation);

        #ifdef USE_GSL
        DELETE(linear_stdev_correlation);
        #endif // USE_GSL

        DELETE(nn_stdev_correlation);

        DELETE(baseline_ave_p_value);
        DELETE(model_ave_p_value);

        #ifdef USE_GSL
        DELETE(linear_ave_p_value);
        #endif // USE_GSL

        DELETE(nn_ave_p_value);

        DELETE(baseline_ave_rmse);
        DELETE(baseline_stdev_rmse);
        DELETE(baseline_ave_R2);
        DELETE(baseline_stdev_R2);
        DELETE(model_ave_rmse);
        DELETE(model_stdev_rmse);
        DELETE(model_ave_R2);
        DELETE(model_stdev_R2);
        DELETE(model_ave_R2_baseline);
        DELETE(model_stdev_R2_baseline);

        DELETE(ave_prediction);
        DELETE(stdev_prediction);
    }
    catch(const char *error){

        cerr << "Caught the error: " << error << endl;
        
        #ifdef USE_MPI
        MPI_Finalize();
        #endif // USE_MPI

        return EXIT_FAILURE;
    }
    catch(...){

        cerr << "Caught an unhandled error" << endl;

        #ifdef USE_MPI
        MPI_Finalize();
        #endif // USE_MPI

        return EXIT_FAILURE;
    }

    #ifdef USE_MPI
    MPI_Finalize();
    #endif // USE_MPI

    return EXIT_SUCCESS;
}

// Report our rank, the signal we caught and the time spent running
void terminate_program(int m_sig)
{
    cerr << "[" << mpi_rank << "] caught signal " << m_sig << endl;
    
    #ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 0);
    #endif // USE_MPI
}

vector< pair<string, Affinity> > parse_affinity(const string &m_affinity_file)
{
    vector< pair<string, Affinity> > ret;

    ifstream fin( m_affinity_file.c_str() );

    if(!fin){
        cerr << "Unable to open " << m_affinity_file << " for reading affinity data" << endl;
        throw __FILE__ ":parse_affinity: Unable to open affinity data file";
    }

    string line;
    size_t line_number = 0;

    // Skip any comments
    while(true){

        ++line_number;

        if( !getline(fin, line) ){
            throw __FILE__ ":parse_affinity: Unable to read header";
        }

        const string::size_type loc = line.find('#');

        if(loc != string::npos){
            line = line.substr(0, loc);
        }

        // Make sure we have a comma in the line we just read
        if(line.find(',') != string::npos){
            break;
        }
    }

    const vector<string> header = split(line, ',');

    const size_t accession_col = get_col(header, "accession");

    if(accession_col == COLUMN_NOT_FOUND){
        throw __FILE__ ":parse_affinity: Unable to find column \"accession\"";
    }

    const size_t affinity_col = get_col(header, "affinity_value");

    if(affinity_col == COLUMN_NOT_FOUND){
        throw __FILE__ ":parse_affinity: Unable to find column \"affinity_value\"";
    }
    
    const size_t affinity_type_col = get_col(header, "affinity_type");

    if(affinity_type_col == COLUMN_NOT_FOUND){
        throw __FILE__ ":parse_affinity: Unable to find column \"affinity_type\"";
    }

    const size_t match_col = get_col(header, "affinity_match");

    if(match_col == COLUMN_NOT_FOUND){
        throw __FILE__ ":parse_affinity: Unable to find column \"affinity_match\"";
    }

    const size_t structure_col = get_col(header, "structure_format");

    if(structure_col == COLUMN_NOT_FOUND){
        throw __FILE__ ":parse_affinity: Unable to find column \"structure_format\"";
    }

    const size_t year_col = get_col(header, "year");

    if(year_col == COLUMN_NOT_FOUND){
        throw __FILE__ ":parse_affinity: Unable to find column \"year\"";
    }

    const size_t structure_ph_col = get_col(header, "structure_pH");
    const size_t structure_temperature_col = get_col(header, "structure_Temperature");

    try{
        while( getline(fin, line) ){

            ++line_number;

            // Skip comments
            const string::size_type loc = line.find('#');

            if(loc != string::npos){
                line = line.substr(0, loc);
            }

            // Make sure we have a comma in the line we just read
            if(line.find(',') == string::npos){
                continue;
            }

            const vector<string> data = split(line, ',');

            if(data.size() != header.size()){
                throw __FILE__ ":parse_affinity: Did not read the expected number of columns";
            }

            string structure_ph;
            string structure_temperature;

            if(structure_ph_col != COLUMN_NOT_FOUND){
                structure_ph = data[structure_ph_col];
            }

            if(structure_temperature_col != COLUMN_NOT_FOUND){
                structure_temperature = data[structure_temperature_col];
            }

            // DEBUG
            //if( (structure_ph == "") || (structure_temperature == "") ){
            //    continue;
            //}

            ret.push_back( make_pair( data[accession_col], 
                Affinity(data[affinity_type_col],
                        data[match_col],
                        data[structure_col],
                        atoi(data[year_col].c_str()),
                        structure_ph, structure_temperature,
                        atof( data[affinity_col].c_str() ) ) ) );
        }
    }
    catch(const char* m_error){
        cerr << "Caught the error: " << m_error << " on line number " << line_number << endl;
        throw __FILE__ ":parse_affinity: Unable to parse affinity file";
    }
    catch(...){
        cerr << "Caught an unhandled error on line number " << line_number << endl;
        throw __FILE__ ":parse_affinity: Unable to parse affinity file";
    }

    // Sort the affinity data by accession for fast lookup
    sort(ret.begin(), ret.end());

    return ret;
}

vector< pair<string, Affinity> >::const_iterator affinity_lookup(const vector< pair<string, Affinity> > &m_data, const string &m_pdb,
    const char *m_file_name, const size_t &m_line_number)
{
    vector< pair<string, Affinity> >::const_iterator ret = lower_bound(m_data.begin(), m_data.end(), m_pdb, find_by_accession());

    if( ( ret == m_data.end() ) || (ret->first != m_pdb) ){
        cerr << "Unable to find accession \"" << m_pdb << "\" at line " << m_line_number << " in file " << m_file_name << endl;
        throw __FILE__ ":affinity_lookup: Unable to lookup PDB accession";
    }

    return ret;
}

float pearson_r(const deque<ClusterPredictions> &m_clusters, const bool m_use_weights)
{
    if(m_use_weights){
        return weighted_pearson_r(m_clusters);
    }

    float norm = 0.0;
    float ave_x = 0.0;
    float ave_y = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_clusters.begin();i != m_clusters.end();++i){

        if( i->empty() ){
            throw __FILE__ ":pearson_r: Empty clusters are not allowed!";
        }

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            norm += 1.0;
            ave_x += j->first;
            ave_y += j->second;
        }
    }

    ave_x /= norm;
    ave_y /= norm;

    float ave_xy = 0.0;
    float ave_xx = 0.0;
    float ave_yy = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_clusters.begin();i != m_clusters.end();++i){

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            const float delta_x = j->first - ave_x;
            const float delta_y = j->second - ave_y;

            ave_xy += delta_x*delta_y;
            ave_xx += delta_x*delta_x;
            ave_yy += delta_y*delta_y;
        }
    }

    return ave_xy/sqrt(ave_xx*ave_yy);
}

// From https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Weighted_correlation_coefficient
float weighted_pearson_r(const deque<ClusterPredictions> &m_clusters)
{
    float total_weight = 0.0;
    float ave_x = 0.0;
    float ave_y = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_clusters.begin();i != m_clusters.end();++i){

        if( i->empty() ){
            throw __FILE__ ":weighted_pearson_r: Empty clusters are not allowed!";
        }

        const float weight = 1.0/i->size();

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            total_weight += weight;

            ave_x += j->first*weight;
            ave_y += j->second*weight;
        }
    }

    ave_x /= total_weight;
    ave_y /= total_weight;

    float ave_xy = 0.0;
    float ave_xx = 0.0;
    float ave_yy = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_clusters.begin();i != m_clusters.end();++i){

        const float weight = 1.0/i->size();

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            const float delta_x = j->first - ave_x;
            const float delta_y = j->second - ave_y;

            ave_xy += delta_x*delta_y*weight;
            ave_xx += delta_x*delta_x*weight;
            ave_yy += delta_y*delta_y*weight;
        }
    }

    return ave_xy/sqrt(ave_xx*ave_yy);
}

// Compute the two-tailed p-value by permutation testing
float pearson_pvalue(float m_r /*copy*/, const deque<ClusterPredictions> &m_clusters, 
    const size_t &m_num_permute, struct drand48_data *m_rand_ptr, const bool m_use_weights)
{
    if(m_num_permute <= 1){
        throw __FILE__ ":weighted_pearson_pvalue: Invalid number of permutations";
    }

    const size_t num_cluster = m_clusters.size();

    if(num_cluster <= 1){
        throw __FILE__ ":weighted_pearson_pvalue: Must have at least two clusters to permute";
    }

    // For a two-tailed p-value, we want the absolute value of the observed correlation
    m_r = fabs(m_r);

    size_t num_more_correlated = 0;

    // DEBUG
    //cerr << "num_cluster = " << num_cluster << endl;
    //cerr << "m_num_permute = " << m_num_permute << endl;

    vector<size_t> index(num_cluster);

    for(size_t i = 0;i < num_cluster;++i){
        index[i] = i;
    }

    for(size_t i = 0;i < m_num_permute;++i){

        // DEBUG
        //cerr << "i = " << i << endl;

        randomize(index.begin(), index.end(), m_rand_ptr);

        deque<ClusterPredictions> clusters(num_cluster);

        for(size_t j = 0;j < num_cluster;++j){

            int min_size = min( m_clusters[ index[j] ].size(), m_clusters[j].size() );
            int max_size = max( m_clusters[ index[j] ].size(), m_clusters[j].size() );

            long int cluster_size;
            
            lrand48_r(m_rand_ptr, &cluster_size);

            cluster_size = min_size + cluster_size%(max_size - min_size + 1);

            clusters[j].resize(cluster_size);

            // DEBUG
            //cerr << "cluster_size = " << cluster_size << endl;
            //cerr << "num_cluster = " << num_cluster << endl;
            //cerr << "index[" << i << "] = " << index[i] << endl;

            for(long int k = 0;k < cluster_size;++k){

                // Combine the x-values from cluster index[i] with the
                // y-values from cluster i. Sample from each set of cluster values
                // with replacement. 
                long int x_index;
                long int y_index;

                lrand48_r(m_rand_ptr, &x_index);
                lrand48_r(m_rand_ptr, &y_index);

                x_index = x_index%m_clusters[ index[j] ].size();
                y_index = y_index%m_clusters[j].size();
                
                clusters[j][k] = make_pair(m_clusters[ index[j] ][x_index].first, m_clusters[j][y_index].second);

                // DEBUG
                //cerr << "(" << clusters[j][k].first << ',' << clusters[j][k].second << ")" << endl;
            }
        }

        const float r = pearson_r(clusters, m_use_weights);

        // DEBUG
        //cerr << "r = " << r << endl;

        // Two-tailed p-value
        if( (r <= -m_r) || (r >= m_r) ){
            ++num_more_correlated;
        }

        // DEBUG
        //cerr << "num_more_correlated = " << num_more_correlated << endl;
    }

    return ( (float)num_more_correlated )/m_num_permute;
}

vector< pair<string, Affinity> > random_affinity(const vector< pair< float, vector<Cluster> > > &m_clusters, 
    struct drand48_data *m_rand_ptr)
{
    unordered_map<string, Affinity> affinity; // <-- individual, per-complex affinity values
    unordered_set<string> assigned_clusters; // <-- cluster that have already been assigned affinity values

    for(vector< pair<float, vector<Cluster> > >::const_iterator i = m_clusters.begin();i != m_clusters.end();++i){

        // The random number amplitude could be constant or a function of cluster size, i.e.:
        //  weight = i->size()
        //  weight = 1.0
        //  weight = 1.0/i->size()
        //  weight = ??
        // Could try fitting the amplitude to obtain the best match to the real affinity data
        const float weight = 1.0;

        for(vector<Cluster>::const_iterator j = i->second.begin();j != i->second.end();++j){
            
            const string id = cluster_id(*j);

            if( assigned_clusters.find(id) != assigned_clusters.end() ){

                // We have already assigned an affinity to this cluster
                continue;
            }

            assigned_clusters.insert(id);

            // Assign a random number to each cluster
            double r = 0.0;

            drand48_r(m_rand_ptr, &r);

            r *= weight;

            for(Cluster::const_iterator k = j->begin();k != j->end();++k){

                // Each member of this cluster adds the same, cluster-specific random number
                affinity[*k].value += r;
            }
        }
    }

    vector< pair<string, Affinity> > ret;

    ret.reserve( affinity.size() );

    for(unordered_map<string, Affinity>::const_iterator i = affinity.begin();i != affinity.end();++i){

        // DEBUG
        //cout << i->first << '\t' << i->second << endl;

        //double r = 0.0;

        //drand48_r(m_rand_ptr, &r);

        ret.push_back( make_pair(i->first, i->second) );
        //ret.push_back( make_pair(i->first, r) );
    }

    // The affinity data structure must be sorted by id to enable the use of lower_bound
    sort(ret.begin(), ret.end());

    return ret;
}

// Convert the set of PDB accessions into a unique cluster id
string cluster_id(Cluster m_cluster /*copy*/)
{
    sort(m_cluster.begin(), m_cluster.end());

    string ret;

    for(Cluster::const_iterator i = m_cluster.begin();i != m_cluster.end();++i){
        
        // Create a unique string by concatinating the indvidual strings
        ret += *i + ':';
    }

    return ret;
}

void remove_atoms(vector<AminoAcidData> &m_chain, const unsigned int &m_atom_type)
{
    for(vector<AminoAcidData>::iterator i = m_chain.begin();i != m_chain.end();++i){

        vector<AtomData> new_atoms;

        new_atoms.reserve( i->atoms.size() );
        
        for(vector<AtomData>::const_iterator j = i->atoms.begin();j != i->atoms.end();++j){
            if(j->type != m_atom_type){
                new_atoms.push_back(*j);
            }
        }

        i->atoms.swap(new_atoms);
    }
}

SymmetricMatrix<float> surrogate_affinity_param(const unordered_map<string, unsigned int> &m_atom_table, struct drand48_data *m_rand_ptr)
{
    // Collect the atom types
    unsigned int max_atom_type = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        max_atom_type = max(max_atom_type, i->second);
    }

    ++max_atom_type; // Increment by one to serve as the number of atom types

    // The weighting of the pairwise distances depends on the atom types
    SymmetricMatrix<float> weight(max_atom_type);

    // Randomly assign atom weights
    for(unsigned int i = 0;i < max_atom_type;++i){
        for(unsigned int j = i;j < max_atom_type;++j){

             // Assign a random number to each cluster
            double r = 0.0;

            drand48_r(m_rand_ptr, &r);

            weight(i, j) = r;
        }
    }

    return weight;
}

// Compute a pair-wise distance based scoring function to serve as an affinity surrogate and positive control
Affinity surrogate_affinity_score(const PDBComplex &m_mol, const SymmetricMatrix<float> &m_param)
{
    Affinity ret;

    ret.affinity_type = Affinity::AFFINITY_SCORE;
    ret.measurement_type = Affinity::MEASURE_EQUAL;
    ret.structure_type = Affinity::STRUCTURE_PREDICTED;
    ret.structure_resolution = 0.0;

    // Accumulate the score as a double
    double score = 0.0;

    // Only include pair-wise interactions *between* proteins and ignore interactions *within* proteins.
    for(vector<AminoAcidData>::const_iterator i = m_mol.first.begin();i != m_mol.first.end();++i){
        for(vector<AtomData>::const_iterator j = i->atoms.begin();j != i->atoms.end();++j){
            for(vector<AminoAcidData>::const_iterator k = m_mol.second.begin();k != m_mol.second.end();++k){
                for(vector<AtomData>::const_iterator m = k->atoms.begin();m != k->atoms.end();++m){
                    
                    const float d2 = j->distance2(*m);

                    // DEBUG
                    //if(d2 > 10.0*10.0){
                    //    continue;
                    //}
                    
                    score += m_param(j->type, m->type)/d2;
                    //score += m_param(j->type, m->type)*d2;
                }
            }
        }
    }

    ret.value = score;
    
    return ret;
}

float feature_similarity(const Features &m_a, const Features &m_b)
{
    const size_t len = m_a.size();

    if( len != m_b.size() ){
        throw __FILE__ ":feature_similarity: Can't compare feature vectors with different lengths";
    }

    float ab = 0.0;
    float aa = 0.0;
    float bb = 0.0;

    for(size_t i = 0;i < len;++i){

        ab += m_a[i] * m_b[i];
        aa += m_a[i] * m_a[i];
        bb += m_b[i] * m_b[i];
    }

    if(ab == 0.0){
        return 0.0;
    }

    return ab/sqrt(aa*bb);
}

float compute_rmse(const deque<ClusterPredictions> &m_clusters, const bool m_use_weights)
{
    float ret = 0.0;
    float norm = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_clusters.begin();i != m_clusters.end();++i){
        
        if( i->empty() ){
            continue;
        }

        const float w = m_use_weights ? 1.0/i->size() : 1.0;

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            const float delta = j->first - j->second;

            ret += w*delta*delta;
        }

        // Each cluster receives a weight of 1.0 and every cluster element receives a weight of 1/N (where N is the
        // number of elements in the given cluster)
        norm += m_use_weights ? 1.0 : i->size();
    }

    if(norm > 0.0){
        ret = sqrt(ret/norm);
    }

    return ret;
}

float compute_R2(const deque<ClusterPredictions> &m_predict, const bool m_use_weights)
{
    float ave = 0.0;
    float norm = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_predict.begin();i != m_predict.end();++i){

        if( i->empty() ){
            continue;
        }

        const float w = m_use_weights ? 1.0/i->size() : 1.0;

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            // ClusterPredictions -> (observed, predicted)
            ave += w*j->first;
        }

        norm += m_use_weights ? 1.0 : i->size();
    }

    if(norm > 0.0){
        ave /= norm;
    }

    float pred_ss = 0.0;
    float ave_ss = 0.0;

    for(deque<ClusterPredictions>::const_iterator i = m_predict.begin();i != m_predict.end();++i){

        if( i->empty() ){
            continue;
        }

        const float w = m_use_weights ? 1.0/i->size() : 1.0;

        for(ClusterPredictions::const_iterator j = i->begin();j != i->end();++j){

            // ClusterPredictions -> (observed, predicted)
            float delta = j->first - j->second;
            pred_ss += w*delta*delta;

            delta = j->first - ave;
            ave_ss += w*delta*delta;
        }
    }

    if(ave_ss <= 0.0){
        throw __FILE__ ":compute_R2: Unable to compute R2 (zero in the denominator)";
    }

    return 1.0 - pred_ss/ave_ss;
}

// Compute R2 from m_predict relative to the prediction in m_reference
float compute_R2(const deque<ClusterPredictions> &m_predict, const deque<ClusterPredictions> &m_reference, 
    const bool m_use_weights)
{
    if(m_predict.size() != m_reference.size()){
        throw __FILE__ ":compute_R2: Can't compute R2 when |m_predict| != |m_reference|";
    }

    const size_t num_clusters = m_predict.size();

    float pred_ss = 0.0;
    float reference_ss = 0.0;

    for(size_t i = 0;i < num_clusters;++i){

        if(m_predict[i].size() != m_reference[i].size()){
            throw __FILE__ ":compute_R2: Can't compute R2 when |m_predict[i]| != |m_reference[i]|";
        }

        const size_t cluster_size = m_predict[i].size();

        if(cluster_size == 0){
            continue;
        }

        const float w = m_use_weights ? 1.0/cluster_size : 1.0;

        for(size_t j = 0;j < cluster_size;++j){

            // Make sure that both the prediction and reference agree on the same observed value
            // ClusterPredictions -> (observed, predicted)
            if(m_predict[i][j].first != m_reference[i][j].first){
                throw __FILE__ ":compute_R2: Observed values must be identitical for prediction and reference";
            }

            const float obs = m_predict[i][j].first;

            float delta = obs - m_predict[i][j].second;
            pred_ss += w*delta*delta;

            delta = obs - m_reference[i][j].second;
            reference_ss += w*delta*delta;
        }
    }

    if(reference_ss <= 0.0){
        throw __FILE__ ":compute_R2: Unable to compute R2 (zero in the denominator)";
    }

    return 1.0 - pred_ss/reference_ss;
}