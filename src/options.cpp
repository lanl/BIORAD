#include "options.h"
#include "parse_util.h"
#include <iostream>
#include <limits>
#include <algorithm>
#include <getopt.h>
#include <sys/stat.h>

using namespace std;

#define     DEFAULT_NUM_FOLD                5
#define     DEFAULT_SEED                    0 // Triggers a time-based random number seed
#define     DEFAULT_NUM_TRIAL               1

// Note that the Random Forest implementation from scikitlearn:
//  https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
// Randomly samples *features* (not data) to generate tree diversity!

#define     DEFAULT_FOREST_SIZE             100
#define     DEFAULT_FOREST_LEAF             1
#define     DEFAULT_FOREST_BAG              0.50     

// Number of permutations for computing correlation p-values
#define     DEFAULT_NUM_PERMUTE             1000

#define     DEFAULT_HIST_MIN_DIST           1.0
#define     DEFAULT_HIST_MAX_DIST           10.0
#define     DEFAULT_HIST_NUM_BIN            9 // (10.0 - 1.0)/9 = 1 angstrom bins

#define     DEFAULT_GAUSS_MAX_DIST          10.0
#define     DEFAULT_GAUSS_NUM_SAMPLE        10
#define     DEFAULT_GAUSS_NUM_SAMPLE_ITER   4
#define     DEFAULT_GAUSS_SAMPLE_FRACTION   0.25

#define     DEFAULT_SELF_MAX_COUNT          3
#define     DEFAULT_HETERO_MAX_COUNT        10
#define     DEFAULT_SELF_DISTANCE           2
#define     DEFAULT_HETERO_DISTANCE         4

// Sort by string length in *descending* order
struct sort_by_string_length
{
    inline bool operator()(const string &m_a, const string &m_b) const
    {
        return m_b.size() < m_a.size();
    };
};

bool file_exists(const std::string &m_filename);

void Options::load(int argc, char* argv[])
{
    const char* options = "o:t:p:h";
    int config_opt = 0;
    int long_index = 0;

    struct option long_opts[] = {
        {"affinity", true, &config_opt, 0},
        {"fold", true, &config_opt, 1},
        {"seed", true, &config_opt, 2},
        {"pdb", true, &config_opt, 3},
        {"cluster", true, &config_opt, 5},
        {"cluster.align", false, &config_opt, 6},
        {"cluster.composition", false, &config_opt, 7},
        {"pdb.ignore-hetatom", false, &config_opt, 8},
        {"pdb.force-dimer", false, &config_opt, 9},
        {"forest.size", true, &config_opt, 10},
        {"forest.leaf", true, &config_opt, 11},
        {"forest.bag.data", true, &config_opt, 12},
        {"greedy", false, &config_opt, 13},
        {"hist.min", true, &config_opt, 14},
        {"hist.max", true, &config_opt, 15},
        {"hist.bins", true, &config_opt, 16},
        {"hist.dir", true, &config_opt, 17},
        {"rosetta", true, &config_opt, 22}, // Append Rosetta scores to structure-based features
        {"affinity.random", false, &config_opt, 23},
        {"forest.bag.feature", true, &config_opt, 24},
        {"ROSETTA", true, &config_opt, 25}, // Only use Rosetta scores as features
        {"score", true, &config_opt, 26}, // Append scores to structure-based features
        {"SCORE", true, &config_opt, 27}, // Only use scores as features
        {"ignore-kd", false, &config_opt, 28},
        {"ignore-ki", false, &config_opt, 29},
        {"ignore-ic50", false, &config_opt, 30},
        {"self.max-count", true, &config_opt, 31},
        {"self.dist", true, &config_opt, 32},
        {"hetero.max-count", true, &config_opt, 33},
        {"hetero.dist", true, &config_opt, 34},
        {"affinity.surrogate", false, &config_opt, 35},
        {"linear", false, &config_opt, 36},
        {"forest.weight", false, &config_opt, 37},
        {"nn", true, &config_opt, 38},
        {"affinity.year.min", true, &config_opt, 39},
        {"affinity.year.max", true, &config_opt, 40},
        {"affinity.min", true, &config_opt, 41},
        {"affinity.max", true, &config_opt, 42},
        {"pdb.ignore-charge", false, &config_opt, 43},
        {"weighted-stats", false, &config_opt, 44},
        {"no-weighted-stats", false, &config_opt, 45},
        {0,0,0,0} // Terminate options list
    };

    int opt_code;
    opterr = 0;

    quit = (argc == 1);
    
    num_fold = DEFAULT_NUM_FOLD;
    seed = DEFAULT_SEED;
    cluster_by = CLUSTER_BY_FILE;
    num_trial = DEFAULT_NUM_TRIAL;
    ignore_hetatom = false;
    ignore_charge = false;
    ignore_kd = false;
    ignore_ki = false;
    ignore_ic50 = false;
    use_linear_regression = false;
    force_dimer = false;
    random_affinity = false;
    surrogate_affinity = false;
    affinity_min_year = std::numeric_limits<unsigned int>::lowest();
    affinity_max_year = std::numeric_limits<unsigned int>::max();

    affinity_min = 0.0;
    affinity_max = 1.0;

    append_rosetta_scores = true; // By default, append Rosetta scores to the structure-based feature vector
    append_scores = true; // By default, append scores to the structure-based feature vector
    forest_size = DEFAULT_FOREST_SIZE;
    forest_leaf = DEFAULT_FOREST_LEAF;
    forest_data_bag = DEFAULT_FOREST_BAG;
    forest_feature_bag = 1.0; // By default, use all features
    use_forest_weight = false; // Don't weight training data

    use_weighted_stats = false; // Don't report cluster weighted statistics

    k_nearest_neighbors = 0; // By default, nearest neighbor regression is disabled

    feature_self_max_count = DEFAULT_SELF_MAX_COUNT;
    feature_hetero_max_count = DEFAULT_HETERO_MAX_COUNT;
    feature_self_distance = DEFAULT_SELF_DISTANCE;
    feature_hetero_distance = DEFAULT_HETERO_DISTANCE;

    num_permute = DEFAULT_NUM_PERMUTE;

    feature_type = FEATURE_DEFAULT;
    
    feature_hist_min = DEFAULT_HIST_MIN_DIST;
    feature_hist_max = DEFAULT_HIST_MAX_DIST;
    feature_hist_num_bin = DEFAULT_HIST_NUM_BIN;

    while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

        switch( opt_code ){
            case 0:
                
                if(config_opt == 0){ // affinity
            
                    affinity_file = optarg;
                    break;
                }

                if(config_opt == 1){ // fold
            
                    num_fold = stoul(optarg);
                    break;
                }

                if(config_opt == 2){ // seed
            
                    seed = stoul(optarg);
                    break;
                }

                if(config_opt == 3){ // pdb
            
                    pdb_dir = optarg;
                    break;
                }

                if(config_opt == 5){ // cluster

                    cluster_file = optarg;
                    break;
                }

                if(config_opt == 6){ // cluster.align

                    cluster_by = CLUSTER_BY_SEQ_ALIGN;
                    break;
                }

                if(config_opt == 7){ // cluster.composition

                    cluster_by = CLUSTER_BY_COMPOSITION;
                    break;
                }

                if(config_opt == 8){ // pdb.ignore-hetatom

                    ignore_hetatom = true;
                    break;
                }

                if(config_opt == 9){ // pdb.force-dimer

                    force_dimer = true;
                    break;
                }

                if(config_opt == 10){ // forest.size

                    forest_size = stoul(optarg);
                    break;
                }

                if(config_opt == 11){ // forest.leaf

                    forest_leaf = stoul(optarg);
                    break;
                }

                if(config_opt == 12){ // forest.bag.data

                    forest_data_bag = atof(optarg);
                    break;
                }

                if(config_opt == 13){ // greedy

                    feature_type = FEATURE_GREEDY;
                    break;
                }

                if(config_opt == 14){ // hist.min

                    feature_hist_min = atof(optarg);
                    break;
                }

                if(config_opt == 15){ // hist.max

                    feature_hist_max = atof(optarg);
                    break;
                }

                if(config_opt == 16){ // hist.bins

                    feature_hist_num_bin = stoul(optarg);
                    break;
                }

                if(config_opt == 17){ // hist.dir

                    hist_dir = optarg;
                    break;
                }

                if(config_opt == 22){ // rosetta

                    rosetta_score_file = optarg;
                    break;
                }

                if(config_opt == 23){ // affinity.random

                    random_affinity = true;
                    break;
                }

                if(config_opt == 24){ // forest.bag.feature

                    forest_feature_bag = atof(optarg);
                    break;
                }

                if(config_opt == 25){ // ROSETTA

                    rosetta_score_file = optarg;
                    append_rosetta_scores = false;
                    break;
                }

                if(config_opt == 26){ // score

                    score_file = optarg;
                    append_scores = true;
                    break;
                }

                if(config_opt == 27){ // SCORE

                    score_file = optarg;
                    append_scores = false;
                    break;
                }

                if(config_opt == 28){ // ignore-kd

                    ignore_kd = true;
                    break;
                }

                if(config_opt == 29){ // ignore-ki

                    ignore_ki = true;
                    break;
                }

                if(config_opt == 30){ // ignore-ic50

                    ignore_ic50 = true;
                    break;
                }

                if(config_opt == 31){ // self.max-count

                    feature_self_max_count = atoi(optarg);
                    break;
                }

                if(config_opt == 32){ // self.dist

                    feature_self_distance = atof(optarg);
                    break;
                }

                if(config_opt == 33){ // hetero.max-count

                    feature_hetero_max_count = atoi(optarg);
                    break;
                }

                if(config_opt == 34){ // hetero.dist

                    feature_hetero_distance = atof(optarg);
                    break;
                }

                if(config_opt == 35){ // affinity.surrogate

                    surrogate_affinity = true;
                    break;
                }

                if(config_opt == 36){ // linear

                    use_linear_regression = true;
                    break;
                }

                if(config_opt == 37){ // forest.weight

                    use_forest_weight = true;
                    break;
                }

                if(config_opt == 38){ // nn

                    k_nearest_neighbors = abs( atoi(optarg) );
                    break;
                }

                if(config_opt == 39){ // affinity.year.min

                    affinity_min_year = abs( atoi(optarg) );
                    break;
                }

                if(config_opt == 40){ // affinity.year.max

                    affinity_max_year = abs( atoi(optarg) );
                    break;
                }

                if(config_opt == 41){ // affinity.min

                    affinity_min = atof(optarg);
                    break;
                }

                if(config_opt == 42){ // affinity.max

                    affinity_max = atof(optarg);
                    break;
                }

                if(config_opt == 43){ // pdb.ignore-charge

                    ignore_charge = true;
                    break;
                }

                if(config_opt == 44){ // weighted-stats

                    use_weighted_stats = true;
                    break;
                }

                if(config_opt == 45){ // no-weighted-stats

                    use_weighted_stats = false;
                    break;
                }

                cerr << "Unknown command line flag!" << endl;
                quit = true;
                return;
            case 'o':
                output_file = optarg;
                break;
            case 't':
                num_trial = stoul(optarg);
                break;
            case 'p':
                num_permute = stoul(optarg);
                break;
            case 'h':
                quit = true;
                break;
            case '?':
            default:
                cerr << '\"' << argv[optind - 1] << "\" is not a valid option!" << endl;
                quit = true;
                return;
        };
    }

    if(quit){

        cerr << "Usage for biorad (v. " << BIORAD_VERSION << ")" << endl;
        cerr << "\t--affinity <CSV affinity data>" << endl;
        cerr << "\t[--weighted-stats (report cluster size-weighted statistics)]" << endl;
        cerr << "\t[--no-weighted-stats (report unweighted statistics)] (default)" << endl;
        cerr << "\t[--affinity.random (use random affinity values as a negative control)]" << endl;
        cerr << "\t[--affinity.surrogate (use a score function as an affinity surrogate for a positive control)]" << endl;
        cerr << "\t[--affinity.year.min <earliest year to include affinity data from>]" << endl;
        cerr << "\t[--affinity.year.max <latest year to include affinity data from>]" << endl;
        cerr << "\t[--affinity.min <smallest allowed affinity (in M)>]" << endl;
        cerr << "\t[--affinity.max <largest allowed affinity (in M)>]" << endl;
        cerr << "\tRosetta features (choose at most one):" << endl;
        cerr << "\t\t[--rosetta <score file>] (optional file of Rosetta scores to *append* to feature vector)" << endl;
        cerr << "\t\t[--ROSETTA <score file>] (optional file of Rosetta scores to use as the *only* feature vector)" << endl;
        cerr << "\tCSV features (choose at most one):" << endl;
        cerr << "\t\t[--score <score file>] (optional CSV file to *append* to feature vector)" << endl;
        cerr << "\t\t[--SCORE <score file>] (optional CSV file to use as the *only* feature vector)" << endl;
        cerr << "\t--pdb <directory of PDB files>" << endl;
        cerr << "\t[--pdb.ignore-hetatom] (ignore hetatoms connected to protein atoms)" << endl;
        cerr << "\t[--pdb.ignore-charge] (ignore atom charge states when labeling atom type)" << endl;
        cerr << "\t[--pdb.force-dimer] (Include all two-chain structures, ignoring the provided structure/unit information)" << endl;
        cerr << "\t[-o <output file>]" << endl;
        cerr << "\t--cluster <file of PDB clusters>" << endl;
        cerr << "\t[--cluster.align] (cluster PDBs by pairwise sequence alignment and save to cluster file)" << endl;
        cerr << "\t[--cluster.composition] (cluster PDBs by sequence composition and save to cluster file)" << endl;
        cerr << "\t[--fold <number of cross validation folds>] (default is " << DEFAULT_NUM_FOLD << ")" << endl;
        cerr << "\t[-t <number of cross validation trials>]" << endl;
        cerr << "\t[--seed <random number seed>] (default is time-based)" << endl;
        cerr << "\t[--forest.size <number of trees in the random forest>] (default is " << DEFAULT_FOREST_SIZE << ")" << endl;
        cerr << "\t[--forest.leaf <minimum number of data point in a leaf>] (default is " << DEFAULT_FOREST_LEAF << ")" << endl;
        cerr << "\t[--forest.bag.data <fraction of data sampled to generate tree diversity>] (default is " << DEFAULT_FOREST_BAG << ")" << endl;
        cerr << "\t[--forest.bag.feature <fraction of features sampled to generate tree diversity>] (default is 1.0)" << endl;
        cerr << "\t[--forest.weight (apply inverse cluster-size training data weights)]" << endl;
        cerr << "\t[-p <number of permutations for computing correlation p-values>] (default is " << DEFAULT_NUM_PERMUTE << ")" << endl;
        cerr << "\t[--greedy] (greedy feature selection)" << endl;
        cerr << "\t[--hist.min <minimum pairwise distance in angstroms>] (default is " << DEFAULT_HIST_MIN_DIST<< ")" << endl;
        cerr << "\t[--hist.max <maximum pairwise distance in angstroms>] (default is " << DEFAULT_HIST_MAX_DIST<< ")" << endl;
        cerr << "\t[--hist.bins <maximum pairwise distance in angstroms>] (default is " << DEFAULT_HIST_NUM_BIN<< ")" << endl;
        cerr << "\t[--hist.dir <directory for writing CSV histogram data>]" << endl;
        cerr << "\t[--self.max-count <maximum self atom count>] (default is " << DEFAULT_SELF_MAX_COUNT << ")" << endl;
        cerr << "\t[--self.dist <self adjacency distance>] (default is " << DEFAULT_SELF_DISTANCE << ")" << endl;
        cerr << "\t[--hetero.max-count <maximum hetero atom count>] (default is " << DEFAULT_HETERO_MAX_COUNT << ")" << endl;
        cerr << "\t[--hetero.dist <hetero adjacency distance>] (default is " << DEFAULT_HETERO_DISTANCE << ")" << endl;
        cerr << "\t[--ignore-kd] (Do not include Kd-based affinity measurements)" << endl;
        cerr << "\t[--ignore-ki] (Do not include Ki-based affinity measurements)" << endl;
        cerr << "\t[--ignore-ic50] (Do not include IC50-based affinity measurements)" << endl;

        #ifdef USE_GSL
        cerr << "\t[--linear] (Include a linear-regression model)" << endl;
        #endif // USE_GSL

        cerr << "\t[--nn <number of nearest neighbors for k-nn regression>]" << endl;
        return;
    }

    // A seed value of 0 triggers the use of a time-based random number seed
    if(seed == 0){
        seed = time(NULL);
    }

    if(num_trial == 0){
        
        cerr << "Please specify at least one trial of cross-validation" << endl;
        quit = true;
        return;
    }

    if( pdb_dir.empty() ){

        cerr << "Please specify a directoy of PDB files" << endl;
        quit = true;
        return;
    }

    if( cluster_file.empty() ){

        cerr << "Please specify a cluster file (--cluster)" << endl;
        quit = true;
        return;
    }

    if(forest_size < 1){

        cerr << "Please specify at least one tree in the random forest (--forest.size)" << endl;
        quit = true;
        return;
    }

    if(forest_leaf < 1){

        cerr << "Please specify at least one data point per leaf in the random forest (--forest.leaf)" << endl;
        quit = true;
        return;
    }

    if( (forest_data_bag <= 0.0) || (forest_data_bag > 1.0) ){

        cerr << "Please specify a random forest data sampling fraction > 0 and <= 1 (--forest.bag.data)" << endl;
        quit = true;
        return;
    }

    if( (forest_feature_bag <= 0.0) || (forest_feature_bag > 1.0) ){

        cerr << "Please specify a random forest feature sampling fraction > 0 and <= 1 (--forest.bag.feature)" << endl;
        quit = true;
        return;
    }

    if(feature_type == FEATURE_NOT_SET){

        cerr << "Please specify a feature type for ML (--hist | --gauss)" << endl;
        quit = true;
        return;
    }
        
    if(feature_hist_min < 0.0){

        cerr << "Please specify a minimum histogram distance >= 0.0" << endl;
        quit = true;
        return;
    }

    if(feature_hist_min > feature_hist_max){

        cerr << "Please specify a maximum histogram distance >= minimum histogram distance" << endl;
        quit = true;
        return;
    }

    if(feature_hist_num_bin == 0){

        cerr << "Please specify a number of histogram bins >= 1" << endl;
        quit = true;
        return;
    }

    if(feature_self_max_count == 0){

        cerr << "Please specify a maximum number of self adjacent atoms >= 1" << endl;
        quit = true;
        return;
    }

    if(feature_self_distance <= 0.0){

        cerr << "Please specify a self adjacency distance > 0.0" << endl;
        quit = true;
        return;
    }

    if(feature_hetero_max_count == 0){

        cerr << "Please specify a maximum number of hetero adjacent atoms >= 1" << endl;
        quit = true;
        return;
    }

    if(feature_hetero_distance <= 0.0){

        cerr << "Please specify a hetero adjacency distance > 0.0" << endl;
        quit = true;
        return;
    }

    if( !random_affinity && affinity_file.empty() ){
        
        cerr << "Please specify a file of binding affinity values (--affinity)" << endl;
        quit = true;
        return;
    }

    if(!score_file.empty() && !rosetta_score_file.empty()){

        if(!append_rosetta_scores || !append_scores){

            cerr << "Rosetta and CSV scores cannot *both* be used *exclusively*" << endl;
            quit = true;
            return;
        }
    }

    // Make sure the we will not attempt to overwrite an existing cluster file
    if( file_exists(cluster_file) && (cluster_by != CLUSTER_BY_FILE) ){

        cerr << "Overwriting an existing cluster file is not allowed" << endl;
        quit = true;
        return;
    }

    if(affinity_min_year > affinity_max_year){

        cerr << "Minimum affinity year must be less than, or equal, to the maximum affinity year" << endl;
        quit = true;
        return;
    }

    if(affinity_min > affinity_max){

        cerr << "Minimum affinity must be less than, or equal, to the maximum affinity (in M, not log(M))" << endl;
        quit = true;
        return;
    }
}

bool file_exists(const std::string &m_filename)
{
        struct stat file_info;

        if(stat(m_filename.c_str(), &file_info) != 0){
                return false;
        }

        // The path exists, make sure it is a file
        return S_ISREG(file_info.st_mode);
}
