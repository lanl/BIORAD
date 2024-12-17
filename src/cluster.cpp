#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <deque>
#include <unordered_set>

#include <math.h>

#include "cluster.h"
#include "parse_util.h"
#include "symmetric_matrix.h"
#include "seq_overlap.h"

using namespace std;

#define     INDEX(I, J, N)  ((I)*(N) + (J))

Cluster parse_cluster(const string &m_str);
vector< pair<float/*min cluster distance*/, vector<Cluster> > > cluster_data(const vector<string> &m_id, 
    const SymmetricMatrix<float> &m_dist, const size_t &m_min_num_cluster);
float cluster_distance(const vector<string> &m_id, const SymmetricMatrix<float> &m_dist, 
    const Cluster &m_a, const Cluster &m_b);
int perfect_BLOSUM62_score(const string &m_seq);

void cluster_by_composition(const string &m_filename, const unordered_map<string, PDBComplex> &m_data,
    const Options &m_opt)
{
    cerr << "Clustering protein complexes by amino acid composition" << endl;

    const size_t num_data = m_data.size();

    // Precompute the pairwise distances
    SymmetricMatrix<float> dist(num_data);

    // Compute the amino acid composition of each complex
    vector< vector<unsigned int> > features;

    // Record the PDB accessions in a fixed order that corresponds to the rows and columns in
    // the distance matrix
    vector<string> id;

    features.reserve(num_data);
    id.reserve(num_data);

    for(unordered_map<string, PDBComplex>::const_iterator i = m_data.begin();i != m_data.end();++i){
        id.push_back(i->first);
    }

    // Allow fast id lookups with a well defined ordering
    sort(id.begin(), id.end());

    for(vector<string>::const_iterator i = id.begin();i != id.end();++i){

        unordered_map<string, PDBComplex>::const_iterator iter = m_data.find(*i);

        if(iter == m_data.end()){
            throw __FILE__ ":cluster_by_composition: Unable to lookup protein complex!";
        }

        vector<unsigned int> f(AminoAcidData::NUM_REAL_AMINO_ACID_TYPES);

        // Count the amino acids in the first protein ...
        for(vector<AminoAcidData>::const_iterator j = iter->second.first.begin();j != iter->second.first.end();++j){

            if(j->type < AminoAcidData::NUM_REAL_AMINO_ACID_TYPES){
                ++f[j->type];
            }
        }

        // ... and the second protein
        for(vector<AminoAcidData>::const_iterator j = iter->second.second.begin();j != iter->second.second.end();++j){

            if(j->type < AminoAcidData::NUM_REAL_AMINO_ACID_TYPES){
                ++f[j->type];
            }
        }

        features.push_back(f);
    }

    cerr << "\tComputing pairwise distances" << endl;

    // Compute the pairwise distances
    for(size_t i = 0;i < num_data;++i){

        for(size_t j = i + 1;j < num_data;++j){

            float d = 0.0;

            for(size_t k = 0;k < AminoAcidData::NUM_REAL_AMINO_ACID_TYPES;++k){

                const float delta = float(features[i][k]) - float(features[j][k]);

                d += delta*delta;
            }

            dist(i, j) = sqrt(d);
        }
    }

    cerr << "\tPerforming single-linkage clustering" << endl;

    write_clusters( m_filename, cluster_data(id, dist, m_opt.num_fold), m_opt);
}

void cluster_by_seq_align(const string &m_filename, const unordered_map<string, PDBComplex> &m_data,
    const Options &m_opt)
{
    cerr << "Clustering protein complexes by sequence alignment" << endl;

    const size_t num_data = m_data.size();

    // Precompute the pairwise distances
    SymmetricMatrix<float> dist(num_data);

    // Record the PDB accessions in a fixed order that corresponds to the rows and columns in
    // the distance matrix
    vector<string> id;

    id.reserve(num_data);

    for(unordered_map<string, PDBComplex>::const_iterator i = m_data.begin();i != m_data.end();++i){
        id.push_back(i->first);
    }

    // Allow fast id lookups with a well defined ordering
    sort(id.begin(), id.end());

    cerr << "\tComputing pairwise distances" << endl;

    // Compute the pairwise distances
    #pragma omp parallel
    {
        SO::SeqOverlap align1(SO::SeqOverlap::SmithWaterman, false /*protein alignment*/);
        SO::SeqOverlap align2(SO::SeqOverlap::SmithWaterman, false /*protein alignment*/);

        #pragma omp for
        for(size_t i = 0;i < num_data;++i){

            unordered_map<string, PDBComplex>::const_iterator iter_i = m_data.find(id[i]);

            if( iter_i == m_data.end() ){
                throw __FILE__ ":cluster_by_seq_align: Unable to look up protein complex (i)";
            }

            const string seq_i1 = protein_sequence(iter_i->second.first);
            const string seq_i2 = protein_sequence(iter_i->second.second);

            const SO::SO_Score perfect_i1 = perfect_BLOSUM62_score(seq_i1);
            const SO::SO_Score perfect_i2 = perfect_BLOSUM62_score(seq_i2);

            align1.set_query(seq_i1);
            align2.set_query(seq_i2);

            for(size_t j = i + 1;j < num_data;++j){

                unordered_map<string, PDBComplex>::const_iterator iter_j = m_data.find(id[j]);

                if( iter_j == m_data.end() ){
                    throw __FILE__ ":cluster_by_seq_align: Unable to look up protein complex (j)";
                }
                
                const string seq_j1 = protein_sequence(iter_j->second.first);
                const string seq_j2 = protein_sequence(iter_j->second.second);

                const SO::SO_Score perfect_j1 = perfect_BLOSUM62_score(seq_j1);
                const SO::SO_Score perfect_j2 = perfect_BLOSUM62_score(seq_j2);

                align1.set_target(seq_j1);

                align1.align(); // First protein in i vs first protein in j

                const SO::SO_Score s_11 = align1.score();
                // Old pairwise protein distance
                //const float max_possible_score_11 = min(perfect_i1, perfect_j1);
                //const float d_11 = float(max_possible_score_11 - s_11)/max_possible_score_11;
                // New pairwise protein distance from Li et al, Cell Systems 10, 308-322.e1-e11, 2020
                const float d_11 = 1.0 - float(s_11)/sqrt(perfect_i1*perfect_j1);

                align1.set_target(seq_j2);

                align1.align(); // First protein in i vs second protein in j

                const SO::SO_Score s_12 = align1.score();
                // Old pairwise protein distance
                //const float max_possible_score_12 = min(perfect_i1, perfect_j2);
                //const float d_12 = float(max_possible_score_12 - s_12)/max_possible_score_12;
                // New pairwise protein distance from Li et al, Cell Systems 10, 308-322.e1-e11, 2020
                const float d_12 = 1.0 - float(s_12)/sqrt(perfect_i1*perfect_j2);

                align2.set_target(seq_j1);

                align2.align(); // Second protein in i vs first protein in j

                const SO::SO_Score s_21 = align2.score();
                // Old pairwise protein distance
                //const float max_possible_score_21 = min(perfect_i2, perfect_j1);
                //const float d_21 = float(max_possible_score_21 - s_21)/max_possible_score_21;
                // New pairwise protein distance from Li et al, Cell Systems 10, 308-322.e1-e11, 2020
                const float d_21 = 1.0 - float(s_21)/sqrt(perfect_i2*perfect_j1);

                align2.set_target(seq_j2);

                align2.align(); // Second protein in i vs second protein in j

                const SO::SO_Score s_22 = align2.score();
                // Old pairwise protein distance
                //const float max_possible_score_22 = min(perfect_i2, perfect_j2);
                //const float d_22 = float(max_possible_score_22 - s_22)/max_possible_score_22;
                // New pairwise protein distance from Li et al, Cell Systems 10, 308-322.e1-e11, 2020
                const float d_22 = 1.0 - float(s_22)/sqrt(perfect_i2*perfect_j2);

                if(d_11 < 0.0){
                    throw __FILE__ ":cluster_by_seq_align: d_11 < 0.0";
                }

                if(d_12 < 0.0){
                    throw __FILE__ ":cluster_by_seq_align: d_12 < 0.0";
                }

                if(d_21 < 0.0){
                    throw __FILE__ ":cluster_by_seq_align: d_21 < 0.0";
                }

                if(d_22 < 0.0){
                    throw __FILE__ ":cluster_by_seq_align: d_22 < 0.0";
                }

                // Use the best average distance between the two possible protein pairs 
                dist(i, j) = 0.5*min(d_11 + d_22, d_12 + d_21);
            }
        }
    }

    cerr << "\tPerforming single-linkage clustering" << endl;

    write_clusters( m_filename, cluster_data(id, dist, m_opt.num_fold), m_opt);
}

// Use single linkage clustering to perform agglomerative clustering. The vector of ids, m_id, must
// be sorted in ascending order
vector< pair<float, vector<Cluster> > > cluster_data(const vector<string> &m_id, const SymmetricMatrix<float> &m_dist,
    const size_t &m_min_num_cluster)
{
    vector< pair<float/*minimum distance between clusters*/, vector<Cluster> > > ret;

    // Initially, each id is placed in a separate cluster
    const size_t num_id = m_id.size();

    vector<Cluster> prev(num_id);

    for(size_t i = 0;i < num_id;++i){
        prev[i].push_back(m_id[i]);
    }

    float prev_min_dist = std::numeric_limits<float>::max();
    
    for(size_t i = 0;i < num_id;++i){
        for(size_t j = i + 1;j < num_id;++j){

            const float d = cluster_distance(m_id, m_dist, prev[i], prev[j]);
            prev_min_dist = min(prev_min_dist,  d);
        }
    }

    ret.push_back( make_pair(prev_min_dist, prev) );

    cerr << "\tNumber of clusters: ";
    string info;
    size_t last_num_clusterings = ret.size();

    // Iteratively merge clusters until we have min_num_cluster or fewer clusters
    while(prev.size() > m_min_num_cluster){

        // Find the closest distance between clusters
        const size_t num_cluster = prev.size();

        // Update the progress
        if( last_num_clusterings != ret.size() ){

            last_num_clusterings = ret.size();

            stringstream ssout;

            ssout << num_cluster << ": min distance = " << prev_min_dist;

            for(string::const_iterator i = info.begin();i != info.end();++i){
                cerr << '\b';
            }

            for(string::const_iterator i = info.begin();i != info.end();++i){
                cerr << ' ';
            }

            for(string::const_iterator i = info.begin();i != info.end();++i){
                cerr << '\b';
            }

            info = ssout.str();

            cerr << info;
        }

        // The minimum distance between clusters
        float min_dist = std::numeric_limits<float>::max();

        // The minimum distance after we merge the two closest clusters
        float next_min_dist = std::numeric_limits<float>::max();

        pair<size_t, size_t> closest_clusters;

        #pragma omp parallel
        {
            float local_min_dist = std::numeric_limits<float>::max();
            float local_next_min_dist = std::numeric_limits<float>::max();
            
            pair<size_t, size_t> local_closest_clusters;

            #pragma omp for
            for(size_t i = 0;i < num_cluster;++i){

                for(size_t j = i + 1;j < num_cluster;++j){
                    
                    const float d = cluster_distance(m_id, m_dist, prev[i], prev[j]);

                    if(d < local_min_dist){

                        local_next_min_dist = local_min_dist;
                        local_min_dist = d;
                        local_closest_clusters = make_pair(i, j);
                    }
                    else if(d < local_next_min_dist){
                        local_next_min_dist = d;
                    }
                }
            }

            #pragma omp critical
            {
                if(local_min_dist < min_dist){

                    next_min_dist = min_dist;
                    min_dist = local_min_dist;
                    closest_clusters = local_closest_clusters;
                }
                else if(local_min_dist < next_min_dist){
                    next_min_dist = local_min_dist;
                }
                else if(local_next_min_dist < next_min_dist){
                    next_min_dist = local_next_min_dist;
                }
            }
        }

        // Merge the closest pair of clusters and copy the unmerged clusters    
        Cluster c;

        for(Cluster::const_iterator j = prev[closest_clusters.first].begin();j != prev[closest_clusters.first].end();++j){
            c.push_back(*j);
        }

        for(Cluster::const_iterator j = prev[closest_clusters.second].begin();j != prev[closest_clusters.second].end();++j){
            c.push_back(*j);
        }

        vector<Cluster> curr;

        curr.push_back(c);

        for(size_t i = 0;i < num_cluster;++i){

            if( (i == closest_clusters.first) || (i == closest_clusters.second) ){
                continue;
            }

            // Since this previous cluster was not merged, it is safe to copy
            curr.push_back(prev[i]);
        }

        // If we have a new minimum distance between clusters that is greater (but not equal to)
        // the previous minimum distance, then save this set of clusters
        if(next_min_dist > prev_min_dist){
            ret.push_back( make_pair(next_min_dist, curr) );
        }

        // Update the previous cluster
        prev.swap(curr);
        prev_min_dist = next_min_dist;
    }

    cerr << endl;

    // Impose a well defined order on the clusters to make it easier to compare clustering results
    for(vector< pair<float, vector<Cluster> > >::iterator i = ret.begin();i != ret.end();++i){
        sort( i->second.begin(), i->second.end() );
    }

    return ret;
}

// Find the single linkage distance (i.e., the closest distance between any element of m_a and any element of m_b)
float cluster_distance(const vector<string> &m_id, const SymmetricMatrix<float> &m_dist, const Cluster &m_a, const Cluster &m_b)
{
    float min_dist = std::numeric_limits<float>::max();

    for(Cluster::const_iterator i = m_a.begin();i != m_a.end();++i){

        vector<string>::const_iterator iter = lower_bound(m_id.begin(), m_id.end(), *i);

        if( (iter == m_id.end()) || (*iter != *i) ){
            throw __FILE__ ":cluster_distance: Unable to lookup index for cluster a";
        }

        const size_t index_a = iter - m_id.begin();

        for(Cluster::const_iterator j = m_b.begin();j != m_b.end();++j){

            iter = lower_bound(m_id.begin(), m_id.end(), *j);

            if( (iter == m_id.end()) || (*iter != *j) ){
                throw __FILE__ ":cluster_distance: Unable to lookup index for cluster b";
            }

            const size_t index_b = iter - m_id.begin();

            min_dist = min(min_dist, m_dist(index_a, index_b));
        }
    }

    return min_dist;
}

void write_clusters(const string &m_filename, const vector< pair<float/*min cluster distance*/, vector<Cluster> > > &m_data, 
    const Options &m_opt)
{
    ofstream fout( m_filename.c_str() );

    if(!fout){
        throw __FILE__ ":write_clusters: Unable to open output file for writing";
    }

    // Write the header information so we know how the clusters were generated
    const time_t timestamp = time(NULL);

    fout << "# Clusters generated on " << ctime(&timestamp); // ctime() output is '\n' terminated!
    fout << "# PDB data directory = " << m_opt.pdb_dir << endl;
    fout << "# Force dimer? = " << (m_opt.force_dimer ? "true" : "false") << endl;
    fout << "# Ignore hetatom? = " << (m_opt.ignore_hetatom ? "true" : "false") << endl;
    
    switch(m_opt.cluster_by){
        case Options::CLUSTER_BY_COMPOSITION:
            fout << "# Clustering by sequence composition" << endl;
            break;
        case Options::CLUSTER_BY_SEQ_ALIGN:
            fout << "# Clustering by pairwise sequence alignment" << endl;
            break;
        default:
            throw __FILE__ ":write_clusters: Unknown clustering algorithm!";
    };

    for(vector< pair<float/*min cluster distance*/, vector<Cluster> > >::const_iterator i = m_data.begin();i != m_data.end();++i){

        // Write the minimum distance between clusters
        fout << i->first << ':';

        for(vector<Cluster>::const_iterator j = i->second.begin();j != i->second.end();++j){

            fout << '(';

            for(Cluster::const_iterator k = j->begin();k != j->end();++k){
                
                fout << *k;

                if( (k + 1) != j->end() ){
                    fout << ',';
                }
            }

            fout << ')';
        }

        fout << endl;
    }
}

vector< pair<float/*min cluster distance*/, vector<Cluster> > > read_clusters(const string &m_filename)
{
    vector< pair<float, vector<Cluster> > > ret;

    ifstream fin(m_filename.c_str());

    if(!fin){

        cerr << "Unable to open " << m_filename << " for reading" << endl;
        throw __FILE__ ":read_clusters: Unable to open cluster file for reading";
    }

    // Parse clusters in the form of:
    // min_d_0:(A) (B) (C) (D)
    // min_d_1:(A, C) (B, D)
    // min_d_2:(A, B, C, D)
    // 
    // where 
    //  - Each line is a separate clustering
    //  - 'A', 'B', 'C', 'D' are strings representing PDB accession
    //  - The symbols '(', ')' and ',' are protected (and must not appear in PDB accession)
    //  - The set of clustered strings must be the same on each line
    //  - The sets that appear on a line must not overlap
    string line;

    size_t line_number = 0;

    while( getline(fin, line) ){

        ++line_number;

        // Skip text that follows a '#'
        string::size_type loc = line.find('#');

        if(loc != string::npos){
            line = line.substr(0, loc);
        }

        // Skip empty lines
        if( line.empty() ){
            continue;
        }

        // Read the minimum distance between clusters
        loc = line.find(':');

        float min_dist = -1.0f;

        if(loc != string::npos){

            min_dist = atof( line.substr(0, loc).c_str() );
            line = line.substr(loc + 1, line.size() - (loc + 1) );
        }

        vector<Cluster> curr;

        string::const_iterator begin = line.end();

        for(string::const_iterator i = line.begin();i != line.end();++i){

            if( begin == line.end() ){

                // Looking for the opening cluster symbol '('
                if(*i == '('){
                    begin = i + 1;
                }
                else if( !isspace(*i) ){

                    cerr << "Error reading clusters from " << m_filename << " on line " << line_number << endl;
                    throw __FILE__ ":read_clusters: Unable to find a valid cluster open symbol";
                }
            }
            else{

                // Looking for the closing cluster symbol ')'
                if(*i == ')'){
                    
                    // Found a matching cluster closing symbol
                    curr.push_back( parse_cluster( string(begin, i) ) );

                    begin = line.end(); 
                }
            }
        }

        if(begin != line.end()){
            throw __FILE__ ":read_clusters: Did not find a cluser closing symbol";
        }

        ret.push_back( make_pair(min_dist, curr) );
    }

    // Validate the clusters
    if( !ret.empty() ){
        
        unordered_set<string> total;

        for(vector< pair<float, vector<Cluster> > >::const_iterator i = ret.begin();i != ret.end();++i){

            if( total.empty() ){

                for(vector<Cluster>::const_iterator j = i->second.begin();j != i->second.end();++j){
                    for(Cluster::const_iterator k = j->begin();k != j->end();++k){

                        if( total.find(*k) != total.end() ){
                            throw __FILE__ ":read_clusters: Overlapping cluster found";
                        }

                        total.insert(*k);
                    }
                }
            }
            else{

                unordered_set<string> curr;

                for(vector<Cluster>::const_iterator j = i->second.begin();j != i->second.end();++j){
                    for(Cluster::const_iterator k = j->begin();k != j->end();++k){

                        if( curr.find(*k) != curr.end() ){
                            throw __FILE__ ":read_clusters: Overlapping cluster found";
                        }

                        curr.insert(*k);
                    }
                }

                if(curr != total){
                    throw __FILE__ ":read_clusters: Inconsistent set elements between clusterings";
                }

            }
        }
    }

    return ret;
}

Cluster parse_cluster(const string &m_str)
{
    Cluster ret;

    const vector<string> elem = split(m_str, ',');

    if( elem.empty() ){
        throw __FILE__ ":parse_cluster: Did not find any set elements";
    }

    for(vector<string>::const_iterator i = elem.begin();i != elem.end();++i){
        
        ret.push_back( strip(*i) );

        if( ret.back().empty() ){
            throw __FILE__ ":parse_cluster: Missing set element";
        }
    }

    return ret;
}

int perfect_BLOSUM62_score(const string &m_seq)
{
    int ret = 0;

    for(string::const_iterator i = m_seq.begin();i != m_seq.end();++i){

        // Scores are from the diagonal elements of the NCBI BLOSUM62 matrix
        switch(*i){
            case 'A':
                ret += 4;
                break;
            case 'R':
                ret += 5;
                break;
            case 'N':
                ret += 6;
                break;
            case 'D':
                ret += 6;
                break;
            case 'C':
                ret += 9;
                break;
            case 'Q':
                ret += 5;
                break;
            case 'E':
                ret += 5;
                break;
            case 'G':
                ret += 6;
                break;
            case 'H':
                ret += 8;
                break;
            case 'I':
                ret += 4;
                break;
            case 'L':
                ret += 4;
                break;
            case 'K':
                ret += 5;
                break;
            case 'M':
                ret += 5;
                break;
            case 'F':
                ret += 6;
                break;
            case 'P':
                ret += 7;
                break;
            case 'S':
                ret += 4;
                break;
            case 'T':
                ret += 5;
                break;
            case 'W':
                ret += 11;
                break;
            case 'Y':
                ret += 7;
                break;
            case 'V':
                ret += 4;
                break;
            default:
                throw __FILE__ ":perfect_BLOSUM62_score: Unknown amino acid";
        };
    }

    return ret;
}