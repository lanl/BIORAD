#include "rosetta.h"
#include "parse_util.h"
#include <fstream>
#include <sstream>
#include <deque>

using namespace std;

void parse_rosetta_scores(const string &m_filename, unordered_map<string, vector<float>> &m_scores)
{
    ifstream fin( m_filename.c_str() );

    if(!fin){
        throw __FILE__ ":parse_rosetta_scores: Unable to open score file for reading";
    }

    string line;

    // Read the header
    if( !getline(fin, line) ){
        throw __FILE__ ":parse_rosetta_scores: Unable to read the first header line";
    }

    if(line.find("SEQUENCE:") != 0){
        throw __FILE__ ":parse_rosetta_scores: Did not read the expected \"SEQUENCE:\" header";
    }

    if( !getline(fin, line) ){
        throw __FILE__ ":parse_rosetta_scores: Unable to read the second header line";
    }

    if(line.find("SCORE:") != 0){
        throw __FILE__ ":parse_rosetta_scores: Did not read the expected \"SCORE:\" header";
    }

    size_t num_features = 0;

    while( getline(fin, line) ){

        if( line.empty() ){
            continue;
        }

        stringstream ssin(line);
        deque<string> data;
        string tmp;

        while(ssin >> tmp){
            data.push_back(tmp);
        }

        if(data.size() < 3){
            throw __FILE__ ":parse_rosetta_scores: Did not read enough columns!";
        }

        if(num_features == 0){

            // Do not include the first or last column
            num_features = data.size() - 2;
        }
        else{
            if( num_features != (data.size() - 2) ){
                throw __FILE__ ":parse_rosetta_scores: Encountered a variable number of columns";
            }
        }

        // The first column must be "SCORE:"
        if(data[0] != "SCORE:"){
            throw __FILE__ ":parse_rosetta_scores: Did not read \"SCORE:\" in the first column";
        }

        vector<float> features(num_features);

        for(size_t i = 1;i <= num_features;++i){
            features[i - 1] = atof( data[i].c_str() );
        }

        string::size_type loc = data.back().find("_0001");

        if(loc == string::npos){
            throw __FILE__ ":parse_rosetta_scores: Unable to find \"_0001\" suffix in PDB accession";
        }

        const string id = data.back().substr(0, loc);

        if( m_scores.find(id) != m_scores.end() ){
            throw __FILE__ ":parse_rosetta_scores: Duplicate PDB accession in score file";
        }

        m_scores[id] = features;
    }
}