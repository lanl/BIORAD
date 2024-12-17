#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits.h>
#include <stdlib.h> // atof
#include "parse_util.h"

using namespace std;

vector<string> split(const string &m_line, const char &m_delim)
{
	// Delimiters are litteral if they are contained within matching protect characters
	const char protect = '"';
	
	vector<string> ret;
	
	size_t len = 1;
	
	size_t protect_count = 0;
	
	for(string::const_iterator i = m_line.begin();i != m_line.end();++i){
		
		// A '\n' or DOS/Windows '\r' symbol forces the end of the line
		if( (*i == '\r') || (*i == '\n') ){
			break;
		}
		
		protect_count += (*i == protect);
		
		len += (*i == m_delim) && (protect_count%2 == 0);
	}
	
	if(protect_count%2 == 1){
		throw __FILE__ ":split: Unmatched protection symbol";
	}
	
	ret.resize(len);
	
	for(size_t i = 0;i < len;++i){
		ret[i].reserve(32);
	}
	
	size_t index = 0;
	
	protect_count = 0;
	
	for(string::const_iterator i = m_line.begin();i != m_line.end();++i){	
		
		// A '\n' or DOS/Windows '\r' symbol forces the end of the line
		if( (*i == '\r') || (*i == '\n') ){
			break;
		}
		
		protect_count += (*i == protect);
		
		if( (*i == m_delim) && (protect_count%2 == 0) ){
			++index;
		}
		else{
			ret[index].push_back(*i);
		}
	}
	
	return ret;
}

size_t get_col(const vector<string> &m_header, const string &m_key)
{
	vector<string>::const_iterator iter = find(m_header.begin(), m_header.end(), m_key);

	if( ( iter == m_header.end() ) || (*iter != m_key) ){
		return COLUMN_NOT_FOUND;
	}
	
	return ( iter - m_header.begin() );
}

size_t string_to_size_t(const string &m_str)
{
	size_t ret = 0;
	
	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){
		
		switch(*i){
			case '0':
				ret = ret*10;
				break;
			case '1':
				ret = ret*10 + 1;
				break;
			case '2':
				ret = ret*10 + 2;
				break;
			case '3':
				ret = ret*10 + 3;
				break;
			case '4':
				ret = ret*10 + 4;
				break;
			case '5':
				ret = ret*10 + 5;
				break;
			case '6':
				ret = ret*10 + 6;
				break;
			case '7':
				ret = ret*10 + 7;
				break;
			case '8':
				ret = ret*10 + 8;
				break;
			case '9':
				ret = ret*10 + 9;
				break;
			default:
				cerr << "Unable to parse: \"" << m_str << "\"" << endl;
				throw __FILE__ ":string_to_size_t: Illegal symbol!";
				break;
		};
	}
	
	return ret;
}

unsigned int string_to_uint(const string &m_str)
{
	const size_t ret = string_to_size_t(m_str);
	
	if(ret > UINT_MAX){
		throw __FILE__ ":string_to_uint: Overflow!";
	}
	
	return (unsigned int)ret;
}

float string_to_float(const string &m_str)
{
	return atof( m_str.c_str() );
}

bool has_digit(const string &m_str)
{
	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){
		if( isdigit(*i) ){
			return true;
		}
	}
	
	return false;
}

string toupper(const string &m_str)
{
    string ret(m_str);

    for(string::iterator i = ret.begin();i != ret.end();++i){
        *i = ::toupper(*i);
    }

    return ret;
}

// Remove leading and trailing white space
string strip(const string &m_str)
{
    string::const_iterator begin = m_str.begin();
    string::const_iterator end = m_str.end();

    while( (begin != end) && isspace(*begin) ){
        ++begin;
    }

    while( (begin != end) && isspace(*(end - 1) ) ){
        --end;
    }

    return string(begin, end);
}

void parse_csv_scores(const std::string &m_filename, std::unordered_map<std::string, std::vector<float> > &m_scores)
{
	ifstream fin( m_filename.c_str() );

	if(!fin){
		throw __FILE__ ":parse_csv_scores: Unable to open input CSV for reading";
	}

	string line;

	int num_col = -1;

	while( getline(fin, line) ){

		// Skip commented lines
		const string::size_type pos = line.find('#');

		if(pos == 0){
			continue;
		}

		const vector<string> data = split(line, ',');

		if(data.size() == 0){
			// Skip blank lines
			continue;
		}

		if(num_col < 0){
			num_col = data.size();
		}
		else{
			if(num_col != int(data.size())){
				throw __FILE__ ":parse_csv_scores: Variable number of columns not allowed";
			}
		}

		if(num_col == 1){
			throw __FILE__ ":parse_csv_scores: Did not find any feature values!";
		}

		if( m_scores.find(data[0]) != m_scores.end() ){
			throw __FILE__ ":parse_csv_scores: Duplicate PDB accession";
		}

		vector<float> features(num_col - 1);

		for(int i = 1;i < num_col;++i){
			features[i - 1] = stof(data[i]);
		}

		m_scores[data[0]] = features;
	}

	fin.close();
}