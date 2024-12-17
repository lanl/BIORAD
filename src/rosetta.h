#ifndef __ROSETTA_UTIL
#define __ROSETTA_UTIL

#include <unordered_map>
#include <vector>
#include <string>

void parse_rosetta_scores(const std::string &m_filename, std::unordered_map< std::string, std::vector<float> > &m_scores);

#endif // __ROSETTA_UTIL