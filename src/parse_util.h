#ifndef __PARSE_CSV
#define __PARSE_CSV

#include <vector>
#include <string>
#include <unordered_map>

#define	COLUMN_NOT_FOUND	0xFFFFFFFFFFFFFFFF

std::vector<std::string> split(const std::string &m_buffer, const char &m_delim);
size_t get_col(const std::vector<std::string> &m_header, const std::string &m_key);
unsigned int string_to_uint(const std::string &m_str);
float string_to_float(const std::string &m_str);
size_t string_to_size_t(const std::string &m_str);
bool has_digit(const std::string &m_str);
std::string toupper(const std::string &m_str);
std::string strip(const std::string &m_str);
void parse_csv_scores(const std::string &m_filename, std::unordered_map<std::string, std::vector<float> > &m_scores);

#endif // __PARSE_CSV
