#ifndef __CLUSTER
#define __CLUSTER

#include <vector>
#include "pdb.h"
#include "options.h"
typedef std::vector<std::string> Cluster;

void cluster_by_composition(const std::string &m_filename, const std::unordered_map<std::string, PDBComplex> &m_data,
    const Options &m_opt);
void cluster_by_seq_align(const std::string &m_filename, const std::unordered_map<std::string, PDBComplex> &m_data,
    const Options &m_opt);
std::vector< std::pair<float/*min cluster distance*/, std::vector<Cluster> > > read_clusters(const std::string &m_filename);
void write_clusters(const std::string &m_filename, 
    const std::vector< std::pair<float/*min cluster distance*/, std::vector<Cluster> > > &m_data, 
    const Options &m_opt);

#endif // __CLUSTER