#include "pdb.h"
#include "tensor.h"
#include <unordered_set>
#include <deque>

#include <iostream> // DEBUG

using namespace std;

#define     MAX_GRID_DIM            20
#define     DONOR_THRESHOLD         3.0 // Max distance (in angstroms) from the donor to the hydrogen atom
#define     ACCEPTOR_THRESHOLD      3.0 // Max distance (in angstroms) from the acceptor to the hydrogen atom

void inventory_hydrogen_bond_donors(vector<AminoAcidData> &m_mol, const unordered_map<string, unsigned int> &m_atom_table, 
    const string &m_pdb /*for debugging*/)
{
    // Identify all of the atom types classified as "hydrogen"
    unordered_set<unsigned int> hydrogen;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){

        // Assume that only hydrogen atoms have a name that starts with the letter 'H'
        if(i->first.find('H') == 0){
            hydrogen.insert(i->second);
        }
    }

    if( hydrogen.empty() ){
        return;
    }

    // For each hydrogen atom, find the closest atom (within the same protein). This atom is assumed to be the donor atom.
    // Accelerate this search by binning the atoms by 3D coordinates.
    float x_lower = std::numeric_limits<float>::max();
    float x_upper = std::numeric_limits<float>::min();
    float y_lower = std::numeric_limits<float>::max();
    float y_upper = std::numeric_limits<float>::min();
    float z_lower = std::numeric_limits<float>::max();
    float z_upper = std::numeric_limits<float>::min();

    // Compute the bounds of the protein atoms in 3D space
    for(vector<AminoAcidData>::const_iterator i = m_mol.begin();i != m_mol.end();++i){
    
        for(vector<AtomData>::const_iterator j = i->atoms.begin();j != i->atoms.end();++j){

            x_lower = min(x_lower, j->x);
            x_upper = max(x_upper, j->x);

            y_lower = min(y_lower, j->y);
            y_upper = max(y_upper, j->y);

            z_lower = min(z_lower, j->z);
            z_upper = max(z_upper, j->z);
        }
    }

    const float dx = x_upper - x_lower;
    const float dy = y_upper - y_lower;
    const float dz = z_upper - z_lower;

    if( (dx <= 0.0) || (dy <= 0.0) || (dz <= 0.0) ){
        throw "inventory_hydrogen_bond_donors: Invalid grid boundary";
    }

    const int GRID_DIM_X = min( MAX_GRID_DIM, int( dx/DONOR_THRESHOLD ) );
    const int GRID_DIM_Y = min( MAX_GRID_DIM, int( dy/DONOR_THRESHOLD ) );
    const int GRID_DIM_Z = min( MAX_GRID_DIM, int( dz/DONOR_THRESHOLD ) );

    // Bin the atoms using a uniform 3D grid to accelerate finding the minimum distance between two atoms.
    // While tedious to code, this approach is O(N), while the brute-force approach of testing all possible
    // atoms pairs is O(N^2). This same approach is used to 3D molecular viewing software to identify covalent bonds
    //vector< vector< vector< deque< std::pair<int, int> > > > > atom_grid(GRID_DIM_X,
    //    vector< vector< deque< std::pair<int, int> > > >(GRID_DIM_Y,
    //        vector< deque< std::pair<int, int> > >(GRID_DIM_Z) ) );
    Tensor< deque< std::pair<int, int> > > atom_grid(GRID_DIM_X, GRID_DIM_Y, GRID_DIM_Z);

    const int num_amino_acid = m_mol.size();

    for(int i = 0;i < num_amino_acid;++i){
    
        const vector<AtomData> &aa = m_mol[i].atoms;

        const int num_atom = aa.size();

        for(int j = 0;j < num_atom;++j){

            const int x_index = min(GRID_DIM_X - 1, int(GRID_DIM_X*(aa[j].x - x_lower)/dx) );
            const int y_index = min(GRID_DIM_Y - 1, int(GRID_DIM_Y*(aa[j].y - y_lower)/dy) );
            const int z_index = min(GRID_DIM_Z - 1, int(GRID_DIM_Z*(aa[j].z - z_lower)/dz) );

            // In the current scheme, atoms must be indexed by amino acid and then by atom
            atom_grid(x_index, y_index, z_index).push_back( make_pair(i, j) );
        }
    }

    for(int i = 0;i < GRID_DIM_X;++i){
        for(int j = 0;j < GRID_DIM_Y;++j){
            for(int k = 0;k < GRID_DIM_Z;++k){

                for(deque< pair<int, int> >::const_iterator h_index = atom_grid(i, j, k).begin();
                    h_index != atom_grid(i, j, k).end();++h_index){

                    // h_index must refer to a hydrogen atom
                    if(hydrogen.find(m_mol[h_index->first].atoms[h_index->second].type) == hydrogen.end()){
                        continue;
                    }

                    pair<int, int> best_donor_index(-1, -1);
                    float best_donor_distance = std::numeric_limits<float>::max();

                    for(int n_i = i - 1;n_i <= i + 1;++n_i){
                        
                        if( (n_i < 0) || (n_i >= GRID_DIM_X) ) {
                            continue;
                        }

                        for(int n_j  = j - 1;n_j <= j + 1;++n_j){
                        
                            if( (n_j < 0) || (n_j >= GRID_DIM_Y) ){
                                continue;
                            }

                            for(int n_k = k - 1;n_k <= k + 1;++n_k){

                                if( (n_k < 0) || (n_k >= GRID_DIM_Z) ){
                                    continue;
                                }

                                for(deque< pair<int, int> >::const_iterator d_index = atom_grid(n_i, n_j, n_k).begin();
                                    d_index != atom_grid(n_i, n_j, n_k).end();++d_index){

                                    // The donor atom *cannot* be a hydrogen
                                    if(hydrogen.find(m_mol[d_index->first].atoms[d_index->second].type) != hydrogen.end()){
                                        continue;
                                    }

                                    const float d = m_mol[h_index->first].atoms[h_index->second].distance(m_mol[d_index->first].atoms[d_index->second]);

                                    if(d < best_donor_distance){

                                        best_donor_distance = d;
                                        best_donor_index = *d_index;
                                    }
                                }
                            }
                        }
                    }

                    if(best_donor_distance <= DONOR_THRESHOLD){
                        m_mol[h_index->first].atoms[h_index->second].h_donor_index = best_donor_index;
                    }
                    // DEBUG -- hydrogen atoms can be "orphaned" by missing donor atoms in the input 3D structures ...
                    //else{
                    //
                    //    cerr << "Unable to find hydrogen bond donor for " << m_pdb << " (best donor distance is " << best_donor_distance 
                    //        << " angstroms): (" 
                    //        << m_mol[h_index->first].atoms[h_index->second].x << ", "
                    //        << m_mol[h_index->first].atoms[h_index->second].y << ", "
                    //        << m_mol[h_index->first].atoms[h_index->second].z << ")" << endl;
                    //}
                }
            }
        }
    }
}

void inventory_interface_hydrogen_bond_acceptors(PDBComplex &m_dimer, const unordered_map<string, unsigned int> &m_atom_table, 
    const string &m_pdb /*for debugging*/)
{
    // Identify all of the atom types classified as "hydrogen"
    unordered_set<unsigned int> hydrogen;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){

        // Assume that only hydrogen atoms have a name that starts with the letter 'H'
        if(i->first.find('H') == 0){
            hydrogen.insert(i->second);
        }
    }

    if( hydrogen.empty() ){
        return;
    }

    // For each hydrogen atom, find the closest non-hydrogen atom within the *other* protein. 
    // If the distance is less than ACCEPTOR_THRESHOLD then this atom is assumed to be an acceptor atom.
    // Accelerate this search by binning the atoms by 3D coordinates.
    float x_lower = std::numeric_limits<float>::max();
    float x_upper = std::numeric_limits<float>::min();
    float y_lower = std::numeric_limits<float>::max();
    float y_upper = std::numeric_limits<float>::min();
    float z_lower = std::numeric_limits<float>::max();
    float z_upper = std::numeric_limits<float>::min();

    // Compute the bounds of the protein atoms in 3D space
    #define EXTRACT_BOUNDS(MOL) \
        for(vector<AminoAcidData>::const_iterator i = MOL.begin();i != MOL.end();++i){ \
            \
            for(vector<AtomData>::const_iterator j = i->atoms.begin();j != i->atoms.end();++j){ \
                \
                x_lower = min(x_lower, j->x); \
                x_upper = max(x_upper, j->x); \
                \
                y_lower = min(y_lower, j->y); \
                y_upper = max(y_upper, j->y); \
                \
                z_lower = min(z_lower, j->z); \
                z_upper = max(z_upper, j->z); \
            } \
        }

    EXTRACT_BOUNDS(m_dimer.first);
    EXTRACT_BOUNDS(m_dimer.second);

    const float dx = x_upper - x_lower;
    const float dy = y_upper - y_lower;
    const float dz = z_upper - z_lower;

    if( (dx <= 0.0) || (dy <= 0.0) || (dz <= 0.0) ){
        throw "inventory_interface_hydrogen_bond_acceptors: Invalid grid boundary";
    }

    const int GRID_DIM_X = min( MAX_GRID_DIM, int( dx/ACCEPTOR_THRESHOLD ) );
    const int GRID_DIM_Y = min( MAX_GRID_DIM, int( dy/ACCEPTOR_THRESHOLD ) );
    const int GRID_DIM_Z = min( MAX_GRID_DIM, int( dz/ACCEPTOR_THRESHOLD ) );

    // Bin the atoms using a uniform 3D grid to accelerate finding the minimum distance between two atoms.
    // While tedious to code, this approach is O(N), while the brute-force approach of testing all possible
    // atoms pairs is O(N^2). This same approach is used to 3D molecular viewing software to identify covalent bonds
    Tensor< deque< std::pair<int, int> > > first_grid(GRID_DIM_X, GRID_DIM_Y, GRID_DIM_Z);
    Tensor< deque< std::pair<int, int> > > second_grid(GRID_DIM_X, GRID_DIM_Y, GRID_DIM_Z);

    #define BIN_ATOMS(GRID, MOL) \
        for(int i = 0;i < int( MOL.size() );++i){ \
            \
            const vector<AtomData> &aa = MOL[i].atoms; \
            \
            const int num_atom = aa.size(); \
            \
            for(int j = 0;j < num_atom;++j){ \
                \
                const int x_index = min(GRID_DIM_X - 1, int(GRID_DIM_X*(aa[j].x - x_lower)/dx) ); \
                const int y_index = min(GRID_DIM_Y - 1, int(GRID_DIM_Y*(aa[j].y - y_lower)/dy) ); \
                const int z_index = min(GRID_DIM_Z - 1, int(GRID_DIM_Z*(aa[j].z - z_lower)/dz) ); \
                \
                /* In the current scheme, atoms must be indexed by amino acid and then by atom */ \
                GRID(x_index, y_index, z_index).push_back( make_pair(i, j) ); \
            } \
        }

    BIN_ATOMS(first_grid, m_dimer.first);
    BIN_ATOMS(second_grid, m_dimer.second);

    for(int i = 0;i < GRID_DIM_X;++i){
        for(int j = 0;j < GRID_DIM_Y;++j){
            for(int k = 0;k < GRID_DIM_Z;++k){

                #define FIND_ACCEPTORS(FIRST_GRID, SECOND_GRID, FIRST_MOL, SECOND_MOL) \
                    for(deque< pair<int, int> >::const_iterator h_index = FIRST_GRID(i, j, k).begin(); \
                        h_index != FIRST_GRID(i, j, k).end();++h_index){ \
                        \
                        /* h_index must refer to a hydrogen atom */ \
                        if(hydrogen.find(FIRST_MOL[h_index->first].atoms[h_index->second].type) == hydrogen.end()){ \
                            continue; \
                        } \
                        \
                        pair<int, int> best_acceptor_index(-1, -1); \
                        float best_acceptor_distance = std::numeric_limits<float>::max(); \
                        \
                        for(int n_i = i - 1;n_i <= i + 1;++n_i){ \
                            \
                            if( (n_i < 0) || (n_i >= GRID_DIM_X) ) { \
                                continue; \
                            } \
                            \
                            for(int n_j  = j - 1;n_j <= j + 1;++n_j){ \
                                \
                                if( (n_j < 0) || (n_j >= GRID_DIM_Y) ){ \
                                    continue; \
                                } \
                                \
                                for(int n_k = k - 1;n_k <= k + 1;++n_k){ \
                                    \
                                    if( (n_k < 0) || (n_k >= GRID_DIM_Z) ){ \
                                        continue; \
                                    } \
                                    \
                                    for(deque< pair<int, int> >::const_iterator a_index = SECOND_GRID(n_i, n_j, n_k).begin(); \
                                        a_index != SECOND_GRID(n_i, n_j, n_k).end();++a_index){ \
                                        \
                                        /* The acceptor atom *cannot* be a hydrogen */ \
                                        if(hydrogen.find(SECOND_MOL[a_index->first].atoms[a_index->second].type) != hydrogen.end()){ \
                                            continue; \
                                        } \
                                        \
                                        const float d = FIRST_MOL[h_index->first].atoms[h_index->second].distance(SECOND_MOL[a_index->first].atoms[a_index->second]); \
                                        \
                                        if(d < best_acceptor_distance){ \
                                            \
                                            best_acceptor_distance = d; \
                                            best_acceptor_index = *a_index; \
                                        } \
                                    } \
                                } \
                            } \
                        } \
                        \
                        if(best_acceptor_distance <= ACCEPTOR_THRESHOLD){ \
                            FIRST_MOL[h_index->first].atoms[h_index->second].interface_h_acceptor_index = best_acceptor_index; \
                        } \
                    }

                FIND_ACCEPTORS(first_grid, second_grid, m_dimer.first, m_dimer.second);
                FIND_ACCEPTORS(second_grid, first_grid, m_dimer.second, m_dimer.first);
            }
        }
    }
}