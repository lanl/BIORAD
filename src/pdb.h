#ifndef __PDB
#define __PDB

#include <math.h>
#include <vector>
#include <unordered_map>
#include <ostream>

#define     INVALID_ATOM_TYPE   0xFFFFFFFF

// Forward declaration of mpi helper functions to keep the compiler happy
template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to
// ensure that enumerated variables and the strings representing their
// names are always synchronized
#define AMINO_ACID_TYPES \
    AMINO_ACID_NAME(ALA) \
    AMINO_ACID_NAME(ARG) \
    AMINO_ACID_NAME(ASN) \
    AMINO_ACID_NAME(ASP) \
    AMINO_ACID_NAME(CYS) \
    AMINO_ACID_NAME(GLN) \
    AMINO_ACID_NAME(GLU) \
    AMINO_ACID_NAME(GLY) \
    AMINO_ACID_NAME(HIS) \
    AMINO_ACID_NAME(ILE) \
    AMINO_ACID_NAME(LEU) \
    AMINO_ACID_NAME(LYS) \
    AMINO_ACID_NAME(MET) \
    AMINO_ACID_NAME(PHE) \
    AMINO_ACID_NAME(PRO) \
    AMINO_ACID_NAME(SER) \
    AMINO_ACID_NAME(THR) \
    AMINO_ACID_NAME(TRP) \
    AMINO_ACID_NAME(TYR) \
    AMINO_ACID_NAME(VAL) \
    AMINO_ACID_NAME(UNKNOWN)

// Enumerate the different groupings of amino acids
#define CHEMICAL_PROPERTIES \
    AMINO_ACID_PROPERTY(AA_POSITIVE) \
    AMINO_ACID_PROPERTY(AA_NEGATIVE) \
    AMINO_ACID_PROPERTY(AA_POLAR) \
    AMINO_ACID_PROPERTY(AA_HYDROPHOBIC) \
    AMINO_ACID_PROPERTY(UNKNOWN_PROPERTY)

extern const char* amino_acid_names[];

#define SINGLE_ARG(...) __VA_ARGS__

class Transform
{
    private:
        std::vector<float> rotation;
        std::vector<float> translation;

    public:
        Transform()
        {
            reset();
        };

        inline void reset()
        {
            rotation.assign(9, 0.0f); // Rotation in 3D
            translation.assign(3, 0.0f); // Translation in 3D

            // Initialize the rotation to the idenitity matrix
            (*this)(0, 0) = 1.0;
            (*this)(1, 1) = 1.0;
            (*this)(2, 2) = 1.0;
        };

        inline float& operator()(const unsigned int m_row, const unsigned int m_col)
        {
            if( (m_row > 2) || (m_col > 2) ){
                throw __FILE__ ":Transform::(): row and/or column index out of bounds";
            }

            return rotation[m_row*3 + m_col];
        };

        inline float operator()(const unsigned int m_row, const unsigned int m_col) const
        {
            if( (m_row > 2) || (m_col > 2) ){
                throw __FILE__ ":Transform::(): row and/or column index out of bounds";
            }

            return rotation[m_row*3 + m_col];
        };

        inline float& operator[](const unsigned int m_index)
        {
            if(m_index > 2){
                throw __FILE__ ":Transform::[]: index out of bounds";
            }

            return translation[m_index];
        };

        inline float operator[](const unsigned int m_index) const
        {
            if(m_index > 2){
                throw __FILE__ ":Transform::[]: index out of bounds";
            }

            return translation[m_index];
        };

        inline bool is_identity() const
        {
            return (translation[0] == 0.0f) && (translation[1] == 0.0f) && (translation[2] == 0.0f) &&
                (rotation[0] == 1.0f) && (rotation[1] == 0.0f) && (rotation[2] == 0.0f) &&
                (rotation[3] == 0.0f) && (rotation[4] == 1.0f) && (rotation[5] == 0.0f) &&
                (rotation[6] == 0.0f) && (rotation[7] == 0.0f) && (rotation[8] == 1.0f);
        };
};

struct AtomData
{
    #define ATOM_MEMBERS \
        VARIABLE(unsigned int, type) \
        VARIABLE(int, index) \
		VARIABLE(float, x) \
        VARIABLE(float, y) \
        VARIABLE(float, z) \
        VARIABLE(float, temperature) \
        VARIABLE(SINGLE_ARG(std::pair<int, int>), h_donor_index) \
        VARIABLE(SINGLE_ARG(std::pair<int, int>), interface_h_acceptor_index)

	#define VARIABLE(A, B) A B;
		ATOM_MEMBERS
	#undef VARIABLE

    AtomData(){};

    AtomData(const unsigned int &m_type, const int &m_index, 
        const float &m_x, const float &m_y, const float &m_z, const float &m_temperature) :
        type(m_type), index(m_index), x(m_x), y(m_y), z(m_z), temperature(m_temperature)
    {
        // If this atom is a hydrogen, the index of the donor atom in the same protein will be stored in interface_h_donor_index.
        // If this atom is a hydrogen and participates in an interface hydrogen bond, the index of acceptor atom in the *other* protein 
        // will be stored in interface_h_acceptor_index.
        h_donor_index = std::make_pair(-1, -1); // amino acid index, atom index
        interface_h_acceptor_index = std::make_pair(-1, -1); // amino acid index, atom index
    };

    inline float distance(const AtomData &m_rhs) const
    {
        return sqrt( fabs(distance2(m_rhs)) );
    };

    inline float distance2(const AtomData &m_rhs) const
    {
        const float dx = x - m_rhs.x;
        const float dy = y - m_rhs.y;
        const float dz = z - m_rhs.z;

        return (dx*dx + dy*dy + dz*dz);
    };

    inline void apply_transformation(const Transform &m_t)
    {
        // Rotate
        float x_new = m_t(0, 0)*x + m_t(0, 1)*y + m_t(0, 2)*z;
        float y_new = m_t(1, 0)*x + m_t(1, 1)*y + m_t(1, 2)*z;
        float z_new = m_t(2, 0)*x + m_t(2, 1)*y + m_t(2, 2)*z;

        // Translate
        x_new += m_t[0];
        y_new += m_t[1];
        z_new += m_t[2];

        x = x_new;
        y = y_new;
        z = z_new;
    };

    inline bool has_interface_hydrogen_bond() const
    {
        if( (h_donor_index.first < 0) || (h_donor_index.second < 0) ){
            return false;
        }

        if( (interface_h_acceptor_index.first < 0) || (interface_h_acceptor_index.second < 0) ){
            return false;
        }

        return true;
    };
};

template<> size_t mpi_size(const AtomData &m_opt);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const AtomData &m_opt);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, AtomData &m_opt);

struct AminoAcidData
{
    // The allowed amino acid types

    #define AMINO_ACID_NAME(X)  X,
    enum{
        AMINO_ACID_TYPES
        NUM_AMINO_ACID_TYPES = UNKNOWN + 1,
        NUM_REAL_AMINO_ACID_TYPES
    };
    #undef AMINO_ACID_NAME

    #define AMINO_ACID_PROPERTY(X)  X,
    enum{
        CHEMICAL_PROPERTIES
        NUM_AMINO_ACID_PROPERTY = UNKNOWN_PROPERTY + 1
    };
    #undef AMINO_ACID_PROPERTY

    #define AMINO_ACID_MEMBERS \
        VARIABLE(unsigned int, type) \
		VARIABLE(std::vector<AtomData>, atoms)
	
	#define VARIABLE(A, B) A B;
		AMINO_ACID_MEMBERS
	#undef VARIABLE

    AminoAcidData(){};

    AminoAcidData(const unsigned int &m_type) :
        type(m_type)
    {
    };

    inline unsigned int property() const
    {
        switch(type){
            case AminoAcidData::ARG: case AminoAcidData::HIS: case AminoAcidData::LYS:
                return AminoAcidData::AA_POSITIVE;
            case AminoAcidData::ASP: case AminoAcidData::GLU:
                return AminoAcidData::AA_NEGATIVE;
            case AminoAcidData::SER: case AminoAcidData::THR: case AminoAcidData::ASN: case AminoAcidData::GLN:
                return AminoAcidData::AA_POLAR;
            case AminoAcidData::ALA: case AminoAcidData::VAL: case AminoAcidData::ILE: case AminoAcidData::LEU:
            case AminoAcidData::MET: case AminoAcidData::PHE: case AminoAcidData::TYR: case AminoAcidData::TRP:
                return AminoAcidData::AA_HYDROPHOBIC;  
        };

        return AminoAcidData::UNKNOWN_PROPERTY;
    };

    inline void apply_transformation(const Transform &m_t)
    {
        for(std::vector<AtomData>::iterator i = atoms.begin();i != atoms.end();++i){
            i->apply_transformation(m_t);
        }
    };
};

template<> size_t mpi_size(const AminoAcidData &m_opt);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const AminoAcidData &m_opt);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, AminoAcidData &m_opt);

// A pair of structures is a PDBComplex. A typedef is used (rather than a wrapper struct or class that
// has the pair of vectors as a base class) so that we do not need to write explicit mpi_pack/unpack/size
// helper function.
typedef std::pair< std::vector<AminoAcidData>, std::vector<AminoAcidData> > PDBComplex;

void parse_pdb_data(std::unordered_map<std::string, PDBComplex> &m_pdb_data, const std::string &m_dir,
    std::unordered_map<std::string, unsigned int> &m_atom_table, 
    const bool &m_ignore_hetatom_chain, const bool &m_ignore_charge, const bool &m_force_dimer, 
    std::ostream &m_out);

std::string protein_sequence(const std::vector<AminoAcidData> &m_seq);

void inventory_hydrogen_bond_donors(std::vector<AminoAcidData> &m_mol, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const std::string &m_pdb);
void inventory_interface_hydrogen_bond_acceptors(PDBComplex &m_dimer, 
    const std::unordered_map<std::string, unsigned int> &m_atom_table, const std::string &m_pdb);

#endif // __PDB