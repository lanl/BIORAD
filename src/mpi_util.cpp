#include "mpi_util.h"
#include "options.h"
#include "pdb.h"
#include "regression.h"
#include "affinity.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<__uint128_t>(const __uint128_t &m_obj)
{
	return sizeof(__uint128_t);
}

template<>
unsigned char* mpi_pack<__uint128_t>(unsigned char* m_ptr, const __uint128_t &m_obj)
{
	memcpy( m_ptr, &m_obj, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<__uint128_t>(unsigned char* m_ptr, __uint128_t &m_obj)
{
	memcpy( &m_obj, m_ptr, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Options
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<Options>(const Options &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		BIORAD_OPTIONS_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<Options>(unsigned char* m_ptr, const Options &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	BIORAD_OPTIONS_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<Options>(unsigned char* m_ptr, Options &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	BIORAD_OPTIONS_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Affinity
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<Affinity>(const Affinity &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		BIORAD_AFFINITY_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<Affinity>(unsigned char* m_ptr, const Affinity &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	BIORAD_AFFINITY_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<Affinity>(unsigned char* m_ptr, Affinity &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	BIORAD_AFFINITY_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for AtomData
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<AtomData>(const AtomData &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		ATOM_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<AtomData>(unsigned char* m_ptr, const AtomData &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	ATOM_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<AtomData>(unsigned char* m_ptr, AtomData &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	ATOM_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for AminoAcidData
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<AminoAcidData>(const AminoAcidData &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		AMINO_ACID_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<AminoAcidData>(unsigned char* m_ptr, const AminoAcidData &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	AMINO_ACID_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<AminoAcidData>(unsigned char* m_ptr, AminoAcidData &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	AMINO_ACID_MEMBERS
	#undef VARIABLE

    return m_ptr;
}
