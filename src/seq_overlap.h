#ifndef __SEQ_OVERLAP
#define __SEQ_OVERLAP

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <deque>

namespace SO {

// Older versions of the Visual Studio C++ compiler do
// not povide std::min or std::max. _MSC_VER < 1300 should
// capture all of these older version (like 6.0).
#if defined(WIN32) && (_MSC_VER < 1300)

#ifdef min
#undef min
#endif 

template <class A>
inline A min(const A &m_a, const A &m_b)
{
	return (m_a < m_b) ? m_a : m_b;
}

#ifdef max
#undef max
#endif

template <class A>
inline A max(const A &m_a, const A &m_b)
{
	return (m_a > m_b) ? m_a : m_b;
}

#endif // WIN32

/////////////////////////////////////////////////////////
// By default, use an integer based alignment
#ifdef FLOATING_POINT_ALIGNMENT
	typedef float SO_Score;
#else // Integer alignment
	typedef int SO_Score;
#endif // FLOATING_POINT_ALIGNMENT

// The definition of the amino acids comes *before* the definition of the nucleic acids
// to allow the nucleic acid gap value to be the same as the amino acid gap value.
namespace AA {
	
	// Allowed amino acid bases
	// The actual values must match the score matrix arrangement
	typedef enum {
		A = 0, // Alanine, Ala
		R, // Arginine, Arg
		N, // Asparagine, Asn
		D, // Aspartic acid, Asp
		C, // Cysteine, Cys
		Q, // Glutamine, Gln
		E, // Glutamic acid, Glu
		G, // Glycine, Gly
		H, // Histidine, His
		I, // Isoleucine, Ile
		L, // Leucine, Leu
		K, // Lysine, Lys
		M, // Methionine, Met
		F, // Phenylalanine, Phe
		P, // Proline, Pro
		S, // Serine, Ser
		T, // Threonine, Thr
		W, // Tryptophan, Trp
		Y, // Tyrosine, Tyr
		V, // Valine, Val
		B, // Aspartic acid or Asparagine, Asx
		Z, // Glutamine or Glutamic acid, Glx
		X, // Any amino acid, Xaa
		GAP = (1 << 5) // Make sure that no bits overlap with any nucleic acid

	} amino_acid;
	
	const size_t NUM_AA = (X + 1); // (20 aa + B, Z and X)
	const size_t MATRIX_SIZE = NUM_AA*NUM_AA; // (20 aa + B, Z and X) by (20 aa + B, Z and X)
}

namespace NA {
	
	// Allowed nucleic acid bases
	//	A = adenine
	//	C = cytosine
	//	G = guanine
	//	T = thymine
	//	M = A or C
	//	R = G or A
	//	S = G or C
	//	V = G or C or A
	//	W = A or T
	//	Y = T or C
	//	H = A or C or T
	//	K = G or T
	//	D = G or A or T
	//	B = G or T or C
	//	N = A or T or C or G
	typedef enum {
		A = (1 << 0), 
		C = (1 << 1), 
		G = (1 << 2), 
		T = (1 << 3), 
		M = (A | C),
		R = (G | A),
		S = (G | C),
		V = (G | C | A),
		W = (A | T),
		Y = (T | C),
		H = (A | C | T),
		K = (G | T),
		D = (G | A | T),
		B = (G | T | C),
		N = (A | T | C | G),
		// It make life easier to share the same value for a GAP between
		// nucleic acids and amino acids
		GAP = AA::GAP
	} nucleic_acid;
}

// Both base types (nucleic acid and amino acid) must fit
// into a variable of size base_type
typedef unsigned char base_type;

class SeqOverlap{

	public:
	
	typedef enum {SmithWaterman, Overlap} AlignmentMode;
	
	private:
		
		struct SO_Elem{

			SO_Elem()
			{
				// Do nothing
			};

			~SO_Elem()
			{
				// Do nothing
			};

			inline void clear()
			{
				// Pick a default score that is lower than any score we will encounter
				M = I_query = I_target = SO_Score(-99999999);
				M_start = I_query_start = I_target_start = std::make_pair(-1, -1);
			};
			
			// Use the notation of "Biological sequence analysis" by Durbin, Eddy, Krogh and Mitchison
			SO_Score M;
			SO_Score I_query;  // insertion in query
			SO_Score I_target; // insertion in target
			
			#ifdef _COMPUTE_NUM_MATCH
			// The number of base matches in the alignment
			unsigned int M_num_match; 	
			unsigned int I_query_num_match;
			unsigned int I_target_num_match;
			#endif // _COMPUTE_NUM_MATCH
						
			// The actual alignment is not needed, just the score and the
			// coordinates of the begining and end of the alignment. For this purpose, propagate
			// the alignment start during the dynamic programming.
			
			// The starting row and column of the alignment
			std::pair<int, int> M_start;
			std::pair<int, int> I_query_start;
			std::pair<int, int> I_target_start;
			
			inline std::pair<int, int> start() const
			{
				// DEBUG
				return M_start;
			};
		};
		
		AlignmentMode mode;
		
		SO_Elem *last_row;
		SO_Elem *curr_row;
		SO_Elem max_elem;
		size_t dp_num_elem;
		
		std::pair<int, int> stop; // The location of the max element
		
		SO_Score match;
		SO_Score mask;
		SO_Score mismatch;
		SO_Score gap_existance;
		SO_Score gap_extension;

		// The input sequences in 5'-3' orientation
		std::vector<base_type> query;
		std::vector<base_type> target;
		
		bool is_na; // If is_na == false, then we are using amino acids
		
		// If mask_N_na == true, then treat 'N' as a masking character that gets a score of mask
		bool mask_N_na;
		
		std::vector<SO_Score> aa_score_matrix;
		
		void init_BLOSUM62(std::vector<SO_Score> &m_matrix);
		
		inline char bits_to_na(const base_type &m_bits) const
		{
			switch(m_bits){
				case NA::A:
					return 'A';
				case NA::C:
					return 'C';
				case NA::G:
					return 'G';
				case NA::T:
					return 'T';
				case NA::M:
					return 'M';
				case NA::R:
					return 'R';
				case NA::S:
					return 'S';
				case NA::V:
					return 'V';
				case NA::W:
					return 'W';
				case NA::Y:
					return 'Y';
				case NA::H:
					return 'H';
				case NA::K:
					return 'K';
				case NA::D:
					return 'D';
				case NA::B:
					return 'B';
				case NA::N:
					return 'N';
				case NA::GAP:
					return '-';
			};

			throw __FILE__ ":bit_to_na: Unknown base!";
			return 'X'; // Keep the compiler happy
		};
		
		inline char bits_to_aa(const base_type &m_bits) const
		{
			switch(m_bits){
				case AA::A:
					return 'A';
				case AA::R:
					return 'R';
				case AA::N:
					return 'N';
				case AA::D:
					return 'D';
				case AA::C:
					return 'C';
				case AA::Q:
					return 'Q';
				case AA::E:
					return 'E';
				case AA::G:
					return 'G';
				case AA::H:
					return 'H';
				case AA::I:
					return 'I';
				case AA::L:
					return 'L';
				case AA::K:
					return 'K';
				case AA::M:
					return 'M';
				case AA::F:
					return 'F';
				case AA::P:
					return 'P';
				case AA::S:
					return 'S';
				case AA::T:
					return 'T';
				case AA::W:
					return 'W';
				case AA::Y:
					return 'Y';
				case AA::V:
					return 'V';
				case AA::B:
					return 'B';
				case AA::Z:
					return 'Z';
				case AA::X:
					return 'X';
				case AA::GAP:
					return '-';
			};

			throw __FILE__ ":bit_to_aa: Unknown base!";
			return '?'; // Keep the compiler happy
		};

		void translate_na_seq(const std::string &m_input, std::vector<base_type> &m_output)
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":translate_na_seq: Only valid for NA sequences";
			}
			
			if( m_input.size() != m_output.size() ){
				throw __FILE__ ":translate_na_seq: input/output size mismatch";
			}
			#endif // _DEBUG
			
			std::string::const_iterator i = m_input.begin();
			std::vector<base_type>::iterator o = m_output.begin();

			for(;i != m_input.end();i++, o++){
			
				switch(*i){
					case 'A': case 'a':
						*o = NA::A;
						break;
					case 'T': case 't':
						*o = NA::T;
						break;
					case 'G': case 'g':
						*o = NA::G;
						break;
					case 'C': case 'c':
						*o = NA::C;
						break;
					case 'M': case 'm':
						*o = NA::M;
						break;
					case 'R': case 'r':
						*o = NA::R;
						break;
					case 'S': case 's':
						*o = NA::S;
						break;
					case 'V': case 'v':
						*o = NA::V;
						break;
					case 'W': case 'w':
						*o = NA::W;
						break;
					case 'Y': case 'y':
						*o = NA::Y;
						break;
					case 'H': case 'h':
						*o = NA::H;
						break;
					case 'K': case 'k':
						*o = NA::K;
						break;
					case 'D': case 'd':
						*o = NA::D;
						break;
					case 'B': case 'b':
						*o = NA::B;
						break;
					case 'N': case 'n':
					case 'I': case 'i': // For now, treat inosine as an 'N'
						*o = NA::N;
						break;
					default:
						throw __FILE__ ":translate_na_seq: Illegal base";
						break;
				};
			}
		};

		void translate_aa_seq(const std::string &m_input, std::vector<base_type> &m_output)
		{
			#ifdef _DEBUG
			if(is_na){
				throw __FILE__ ":translate_aa_seq: Only valid for AA sequences";
			}
			
			if( m_input.size() != m_output.size() ){
				throw __FILE__ ":translate_aa_seq: input/output size mismatch";
			}
			#endif // _DEBUG
			
			std::string::const_iterator i = m_input.begin();
			std::vector<base_type>::iterator o = m_output.begin();

			for(;i != m_input.end();i++, o++){
			
				switch(*i){
					case 'A': case 'a':
						*o = AA::A;
						break;
					case 'R': case 'r':
						*o = AA::R;
						break;
					case 'N': case 'n':
						*o = AA::N;
						break;
					case 'D': case 'd':
						*o = AA::D;
						break;
					case 'C': case 'c':
						*o = AA::C;
						break;
					case 'Q': case 'q':
						*o = AA::Q;
						break;
					case 'E': case 'e':
						*o = AA::E;
						break;
					case 'G': case 'g':
						*o = AA::G;
						break;
					case 'H': case 'h':
						*o = AA::H;
						break;
					// J = I or L. Since there is no 'J' in the
					// BLOSUM62 matrix, approximate J = I
					case 'J': case 'j':
					case 'I': case 'i':
						*o = AA::I;
						break;
					case 'L': case 'l':
						*o = AA::L;
						break;
					case 'K': case 'k':
						*o = AA::K;
						break;
					case 'M': case 'm':
						*o = AA::M;
						break;
					case 'F': case 'f':
						*o = AA::F;
						break;
					case 'P': case 'p':
						*o = AA::P;
						break;
					case 'S': case 's':
						*o = AA::S;
						break;
					case 'T': case 't':
						*o = AA::T;
						break;
					case 'W': case 'w':
						*o = AA::W;
						break;
					case 'Y': case 'y':
						*o = AA::Y;
						break;
					case 'V': case 'v':
						*o = AA::V;
						break;
					case 'B': case 'b':
						*o = AA::B;
						break;
					case 'Z': case 'z':
						*o = AA::Z;
						break;
					case 'X': case 'x':
					case 'U': case 'u': // Treat selenocysteine as 'X' (same as BLAST)
						*o = AA::X;
						break;
					default:
						throw __FILE__ ":translate_aa_seq: Illegal base";
						break;
				};
			}
		};
		
		void reverse_complement(std::vector<base_type> &m_seq) const
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":reverse_complement: Only valid for NA sequences";
			}
			#endif // _DEBUG
			
			size_t len = m_seq.size();

			len = len/2 + (len%2);

			std::vector<base_type>::iterator f = m_seq.begin();
			std::vector<base_type>::reverse_iterator r = m_seq.rbegin();

			for(size_t i = 0;i < len;i++, f++, r++){
				
				const base_type tmp = complement(*f);
				*f = complement(*r);
				*r = tmp;
			}
		};

		inline base_type complement(const base_type &m_base) const
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":complement: Only valid for NA sequences";
			}
			#endif // _DEBUG
			
			switch(m_base){
				case NA::A:
					return NA::T;
				case NA::C:
					return NA::G;
				case NA::G:
					return NA::C;
				case NA::T:
					return NA::A;
				case NA::M:
					return NA::K;
				case NA::R:
					return NA::Y;
				case NA::S:
					return NA::S;
				case NA::V:
					return NA::B;
				case NA::W:
					return NA::W;
				case NA::Y:
					return NA::R;
				case NA::H:
					return NA::D;
				case NA::K:
					return NA::M;
				case NA::D:
					return NA::H;
				case NA::B:
					return NA::V;
				case NA::N:
					return NA::N;
				case NA::GAP:
					return NA::GAP;
			};

			throw __FILE__ ":complement: Unknown base";
			return NA::GAP; // Keep the compiler happy
		};

		void align_overlap();
		void align_smith_waterman();
		
		SO_Score trace_back();

	public:
		
		SeqOverlap(const AlignmentMode &m_mode, const bool &m_is_na);
			
		~SeqOverlap()
		{
			if(last_row != NULL){
				delete [] last_row;
			}
			
			if(curr_row != NULL){
				delete [] curr_row;
			}
		};

		SO_Score align()
		{
			switch(mode){
				case Overlap:
					align_overlap();
					break;
				case SmithWaterman:
					align_smith_waterman();
					break;
				default:
					throw __FILE__ ":align: Unknown mode";
			};
			
			return score();
		}

		inline void clear()
		{
			query.clear();
			target.clear();

			if(last_row != NULL){
				delete [] last_row;
			}
			
			if(curr_row != NULL){
				delete [] curr_row;
			}
		};
		
		inline void clear_query()
		{
			query.clear();
		};
		
		inline void clear_target()
		{
			target.clear();
		};
		
		inline void set_query(const std::string &m_query)
		{
			std::vector<base_type>( m_query.size() ).swap(query);
			
			if(is_na){
				translate_na_seq(m_query, query);
			}
			else{
				translate_aa_seq(m_query, query);
			}
		};
		
		inline void set_query_reverse_complement(const std::string &m_query)
		{
			if(!is_na){
				throw __FILE__ ":set_query_reverse_complement: not defined for AA sequences";
			}
			
			std::vector<base_type>( m_query.size() ).swap(query);
			
			translate_na_seq(m_query, query);

			reverse_complement(query);
		};
		
		inline void set_target(const std::string &m_target)
		{
			std::vector<base_type>( m_target.size() ).swap(target);
			
			if(is_na){
				translate_na_seq(m_target, target);
			}
			else{
				translate_aa_seq(m_target, target);
			}
		};
		
		inline void set_target_reverse_complement(const std::string &m_target)
		{
			if(!is_na){
				throw __FILE__ ":set_target_reverse_complement: not defined for AA sequences";
			}
			
			std::vector<base_type>( m_target.size() ).swap(target);
			
			translate_na_seq(m_target, target);

			reverse_complement(target);
		};
		
		inline size_t size_query() const
		{
			return query.size();
		};
		
		inline size_t size_target() const
		{
			return target.size();
		};
			
		// The coordinates of the first and last aligned base in the query
		inline std::pair<int, int> alignment_range_query() const
		{
			return std::make_pair(max_elem.start().first, stop.first);
		};
		
		// The coordinates of the first and last aligned base in the query
		inline std::pair<int, int> alignment_range_target() const
		{
			return std::make_pair(max_elem.start().second, stop.second);
		};
		
		inline void alignment_range(std::pair<int, int> &m_query_range,
			std::pair<int, int> &m_target_range) const
		{
		
			m_query_range = alignment_range_query();
			m_target_range = alignment_range_target();
		};
		
		// Return the query and target sequences
		std::string query_seq() const;
		std::string target_seq() const;
		
		inline SO_Score score() const
		{
			return max_elem.M;
		};
		
		#ifdef _COMPUTE_NUM_MATCH
		inline unsigned int num_match() const
		{
			return max_elem.M_num_match;
		};
		#endif // _COMPUTE_NUM_MATCH
		
		inline void set_mode(const AlignmentMode &m_mode)
		{
			mode = m_mode;
		};
		
		inline void enable_mask_N_na()
		{
			mask_N_na = true;
		};
		
		inline void disable_mask_N_na()
		{
			mask_N_na = false;
		};
};

} // namespace::SO

#endif // __SEQ_OVERLAP
