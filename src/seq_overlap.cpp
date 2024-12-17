#include "seq_overlap.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace SO;

SeqOverlap::SeqOverlap(const AlignmentMode &m_mode, const bool &m_is_na) :
	is_na(m_is_na)
{	
	mode = m_mode;
	
	mask_N_na = false;
	
	last_row = curr_row = NULL;
	dp_num_elem = 0;
	
	max_elem.clear();
	
	if(is_na){
	
		// These are blastn defaults, note that other nucleotide alignment algoirhtms (like FASTA, megablast, etc.) 
		// use difference paramters
		match = 2;
		mask = 0;
		mismatch = -3;
		gap_existance = 5;
		gap_extension = 2;
	}
	else{
		// Gaps costs from BLASTP default parameters
		gap_existance = 11;
		gap_extension = 1;
		
		// There are number of possible protein score matricies: 
		// PAM30, PAM70,  BLOSUM45, BLOSUM62, BLOSUM80
		// By default, BLASTP used BLOSUM62
		init_BLOSUM62(aa_score_matrix);
	}
};

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqOverlap::align_overlap()
{
	// Add one to allow for an initial row and column of zeros
	const size_t query_size = query.size();
	const size_t target_size = target.size();
	
	const size_t num_elem = target_size + 1;

	// Resize the dynamic programming matrix if needed
	if( !last_row || !curr_row || (num_elem > dp_num_elem) ){

		if(last_row){
			delete [] last_row;
		}
		
		if(curr_row){
			delete [] curr_row;
		}
		
		dp_num_elem = num_elem;
		
		last_row = new SO_Elem [dp_num_elem];

		if(!last_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for last_row";
		}
		
		curr_row = new SO_Elem [dp_num_elem];

		if(!curr_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for curr_row";
		}
	}

	// Initialize the dynamic programming matrix to have zeros along the first row 
	for(size_t j = 0;j <= target_size;j++){

		SO_Elem &elem_ref = last_row[j];

		elem_ref.M = SO_Score(0);
		elem_ref.I_query = elem_ref.I_target = -gap_existance;// - int(j)*gap_extension;
		
		elem_ref.M_start = make_pair(0, j);
		
		#ifdef _COMPUTE_NUM_MATCH
		elem_ref.M_num_match = elem_ref.I_query_num_match = elem_ref.I_target_num_match = 0;
		#endif // _COMPUTE_NUM_MATCH
	}

	// Reset the maximum element
	max_elem.clear();
	
	vector<base_type>::const_iterator query_iter = query.begin();
	
	for(size_t i = 0;i < query_size;i++, query_iter++){
		
		// Initialize the dynamic programming matrix to have zeros along the first column
		SO_Elem &elem_ref = curr_row[0];

		elem_ref.M = SO_Score(0);
		elem_ref.I_query = elem_ref.I_target = -gap_existance;// - int(i)*gap_extension;
		
		elem_ref.M_start = make_pair(i + 1, 0);
		
		#ifdef _COMPUTE_NUM_MATCH
		elem_ref.M_num_match = elem_ref.I_query_num_match = elem_ref.I_target_num_match = 0;
		#endif // _COMPUTE_NUM_MATCH
		
		// The dp matrix has query_size rows and target_size columns
		// A B
		// C X <-- dp[i][j]
		SO_Elem *A_ptr = last_row;
		SO_Elem *B_ptr = A_ptr + 1;
		SO_Elem *C_ptr = curr_row;
		SO_Elem *X_ptr = C_ptr + 1;

		vector<base_type>::const_iterator target_iter = target.begin();
		
		for(size_t j = 0;j < target_size;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++, target_iter++){
			
			// Match or mismatch
			const SO_Score &s1 = A_ptr->M;

			// gap the query
			const SO_Score &s2 = A_ptr->I_query;

			// gap the target
			const SO_Score &s3 = A_ptr->I_target;
			
			// All three possible outcomes will be added to same match or mismatch score
			if(is_na){
			
				if(*query_iter & *target_iter){
				
					if(mask_N_na && ( (*query_iter == NA::N) || (*target_iter == NA::N) ) ){
						
						X_ptr->M = mask;
						
						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match = 0;
						#endif // _COMPUTE_NUM_MATCH
					}
					else{
						X_ptr->M = match;

						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match = 1;
						#endif // _COMPUTE_NUM_MATCH
					}
				}
				else{
					X_ptr->M = mismatch;

					#ifdef _COMPUTE_NUM_MATCH
					X_ptr->M_num_match = 0;
					#endif // _COMPUTE_NUM_MATCH
				}
			}
			else{
				X_ptr->M = aa_score_matrix[*query_iter*AA::NUM_AA + *target_iter];
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->M_num_match = (*query_iter == *target_iter) ? 1 : 0;
				#endif // _COMPUTE_NUM_MATCH
			}
			
			if(s1 >= s2){
				
				if(s1 >= s3){
					
					X_ptr->M += s1;
					X_ptr->M_start = A_ptr->M_start;
					
					#ifdef _COMPUTE_NUM_MATCH
					X_ptr->M_num_match += A_ptr->M_num_match;
					#endif // _COMPUTE_NUM_MATCH
				}
				else{ // s1 < s3
				
					X_ptr->M += s3;
					X_ptr->M_start = A_ptr->I_target_start;
					
					#ifdef _COMPUTE_NUM_MATCH
					X_ptr->M_num_match += A_ptr->I_target_num_match;
					#endif // _COMPUTE_NUM_MATCH
				}
			}
			else{ // s1 < s2
				
				if(s2 >= s3){
				
					X_ptr->M += s2;
					X_ptr->M_start = A_ptr->I_query_start;
					
					#ifdef _COMPUTE_NUM_MATCH
					X_ptr->M_num_match += A_ptr->I_query_num_match;
					#endif // _COMPUTE_NUM_MATCH
				}
				else{ // s2 < s3
				
					X_ptr->M += s3;
					X_ptr->M_start = A_ptr->I_target_start;
					
					#ifdef _COMPUTE_NUM_MATCH
					X_ptr->M_num_match += A_ptr->I_target_num_match;
					#endif // _COMPUTE_NUM_MATCH
				}
			}
			
			SO_Score insert_gap = C_ptr->M - gap_existance;
			SO_Score extend_gap = C_ptr->I_query - gap_extension;

			if(insert_gap >= extend_gap){
			
				X_ptr->I_query = insert_gap;
				X_ptr->I_query_start = C_ptr->M_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_query_num_match = C_ptr->M_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			else{
			
				X_ptr->I_query = extend_gap;
				X_ptr->I_query_start = C_ptr->I_query_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_query_num_match = C_ptr->I_query_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			
			insert_gap = B_ptr->M - gap_existance;
			extend_gap = B_ptr->I_target - gap_extension;
						
			if(insert_gap >= extend_gap){

				X_ptr->I_target = insert_gap;
				X_ptr->I_target_start = B_ptr->M_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_target_num_match = B_ptr->M_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			else{
			
				X_ptr->I_target = extend_gap;
				X_ptr->I_target_start = B_ptr->I_target_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_target_num_match = B_ptr->I_target_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			
			// Is this the last row or last column?
			if( ( i == (query_size - 1) ) || (j == (target_size - 1) ) ){
				
				if(X_ptr->M >= max_elem.M){
				
					max_elem = *X_ptr;
					stop = make_pair(i, j); // Save the location of the maximum element
				}
			}
		}
		
		// Swap the last_row and curr_row
		swap(last_row, curr_row);
	}
}

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqOverlap::align_smith_waterman()
{
	// Add one to allow for an initial row and column of zeros
	const size_t query_size = query.size();
	const size_t target_size = target.size();
	
	const size_t num_elem = target_size + 1;

	// Resize the dynamic programming matrix if needed
	if( !last_row || !curr_row || (num_elem > dp_num_elem) ){

		if(last_row){
			delete [] last_row;
		}
		
		if(curr_row){
			delete [] curr_row;
		}
		
		dp_num_elem = num_elem;
		
		last_row = new SO_Elem [dp_num_elem];

		if(!last_row){
			throw __FILE__ ":SeqOverlap::align_smith_waterman: Unable to allocate memory for last_row";
		}
		
		curr_row = new SO_Elem [dp_num_elem];

		if(!curr_row){
			throw __FILE__ ":SeqOverlap::align_smith_waterman: Unable to allocate memory for curr_row";
		}
	}

	// Initialize the dynamic programming matrix to have zeros along the first row 
	for(size_t j = 0;j <= target_size;j++){

		SO_Elem &elem_ref = last_row[j];

		elem_ref.M = SO_Score(0);
		elem_ref.I_query = elem_ref.I_target = -gap_existance;// - int(j)*gap_extension;
		
		elem_ref.M_start = make_pair(0, j);
		
		#ifdef _COMPUTE_NUM_MATCH
		elem_ref.M_num_match = elem_ref.I_query_num_match = elem_ref.I_target_num_match = 0;
		#endif // _COMPUTE_NUM_MATCH
	}

	// Reset the maximum element
	max_elem.clear();
	
	vector<base_type>::const_iterator query_iter = query.begin();
	
	for(size_t i = 0;i < query_size;i++, query_iter++){
		
		// Initialize the dynamic programming matrix to have zeros along the first column
		SO_Elem &elem_ref = curr_row[0];

		elem_ref.M = SO_Score(0);
		elem_ref.I_query = elem_ref.I_target = -gap_existance;// - int(i)*gap_extension;
		
		elem_ref.M_start = make_pair(i + 1, 0);
		
		#ifdef _COMPUTE_NUM_MATCH
		elem_ref.M_num_match = elem_ref.I_query_num_match = elem_ref.I_target_num_match = 0;
		#endif // _COMPUTE_NUM_MATCH
		
		// The dp matrix has query_size rows and target_size columns
		// A B
		// C X <-- dp[i][j]
		SO_Elem *A_ptr = last_row;
		SO_Elem *B_ptr = A_ptr + 1;
		SO_Elem *C_ptr = curr_row;
		SO_Elem *X_ptr = C_ptr + 1;

		vector<base_type>::const_iterator target_iter = target.begin();
		
		for(size_t j = 0;j < target_size;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++, target_iter++){
			
			SO_Score s1 = A_ptr->M;
			SO_Score s2 = A_ptr->I_query;
			SO_Score s3 = A_ptr->I_target;
			
			bool s1_new_start = false;
			bool s2_new_start = false;
			bool s3_new_start = false;
			
			// Match or mismatch. Note that we only re-start the alignment when the score is *less* than
			// zero. Scores that are *equal* to zero keep their original start (i.e. always save the
			// longest sequence alignment when there are multiple, equal scoring alignments to choose
			// from).
			if(SO_Score(0) > A_ptr->M){
			
				s1 = SO_Score(0);
				s1_new_start = true;
			}
			
			if(SO_Score(0) > A_ptr->I_query){
			
				s2 = SO_Score(0);
				s2_new_start = true;
			}
			
			if(SO_Score(0) > A_ptr->I_target){
			
				s3 = SO_Score(0);
				s3_new_start = true;
			}
				
			// All three possible outcomes will be added to same match or mismatch score
			if(is_na){
			
				if(*query_iter & *target_iter){
				
					if(mask_N_na && ( (*query_iter == NA::N) || (*target_iter == NA::N) ) ){
					
						X_ptr->M = mask;

						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match = 0;
						#endif // _COMPUTE_NUM_MATCH
					}
					else{
						X_ptr->M = match;

						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match = 1;
						#endif // _COMPUTE_NUM_MATCH
					}
				}
				else{
					X_ptr->M = mismatch;

					#ifdef _COMPUTE_NUM_MATCH
					X_ptr->M_num_match = 0;
					#endif // _COMPUTE_NUM_MATCH
				}
			}
			else{
				X_ptr->M = aa_score_matrix[*query_iter*AA::NUM_AA + *target_iter];
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->M_num_match = (*query_iter == *target_iter) ? 1 : 0;
				#endif // _COMPUTE_NUM_MATCH
			}
			
			if(s1 >= s2){
				
				if(s1 >= s3){
				
					X_ptr->M += s1;
					
					if(s1_new_start){
						X_ptr->M_start = make_pair(i, j);
					}
					else{
						X_ptr->M_start = A_ptr->M_start;
						
						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match += A_ptr->M_num_match;
						#endif // _COMPUTE_NUM_MATCH
					}
				}
				else{ // s1 < s3
				
					X_ptr->M += s3;
						
					if(s3_new_start){
						X_ptr->M_start = make_pair(i, j);
					}
					else{
						X_ptr->M_start = A_ptr->I_target_start;
						
						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match += A_ptr->I_target_num_match;
						#endif // _COMPUTE_NUM_MATCH
					}
				}
			}
			else{ // s1 < s2
				
				if(s2 >= s3){
				
					X_ptr->M += s2;
					
					if(s2_new_start){
						X_ptr->M_start = make_pair(i, j);
					}
					else{
						X_ptr->M_start = A_ptr->I_query_start;
						
						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match += A_ptr->I_query_num_match;
						#endif // _COMPUTE_NUM_MATCH
					}
				}
				else{ // s2 < s3
				
					X_ptr->M += s3;
						
					if(s3_new_start){
						X_ptr->M_start = make_pair(i, j);
					}
					else{
						X_ptr->M_start = A_ptr->I_target_start;
						
						#ifdef _COMPUTE_NUM_MATCH
						X_ptr->M_num_match += A_ptr->I_target_num_match;
						#endif // _COMPUTE_NUM_MATCH
					}
				}
			}
			
			SO_Score insert_gap = (SO_Score(0) < C_ptr->M) ? C_ptr->M - gap_existance :
				- gap_existance;
			SO_Score extend_gap = (SO_Score(0) < C_ptr->I_query) ? C_ptr->I_query - gap_extension :
				- gap_extension;
				
			if(insert_gap >= extend_gap){
			
				X_ptr->I_query = insert_gap;
				X_ptr->I_query_start = C_ptr->M_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_query_num_match = C_ptr->M_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			else{
			
				X_ptr->I_query = extend_gap;
				X_ptr->I_query_start = C_ptr->I_query_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_query_num_match = C_ptr->I_query_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			
			insert_gap = (SO_Score(0) < B_ptr->M) ? B_ptr->M - gap_existance :
				- gap_existance;
			extend_gap = (SO_Score(0) < B_ptr->I_target) ? B_ptr->I_target - gap_extension :
				- gap_extension;
						
			if(insert_gap >= extend_gap){

				X_ptr->I_target = insert_gap;
				X_ptr->I_target_start = B_ptr->M_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_target_num_match = B_ptr->M_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			else{
			
				X_ptr->I_target = extend_gap;
				X_ptr->I_target_start = B_ptr->I_target_start;
				
				#ifdef _COMPUTE_NUM_MATCH
				X_ptr->I_target_num_match = B_ptr->I_target_num_match;
				#endif // _COMPUTE_NUM_MATCH
			}
			
			
			// Is this the largest element we've found so far?
			if(X_ptr->M >= max_elem.M){

				max_elem = *X_ptr;
				stop = make_pair(i, j); // Save the location of the maximum element
			}
		}
		
		// Swap the last_row and curr_row
		swap(last_row, curr_row);
	}
}

string SeqOverlap::query_seq() const
{
	const size_t len = query.size();
	string tmp(len, '-');

	if(is_na){
		for(size_t i = 0;i < len;i++){
			tmp[i] = bits_to_na(query[i]);
		}
	}
	else{
		for(size_t i = 0;i < len;i++){
			tmp[i] = bits_to_aa(query[i]);
		}
	}
	
	return tmp;
};

string SeqOverlap::target_seq() const
{
	const size_t len = target.size();
	string tmp(len, '-');

	if(is_na){
		for(size_t i = 0;i < len;i++){
			tmp[i] = bits_to_na(target[i]);
		}
	}
	else{
		for(size_t i = 0;i < len;i++){
			tmp[i] = bits_to_aa(target[i]);
		}
	}
	
	return tmp;
};

void SeqOverlap::init_BLOSUM62(vector<SO_Score> &m_matrix)
{
	// From the NCBI blast source distribution: ~/ncbi/data/BLOSUM62
	const SO_Score matrix[AA::MATRIX_SIZE] = {
		4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1,
		-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1,
		-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1,
		-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1,
		0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -1,
		-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1,
		-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
		0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1,
		-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1,
		-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1,
		-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1,
		-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1,
		-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1,
		-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1,
		-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -1,
		1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, -1,
		0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1,
		-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -1,
		-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1,
		0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1,
		-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1,
		-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};

	m_matrix = vector<SO_Score>(matrix, matrix + AA::MATRIX_SIZE);
}

