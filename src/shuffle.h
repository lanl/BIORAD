#ifndef __RANDOM_SHUFFLE
#define __RANDOM_SHUFFLE

#include <stdlib.h>

// A random_shuffle-like function that uses drand48() -- which is NOT a thread safe random number generator!!
template <class T>
void randomize(const T &m_begin, const T &m_end)
{
	const size_t len = m_end - m_begin;
	
	for(size_t i = 0;i < len;++i){
	
		// Generate a [non-thread safe] random number between [0, len)
		size_t index = size_t( drand48()*len );
		
		while(index == len){
			index = size_t( drand48()*len );
		}
		
		std::swap( *(m_begin + i), *(m_begin + index) );
	}
}

#endif // __RANDOM_SHUFFLE
