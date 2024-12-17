#ifndef __RANDOM_SHUFFLE
#define __RANDOM_SHUFFLE

#include <stdlib.h>

// A random_shuffle-like function that uses the Linux-specific, re-entrant random number generator
template <class T>
void randomize(const T &m_begin, const T &m_end, struct drand48_data *m_rand_ptr)
{
	if(m_rand_ptr == NULL){
		throw __FILE__ ":randomize: Invalid state pointer";
	}

	const size_t len = m_end - m_begin;
	
	for(size_t i = 0;i < len;++i){
	
		double r;
		
		// Note that we're using drand48_r() and not drand48() to make the ranomdize() function
		// thread safe. The previous version used "double erand48(unsigned short xsubi[3])" as I was
		// under the mistaken impression that the xsubi[3] variable contained all of the state information
		// required to make erand48 thread-safe. However, according to the man page for erand48(), this
		// function is *not* thread safe!
		drand48_r(m_rand_ptr, &r);
		
		// Generate a random number between [0, len)
		size_t index = size_t( r*len );
		
		while(index == len){
			
			drand48_r(m_rand_ptr, &r);
		
			index = size_t( r*len );
		}
		
		std::swap( *(m_begin + i), *(m_begin + index) );
	}
}

#endif // __RANDOM_SHUFFLE
