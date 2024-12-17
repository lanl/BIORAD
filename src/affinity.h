#ifndef __AFFINITY
#define __AFFINITY

// Forward declaration of mpi helper functions to keep the compiler happy
template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

struct Affinity
{
    typedef enum {AFFINITY_KD, AFFINITY_KI, AFFINITY_IC50, AFFINITY_SCORE, AFFINITY_UNKNOWN, NUM_AFFINITY} AffinityType;
    typedef enum {MEASURE_EQUAL, MEASURE_GREATER, MEASURE_LESS, MEASURE_APPROXIMATE, MEASURE_UNKNOWN, NUM_MEASURE} MeasurementType;
    typedef enum {STRUCTURE_X_RAY, STRUCTURE_NMR, STRUCTURE_PREDICTED, STRUCTURE_UNKNOWN, NUM_STRUCTURE} StructureType;

    // Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define BIORAD_AFFINITY_MEMBERS \
        VARIABLE(AffinityType, affinity_type) \
        VARIABLE(MeasurementType, measurement_type) \
        VARIABLE(StructureType, structure_type) \
        VARIABLE(unsigned int, year) \
        VARIABLE(float, structure_resolution) \
        VARIABLE(float, structure_ph) \
        VARIABLE(float, structure_temperature) \
        VARIABLE(float, value)

    #define VARIABLE(A, B) A B;
		BIORAD_AFFINITY_MEMBERS
	#undef VARIABLE

    Affinity()
    {
        affinity_type = AFFINITY_UNKNOWN;
        measurement_type = MEASURE_UNKNOWN;
        structure_type = STRUCTURE_UNKNOWN;
        year = 0;
        structure_resolution = -1.0;
        structure_ph = -1.0;
        structure_temperature = -1.0;
        value = 0.0;
    };

    Affinity(const std::string &m_affinity_type, const std::string &m_measure_type, 
        const std::string &m_structure_type, const unsigned int &m_year, 
        const std::string &m_structure_ph, const std::string &m_structure_temperature, 
        const float &m_value)
    {
        if(m_affinity_type == "Kd"){
            affinity_type = AFFINITY_KD;
        }
        else if(m_affinity_type == "Ki"){
            affinity_type = AFFINITY_KI;
        }
        else if(m_affinity_type == "IC50"){
            affinity_type = AFFINITY_IC50;
        }
        else if(m_affinity_type == "score"){
            affinity_type = AFFINITY_SCORE;
        }
        else{
            throw __FILE__ ":Affinity: Unknown affinity type";
        }

        if(m_measure_type == "="){
            measurement_type = MEASURE_EQUAL;
        }
        else if(m_measure_type == ">"){
            measurement_type = MEASURE_GREATER;
        }
        else if(m_measure_type == "<"){
            measurement_type = MEASURE_LESS;
        }
        else if(m_measure_type == "~"){
            measurement_type = MEASURE_APPROXIMATE;
        }
        else{
            throw __FILE__ ":Affinity: Unknown measure type";
        }

        if(m_structure_type.find("X-RAY") != std::string::npos){
            
            structure_type = STRUCTURE_X_RAY;

            std::string::size_type loc = m_structure_type.find('=');

            if(loc == std::string::npos){
                throw __FILE__ ":Affinity: Unable to find x-ray resolution delimeter";
            }

            structure_resolution = atof( m_structure_type.substr(loc + 1, m_structure_type.size() - (loc + 1) ).c_str() );
        }
        else if(m_structure_type.find("NMR") != std::string::npos){
            structure_type = STRUCTURE_NMR;
        }
        else{
            throw __FILE__ ":Affinity: Unknown structure type";
        }

        year = m_year;

        if(m_structure_ph == ""){
            structure_ph = -1.0;
        }
        else{
            structure_ph = atof( m_structure_ph.c_str() );
        }

        if(m_structure_temperature == ""){
            structure_temperature = -1.0;
        }
        else{
            structure_temperature = atof( m_structure_temperature.c_str() );

            if(structure_temperature < 200){ // Centigrade to Kelvin
                structure_temperature += 273.15;
            }
        }

        value = m_value;
    };

    inline bool operator<(const Affinity &m_rhs) const
    {
        return (value < m_rhs.value);
    };
};

template<> size_t mpi_size(const Affinity &m_opt);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Affinity &m_opt);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Affinity &m_opt);

#endif // __AFFINITY