#include "pdb.h"
#include "parse_util.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>

#include <fstream>
#include <regex>
#include <unordered_set>
#include <iostream> // debug

using namespace std;

#define     PATH_SEPARATOR      "/"

#define     PDB_EXTENSION       ".pdb"

#define AMINO_ACID_NAME(X)      #X,
const char* amino_acid_names[] = {
    AMINO_ACID_TYPES
};
#undef AMINO_ACID_NAME

string extract_pdb_accession(const string &m_filename);
bool parse_pdb_complex(const string &m_filename, PDBComplex &m_pdb, unordered_map<string, unsigned int> &m_atom_table,
    const bool &m_ignore_hetatom_chain, const bool &m_ignore_charge, const bool &m_force_dimer, ostream &m_out);
unsigned int parse_atom_type(const string &m_type, unordered_map<string, unsigned int> &m_atom_table, 
    const bool &m_ignore_charge);
unsigned int parse_amino_acid_type(const string &m_type);
vector<float> parse_biomt_col(const string &m_buffer);
bool check_atomic_distances(const string &m_accession, const PDBComplex &m_pdb, const float &m_threshold);
vector<int> split_conect(const string &m_buffer);

// Attempt to read all of the PDB files in m_dir
void parse_pdb_data(unordered_map<string, PDBComplex> &m_pdb_data, const string &m_dir,
    unordered_map<string, unsigned int> &m_atom_table, const bool &m_ignore_hetatom_chain, 
    const bool &m_ignore_charge, const bool &m_force_dimer, ostream &m_out)
{
    // Is the input path a directory?
    struct stat info;

    if( stat(m_dir.c_str(), &info) != 0){
        throw __FILE__ ":parse_pdb_data: Unable to stat PDB directory";
    }

    // Is this a directory
    if( !S_ISDIR(info.st_mode) ){
        throw __FILE__ ":parse_pdb_data: PDB directory is not valid";
    }

    // Read the contents of the directory
    DIR *dp = opendir( m_dir.c_str() );

    if(dp == NULL){
        throw __FILE__ ":parse_pdb_data: Unable to open PDB directory for reading";
    }

    struct dirent *d = NULL;

    while( ( d = readdir(dp) ) ){

        // Skip any removed files or diretories
        if(d->d_ino == 0){
                continue;
        }

        // Skip the special directories "." and ".."
        if( (strcmp(d->d_name, ".") == 0) ||
            (strcmp(d->d_name, "..") == 0) ){
                continue;
        }

        const string name = d->d_name;

        if( ( name.find(PDB_EXTENSION) + strlen(PDB_EXTENSION) )!= name.size() ){
            continue;
        }

        const string full_path = m_dir + PATH_SEPARATOR + name;

        if( stat(full_path.c_str(), &info) != 0){
            throw __FILE__ ":parse_pdb_data: Unable to stat PDB file";
        }

        if( S_ISREG(info.st_mode) ){
            
            const string accession = extract_pdb_accession(name);

            if( m_pdb_data.find(accession) != m_pdb_data.end() ){
                throw __FILE__ ":parse_pdb_data: DUplicate PDB accession!";
            }

            try{

                PDBComplex pdb;

                if( parse_pdb_complex(full_path, pdb, m_atom_table, m_ignore_hetatom_chain, m_ignore_charge, m_force_dimer, m_out) ){

                    // Check for atoms that are too close ...
                    //check_atomic_distances(accession, pdb, 1.0 /*Minimum allowed interatomic distance*/);

                    m_pdb_data[accession] = pdb;
                }
            }
            catch(const char* error){

                cerr << "Error parsing " << full_path << ": " << error << endl;
                throw __FILE__ ":parse_pdb_data: Unable to parse PDB file";
            }
            catch(...){

                cerr << "Unknown error parsing " << full_path << endl;
                throw __FILE__ ":parse_pdb_data: Unable to parse PDB file";
            }
        }
    }

    closedir(dp);
}

// Parse PDB file names to extract the PDB accession and, optionally, the
// record id number:
//  <accession>.pdb
//  or
//  <accession>.<id>.pdb
//
// If an id is present, it is appended to the accession: <accession>.<id>
string extract_pdb_accession(const string &m_filename)
{

    string accession;
    
    string::const_iterator i = m_filename.begin();

    while( (i != m_filename.end()) && isalnum(*i) ){

        accession.push_back(*i);
        ++i;
    }

    #ifdef ACCESSION_SUFFIX
    if( ( i != m_filename.end() ) && (*i == '.') ){
        
        ++i;

        accession.push_back('.');

        if( ( i != m_filename.end() ) && isdigit(*i) ){

            while( ( i != m_filename.end() ) && isdigit(*i) ){
                
                accession.push_back(*i);
                ++i;
            }
        }
        else{
            // Default id is 0
            accession.push_back('0');
        }
    }
    #endif // ACCESSION_SUFFIX

    return accession;
}

bool parse_pdb_complex(const string &m_filename, PDBComplex &m_pdb, unordered_map<string, unsigned int> &m_atom_table,
    const bool &m_ignore_hetatom_chain, const bool &m_ignore_charge, const bool &m_force_dimer, ostream &m_out)
{
    // Make a local copy of the atom table that will only be return to the calling function if this PDB file
    // is valid
    unordered_map<string, unsigned int> local_atom_table(m_atom_table);

    ifstream fin( m_filename.c_str() );

    if(!fin){
        
        cerr << "Unable to open PDB file: " << m_filename << endl;
        throw __FILE__ ":parse_pdb_complex: Unable to open PDB file";
    }

    string line;
    size_t line_number = 0;

    // While parsing a PDB file, we will store atom data in a nested map structure:
    // chain id (char) -> amino acid id (unsigned int) ->AminoAcidData
    unordered_map<char, unordered_map<unsigned int, AminoAcidData> > target_chains;

    // The amino acid sequence reported by the header. Used to double check the amino
    // acid sequence extracted by parsing the ATOM records
    unordered_map<char /*chain*/, deque<unsigned int> > seqres;

    // Track the connections between atoms to allow HETATM chains to merged with amino acid chains.
    // This is needed for modified amino acids (phosphorylation, glycosylation, acetylation, etc.).
    deque< pair<int, int> > conect; // Store CONECT records as: (min index, max index)
    unordered_map<int/*atom index*/, AtomData> hetatm;
    bool conect_overflow = false;
    size_t num_atoms = 0;

    regex biomolecule_regex("REMARK 350 BIOMOLECULE:\\s+(\\d+)\\s+$");
    regex chains_regex("REMARK 350 APPLY THE FOLLOWING TO CHAINS:\\s+(\\w),\\s+(\\w)\\s+$");

    // The vast majority of PDB files will have the target heterodimer as biomolecule 1. However, a few
    // outliers (e.g., 5ywr and 5dcq) may not ...
    int target_biomolecule = 1;
    bool biomolecule_header = false;
    bool biomolecule_dimer = false;
    bool found_remark_350 = false;

    // The chains that are the target of the current BIOMT rotation and translation
    deque<char> transformation_target_chains;

    // The per-chain biomolecule transformations
    unordered_map< char /*chain*/, deque<Transform> > biomt;

    Transform curr_transform;

    while( getline(fin, line) ){

        ++line_number;

        // Extract representative chain ids from BIOMOLECULE 1
        if(line.find("REMARK 350") == 0){

            found_remark_350 = true;

            smatch match;

            if(regex_search(line, match, biomolecule_regex)){
                
                // match[0] is the input string
                // match[1] is the first group match!
                biomolecule_header = ( std::stoi(match[1]) == target_biomolecule);
            }

            if(biomolecule_header){

                if( (line.find("REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC") == 0) || 
                    (line.find("REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC") == 0) ){

                    biomolecule_dimer = true;
                }

                string::size_type loc = line.find("REMARK 350 APPLY THE FOLLOWING TO CHAINS:");

                if( (loc == 0) && biomolecule_dimer ){

                    const size_t len = 41; // strlen("REMARK 350 APPLY THE FOLLOWING TO CHAINS:")
                    const vector<string> chains = split(line.substr(len, line.size() - len), ',');

                    #ifdef SKIP_SINGLE_CHAIN_TRANSFORMS
                    if(chains.size() == 1){

                        // We're looking for heterodimers, not single chains. A single chain likely indicates
                        // a heterodimer where the individual chains require separate biomolecular coordinate
                        // transformations. To avoid applying transformations (for now), check
                        // the next biomolecule for a valid dimer.
                        target_biomolecule += 1; // Test the next biomolecule
                        biomolecule_header = false;
                        biomolecule_dimer = false;
                        target_chains.clear();
                        continue;
                    }
                    #endif // SKIP_SINGLE_CHAIN_TRANSFORMS

                    // Don't remove existing chain information, as some proteins have biomolecule chains
                    // specified separately (to allow per-chain rotation and translation transformations)
                    //target_chains.clear();

                    transformation_target_chains.clear();

                    for(vector<string>::const_iterator i = chains.begin();i != chains.end();++i){

                        const string c = strip(*i);

                        if(c.size() == 1){

                            target_chains[c[0]] = unordered_map<unsigned int, AminoAcidData>();
                            transformation_target_chains.push_back(c[0]);
                        }
                    }
                }

                loc = line.find("REMARK 350   BIOMT1   ");

                if(loc == 0){

                    curr_transform.reset();

                    const size_t len = 22; // strlen("REMARK 350   BIOMT1   ")
                    const vector<float> data = parse_biomt_col( line.substr(len, line.size() - len) );

                    if(data.size() != 5){
                        throw __FILE__ ":parse_pdb_complex: Did not read the expected number of BIOMT1 columns";
                    }

                    curr_transform(0, 0) = data[1];
                    curr_transform(0, 1) = data[2];
                    curr_transform(0, 2) = data[3];

                    curr_transform[0] = data[4];
                }

                loc = line.find("REMARK 350   BIOMT2   ");

                if(loc == 0){

                    const size_t len = 22; // strlen("REMARK 350   BIOMT2   ")
                    const vector<float> data = parse_biomt_col( line.substr(len, line.size() - len) );

                    if(data.size() != 5){
                        throw __FILE__ ":parse_pdb_complex: Did not read the expected number of BIOMT2 columns";
                    }

                    curr_transform(1, 0) = data[1];
                    curr_transform(1, 1) = data[2];
                    curr_transform(1, 2) = data[3];

                    curr_transform[1] = data[4];
                }

                loc = line.find("REMARK 350   BIOMT3   ");

                if(loc == 0){

                    const size_t len = 22; // strlen("REMARK 350   BIOMT3   ")
                    const vector<float> data = parse_biomt_col( line.substr(len, line.size() - len) );

                    if(data.size() != 5){
                        throw __FILE__ ":parse_pdb_complex: Did not read the expected number of BIOMT3 columns";
                    }

                    curr_transform(2, 0) = data[1];
                    curr_transform(2, 1) = data[2];
                    curr_transform(2, 2) = data[3];

                    curr_transform[2] = data[4];

                    // Store this transformation for the current chains. Don't bother applying the identity transformation
                    // (i.e., no rotation and no translation)
                    if( !curr_transform.is_identity() ){

                        for(deque<char>::const_iterator c = transformation_target_chains.begin();c != transformation_target_chains.end();++c){
                            biomt[*c].push_back(curr_transform);
                        }
                    }
                }
            }
        }

        if(line.find("SEQRES") == 0){

            // Extract the amnio acid sequence of each chain as reported by the PDB file header.
            // This sequence will be used to double check our PDB file parsing!
            // Expected format: SEQRES  line_number chain_id    number_of_aa_in_chain   AA1 AA2 AA3 ...
            stringstream ssin(line);

            string tmp;

            if( !(ssin >> tmp) ){
                throw __FILE__ ":parse_pdb_complex: Error reading SEQRES header information";
            }

            if(tmp != "SEQRES"){
                throw __FILE__ ":parse_pdb_complex: Error reading SEQRES keyword";
            }

            if( !(ssin >> tmp) ){
                throw __FILE__ ":parse_pdb_complex: Error reading SEQRES line number";
            }

            if( !(ssin >> tmp) ){
                throw __FILE__ ":parse_pdb_complex: Error reading SEQRES chain id";
            }

            if(tmp.size() != 1){
                throw __FILE__ ":parse_pdb_complex: Did not read a single character SEQRES chain id";
            }

            const char chain_id = tmp[0];

            if( !(ssin >> tmp) ){
                throw __FILE__ ":parse_pdb_complex: Error reading SEQRES number of amino acids";
            }

            deque<unsigned int > &aa = seqres[chain_id];

            while(ssin >> tmp){
                aa.push_back( parse_amino_acid_type(tmp) );
            }
        }

        // Update: March 5, 2023: Don't terminate the chain when we encounter the "TER" symbol,
        // as some PDB files (i.e. 5mad) contain chains that stop and then restart after a non-standard
        // amino acid or chemical group!
        // Is this the end of a chain?
    
        // For NMR models, only load the first model
        if(line.find("ENDMDL") == 0){

            if( target_chains.size() == 2 ){
                break;
            }

            throw __FILE__ ":parse_pdb_complex: Encountered \"ENDMDL\" (end-of-model) without reading two chains";
        }

        if(line.find("CONECT") == 0){

            const vector<int> indicies = split_conect(line);

            if(indicies.size() <= 1){
                throw __FILE__ ":parse_pdb_complex: Found an empty CONECT record";
            }

            if(indicies.size() > 5){

                if(!conect_overflow){
                    conect_overflow = true;
                }
            }
            else{

                const int primary_index = indicies[0];

                for(vector<int>::const_iterator secondary = indicies.begin() + 1;secondary != indicies.end();++secondary){
                    // Store CONECT records as: (min index, max index)
                    conect.push_back( make_pair( min(primary_index, *secondary), max(primary_index, *secondary) ) );
                }
            }

            continue;
        }

        // Only consider lines that start with "ATOM" or "HETATM"
        const bool is_atom = (line.find("ATOM") == 0);
        const bool is_hetatom = (line.find("HETATM") == 0);

        if(!is_atom && !is_hetatom){
            continue;
        }

        // Some PDB files (like the predicted PDB file output by AlphaFold) do
        // not contain a "REMARK 350" header. In these cases, assume that
        // biomolecule_dimer is true  
        if( found_remark_350 && !m_force_dimer && !biomolecule_dimer ){

           // m_out << m_filename << " does not have a dimeric biomolecule" << endl;
            return false;
        }

        if(line.size() < 54){
            throw __FILE__ ":parse_pdb_complex: Line too short to parse!";
        }

        const char chain = line[21];

        unordered_map<char, unordered_map<unsigned int, AminoAcidData> >::iterator chain_iter = 
            target_chains.find(chain);

        ++num_atoms;

        if(is_atom){

            // There can be multiple copies of the protein complex in a single PDB file.
            // For now, only read the representative complex.
            //if(found_remark_350 && (target_chains.size() != 2) ){

            //    m_out << m_filename << " does not have two chains in BIOMOLECULE" << endl;
            //    return false;
            //}

            if(found_remark_350 && ( chain_iter == target_chains.end() ) ){
                continue;
            }

            if(!found_remark_350){

                // We need to build the target chains as we go
                if( chain_iter == target_chains.end() ){
                    if(target_chains.size() >= 2){
                        throw __FILE__ ":parse_pdb_complex: Found more than two chains!";
                    }
                    else{

                        // Make sure that chain_iter is valid if we just added a new chain
                        chain_iter = target_chains.insert( make_pair(chain, unordered_map<unsigned int, AminoAcidData>() ) ).first;
                    }
                }
            }

            if( chain_iter == target_chains.end() ){
                // This line does not contain information relevant to the chains of interest
                continue;
            }
        }

        // For the definition of PDB columns, see: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
        if(line.size() < 80){
            throw __FILE__ ":parse_pdb_complex: Unexpected short line";
        }

        // Multi-character element name (with or without charge information depending on m_ignore_charge)
        const unsigned int atom_type = parse_atom_type(line.substr(76, line.size() - 76), local_atom_table, m_ignore_charge);

        const int atom_index = atoi(line.substr(6, 5).c_str());
        const unsigned int amino_acid_type = parse_amino_acid_type(line.substr(17, 3));
        const int loc = atoi(line.substr(22, 4).c_str());

        const float x = atof(line.substr(29, 8).c_str());
        const float y = atof(line.substr(38, 8).c_str());
        const float z = atof(line.substr(46, 8).c_str());

        const float temperature_factor = atof(line.substr(60, 7).c_str());
        
        // DEBUG
        #define DUMP_DEBUG \
            cout << "-----------------------------------------" << endl; \
            cout << "pdb = " << m_filename << endl; \
            cout << "line number = " << line_number << endl; \
            cout << "atom_type = " << atom_type << endl; \
            cout << "amino_acid_type = " << amino_acid_type << endl; \
            cout << "chain = " << chain << endl; \
            cout << "loc = " << loc << endl; \
            cout << "x = " << x << endl; \
            cout << "y = " << y << endl; \
            cout << "z = " << z << endl; \
            cout << "|chain| = " << target_chains.size() << endl;

        const AtomData curr_atom(atom_type, atom_index, x, y, z, temperature_factor);

        if(is_atom){

            // Add the current atom to the current, or next, amino acid in the current protein chain
            unordered_map<unsigned int, AminoAcidData>::iterator aa_iter = chain_iter->second.find(loc);

            if( aa_iter == chain_iter->second.end() ){
                aa_iter = chain_iter->second.insert( make_pair( loc, AminoAcidData(amino_acid_type) ) ).first;
            }
            
            aa_iter->second.atoms.push_back(curr_atom);
        }
        else if(is_hetatom && !m_ignore_hetatom_chain){
            
            // Add the current atom to the HETATM table. If these atoms are connected (directly or indirectly) to a protein
            // chain then they will be added later.
            hetatm[atom_index] = curr_atom;
        }
    }

    if(num_atoms > 99999){
        conect_overflow = true;
    }

    if(conect_overflow){

        cerr << "Warning! Overflowed the CONECT records for " << m_filename << "; found " 
            << num_atoms << " atoms" << endl;

        conect.clear();
    }

    // Make the list of CONECT records unique
    sort( conect.begin(), conect.end() );
    conect.erase( unique( conect.begin(), conect.end() ), conect.end() );

    // Add hetam records that are connected to protein chains. To simpify this process, associate every protein atom
    // with the pointer to its parent AminoAcid data
    unordered_map<int, AminoAcidData*> index_to_ptr;
    
    for(unordered_map<char, unordered_map<unsigned int, AminoAcidData> >::iterator i = target_chains.begin();
        i != target_chains.end();++i){
        
        for(unordered_map<unsigned int, AminoAcidData>::iterator j = i->second.begin();j != i->second.end();++j){

            for(vector<AtomData>::iterator k = j->second.atoms.begin();k != j->second.atoms.end();++k){
                index_to_ptr[k->index] = &(j->second);
            }
        }
    }

    // Iteratively add HETATM records to protein chains that they are connected to
    while(true){

        // DEBUG
        //cerr << conect.size() << '\t' << m_filename << endl;

        // Build a new CONECT list
        deque< pair<int, int> > new_conect;

        for(deque< pair<int, int> >::const_iterator i = conect.begin();i != conect.end();++i){

            // Is the first index a protein index?
            unordered_map<int, AminoAcidData*>::iterator aa_iter = index_to_ptr.find(i->first);
            unordered_map<int, AtomData>::iterator hetatm_iter;

            if( aa_iter != index_to_ptr.end() ){ // First is protein index

                hetatm_iter = hetatm.find(i->second); // Look for the second index in the hetatm
            }
            else{ // First is HETATM index

                hetatm_iter = hetatm.find(i->first); // Look for the first in the hetatm index
                aa_iter = index_to_ptr.find(i->second); // Look for the second in the protein index
            }

            // Do we have a connection that involves both a HETATM and a protein?
            if( hetatm_iter != hetatm.end() ){

                // We have a HETATM ...
                if( aa_iter != index_to_ptr.end() ){

                    // ... and a protein atom!
                    // Add the HETATM to the parent amino acid
                    aa_iter->second->atoms.push_back(hetatm_iter->second);

                    // Update pointer table to capture connected HETATMs in future iterations
                    // (i.e. this atom is now a "protein" atom for the sake of recruiting additional
                    // connected HETATMs).
                    index_to_ptr[hetatm_iter->first] = aa_iter->second;

                    // Remove the newly added atom from the source HETATM list to prevent adding the same
                    // atom multiple time when there are multiple connections to different protein atoms.
                    // When multiple connections are present between a single HETATM and multiple
                    // protein atoms, only add the HETAM to the *first* listed protein atom.
                    hetatm.erase(hetatm_iter);
                }
                else{
                    // This is a potential connection between two HETAMs -- save it for the next iteration
                    new_conect.push_back(*i);
                }
            }
        }

        // Did we change the number of CONECT records?
        if( new_conect.size() == conect.size() ){
            break;
        }

        // Update the CONECT records for the next iteration
        conect.swap(new_conect);
    }

    // Remove any chains that do not have any associated atom information (as some chains can be comprised only of
    // HETATOM records)
    deque<unordered_map<char, unordered_map<unsigned int, AminoAcidData> >::iterator> reaper;

    for(unordered_map<char, unordered_map<unsigned int, AminoAcidData> >::iterator i = target_chains.begin();
        i != target_chains.end();++i){

        if( i->second.empty() ){
            reaper.push_back(i);
        }
    }

    while( !reaper.empty() ){
        
        target_chains.erase( reaper.back() );
        reaper.pop_back();
    }

    #ifdef SEQRES_COMPARISON
    // Compare the extracted amnio acid sequence data to the amino acid provided by the header (SEQRES).
    // Update: Unfortunately, this is of only limited utility, as the SEQRES information refers to the original
    // biological source sequence and is typically a *superset* of the actual protein sequence for which there
    // is structural information in the PDB file 
    unordered_map<char /*chain*/, deque<unsigned int> >::const_iterator seq_iter = seqres.find(assigned_chains.first);

    if( seq_iter == seqres.end() ){
        cerr << "Did not find SEQRES information for chain " << assigned_chains.first << " in " << m_filename << endl;
    }
    else{

        if( int(seq_iter->second.size()) - int(m_pdb.first.size()) > 0 ){
            cerr << m_filename << "\t" << assigned_chains.first << '\t' << seq_iter->second.size() << '\t' 
                << m_pdb.first.size() << endl;
        }
    }

    seq_iter = seqres.find(assigned_chains.second);

    if( seq_iter == seqres.end() ){
        cerr << "Did not find SEQRES information for chain " << assigned_chains.second << " in " << m_filename << endl;
    }
    else{

        if( int(seq_iter->second.size()) - int(m_pdb.second.size()) > 0 ){
            cerr << m_filename << "\t" << assigned_chains.second << '\t' << seq_iter->second.size() << '\t' 
                << m_pdb.second.size() << endl;
        }
    }
    #endif // SEQRES_COMPARISON

    if(target_chains.size() != 2){

        m_out << m_filename << " contains " << target_chains.size() << " protein chain(s):";

        for(unordered_map<char, unordered_map<unsigned int, AminoAcidData> >::const_iterator i = target_chains.begin();
            i != target_chains.end();++i){
            m_out << ' ' << i->first;
        }

        m_out << "; DNA or other non-protein chain?" << endl;

        return false;
    }

    unordered_map<char, unordered_map<unsigned int, AminoAcidData> >::iterator chain_iter = target_chains.begin();
    
    // Pack the first protein chain
    m_pdb.first.reserve( chain_iter->second.size() );
    
    for(unordered_map<unsigned int, AminoAcidData>::const_iterator i = chain_iter->second.begin();i != chain_iter->second.end();++i){
        m_pdb.first.push_back(i->second);
    }

    // It appears that BIOMT transformations are *not* used to obtain heterodimer structures ...
    //#define APPLY_BIOMT_TRANSFORMATIONS
    #ifdef APPLY_BIOMT_TRANSFORMATIONS

    unordered_map< char /*chain*/, deque<Transform> >::const_iterator biomt_iter = biomt.find(chain_iter->first);

    if(biomt_iter != biomt.end()){

        m_out << "Applying " << biomt_iter->second.size() << " coordinate transformation(s) to chain " << chain_iter->first << " in " 
            << m_filename << endl;

        // Apply the specified transformations to all atoms in this chain
        for(vector<AminoAcidData>::iterator i = m_pdb.first.begin();i != m_pdb.first.end();++i){

            for(deque<Transform>::const_iterator j = biomt_iter->second.begin();j != biomt_iter->second.end();++j){
                i->apply_transformation(*j);
            }
        }
    }
    #endif // APPLY_BIOMT_TRANSFORMATIONS

    // Pack the second protein chain
    ++chain_iter;
    m_pdb.second.reserve( chain_iter->second.size() );
    
    for(unordered_map<unsigned int, AminoAcidData>::const_iterator i = chain_iter->second.begin();i != chain_iter->second.end();++i){
        m_pdb.second.push_back(i->second);
    }

    #ifdef APPLY_BIOMT_TRANSFORMATIONS

    biomt_iter = biomt.find(chain_iter->first);

    if(biomt_iter != biomt.end()){

        m_out << "Applying " << biomt_iter->second.size() << " coordinate transformation(s) to chain " << chain_iter->first << " in " 
            << m_filename << endl;

        // Apply the specified transformations to all atoms in this chain
        for(vector<AminoAcidData>::iterator i = m_pdb.second.begin();i != m_pdb.second.end();++i){

            for(deque<Transform>::const_iterator j = biomt_iter->second.begin();j != biomt_iter->second.end();++j){
                i->apply_transformation(*j);
            }
        }
    }
    #endif // APPLY_BIOMT_TRANSFORMATIONS

    // Since this is a valid PDB file, we can propagate any changes to the atom table
    m_atom_table = local_atom_table;

    return true;
}

unsigned int parse_atom_type(const string &m_type, unordered_map<string, unsigned int> &m_atom_table,
    const bool &m_ignore_charge)
{
    // Remove leading and trailing white space and make the matching case-insensitive 
    string name = toupper( strip(m_type) );

    if(m_ignore_charge){

        // Remove charge information
        while( !name.empty() ){

            bool removed_charge_info = false;

            switch( name.back() ){
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                case '0':
                case '+':
                case '-':
                    name.pop_back();
                    removed_charge_info = true;
                    break;
            };

            if(!removed_charge_info){
                break;
            }
        }
    }

    if( name.empty() ){
        throw __FILE__ ":parse_atom_type: Missing atom type";
    }

    // The "*" symbol is used in 4q96 to indicate an unknown element type
    if(name == "*"){
        name = "UNKNOWN";
    }

    unordered_map<string, unsigned int>::const_iterator iter = m_atom_table.find(name);
    
    if( iter != m_atom_table.end() ){
        return iter->second;
    }

    // Add this atom type
    if( m_atom_table.empty() ){
        
        // Start numbering atom types from 0
        m_atom_table[name] = 0;

        return 0;
    }

    unsigned int id = 0;

    for(unordered_map<string, unsigned int>::const_iterator i = m_atom_table.begin();i != m_atom_table.end();++i){
        id = max(id, i->second);
    }

    // Increment the id
    ++id;

    m_atom_table[name] = id;

    return id;
}

unsigned int parse_amino_acid_type(const string &m_type)
{
    // Remove leading and trailing white space and make the matching case-insensitive 
    const string local = toupper( strip(m_type) );

    if(local.empty() ){
        return AminoAcidData::UNKNOWN;
    }

    #define MATCH_AA(X) \
        if(local == #X){ \
            return AminoAcidData::X; \
        }

    MATCH_AA(ALA);
    MATCH_AA(ARG);
    MATCH_AA(ASN);
    MATCH_AA(ASP);
    MATCH_AA(CYS);
    MATCH_AA(GLN);
    MATCH_AA(GLU);
    MATCH_AA(GLY);
    MATCH_AA(HIS);
    MATCH_AA(ILE);
    MATCH_AA(LEU);
    MATCH_AA(LYS);
    MATCH_AA(MET);
    MATCH_AA(PHE);
    MATCH_AA(PRO);
    MATCH_AA(SER);
    MATCH_AA(THR);
    MATCH_AA(TRP);
    MATCH_AA(TYR);
    MATCH_AA(VAL);

    // Handle non-standard amino acids
    if(local == "MSE"){
        // Replace Selenomethionine with methionine
        return AminoAcidData::MET;
    }

    if(local == "SEP"){
        // Replace phosphoserine with serine
        return AminoAcidData::SER;
    }

    return AminoAcidData::UNKNOWN;
}

string protein_sequence(const vector<AminoAcidData> &m_seq)
{
    string ret;

    ret.reserve( m_seq.size() );

    for(vector<AminoAcidData>::const_iterator i = m_seq.begin();i != m_seq.end();++i){

        switch(i->type){
            case AminoAcidData::ALA:
                ret.push_back('A');
                break;
            case AminoAcidData::ARG:
                ret.push_back('R');
                break;
            case AminoAcidData::ASN:
                ret.push_back('N');
                break;
            case AminoAcidData::ASP:
                ret.push_back('D');
                break;
            case AminoAcidData::CYS:
                ret.push_back('C');
                break;
            case AminoAcidData::GLN:
                ret.push_back('Q');
                break;
            case AminoAcidData::GLU:
                ret.push_back('E');
                break;
            case AminoAcidData::GLY:
                ret.push_back('G');
                break;
            case AminoAcidData::HIS:
                ret.push_back('H');
                break;
            case AminoAcidData::ILE:
                ret.push_back('I');
                break;
            case AminoAcidData::LEU:
                ret.push_back('L');
                break;
            case AminoAcidData::LYS:
                ret.push_back('K');
                break;
            case AminoAcidData::MET:
                ret.push_back('M');
                break;
            case AminoAcidData::PHE:
                ret.push_back('F');
                break;
            case AminoAcidData::PRO:
                ret.push_back('P');
                break;
            case AminoAcidData::SER:
                ret.push_back('S');
                break;
            case AminoAcidData::THR:
                ret.push_back('T');
                break;
            case AminoAcidData::TRP:
                ret.push_back('W');
                break;
            case AminoAcidData::TYR:
                ret.push_back('Y');
                break;
            case AminoAcidData::VAL:
                ret.push_back('V');
                break;
            case AminoAcidData::UNKNOWN:
                // Skip this amino acid!
                break;
            default:
                throw __FILE__ ":protein_sequence: Unknown amino acid!";
        };
    }

    return ret;
}

vector<float> parse_biomt_col(const string &m_buffer)
{
    stringstream ssout(m_buffer);

    vector<float> ret;

    // We expect five columns
    ret.reserve(5);

    float value;

    while(ssout >> value){
        ret.push_back(value);
    }

    return ret;
}

bool check_atomic_distances(const string &m_accession, const PDBComplex &m_pdb, const float &m_threshold)
{
    bool valid = true;

    for(vector<AminoAcidData>::const_iterator aa_1 = m_pdb.first.begin();aa_1 != m_pdb.first.end();++aa_1){
        for(vector<AtomData>::const_iterator atom_1 = aa_1->atoms.begin();atom_1 != aa_1->atoms.end();++atom_1){
            
            for(vector<AminoAcidData>::const_iterator aa_2 = m_pdb.second.begin();aa_2 != m_pdb.second.end();++aa_2){
                for(vector<AtomData>::const_iterator atom_2 = aa_2->atoms.begin();atom_2 != aa_2->atoms.end();++atom_2){
                    
                    const float d = atom_1->distance(*atom_2);
                    
                    if(d < m_threshold){
                        cout << "atomic clash: " << d << ' ' << m_accession << " " << atom_1->index 
                            << " (first atom); " << atom_2->index << " (second atom)" << endl;
                        
                        valid = false;
                    }
                }
            }
        }
    }

    return valid;
}

// From https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html
// COLUMNS       DATA  TYPE      FIELD        DEFINITION
// -------------------------------------------------------------------------
//  1 -  6        Record name    "CONECT"
//  7 - 11        Integer        serial       Atom  serial number
// 12 - 16        Integer        serial       Serial number of bonded atom
// 17 - 21        Integer        serial       Serial  number of bonded atom
// 22 - 26        Integer        serial       Serial number of bonded atom
// 27 - 31        Integer        serial       Serial number of bonded atom
vector<int> split_conect(const string &m_buffer)
{
    vector<int> ret;

    ret.reserve(5 /*maximum number of atoms in a CONECT record*/);

    if(m_buffer.find("CONECT") != 0){
        throw __FILE__ ":split_conect: Unable to parse CONECT header";
    }

    int index = 0;
    int index_len = 0;

    const size_t buffer_len = m_buffer.size();

    for(size_t i = 6 /*= strlen("CONECT")*/;i < buffer_len;++i){

        if(isdigit(m_buffer[i])){

            index = index*10 + (m_buffer[i] - '0');
            ++index_len;
        }

        if(i%5 == 0){ // separate at positions 10, 15, 20, 25, and 30
            
            if(index_len > 0){
                ret.push_back(index);

                index = 0;
                index_len = 0;
            }
        }
    }

    if(index_len > 0){
        ret.push_back(index);
    }

    // Handle this error in the calling function
    //if(ret.size() > 5){
    //    throw __FILE__ ":split_conect: Found an invalid CONECT record with > 5 entries";
    //}

    return ret;
}