#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include <map>
#include <math.h> 
#include <algorithm>

#include "pw_alignment.hpp"

using namespace std;
//#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
//#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

//ICL library header files: 
//#include <boost/icl/discrete_interval.hpp>
//#include <boost/icl/interval_map.hpp>
//#include <boost/icl/interval_set.hpp>

// Library to read Sam File
//#include "SamFile.h"
//#include "SamValidation.h"
//#include "Cigar.h"
//#include "GenomeSequence.h"
//#include <SamFlag.h>

#include <boost/iostreams/stream.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > btokenizer;
void strsep(std::string str, const char * sep, std::vector<std::string> & parts);

class dnastring {
	public:
	dnastring(std::string str);
	dnastring();
	dnastring(const dnastring &d);
	~dnastring();

	char at(size_t pos) const;
	size_t length() const;
	std::string str() const;

	static bool found_iupac_ambiguity;
	static char complement(char c);
	static size_t base_to_index(char base);
	static char index_to_base(size_t index);


	private:
	std::vector<bool> bits;

	static char base_translate_back(bool bit1, bool bit2, bool bit3);

};



class all_data {

	public:
		all_data();

		void read_fasta_sam(std::string fasta_all_sequences, std::string sam_all_alignments);
		void read_accknown_fasta_sam(std::string fasta_all_sequences, std::string sam_all_alignments);
		void read_fasta_maf(std::string fasta_all_sequences, std::string maf_all_alignments);
		void read_fasta_maf_forward_read_only(std::string fasta_all_sequences, std::string maf_all_alignments);
		void read_fasta(std::string fasta_reads);//It is used to read the reads for mapping.

		// no copy constructor, never copy all data
		~all_data();


		const dnastring & getSequence(size_t index) const;
		const pw_alignment & getAlignment(size_t index) const;
		const located_alignment & getLAlignment(size_t index) const;

		void add_pw_alignment(const pw_alignment& p);
		const std::vector<pw_alignment>& getAlignments()const;
	//	const multistd::map<size_t, size_t> & getAlOnRefMap(size_t seq_idx) const;

		size_t numSequences() const;
		size_t numAlignments() const;
		const size_t numAcc() const;
		bool alignment_fits_ref(const pw_alignment & al) const;
		void print_ref(const pw_alignment * al)const;
		const std::vector<size_t> & getAcc(size_t acc)const;//return the vector of all sequences for a certain acc.
		const size_t get_seq_number_of_acc(size_t & acc)const;
		size_t accNumber(size_t sequence_id) const;
		const std::string get_acc(size_t acc)const; //get the accession number and return its name.
		const size_t get_acc_id(std::string acc)const;
		std::string get_seq_name(size_t s) const;
		size_t get_seq_size(size_t s) const;
		void add_accession(const std::string & acc);
		size_t numOfAcc() const;
		const std::map< std::string, size_t>& getLongname2seqidx()const;
		void compare_seq_with_decoding(std::ifstream &);
		void compare_seq_with_read(std::ifstream& in, size_t seq_id, size_t & from, size_t & to);
		void checkAlignmentRange(const pw_alignment &)const;
		std::vector<std::string> get_reads()const;
		void make_fasta(std::ostream &, std::vector<std::string> & );
		std::string extract_seq_part(unsigned int & ref, size_t & from , size_t & to)const;
		std::string extract_reverse_seq_part(unsigned int & ref, size_t & from , size_t & to)const;
		std::string get_read(unsigned int read_id)const{
                        return reads.at(read_id);
                }
		bool alignmentAlreadyExist(const pw_alignment &a)const; //XXX Added by Marion. 
		int findIdAlignment(std::string name,int startOnRead, int endOnRead,int startOnNode,int endOnNode);//XXX Added by Marion
		int findIdAlignment(std::string name,unsigned int start, unsigned int end);//XXX Added by Marion

	private:
		// data
		std::vector<dnastring> sequences;
		std::vector<std::string> sequence_names;
		std::vector<pw_alignment> alignments; // TODO  remove this (slower, more memory and duplicate with l_alignments)
		std::vector<located_alignment> l_alignments; 
		std::vector<std::string> reads;
		// fast access indices
		std::map< std::string, std::vector< size_t> > acc_sequences; // acc name -> sequences of that acc
		std::vector<size_t> sequence_to_accession; // sequence id -> accession id
		std::map<std::string, size_t> accession_name; // accession name -> accession id
		std::map< std::string, size_t> longname2seqidx; // long sequence name ("Accession:sequence name") -> sequence index

		

		void insert_sequence(const std::string & acc, const std::string & seq_name, const std::string & dna);
		void insert_read(const std::string & read);
		static void name_split(const std::string & longname, std::string & acc, std::string & name);



};

template<typename overlap_type>
class compute_cc;
class overlap{
public:
	overlap(const all_data&);
	overlap(const overlap & o);
	~overlap();
	void insert_without_partial_overlap(const pw_alignment & p);
	// This function removes an alignment with adress identity to remove from overlap, then deletes remove
	void remove_alignment(const pw_alignment & remove);

	void test_all() const;
	void test_all_part()const;// Change it later to check all those pieces with gap in one sample. Put them all in a std::set and check if the only missed parts of coverage are those parts.
	void test_overlap()const;
	void test_no_overlap_between_ccs(const pw_alignment &, std::set<const pw_alignment*, compare_pointer_pw_alignment> &)const;
	void print_all_alignment() const;
	const pw_alignment * get_al_at_left_end(size_t ref1, size_t ref2, size_t left1, size_t left2) const;
	std::multimap<size_t, const pw_alignment &>& get_als_on_reference(size_t sequence) ;
	const std::multimap<size_t, const pw_alignment & >& get_als_on_reference_const(size_t sequence) const ;
//	const boost::icl::interval_map<size_t, std::set<const pw_alignment*> > & get_alignments_interval(size_t sequence)const;
	void test_multimaps()  ;
	bool checkAlignments(const pw_alignment & p)const;

	const std::set<pw_alignment, compare_pw_alignment> & get_all() const;

	void test_partial_overlap() const;
	void test_al_on_ref()const;
	bool check_alignment_address(const pw_alignment & , const pw_alignment * )const;
	static void test_partial_overlap_set(std::set< const pw_alignment *, compare_pw_alignment> & als);
	static void test_partial_overlap_vec(std::vector< const pw_alignment *> & als);
	static bool check_po(size_t l1, size_t r1, size_t l2, size_t r2);

	size_t size() const;
private:
	const all_data & data;
	std::set<pw_alignment, compare_pw_alignment> alignments;
	std::vector< std::multimap< size_t, const pw_alignment &> > als_on_reference; // sequence index -> pos on that sequence -> alignment reference
//	std::vector< boost::icl::interval_map<size_t , std::set<const pw_alignment*> > > alignments_intervals;
	
	
};



class splitpoints {
	public:
	splitpoints(const pw_alignment & , const overlap &, const all_data &);
	~splitpoints();
	void find_initial_split_points(size_t sequence, size_t left, size_t right);
	void find_initial_split_points_nonrecursive(size_t sequence, size_t left, size_t right);
	void recursive_splits();//initial split points + all induced split points everywhere
	void nonrecursive_splits(); // initial split points only 
	void insert_split_point(size_t sequence, size_t position);//recursively find the other split points
	void insert_split_point_nonrecursive(size_t sequence, size_t position);//insert without recursion
	void split_all(std::set<pw_alignment, compare_pw_alignment> & remove_alignments, std::vector<pw_alignment> & insert_alignments);
	void splits(const pw_alignment & p,  std::vector<pw_alignment> & insert_alignments);
	bool onlyGapSample(const pw_alignment & p) const;

	// std::vector<pw_alignment>  get_insert () const;

	private:
	const overlap & overl;
	const pw_alignment & newal;
	const all_data & data;
	std::vector<std::set<size_t> > split_points;
	std::vector<pw_alignment> insert_alignments;	
/*
	initial points in std::sets (method with and without sets)
	compute remove and insert alignments
	cost change level 1, abort if no gain
	if there is gain recurse


*/
	
};

class model{
public:
	model(all_data &);
	~model();
	void acc_base_frequency();
	void alignment_modification();
	void cost_function(pw_alignment& p) const;
	void cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2) const ;
	void cost_function(const pw_alignment& p, size_t & acc1, size_t & acc2, double & c1, double & c2, double & m1, double & m2) const ;
	void cost_function(const pw_alignment & p, double & m1, double & m2, size_t & refacc, size_t & readacc)const;
	void gain_function(const pw_alignment& p, double & g1, double & g2) const;
	void gain_function(const pw_alignment& p, size_t & acc1, size_t & acc2, double & g1, double & g2) const;
	const std::vector<std::vector<double> > get_al_modification(size_t & readacc , size_t & refacc) const;

	void total_information(size_t&);
	void train();

private:
	all_data & data;
//	std::vector<vector<size_t> > transform;
	std::vector<std::vector<double> > cost_on_acc;
	std::vector<std::vector<std::vector<std::vector<double> > > >modification;
};
#define Alignment_level 1
#define Sequence_level 2
#define NUM_DELETE 5
#define NUM_KEEP 10

class abstract_context_functor {
	public:
	abstract_context_functor();
	virtual void see_context(size_t acc1, size_t acc2,const pw_alignment& p, size_t pos, std::string context, char last_char);
	virtual void see_entire_context(size_t acc1, size_t acc2, std::string entireContext);

};


class counting_functor : public abstract_context_functor {
	public:
	counting_functor(all_data &);
	virtual void see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, std::string context, char last_char);
	const std::map<std::string, std::vector<size_t> > & get_context(size_t acc1, size_t acc2)const;
	void total_context();
	void total_context(size_t &, size_t &);
	double get_total(size_t acc1, size_t acc2, std::string context)const;
	void create_context(size_t acc1, size_t acc2, std::string context);
	private:
	all_data & data;
	std::vector<std::vector<std::map<std::string, std::vector<size_t> > > >successive_modification;
	std::vector<std::vector<std::map <std::string, double > > > total;	




};
class mc_model;


class cost_functor : public abstract_context_functor {
	public:
	cost_functor(all_data &, const std::vector<std::vector<std::map<std::string,std::vector<double> > > >&);
	virtual void see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, std::string context, char last_char);
	double get_modify(const pw_alignment & p, size_t acc1, size_t acc2)const;
	private:
	std::vector<std::vector< std::map<std::string, std::vector<double> > > >  modification; 
	all_data & data;
	double modify1;
	double modify2;	

};


class adding_functor : public abstract_context_functor {
	public:

};
/*class encoding_functor : public abstract_context_functor {
	public:
	encoding_functor(all_data& , mc_model* , wrapper &, dlib::entropy_encoder_kernel_1 &);
	virtual void see_context(size_t acc1, size_t acc2,const pw_alignment& p, size_t pos, std::string context, char last_char);	
	virtual void see_entire_context(size_t acc1, size_t acc2, std::string entireContext);
	const std::map<std::string, std::vector<double> > & get_alignment_context()const;
//	std::vector<std::string> & get_alignment_context(pw_alignment& p)const;
	private:
	all_data & data;
	mc_model * model;
	wrapper& wrappers;
	dlib::entropy_encoder_kernel_1 & enc;
	std::string alignment_pattern;//shayadam behtar bashe ye std::vector of string tarif konam
	std::map<std::string, std::vector<double> > alignment_context;

};*/

class clustering_functor : public abstract_context_functor{
	public:
	virtual void see_context(size_t acc1, size_t acc2, const pw_alignment& p, size_t pos, std::string context, char last_char);//computing_modification_oneToTwo is used to fill in the map of modification between center and its associated member.

	private:
	std::map<std::string, std::vector<double> >modification;


};
class decoding_functor : public abstract_context_functor {


};

class mc_model{
	public:
		mc_model(all_data&);
		~mc_model();
		void markov_chain();
		void markov_chain_alignment();
		const std::vector<std::vector<std::map <std::string, vector<double> > > >& get_mod_cost()const;
		void cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2)const ;
		void gain_function(const pw_alignment& p, double & g1, double & g2)const ;
		double get_the_gain(const pw_alignment &, std::string & )const;
		void train(std::ofstream &);
		char modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const;
		void modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep)const;
		void computing_modification_oneToTwo(const pw_alignment & p, abstract_context_functor & functor)const;
		void computing_modification_twoToOne(const pw_alignment & p, size_t & , abstract_context_functor & functor)const;
		void cost_function( pw_alignment& p) const;
		std::string print_modification_character(char enc)const;
		const std::map<std::string, std::vector<double> > & getPattern(size_t acc)const;
		const std::vector<double> & get_create_cost(size_t acc) const;
		const std::vector<std::map<std::string, vector<double> > > & model_parameters()const;
		void write_parameters(std::ofstream &);
		void write_alignments_pattern(std::ofstream&);
		std::vector<unsigned int> get_high_at_position(size_t seq_index, size_t position) const;
		std::vector<unsigned int> get_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const;
		std::vector<unsigned int> get_reverse_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const;
		const std::map<std::string, std::vector<unsigned int> > & get_high(size_t acc)const;
		void make_all_the_patterns();
		void make_all_alignments_patterns();
		void set_patterns(std::ifstream&);
		void set_alignment_pattern(std::ifstream&);
		std::string get_context(size_t position, size_t seq_id)const;
		std::vector<size_t> get_powerOfTwo()const;
		std::string get_firstPattern()const;
		std::string get_firstAlignmentPattern() const;
		const std::map<std::string, std::vector<unsigned int> >& get_highValue(size_t acc1, size_t acc2)const;
		void computing_modification_in_cluster(std::string center, std::string member)const;
		size_t modification_length(char mod)const;
	//	void get_encoded_member(const pw_alignment & al, size_t center_ref, size_t center_left, encoding_functor & functor,std::ofstream&)const;
	//	void get_encoded_member_long_center(const pw_alignment & , unsigned int& , encoding_functor & functor,std::ofstream&)const;
		void get_insertion_high_between_centers(size_t& seq_id ,char & seq_base, char & last_center_base, unsigned int& center_ref,std::string & ,unsigned int& high, unsigned int& low)const;
		void get_a_keep(unsigned int & , unsigned int & , size_t & , std::vector<unsigned int> & , std::vector<unsigned int> & )const;
	private:
	all_data & data;
	std::vector<std::map<std::string, std::vector<double> > >sequence_successive_bases;
	std::vector<std::vector<double> > create_cost;
	std::vector<size_t> powersOfTwo;
	std::vector<std::vector<std::map<std::string, std::vector<double> > > >mod_cost; // alignment modificaton information cost
	std::vector<std::map<std::string, std::vector<unsigned int> > > high;//sequences patterns
	std::map<std::string,std::vector<double> > all_the_patterns; // TODO vector<double> part is wrong (independent of accession) --- it is because at the beginning all the values are set to 0. And then later on in high we set them accession dependent. I need to make a better way to represent all the pattern
	std::set<std::string> all_alignment_patterns; // all possible alignment patterns
	std::vector<std::vector<std::map<std::string , std::vector<unsigned int> > > >highValue;//alignments patterns
	mutable std::string last_context; //last context of an alignment. it is used for arith encoding of long centers

		

};



#endif
