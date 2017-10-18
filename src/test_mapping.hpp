#ifndef TEST_MAPPING_HPP
#define TEST_MAPPING_HPP

#include "iostream"
#include <map>
#include <vector>

#include "pw_alignment.hpp"

class map_check{
	public:
	map_check(){
		error_counter = 0;
	}
	~map_check(){}
	void read_graph_maf(std::ifstream & graph_maf);
	void read_alignments(std::ifstream & alignments);//Read all the intial created alignments before mapping any read against a graph
	void check_output(std::ifstream & mapping_maf,std::string & seqname);
	void check_an_alignment(unsigned int & ref1, const std::string & seqname, size_t & left2 , size_t & right2);
	void read_txt_file_long_center(std::ifstream & , std::string & );
	void read_txt_file(std::ifstream & , std::string &);
	std::multimap<std::string, std::pair<unsigned int, unsigned int> > get_nodes()const{
		return nodes;
	}
	void print_nodes(){
		for(std::multimap<std::string , std::pair<unsigned int ,unsigned int> >::iterator it = nodes.begin() ; it!= nodes.end() ; it++){
			std::vector<std::string> parts;
			strsep(it->first,":",parts);
			std::cout<< "between node " << parts.at(0) << " and "<< parts.at(1) << " : from " << it->second.first << " to " << it->second.second <<std::endl;
		}

	}
	void print_graph(){
		for(std::map<unsigned int , std::vector<std::string> >::iterator it=clusters.begin(); it!= clusters.end(); it++){
			std::cout<< it->first<< std::endl;
			if(it->first==8){
				std::cout<<"if 8 "<<std::endl;
				std::cout<< it->second << std::endl;
			}
		}
	}
	const std::map<unsigned int , std::vector<std::string> > & get_clusters()const;
	const std::multimap<std::string , std::pair<size_t, size_t> > & get_intervals()const;
	const std::map<int,std::pair<int, int> > & get_nonaligned()const;
	const unsigned int get_node_length(unsigned int &)const;
	const std::multimap<size_t , const pw_alignment> & get_als()const{
		return al_from_graph;
	}
	private:
	std::map<unsigned int , std::vector<std::string> > clusters;
	std::multimap<std::string , std::pair<size_t, size_t> > boundries;//Multimap because several part of a sequence can be clustered together in one cluster.
	std::multimap<std::string , std::pair<size_t, size_t> > test_boundries;
	std::multimap<std::string , std::pair<unsigned int ,unsigned int> > nodes;
	std::map<int,std::pair<int, int> > non_aligned; //cluster id, pair(position, length)

	std::map<unsigned int, unsigned int> node_length;
	std::multimap<size_t , const pw_alignment> al_from_graph;
	size_t error_counter;

};


class test_sim_reads_mapping{
	public:
		test_sim_reads_mapping(map_check & mp):mp_check(mp){}
		~test_sim_reads_mapping(){}
		void read_sim_maffile(std::ifstream & sim_maffile, std::map<std::string , std::pair<size_t , size_t> > & onreads);
		void this_part_position_on_seq(bool & , size_t & , size_t & , size_t & , size_t &, size_t &, size_t & );
		void check_output(std::ifstream & mapping_maf, std::map<size_t, std::pair<size_t,size_t> > & nodes_on_ref_graph);
		void check_if_shifted(unsigned int & ref_id, bool & IsMapped);
		void compare_with_my_graph(unsigned int & ref_id, size_t & start_on_seq, size_t & end_on_seq,size_t & begin_on_ref1, size_t & end_on_ref1, size_t & begin_on_ref2, size_t & end_on_ref2,size_t & read_length);
		void compare_with_reveal_graph(std::map<size_t, std::pair<size_t,size_t> > & nodes_on_ref_graph, size_t & ref_id, size_t & begin , size_t & end, size_t & error, bool & GAPONLY);
		void find_position_on_member(unsigned int & center, size_t & begin_on_center , size_t & end_on_center,size_t & begin_on_mem , size_t & end_on_mem,int & position);
	private:
		map_check mp_check;

		std::map<size_t,std::pair<size_t,size_t> > reads_position_on_seq;
		std::map<size_t , const pw_alignment> alignments;

};

class test_reveal{
	public:
		test_reveal(){};
		~test_reveal(){};
		void read_gfa(std::ifstream & gfa);
		void read_the_result(std::ifstream &, std::map<std::string, std::pair<size_t, size_t> > & contigs);
		void compare_with_reveal();
		void compare_with_path(std::ifstream &);//TOO specific!!
		const std::map<size_t, std::pair<size_t,size_t> > get_ref_graph_nodes()const{
			return ref_graph_nodes;
		}
	private:
		std::vector<std::string> nodes;
		std::multimap<size_t , std::pair<size_t , size_t> >from_mapping_output; //for each read(name and direction, from-to on the ref) //TODO direction
		std::map<size_t , std::pair<size_t, size_t> > ref_graph_nodes; //size_t is node name from gfa file and string is its content.
};

#endif
