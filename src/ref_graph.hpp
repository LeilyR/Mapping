#ifndef REF_GRAPH_HPP
#define REF_GRAPH_HPP

#include "data.hpp"
#include "pw_alignment.hpp"
//#include "hash.hpp"
#define MAXGAP 400
#define K 25
//This class is only used to read the dot file includes the reference graph
class ref_graph{
	public:
	ref_graph( const all_data & d):data(d){
		longname2seqidx = data.getLongname2seqidx();
	//	std::map<std::string, size_t>::iterator longname = longname2seqidx.begin();
	//	std::cout << "the begin is "<< longname->first<<std::endl;
	}
	~ref_graph(){}
	void read_dot_file(std::string &, std::string &);
	void add_adjacencies(std::string & , std::string & , std::string & , std::string & );
	void add_adjacencies(std::string & , std::string & );
	void add_predecessor(int & this_node , int & pre_node);
	const std::map<int , std::set<int> > & get_adjacencies()const;
	std::set<std::vector<int> > get_predecessor(unsigned int &, bool , size_t &, size_t &);
	void find_predecessors_bfs(int & startnode, std::string & refacc, size_t & left_on_ref, std::set<std::vector<int> > & pre_paths , size_t &);
	const void look_for_predecessor(int & node , size_t & length , size_t & current_length, std::string & acc_name,std::set<int> & visited,std::vector<int> & this_pre_nodes, std::vector<std::vector<int> > & all_pre_nodes)const;
	const std::vector<std::vector<int> > get_successor(unsigned int & ref_id, bool dir , size_t & length)const;
	const void look_for_successor(int & node, size_t & length, size_t & current_length, std::string & acc_name, std::set<int> & visited, std::vector<int> & this_adjacencies, std::vector<std::vector<int> > & all_adjacencies)const;
	void read_gfa_for_adj(std::string &);
	void read_gfa_for_seq(std::string & gfafile , std::ofstream & fasta , std::ofstream & paths);
	const unsigned get_refid(size_t & , int& )const;

	size_t seq_length(std::string &, std::string &);
	std::string seqname(int & );
	void deep_first_search(int &, std::string & , size_t & );
	void look_for_neighbors(int & node, std::map<int,bool> & visited , std::string & refacc);
	void bfs(int & startnode, std::string & refacc, size_t & right_on_ref);
	void delete_path(std::vector<int> & this_path);
	const std::set<vector<int> > get_paths()const;
	const std::set<int> get_nodes()const;
	const std::set<int> get_nodes_on_paths(const int & next)const;

	const std::set<int> get_subgraph_nodes() const;
	std::set<int> get_adjacencies(int & node)const{
		std::set<int> nodes;
		std::map<int, std::set<int> >::const_iterator it = adjacencies.find(node);
		if(it == adjacencies.end()){
			std::cout << "END"<<std::endl;
		}else{
		//	std::cout << "HERE"<<std::endl;	
			nodes = it->second;		
		}
		return nodes;
	}
	std::set<int> get_predecessor(int& node)const {
		std::set<int> nodes;
		std::cout << "node: "<< node <<std::endl;
		std::map<int, std::set<int> >::const_iterator it = predecessors.find(node);
		if( it != predecessors.end()){
			nodes = it->second;
		}
		return nodes;
	}
	void make_sub_graph(int & , std::string & , size_t & );
	std::set<std::pair<std::vector<int>,size_t> > get_parent_length(int & node)const{ 
		std::set<std::pair<std::vector<int>,size_t> > temp;
		std::map<int , std::set<std::pair<std::vector<int>,size_t> > >::const_iterator it= parent_length.find(node);
		if(it != parent_length.end()){
			return it->second;
		}
		return temp;
	}
	const std::map<int, std::set<std::pair<std::vector<int>,size_t> > > get_parent_length()const{ 
		return parent_length;
	}
	const std::map<int , std::string> get_nodes_content()const{
		return nodes_content;
	}

	const std::string get_content(int & node)const{
		std::map<int, std::string>::const_iterator it = nodes_content.find(node);
		assert(it != nodes_content.end());
		return it->second;

	}

	void merge_nodes(std::string &, size_t &);
	void merge_nodes_bfs(std::string & ref_acc, size_t & length);
	void edge_contraction();
	void find_d1_nodes(int &, std::set<int> & seen);
	void write_gfa(std::ofstream &);
	private:
	const all_data & data;
	std::map<std::string, size_t> longname2seqidx;

	std::map<int, std::set<int> > adjacencies;
	std::map<int, std::set<int> > predecessors;

	std::map<int, std::set<int> > dynamic_adjacencies; //It is used for bfs. It has the same content as adjacencies at the begining but whenever bfs is called they gradually reduced!  //XXX Guess it is not neccessary anymore check and delete if it is not used!!

	std::set<vector<int> > paths;
	std::set<int> nodes_on_paths;
	std::map<int, size_t> path_length;

//NEW:
	std::set<int> sub_graph_nodes;
	std::map<int , std::set<int> > sub_graph;
	std::map<int, std::set<int> > pre_nodes_on_subgraph;

	std::map<int , std::set<std::pair<std::vector<int>,size_t> > > parent_length;
	std::map<int, std::string> nodes_content; //index , string --> the node content
	std::map<int, std::string> new_nodes_content;
	std::map<int , std::set<int> > new_edges; //XXX Add the end dont forget that the negative directions need to be added.
	std::map<std::vector<int>, int > old_edges; //Has merged nodes and from which nodes they came

	std::set<int> begin_nodes;

};


class Kgraph{
	public:
		Kgraph(const ref_graph & rg):rgraph(rg){}
		void make_Kgraph();
		void get_reverse_complement(std::string & sequence_in , std::string & sequence_out){
			for(size_t i = sequence_in.size(); i >0; i--){
				sequence_out += dnastring::complement(sequence_in.at(i-1));
			}
		}
		void bfs(std::map<std::pair<int,int>, std::vector<std::string> > & remainder, int  this_node, int  pre_node , std::map<int,bool> & visited);
		const std::map<std::string, std::set<std::string> > get_edges()const{
			return edges;
		}
		private:
		const ref_graph & rgraph;
		std::map<std::string , std::set<std::string> > edges; // Add hash value of nodes ??!!
			

};
#endif
