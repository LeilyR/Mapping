#ifndef REF_GRAPH_HPP
#define REF_GRAPH_HPP

#include "data.hpp"
#include "pw_alignment.hpp"
#define MAXGAP 400

//This class is only used to read the dot file includes the reference graph
class ref_graph{
	public:
	ref_graph( const all_data & d):data(d){
		longname2seqidx = data.getLongname2seqidx();
		std::map<std::string, size_t>::iterator longname = longname2seqidx.begin();
		std::cout << "the begin is "<< longname->first<<std::endl;
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
	void look_for_neighbors(int & node, std::map<int,bool> & visited , std::string & refacc, std::map<int , std::vector<std::pair<std::vector<int>,size_t> > > & parent_length);
	void bfs(int & startnode, std::string & refacc, size_t & right_on_ref);
	void delete_path(std::vector<int> & this_path);
	const std::set<vector<int> > get_paths()const;
	const std::set<int> get_nodes()const;
	const std::set<int> get_nodes_on_paths(const int & next)const;

	const std::set<int> get_subgraph_nodes() const;
	std::set<int> get_adjacencies(int & node)const{
		std::set<int> nodes;
		std::map<int, std::set<int> >::const_iterator it = adjacencies.find(node);
		if( it != adjacencies.end()){
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
	private:
	const all_data & data;
	std::map<std::string, size_t> longname2seqidx;

	std::map<int, std::set<int> > adjacencies;
	std::map<int, std::set<int> > predecessors;

	std::map<int, std::set<int> > dynamic_adjacencies; //It is used for bfs. It has the same content as adjacencies at the begining but whenever bfs is called they gradually reduced! 

	std::set<vector<int> > paths;
	std::set<int> nodes_on_paths;
	std::map<int, size_t> path_length;

//NEW:
	std::set<int> sub_graph_nodes;
	std::map<int , std::set<int> > sub_graph;
	std::map<int, std::set<int> > pre_nodes_on_subgraph;
};

#endif
