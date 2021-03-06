#ifndef MAP_READ_HPP
#define MAP_READ_HPP

#include "data.hpp"
#include "pw_alignment.hpp"
#include "ref_graph.hpp"
#include "needleman_wunsch.hpp"
#include "simple_NW_algo.hpp"
#include "dijkstra.hpp"
#include "fibonacci_heap.hpp"
#include "fibonacci.hpp"

#include <time.h>
#include <stdio.h>
#include <list>

/*class deep_first{//It runs deep first algorithm on reference graph to retrun parts of graph that the distance between the first and the last node is less than MAXGAP.
	public:
	deep_first(const all_data & d, std::map<int, std::set<int> > & adjacencies):data(d){
		this -> adjacencies = adjacencies;
	//	for(std::map<int, std::set<int> >::iterator it = adjacencies.begin(); it != adjacencies.end(); it++){
	//		std::cout << it->first << std::endl;
	//	}

	}
	~deep_first(){}
	size_t seq_length(std::string &, std::string &);
	std::string seqname(int & );
	void deep_first_search(int &, std::string & , size_t & );
	void look_for_neighbors(int & , std::map<int,bool> &  , std::string &, int & accu_length, std::vector<int> & apath, size_t & );
	void bfs(int & startnode, std::string & refacc, size_t & right_on_ref);
	const vector<vector<int> > get_paths()const;
	const std::set<int> get_nodes()const;
	private:
	const all_data & data;
	std::map<int, std::set<int> > adjacencies;
	vector<vector<int> > paths;
	std::set<int> nodes_on_paths;
	std::map<int, size_t> path_length;
	std::multimap<int,int> parental_relation; //first one shows the child node ,second one shows the parent

};*/

//This class is used to make several components on a read where each component contains aligned regions with less than MAX_GAP base distance
template<typename T>
class als_components{

	public:
	als_components(const all_data & d , ref_graph & refg,  const T & m ,const std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> >  & ccs):data(d), rgraph(refg),model(m){
		index = 2;
		alignments_left.resize(ccs.size());
		alignments_right.resize(ccs.size());

		for(size_t i = 0 ; i < ccs.size(); i++){
		//	std::cout << "cc at " << i << " : "<<std::endl;
			for(std::set<const pw_alignment* , compare_pointer_pw_alignment >::iterator it = ccs.at(i).begin() ; it != ccs.at(i).end() ; it++){
				const pw_alignment * al = *it;
				alignments_left.at(i).insert(*al);
				alignments_right.at(i).insert(*al);

			//	al->print();

			}
			std::multiset<pw_alignment , sort_pw_alignment_by_left >::iterator first = alignments_left.at(i).begin();
			pw_alignment p = *first;
			size_t l2,r2;
			p.get_lr2(l2,r2);
			size_t lower = l2;
			std::cout << "lower " << lower <<std::endl;//The least left
			std::multiset<pw_alignment , sort_pw_alignment_by_right>::iterator last = alignments_right.at(i).end(); 
			--last;
			p = *last;
			p.get_lr2(l2,r2);
			size_t upper = r2;
			std::cout << "upper "<< upper<<std::endl;
			boundries.insert(std::make_pair(lower,std::make_pair(upper,i)));
			
			
		}
	}
	~als_components(){}
	void merg_components(){
		std::multiset<pw_alignment,sort_pw_alignment_by_left > templ;
		std::multiset<pw_alignment,sort_pw_alignment_by_right > tempr;
		size_t previous_right, current_left;
		size_t counter = 0;
		for(std::map<size_t, std::pair<size_t,size_t> >::iterator it = boundries.begin(); it != boundries.end();it++){
			std::cout << it->first <<  " "<< it->second.first<<std::endl;
			if(it == boundries.begin()){
				previous_right = it->second.first;
				templ = alignments_left.at(it->second.second);
				tempr = alignments_right.at(it->second.second);
				counter++;
			}else{
				current_left = it->first;
				if(current_left - previous_right < MAXGAP){
					templ.insert(alignments_left.at(it->second.second).begin(), alignments_left.at(it->second.second).end());
					tempr.insert(alignments_right.at(it->second.second).begin(),alignments_right.at(it->second.second).end());
					previous_right = it->second.first;
					counter++;
				}else{
					std::cout<<"members of a new component"<<std::endl;
					for(std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it1 = templ.begin() ; it1 != templ.end();it1++){
						pw_alignment pw = *it1;
						pw.print();
					}
					alignments.push_back(templ);
					ordered_als_by_right.push_back(tempr);
					assert(tempr.size()==templ.size());
					templ = alignments_left.at(it->second.second);
					tempr = alignments_right.at(it->second.second);
					previous_right = it->second.first;					
					counter++;
				}
			}
		}
		if(templ.size() !=0){//Adding the last one
			assert(tempr.size() != 0);
			assert(tempr.size()==templ.size());
			ordered_als_by_right.push_back(tempr);
			alignments.push_back(templ);
		}
		assert(alignments.size()==ordered_als_by_right.size());
		std::cout << "als size " << alignments.size()<<std::endl;
		adjacencies.resize(alignments.size());
		weight_of_edges.resize(alignments.size());
		indices_nodes.resize(alignments.size());
		node_indices.resize(alignments.size());
	}
	void filter_als(size_t & comp);
	/*Go over the alignments, for each al find all the possible paths on the graph that starts with that alignment and has <=MAXGAP length,
	then find the alignments that their left is closer than MAXGAP bases to the current one.Check them on the graph. If they are on those 
	paths that we already picked then we make new als using needleman wunsch algorithm.
	For those which are not on those paths we check if their left is smaller than beginning of the current alignment component- MAXGAP/2, 
	if so then we try to find the graph paths for it and go further otherwise discard it.*/
	void find_als_on_paths(std::ofstream & output, size_t & refacc, size_t & readacc);
	void finding_successors_on_subref(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors,size_t & refacc, size_t & readacc);
	void find_best_path_to_neighbour(size_t & comp , size_t & from_on_read, size_t & to_on_read, size_t & to_on_ref ,std::map<std::pair<size_t, size_t> , double> & weight, std::map<std::pair<int, int>, std::set<size_t> > & seen_edges, std::map<size_t,size_t> & cur_pre, size_t & cur_index, const pw_alignment & cur_al, const pw_alignment & nex_al,int & last_node, int & current_node,size_t & refacc, size_t & readacc,std::set<size_t> & visited);
	void make_last_als_on_a_path(size_t & comp, unsigned int & ref_id, unsigned int & read_id, size_t & to_on_ref, size_t & from_on_read, size_t & to_on_read, int & cur_node, int & last_node, size_t & cur_index , const pw_alignment & cur_al, const pw_alignment & last_al ,std::map<std::pair<size_t, size_t> , double> & weight, std::map<std::pair<int, int>, std::set<size_t> > & seen_edges, std::map<size_t,size_t> & cur_pre, size_t & readacc , size_t refacc);
	void find_best_path_to_neighbour(size_t comp,int & next_name, int & current_name,size_t & cur_index, const pw_alignment & cur_al, const pw_alignment & neighbour_al , size_t & to_on_read, size_t & last_pos_on_ref , size_t & refacc, size_t & readacc, std::set<int> & nodes_on_paths, std::map<std::pair<int, int> , std::set<size_t> > & seen_path, std::map<std::pair<int, int> , std::pair<size_t, double> > & seen_indices, std::map<size_t , double> & weight, std::map<size_t , size_t> & cur_pre);
	void finding_successors_on_subref(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors);
	void add_als(size_t & comp, size_t & refacc, size_t & readacc, std::vector<int> & path, std::vector<size_t> & refs, std::string & onref , std::string & onread, const pw_alignment & cur_al , const pw_alignment & nex_al, std::string & from_current, std::string & from_next);
	void make_als(size_t & refacc, size_t & readacc, std::vector<int> & path, std::vector<size_t> & refs, std::string & onref , std::string & onread, std::vector<pw_alignment> & als, size_t & pos, size_t & cur_ref, size_t & nex_ref, std::string & from_current, std::string & from_next, size_t & readid, size_t & end_on_read_pos);
	void get_al_samples(std::string & onref, std::string & onread, size_t & node_length , std::string & refout, std::string & readout, size_t & pos_on_read , bool & GAP);
	void add_to_graph(size_t & comp, std::vector<pw_alignment> & als, const pw_alignment & cur_al , const pw_alignment & nex_al,size_t & refacc, size_t & readacc);
	std::multiset<pw_alignment,sort_pw_alignment_by_left> find_successors(const pw_alignment & p , size_t &);
	void get_subgraph(const pw_alignment &);
	void look_for_successors_on_paths(size_t & comp, const pw_alignment& , std::set<int>& , std::multiset<pw_alignment,sort_pw_alignment_by_left> &);
	void get_previous_als(size_t & comp , int & index, std::vector<pw_alignment> & pre_al, std::vector<size_t> & onread_from,size_t & to,size_t & ref_acc, size_t & read_acc);
	void add_to_containers(size_t & comp , const pw_alignment & p1 , const pw_alignment & p2, size_t & refacc , size_t & readacc, int & current_name , int & next_name);
	void add_nodes_to_the_graph(size_t & type, std::string & seq_from_read, std::string & seq_from_ref, size_t & refacc , size_t & readacc, pw_alignment & nw_alignment, unsigned int & read_id , unsigned int & ref_id,size_t & from, size_t & to);
	void make_al_on_a_node(size_t & comp, const pw_alignment&, const pw_alignment& ,bool dir, unsigned int & ref1, size_t & left1, size_t & right1, unsigned int & ref2, size_t & left2, size_t & right2, size_t & refacc, size_t & readacc);
	void get_paths(const std::set<int> & bfs_node , int & node_name , int & current_node_name, std::set<std::vector<int> > & paths);
	void append_nodes(const std::set<int> & bfs_nodes , int & name, int & current_node_name, unsigned int & next_ref_id, unsigned int & current_ref_id, size_t & refacc , std::vector<std::vector<size_t> >& all_refs, std::vector<std::vector<int> >& all_paths, std::vector<std::vector<std::string> >& all_strings);
	void make_alignments(size_t & comp, size_t & begin_on_read, size_t & end_on_read, const pw_alignment & first_al, const pw_alignment & last_al, std::string & read_out, std::string & ref_out, std::vector<size_t> &this_refs, std::vector<int> & this_path, unsigned int & ref2, size_t & first_length, size_t & last_length,size_t &refacc,size_t&readacc);
	void compute_samples(bool tillend, size_t node_length, size_t & current_pos, size_t & read_counter, size_t & ref_counter,std::string & read_in, std::string & read_out, std::string & ref_in, std::string & ref_out);
	void get_reverse_complement(std::string & sequence_in , std::string & sequence_out);
	void add_first_and_last_als(size_t & i, const pw_alignment & p, size_t & , size_t & , size_t & );
	void looking_for_first_al(size_t & ,const pw_alignment & p,size_t &);
	void looking_for_last_al(size_t & ,const pw_alignment & p,size_t &);
	void make_first_nw_als(size_t & comp, const pw_alignment & p ,std::set<std::vector<int> > & refs, std::string & seq_from_ref, unsigned int & read_id, size_t & read_from , size_t & read_to , size_t & read_acc, size_t & ref_acc);
	void make_last_nw_als(size_t & comp,int& node_name, const pw_alignment & current_al, std::set<int> & nodes_on_paths ,std::string & seq_from_read, unsigned int & ref_id, unsigned int & read_id, size_t & read_from , size_t & read_to , size_t & read_acc, size_t & ref_acc);

	double get_cost(const pw_alignment & p, size_t & acc1 , size_t & acc2);
	void make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout, size_t & left_on_current_node,size_t & refacc, size_t & readacc, size_t & read_id,std::vector<pw_alignment> & first_als, size_t &);
	void make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & left_on_current_node,size_t & refacc,size_t & read_id,size_t & current_node, std::vector<pw_alignment> & first_als);
	void make_last_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & right_on_current_node,size_t & refacc, size_t & read_id, size_t & current_node, std::vector<pw_alignment> & last_als);
	const std::vector<pw_alignment> find_the_best_als(std::vector<std::vector<pw_alignment> >&, size_t & refacc, size_t & readacc)const;
	void add_adjacencies(size_t & comp, const pw_alignment &, const pw_alignment &,size_t & ref_acc, size_t & read_acc);
	void add_edge_to_begin(size_t & comp, const pw_alignment &,size_t & read_acc, size_t & ref_acc);
	void add_edge_to_end(size_t & comp, const pw_alignment &);
	void add_expensive_edges(size_t & comp,size_t & refacc, size_t & readacc);
	void add_to_maf(const pw_alignment & al, std::ofstream & output, bool & firstal);
	const pw_alignment & get_alignment(size_t & , size_t &)const;
	void make_components(size_t & comp, std::vector<std::map<std::pair<size_t,size_t>,double> > & components, std::vector<std::set<size_t> > & comp_nodes, bool & FULLCOMP, size_t & no_comp); //XXX Just started thinking of spliting alignment graph into several components andfind the shortest path on each of these components
	void find_shortest_path_to_each_successors(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors,size_t & refacc, size_t & readacc);
	void heap_dijkstra(DotFibonacciHeap & myHeap, const pw_alignment& current, const pw_alignment& neighbor, std::map<int, int> & cur_pre, int & cur_name, int & next_name, size_t & to_on_read, size_t & refacc, size_t & readacc, std::map<int,const pw_alignment> & als_on_nodes);
	void run_nw_algorithm(DotFibonacciHeap & myHeap, unsigned int read_id, int & startnode, const pw_alignment& neighbor, int & node_name, std::set<int> & seen , std::set<int> & frontier,size_t & refacc, size_t & readacc, size_t & from_on_read, size_t & to_on_read, bool LASTNODE,size_t & ref_from, size_t & ref_to, std::set<int> & nodes_on_paths, std::map<int, double> & distance, std::map<double,std::vector<int> > & distance_node, double & parent_dist,std::map<int, int> & cur_pre, std::map<int,const pw_alignment> & als_on_nodes);
	void get_single_shortest_path(size_t & comp, std::map<int,int> & cur_pre , int & lastnode, int & firstnode, size_t & readacc, size_t & refacc,std::map<int, const pw_alignment> & als_on_nodes);
	private:
		const all_data & data;
	//	deep_first & dfirst;
		ref_graph & rgraph;
		const T & model;
		std::vector<std::multiset<pw_alignment , sort_pw_alignment_by_left> >alignments_left;
		std::vector<std::multiset<pw_alignment , sort_pw_alignment_by_right> >alignments_right;
		std::map<size_t,std::pair<size_t,size_t> > boundries;
		std::vector<std::multiset<pw_alignment,sort_pw_alignment_by_left> > alignments;//Includes components of alignments that have more than MAXGAP distance with each other.
		std::vector<std::multiset<pw_alignment , sort_pw_alignment_by_right> >ordered_als_by_right;
	//	std::vector<std::set<const pw_alignment*> >add_to_start;
	//	std::vector<std::set<const pw_alignment*> >add_to_end;
		std::vector<std::vector<size_t> > shortest_path;
		std::vector< std::map<size_t , std::set<size_t> > >adjacencies; //from_al index , to_al index , vector goes over different components (Includes alignment graph)
		std::vector<std::map<pw_alignment, size_t, compare_pw_alignment> >node_indices;//al and its index
		std::vector<std::map<std::pair<size_t,size_t>,double> >weight_of_edges; //Edges and their weight. In the case of exisitng adjacencies, weight is always the modification of the end node of an edge 
		std::vector<std::map<size_t, const pw_alignment> >indices_nodes;
	//	std::map<std::vector<int> , size_t > alignments_on_each_ref_node; //keeps the index of each alignment is created on each node of the reference graph. int-->index from ref graph , size_t al_indices
		std::map<int , std::set<size_t> > alignments_on_each_ref_node;
	//	std::map<int,const pw_alignment> als_on_nodes;
		size_t index;
		size_t previous_right2;
		size_t previous_right1;
		size_t previous_left1;
		size_t previous_ref;
		bool forward;


};
#endif
	
#include "needleman_wunsch.cpp"
