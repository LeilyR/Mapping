#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP
#include"pw_alignment.hpp"
#include "dynamic_mc.hpp"
#include <limits>
#include <algorithm>
#include <map>

#define MAXGAP 400

class dijkstra{//Edges weights are cost of creating them and nodes weight are cost modifying the ref to the read.//XXX Should be careful about als with artificial references that I created to fill in the gaps!

	public:
	dijkstra(const all_data & d, const std::map<std::pair<size_t, size_t>, double> & edges):data(d){
		this->edges =edges;
		size_t srcnode = 0;
		double dist = std::numeric_limits<double>::infinity();
		std::cout<< "adj size "<< edges.size()<<std::endl;
		for(std::map<std::pair<size_t, size_t>, double>::const_iterator it = edges.begin(); it != edges.end(); it++){
			std::pair<size_t , size_t> adj_nodes = it->first;
			adjacencies.insert(it->first);
			std::map<const size_t, double>::iterator distance = al_distance.find(it->first.first);
			std::cout << it->first.first <<std::endl;
			if(distance == al_distance.end()){
				std::cout << "add "<<std::endl;
				al_distance.insert(std::make_pair(it->first.first, dist));
				if(it->first.first != 0){
					distance_al.insert(std::make_pair(dist,it->first.first));
				}
				unvisited.insert(it->first.first);
			}
			std::cout << it->first.second<<std::endl;
			distance = al_distance.find(it->first.second);
			if(distance == al_distance.end()){
				std::cout << "add! "<<std::endl;
				al_distance.insert(std::make_pair(it->first.second, dist));		
				distance_al.insert(std::make_pair(dist,it->first.second));
				unvisited.insert(it->first.second);
			}
		}
		std::map<const size_t, double>::iterator it = al_distance.find(srcnode);
		assert(it != al_distance.end());
		it->second = 0.0;
		distance_al.insert(std::make_pair(0.0,srcnode));
		std::cout << "distance al size " << distance_al.size() << " al distance size "<< al_distance.size() <<std::endl;
		assert(distance_al.size() == al_distance.size());
		
	}
	~dijkstra(){}
//	void compute_cost(const pw_alignment*, double & cost);
	void add_to_map(const size_t & from_node, const size_t & this_node, double & cost);
	void find_shortest_path(size_t & first_left , size_t & last_right , std::vector<size_t> & shortest_path);
	void add_the_path(std::vector<size_t> & shortest_path);
	private:
	const all_data & data;
	std::multimap<const size_t,const size_t> adjacencies;
	std::map<const size_t, double> al_distance;//distance from src to v
	std::multimap<double, const size_t> distance_al;// distance from src to v of unvisited nodes, we remove those are visited.
	std::set<size_t> unvisited; //it is has all the nodes at the beginning
	std::map<const size_t , size_t> current_and_previous; //current visited node, the node we got to the current node from.
	std::map<std::pair<size_t, size_t>, double> edges;


};
#endif
