#include "dijkstra.hpp"

/*void dijkstra::compute_cost(const pw_alignment* p,double & cost){
	double c1,c2,m1,m2;
	model.cost_function(*p,c1,c2,m1,m2);
	cost += m1;
}*/
void dijkstra::add_to_map(const size_t & from_node, const size_t & this_node, double & cost){
	std::map<const size_t, double>::iterator it= al_distance.find(this_node);
	assert(it != al_distance.end());
	std::cout << "for the node "<< this_node <<"current " << cost << "prev cost "<< it->second <<std::endl;
	if(this_node == 50){
		std::cout << "on 50 before: "<<std::endl;
		for(std::multimap<double, const size_t>::iterator it2=distance_al.begin() ; it2!= distance_al.end(); it2++){
			if(it2->second == 50){
				std::cout<< it2->first << " "<< it2->second <<std::endl;
			}
		}
	}

	if(cost<it->second){
		std::cout << "here" << it->second <<std::endl;
		std::multimap<double, const size_t>::iterator it1=distance_al.find(it->second);	
		assert(it1 != distance_al.end());//TODO
	//	if(it1 != distance_al.end()){
		std::pair<std::multimap<double, const size_t>::iterator, std::multimap<double, const size_t>::iterator> findnode = distance_al.equal_range(it->second);
		bool EXISTS = false;
		for(std::multimap<double, const size_t>::iterator it2 = findnode.first; it2 != findnode.second; it2++){
			if(it2->second == this_node){
				distance_al.erase(it2);
				distance_al.insert(std::make_pair(cost,this_node));
				EXISTS = true;
				break;
			}
	
		}
		assert(EXISTS == true);
	//	distance_al.erase(it1);
	//	distance_al.insert(std::make_pair(cost,this_node));
	//	}else{
	//		distance_al.insert(std::make_pair(cost,this_node));			
	//	}
		it->second = cost;
		std::map<const size_t , size_t>::iterator it2 = current_and_previous.find(this_node);
		if(it2 == current_and_previous.end()){
			std::cout << "got to " << this_node << " from " << from_node<<std::endl;
			current_and_previous.insert(std::make_pair(this_node, from_node));
		}else{
			std::cout << "exists! "<<std::endl;
			std::cout << "got to " << this_node << " from " << from_node<<std::endl;			
			it2->second = from_node;
		}
	}
	if(this_node == 50){
		std::cout << "on 50 after: "<<std::endl;
		for(std::multimap<double, const size_t>::iterator it2=distance_al.begin() ; it2!= distance_al.end(); it2++){
			if(it2->second == 50){
				std::cout<< it2->first << " "<< it2->second <<std::endl;
			}
		}
	}

}
void dijkstra::find_shortest_path(size_t & first_left , size_t & last_right , std::vector<size_t> & shortest_path){
	std::cout<< "unvisited size "<<unvisited.size()<<std::endl;
	while(unvisited.size()!=0){
		std::multimap<double,const size_t>::iterator it = distance_al.begin();
		//Delete it from unvisited and distance_al.
		size_t this_node = it->second;
		double current_cost = it->first;
		std::cout << "this node "<< this_node  << "its cost " << current_cost<<std::endl;
		std::set<size_t>::iterator unseen=unvisited.find(this_node);
		if(unseen != unvisited.end()){
			unvisited.erase(unseen);
		}
		std::cout<< "unvisited size1: "<<unvisited.size()<<std::endl;
		distance_al.erase(it);
		if(this_node!=1 && (unseen != unvisited.end())){
			std::cout << "here!!"<<std::endl;
			std::pair<std::multimap<const size_t, const size_t>::iterator , std::multimap<const size_t, const size_t>::iterator> this_adjs=adjacencies.equal_range(this_node);
			for(std::multimap<const size_t, const size_t>::iterator adj = this_adjs.first ; adj != this_adjs.second; adj++){
				double cost = current_cost;
				/*if(this_node==srcnode){//If the previous one was the begin node
					const pw_alignment* al = adj->second;
					unsigned int read_id = al->getreference2();
					size_t l2,r2;
					al->get_lr2(l2,r2);
					al->print();
					std::cout << l2 << " "<< r2 << std::endl;
					std::cout << "num seq "<< data.numSequences() << " "<< std::endl;
					assert((l2<MAXGAP) || (l2-(first_left-MAXGAP/2)<MAXGAP));
					std::string seq;
					size_t from =0;
					size_t to = l2-1;
					if(l2< MAXGAP){
						std::cout <<"read id " << read_id << " from "<< from << " to "<< to <<std::endl;
						seq = data.extract_seq_part(read_id,from,to);
					}else{
						assert(l2-(first_left-MAXGAP/2)<MAXGAP);
						from = first_left-MAXGAP/2;
						seq = data.extract_seq_part(read_id,from,to);
					}
					cost += model.get_cost(seq,read_id);
					compute_cost(adj->second , cost);
				}else if(adj->second==endnode){//if it is connected to the end node
					unsigned int read_id = p->getreference2();
					size_t seqlength = data.getSequence(read_id).length();
					size_t l2,r2;
					p->get_lr2(l2,r2);
					assert((r2 +MAXGAP>seqlength) || (r2 > (MAXGAP/2 + last_right)-MAXGAP));
					std::string seq;
					size_t from =r2+1;
					size_t to = seqlength-1;
					if(r2+MAXGAP > seqlength){
						seq = data.extract_seq_part(read_id,from,to);						
					}else{
						assert(r2 > (MAXGAP/2 + last_right)-MAXGAP);
						to = MAXGAP/2 + last_right;
						seq = data.extract_seq_part(read_id,from,to);
					}
					cost += model.get_cost(seq,read_id);
				}else{*/
			//	std::set<size_t>::iterator seen=unvisited.find(adj->second);
			//	if(seen != unvisited.end()){
					std::pair<size_t,size_t> this_edge(adj->first,adj->second);
					std::map<std::pair<size_t,size_t>,double >::iterator it = edges.find(this_edge);
					std::cout<< adj->second << " this cost "<<it->second<<std::endl;
					cost+=it->second;
			//	compute_cost(this_edge, cost);
				//}
				//Add the cost to the  distance map if it is smaller than the current one it has, otherwise keep it as it is. Also fill in previous map if its value is changed
					add_to_map(adj->first,adj->second,cost);
			//	}
			}
		}else if(this_node==1){
			for(std::map<const size_t , size_t>::iterator it=current_and_previous.begin() ; it != current_and_previous.end() ; it++){
				std::cout<< it->first << " " << it->second << std::endl;
			}
			break;
		}

	}
	//Add the shortest path to shorest_path
	add_the_path(shortest_path);
}

void dijkstra::add_the_path(std::vector<size_t> & shortest_path){
	size_t node = 1; //it is initiallized with the end node.
	shortest_path.push_back(node);
	std::cout<< "node "<<node << " ";
	std::map<const size_t, size_t>::iterator it = current_and_previous.find(node);
	assert(it != current_and_previous.end());
	while(it!= current_and_previous.end()){
		shortest_path.push_back(it->second);
		std::cout << "pre" <<it->second << " ";
		node = it->second;
		it = current_and_previous.find(node);
	}
	assert(shortest_path.back()==0);
	std::reverse(shortest_path.begin(),shortest_path.end());

}
