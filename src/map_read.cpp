#include "map_read.hpp"

#ifndef MAP_READ_CPP
#define MAP_READ_CPP

template<typename T>
void als_components<T>::find_als_on_paths(std::ofstream & output,size_t & refacc, size_t & readacc){
//	add_to_start.resize(alignments.size());
//	add_to_end.resize(alignments.size());
//	shortest_path.resize(alignments.size());//Shortest path is found for each component
	for(size_t i =0; i < alignments.size(); i++){//Over components of alignments between one specific read and the reference graph
		std::cout << "all als at comp "<< i << " are " << alignments.at(i).size() << std::endl;
		size_t first_left;
		std::multiset<pw_alignment,sort_pw_alignment_by_right>::reverse_iterator revit = ordered_als_by_right.at(i).rbegin(); //finds the last alignment of the component
		pw_alignment last_p = *revit;
		size_t lastr, lastl;
		last_p.get_lr2(lastl,lastr);
		size_t counter = 0;
		for(std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it = alignments.at(i).begin() ; it != alignments.at(i).end(); it++){//For each alignment all of its successors are found
			pw_alignment p = *it;
			p.print();
			std::cout << "al number  "<< counter <<std::endl;
			counter++;
			size_t l,r;
			p.get_lr2(l,r);
			if(it==alignments.at(i).begin()){
				first_left = l;
			}
			//All the successor alignments, ordered by left:
			std::multiset<pw_alignment,sort_pw_alignment_by_left> successors = find_successors(p,i);
			std::cout<<"successors size is "<< successors.size()<<std::endl;//Closeer than MAXGAP to the current one
			if(successors.size() != 0){
			//	std::set<std::vector<int> > paths;
				get_subgraph(p);//deep_first_search is called in the get_paths function!
				std::set<int> nodes = rgraph.get_nodes();//nodes from the ref graph that have less than MAXGAP space to the node that this al is aligned on it
			//	std::set<int> nodes = rgraph.get_subgraph_nodes();
			//	std::cout<< "nodes are "<<std::endl;
			//	for(std::set<int>::iterator it = nodes.begin(); it!=nodes.end();it++){
			//		std::cout<< *it<<std::endl;
			//	}
				if(nodes.size()!=0){
				//	look_for_successors_on_paths(i, p, nodes, successors); //TODO
				//	finding_successors_on_subref(i, p, nodes, successors);
					finding_successors_on_subref(i, p, nodes, successors,refacc,readacc);
				}
			}
			std::cout << "find first and last als"<<std::endl;
			size_t successors_size = successors.size();
			add_first_and_last_als(i,p, successors_size,first_left,lastr);//Make alignment to add the current al to the start and end node, use NW with the right type--> it can start at any cell on the last row and ends at any cell on the first row (type 3)
			std::cout << "adj size is "<< adjacencies.size()<<std::endl;
			for(std::map<size_t, std::set<size_t> >::iterator it = adjacencies.at(i).begin(); it != adjacencies.at(i).end(); it++){
				std::cout << it->first << " connects to: " << std::endl;
				for(std::set<size_t>::iterator s = it->second.begin() ; s != it->second.end() ; s++){
					std::cout << *s << " ";
				}
				std::cout << " "<<std::endl;

			}
		//	if(counter == 31) exit(0);
		}
	//	add_expensive_edges(i,refacc,readacc);//i-->component //XXX Temporary commented out //XXX Replaced it by mapping over components 
		std::vector<std::map<std::pair<size_t,size_t>,double> > components;
		std::vector<std::set<size_t> > comp_nodes;
		std::map<std::pair<size_t,size_t>,double> weight_to_use;
		bool FULLCOMP = false;
		size_t no_comp;
		make_components(i,components, comp_nodes, FULLCOMP , no_comp);
		if(components.size() >1){
			std::cout << "components size is "<< components.size() << "bool " << FULLCOMP <<std::endl;
		}
		if(components.size()> 1 && FULLCOMP == false){
			//TODO complicated case
			add_expensive_edges(i,refacc,readacc);
			weight_to_use = weight_of_edges.at(i);//TODO has to be replaced by something ore thoughtful!! this way will be so slow!

		}else if(components.size() > 1 && FULLCOMP == true){
			weight_to_use = components.at(no_comp);

		}else{
			assert(components.size()==1);
			weight_to_use = weight_of_edges.at(i);

		}
	/*	for(size_t i = 0; i < components.size() ; i++){
			std::cout << "at "<<i << std::endl;
			for(std::map<std::pair<size_t, size_t>, double>::iterator it = components.at(i).begin() ; it != components.at(i).end() ;it++){
				std::cout << it->first.first << " " << it->first.second << std::endl;
			}
		}*/
	/*	std::cout << "adj size is "<< adjacencies.size()<<std::endl;
		for(std::map<size_t, std::set<size_t> >::iterator it = adjacencies.at(i).begin(); it != adjacencies.at(i).end(); it++){
			std::cout << it->first << " to: " << std::endl;
			for(std::set<size_t>::iterator s = it->second.begin() ; s != it->second.end() ; s++){
				std::cout << *s << " ";
			}
			std::cout << " "<<std::endl;
		}*/
		for(std::map<const pw_alignment, size_t , compare_pw_alignment>::iterator it = node_indices.at(i).begin(); it != node_indices.at(i).end();it++){
			std::cout<< it->second << " : " <<std::endl;
			it->first.print();
		}
		for(std::map<std::pair<size_t,size_t>,double>::iterator it = weight_of_edges.at(i).begin(); it != weight_of_edges.at(i).end();it++){
			std::pair<size_t , size_t> this_edge = it->first;
			std::cout<< "from " <<this_edge.first << " to "<<this_edge.second << " costs " << it->second << std::endl;
		}//TODO dijkstra on members of components
		dijkstra dij(data,weight_to_use);
		std::vector<size_t> this_short_path;
	//	dij.find_shortest_path(first_left,lastr,shortest_path.at(i));//Input: begin and end of the component, Output: shortest path
		dij.find_shortest_path(first_left,lastr,this_short_path); //TODO first_left and lastr can be removed
		shortest_path.push_back(this_short_path);
		std::cout<< "shortest path is: " << std::endl;
		double sum = 0.0;
		for(size_t j = 0; j < this_short_path.size(); j++){
			if(this_short_path.at(j) != 0 && this_short_path.at(j) !=1){
				std::multimap<size_t,const pw_alignment>::iterator it=indices_nodes.at(i).find(shortest_path.at(i).at(j));
				assert(it != indices_nodes.at(i).end());
				const pw_alignment p = it->second;
				p.print();
				double m1,m2;
				model.cost_function(p, m1, m2, refacc,readacc);
				sum+=m1;
				std::cout << "m1 "<< m1 << std::endl;
			}
		}		
		std::cout<< "sum is "<< sum <<std::endl;
		for(size_t j = 0; j < this_short_path.size(); j++){
			std::cout << this_short_path.at(j) << " ";
			if(this_short_path.at(j) != 0 && this_short_path.at(j) !=1){
				std::multimap<size_t,const pw_alignment>::iterator it=indices_nodes.at(i).find(this_short_path.at(j));
				assert(it != indices_nodes.at(i).end());
				const pw_alignment p = it->second;	

			/*	if(p.getreference1()== data.numSequences()+1 && j != shortest_path.at(i).size()-1){//if it is only gap
					std::multimap<size_t,const pw_alignment>::iterator it1=indices_nodes.at(i).find(shortest_path.at(i).at(j+1));
					assert(it != indices_nodes.at(i).end());
					const pw_alignment post_p = it1->second;

				}else{
					if(GapOnly == true){
						size_t begin, end ;
						std::string sample1 = pre_al.get_al_ref1() + p.get_al_ref1();
						std::string sample2 = pre_al.get_al_ref2() + p.get_al_ref2();
						if(p.getbegin2()<p.getend2()){
							assert(pre_p.getbegin2() < pre_p.getend2());
							begin = pre_p.getbegin2();
							end = p.getned2();
						}else{
							assert(pre_p.getbegin2() > pre_p.getend2());
							assert(p.getbegin2()>p.getend2());
							begin = p.getbegin2();
							end = pre_p.getend2();
								
						}
						const pw_alignmet al(sample1, sample2,p.getbegin1(), begin ,p.getend1(), end ,p.getreferencde1(),p.getreference2());
						add_to_maf(al,output);
						GapOnly = false;
					}else{*/
						//Write it to the output file
						bool firstal = false;
						if(j == 1){
							firstal = true;
							if(p.getreference1()== data.numSequences()+1){
								std::cout << "gap at begin! "<<std::endl;
								assert(j < shortest_path.at(i).size()-2);
								std::multimap<size_t,const pw_alignment>::iterator it1=indices_nodes.at(i).find(shortest_path.at(i).at(j+1));
								assert(it != indices_nodes.at(i).end());
								const pw_alignment post_p = it1->second;
								std::string sample1 = p.get_al_ref1() + post_p.get_al_ref1();
								std::string sample2 = p.get_al_ref2() + post_p.get_al_ref2();
								size_t begin = p.getbegin2();
								size_t end = post_p.getend2();
								assert(p.getreference2() == post_p.getreference2());
								const pw_alignment al(sample1, sample2,post_p.getbegin1(), begin ,post_p.getend1(), end ,post_p.getreference1(),p.getreference2());
								add_to_maf(al,output,firstal);
								j++;

							}else{
								add_to_maf(it->second,output,firstal);
							}
						}else{
							add_to_maf(it->second,output,firstal);
						}
				//	}
			//	}
			}
		}
		index =2;
	}

}

//XXX Creatig the best path rather than all the poosible ones to reduce the running time.
template<typename T>
void als_components<T>::finding_successors_on_subref(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors,size_t & refacc, size_t & readacc){
	//Information of the current alignment:
	size_t cur_l1,cur_r1,cur_l2,cur_r2;
	current_al.get_lr1(cur_l1,cur_r1);	
	current_al.get_lr2(cur_l2,cur_r2);
	unsigned int cur_ref1 = current_al.getreference1();
	unsigned int cur_ref2 = current_al.getreference2();
	std::string ref_accession = data.get_acc(refacc);
	size_t seqlength = data.get_seq_size(cur_ref1);//We use this seq from its r1+1 to length-1. It is part of the ref node that comes after current alignment.
	std::string current_ref_node_name = data.get_seq_name(current_al.getreference1());
	int current_name = std::stoi(current_ref_node_name);
	std::cout << "current node name "<< current_name << std::endl;
	std::string from_current_node;
	size_t from = cur_r1+1;
	size_t to = seqlength-1;
	from_current_node = data.extract_seq_part(cur_ref1,from,to);
	if(current_al.getbegin1() > current_al.getend1()){//If reverse
		std::string temp;
		get_reverse_complement(from_current_node,temp);
		from_current_node = temp;
		current_name = -1*current_name;
	}
	std::cout << "to "<< to << " from "<< from << " seqlength "<< seqlength << std::endl;
	size_t current_remainder = to-from +1;
	std::cout << "current remainder is "<< to-from+1 << std::endl;
	std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator id = node_indices.at(comp).find(current_al);
	if(id == node_indices.at(comp).end()){
		node_indices.at(comp).insert(std::make_pair(current_al,index));
		indices_nodes.at(comp).insert(std::make_pair(index,current_al));
		index++;
		id = node_indices.at(comp).find(current_al);
	}
	size_t cur_index = id->second;
//	std::set<std::vector<int> > all_p = rgraph.get_paths();
	//Make a seen container to not go over the same path for the same part of read again and again
	std::map<std::pair<int, int> , std::set<size_t> > seen_path; //'int's are adjacent node, 'size_t's are positions on the read
	std::map<std::pair<int,int> , std::pair<size_t , double> > seen_indices; //TODO second pair should change to a set of pairs to be able to handle the loops
	std::map<size_t , double>  weight;
	std::map<size_t , size_t> cur_pre;
	std::multimap<double, size_t> al_distance;
	for(std::multiset<pw_alignment, sort_pw_alignment_by_left>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
		//Information of a neighbor alignment:
		size_t nex_l1, nex_r1 ,nex_l2, nex_r2;
		const pw_alignment nex_al = *it;
		nex_al.get_lr1(nex_l1,nex_r1);
		nex_al.get_lr2(nex_l2,nex_r2);
		unsigned int nex_ref1 = nex_al.getreference1();

		std::string nex_seqname = data.get_seq_name(nex_al.getreference1());
		int nex_name = std::stoi(nex_seqname);
		if(nex_al.getbegin1()>nex_al.getend1()){
			nex_name = -1*nex_name;
		}
		std::cout << "neighbor's name is " << nex_name <<std::endl;
		nex_al.print();
		std::set<int>::iterator it1 = nodes.find(nex_name);
		if(it1 != nodes.end() && nex_name != current_name){//If they are neighbor but are on two separate node //XXX check the length of seq on reference graph in between them if > MAXGAP break;
			size_t last_pos_on_read = nex_l2 -1;
			size_t last_pos_on_ref = nex_l1;
			if(nex_l1 >= MAXGAP) continue;
			assert(nex_l1 < MAXGAP);
			std::set<std::vector<int> > all_p = rgraph.get_paths();
			std::set<std::vector<int> > PATH;
			std::set<int> nodes_on_paths;
			std::cout << "all_p size: "<< all_p.size()<<std::endl;
			for(std::set<std::vector<int> >::iterator it = all_p.begin(); it!= all_p.end() ; it++){
			//	size_t length = from_current_node.length();
			//	std::cout << "l "<< length << std::endl;
				std::vector<int> temp = *it;
				std::vector<int> this_p;
				for(size_t j = 0; j < temp.size(); j++){
					this_p.push_back(temp.at(j));
				//	std::string seq_name = rgraph.seqname(temp.at(j)); //XXX Good to have it for less complex graphs
				//	if(j != 0){
				//		length += rgraph.seq_length(seq_name, ref_accession);
				//	}
				//	std::cout << "+l "<< length << std::endl;
					if(temp.at(j)==nex_name){
				//		length -= rgraph.seq_length(seq_name, ref_accession);
				//		length += nex_l1;
				//		std::cout << "-+l "<< length << std::endl;
				//		if(length <= MAXGAP){
						//	if(nex_name == 1068){
						//		std::cout << "this_p: " <<this_p << std::endl;
						//		std::cout << "---------------------"<<std::endl;
						//	}
							PATH.insert(this_p);
							for(size_t i = 0; i < this_p.size();i++){
								nodes_on_paths.insert(this_p.at(i));
							}
				//		}
						break;
					}
				}
			} 
			//Make an for the remainder,add it to the current_al and use it as a previous alignment for the neigbors
			pw_alignment this_al = current_al;
			if(from_current_node.length() !=0){
				size_t from_on_read = cur_r2+1;
				size_t to_on_read = nex_l2-1;
				std::string readin = data.extract_seq_part(cur_ref2,from_on_read,to_on_read); 
				needleman<T> nl(data, model, readin, from_current_node);
			//XXX	simpleNW nl(model, readin, from_current_node);

				size_t type = 5;
				std::string readout, refout;
				nl.run_needleman(readacc,refacc,type,readout,refout);
				size_t readcounter = 0;
				for(size_t i = 0; i < readout.size();i++){
					if(readout.at(i)!= '-') readcounter++;
				}
				if(readcounter>0 ){
					to_on_read = from_on_read + readcounter -1;
				}else{
					to_on_read = data.getSequence(cur_ref2).length();
				}
				if(current_al.getbegin1() < current_al.getend1()){
					pw_alignment nwal(refout,readout,from,from_on_read,to,to_on_read,cur_ref1,cur_ref2);
					this_al = nwal;
				}else{
					pw_alignment nwal(refout,readout,to,from_on_read,from,to_on_read,cur_ref1,cur_ref2);
					this_al = nwal;					
				}
				add_adjacencies(comp, current_al , this_al , refacc, readacc);
				std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator it = node_indices.at(comp).find(this_al);
				assert(it != node_indices.at(comp).end());
				cur_index = it->second;
				if(to_on_read == nex_l2-1){
					std::cout << "end of read is reached! "	<<std::endl;
				} //We can't stop here because we want to reach to the nex_al!!!
			}
		//	std::cout << "PATH size "<< PATH.size()<<std::endl; 
		//	if(PATH.size()==1){
		//		std::set<std::vector<int> >::iterator test_it = PATH.begin();
		//		std::cout << *test_it <<std::endl;
		//	}
		//	std::set<int> nodes_on_paths = rgraph.get_nodes_on_paths(nex_name);
			if(nodes_on_paths.size() != 0){
		//	if(nodes.size() != 0){
				find_best_path_to_neighbour(comp,nex_name,current_name, cur_index , this_al,nex_al,last_pos_on_read, last_pos_on_ref,refacc ,readacc, nodes_on_paths,seen_path,seen_indices,weight, cur_pre);
			}
		}
		if(it1 != nodes.end() && nex_name == current_name && cur_r1<nex_l1 && nex_l1-cur_r1<=MAXGAP & nex_name > 0){
			std::cout << "on the same node "<< std::endl;
			//if both are forwards or both are backwards
			size_t from = cur_r1+1;
			size_t to = nex_l1-1;
			size_t onread_from = cur_r2+1;
			size_t onread_to = nex_l2 -1;
			std::cout<<"when they share a node forward"<<std::endl;
			current_al.print();
			nex_al.print();
			if((cur_r2 == nex_l2-1) && (cur_r1 != nex_l1-1)){
				std::cout<<"only gap on read "<<std::endl;
				//only gap on read 
				onread_from = cur_r2;
				onread_to = data.get_seq_size(cur_ref2);
				make_al_on_a_node(comp,current_al,nex_al,true,cur_ref1,from,to,cur_ref2,onread_from,onread_to,refacc,readacc);
			}
			else if((cur_r2 != nex_l2-1) && (cur_r1 == nex_l1-1)){
				unsigned int temp = data.numSequences()+1;//only gap on ref
				from = 0;
				to = 0;
				make_al_on_a_node(comp,current_al,nex_al,true,temp,from, to,cur_ref2,onread_from,onread_to,refacc,readacc);
			}
			else if((cur_r2 != nex_l2-1) && (cur_r1 != nex_l1-1)){
				make_al_on_a_node(comp,current_al,nex_al,true,cur_ref1,from,to,cur_ref2,onread_from,onread_to,refacc,readacc);
			}else{
				assert((cur_r2 == nex_l2-1) && (cur_r1 == nex_l1-1));
				add_adjacencies(comp,current_al,nex_al,refacc,readacc);
			}
		}
		if(it1 != nodes.end() && nex_name == current_name && nex_r1<cur_l1 && cur_l1-nex_r1<=MAXGAP & nex_name < 0){
			std::cout << "negatively share a node! "<<std::endl;
		}

	}
/*	for(std::map<size_t , std::set<size_t> >::iterator it =adjacencies.at(comp).begin() ; it != adjacencies.at(comp).end(); it++){
		std::cout << "from "<< it-> first << " to " <<std::endl;
		for(std::set<size_t>::iterator it1 = it->second.begin() ; it1 != it->second.end() ;it1++){
			std::cout << *it1 << " ";
		}
	}
	exit(0);*/
}
template<typename T>
void als_components<T>::find_best_path_to_neighbour(size_t comp, int & next_name, int & current_name, size_t & cur_index, const pw_alignment & cur_al, const pw_alignment & neighbour_al , size_t & to_on_read, size_t & last_pos_on_ref , size_t & refacc, size_t & readacc, std::set<int> & nodes_on_paths, std::map<std::pair<int, int> , std::set<size_t> > & seen_path, std::map<std::pair<int, int> , std::pair<size_t, double> > & seen_indices, std::map<size_t , double> & weight, std::map<size_t , size_t>  & cur_pre){
	std::set<size_t> seen_als; //It includes the al indices, if they are seen
//	std::map<size_t , size_t> cur_pre;
//	std::map<size_t , double> weight; //al index and it is accumulative weight from the beginning of the path
	weight.insert(std::make_pair(cur_index, 0));
	std::multimap<double, size_t> al_distance; //Double --> accumulative gain, size_t --> node index
	pw_alignment this_al = cur_al;
	unsigned int read_id = this_al.getreference2();
	int startnode = current_name;
	size_t from_on_read;
	size_t pre_from_on_read = data.getSequence(read_id).length();
	std::set<int> visited; //For loops, may delete it later if seen_als handels all the situations
	//Create the first al from the remaining part of the current node!

	while(startnode!=next_name){
		std::cout << "start node "<< startnode <<std::endl;
		std::cout << "al distance size "<< al_distance.size()<<std::endl;
		std::set<int> adjs = rgraph.get_adjacencies(startnode); 
		assert(adjs.size() != 0); //XXX if we dynamically erase the path then it should be removed!
		std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator id = node_indices.at(comp).find(this_al);
		assert(id != node_indices.at(comp).end());
		std::set<size_t>::iterator find_seens = seen_als.find(id->second);
		if(find_seens==seen_als.end()){
		//	startnode = al_distance.begin()->second;
		//	al_distance.erase(al_distance.begin());
		//	continue;
		//}
		seen_als.insert(id->second);
		std::map<size_t , double>::iterator pre_wei = weight.find(id->second);
	//	if(pre_wei == weight.end()){
	//		std::cout << id->second;
	//		this_al.print();
	//		std::cout << "-----"<<std::endl;
	//	}
		assert(pre_wei != weight.end());
		size_t l2,r2;
		this_al.get_lr2(l2,r2);
		if(r2 != data.getSequence(read_id).length()){
			from_on_read = r2 + 1;
		}else{
			from_on_read = l2;
		}
	/*	if(from_on_read == to_on_read +1){
			END_OF_READ_IS_REACHED = true;
			visited.insert(startnode);//TODO
		}*/
		std::cout << "from on read "<< from_on_read << std::endl;
		std::cout << "to on read "<<to_on_read <<std::endl;
	/*	bool END_REACHED = false;
		for(std::set<int>::iterator it = adjs.begin() ; it != adjs.end() ; it++){ //TODO if *it == next_name
			int this_name = *it;
			if(this_name == next_name){
				std::cout << "end reached! "<< startnode << " "<< this_name <<std::endl;
				startnode = next_name;
				END_REACHED = true;
				break;				
			}
		}
		if(END_REACHED == false){*/
		for(std::set<int>::iterator it = adjs.begin() ; it != adjs.end() ; it++){ //TODO if *it == next_name
			std::cout << "this adj " << *it <<std::endl;
			std::set<int>::iterator vis = visited.find(*it);
			if(vis != visited.end()){//Is used for handeling loops
				std::cout << "visited! "<<std::endl;
				continue;
			}
		//	std::vector<int> this_path;
		//	this_path.push_back(startnode);
		//	this_path.push_back(*it);
		//	rgraph.delete_path(this_path); XXX temp
			bool REMAINDER = true;
			int this_name = *it; //If this adjacent node is on a path that doesn't end to next_name we skip it! 
		/*	if(this_name == 31890 || this_name == 31891){
				std::cout << "nodes on path contains: "<<std::endl;
				for(std::set<int>::iterator it = nodes_on_paths.begin(); it != nodes_on_paths.end() ; it++){
					std::cout << *it << std::endl;
				}
			}*/
			std::set<int>::iterator onPaths = nodes_on_paths.find(this_name);
			if(onPaths == nodes_on_paths.end()){
				std::cout << "not on path! "<<std::endl;
				continue;
			}
			bool DIFFERENT_POS = false;
			std::map<std::pair<int, int>, std::set<size_t> >::iterator seen = seen_path.find(std::make_pair(startnode,*it));
		 	if(seen != seen_path.end()){ //position on read is checked
				std::set<size_t>::iterator pos = seen->second.find(from_on_read);
				if(pos != seen->second.end()){ 
					if(from_on_read == to_on_read +1 || (r2 == data.getSequence(read_id).length())){//To deal with loops
				//	if(from_on_read == to_on_read +1){//To deal with loops
						std::cout << "add to visited "<< *it <<std::endl;
						visited.insert(*it);
					}
					std::cout << "has been seen " << startnode << " to " << *it <<std::endl;
					std::map<std::pair<int, int> , std::pair<size_t, double> >::iterator seen_weight = seen_indices.find(seen->first);
					if(seen_weight == seen_indices.end()){
						std::cout << "is not seen "<< seen->first.first << " " << seen->first.second << std::endl;
						continue;
					}
					assert(seen_weight != seen_indices.end());
					al_distance.insert(std::make_pair(seen_weight->second.second,seen_weight->second.first)); //if equal to end!
					cur_pre.insert(std::make_pair(seen_weight->second.first,id->second));
				//	weight.insert(std::make_pair(seen_weight->second.second,seen_weight->second.first));
				}else{
					DIFFERENT_POS = true;
					std::cout << "DIF_POS "<< from_on_read<<std::endl;
				}
			}
			if(seen == seen_path.end() || DIFFERENT_POS ==true){
				std::cout << "HERE! " << from_on_read<<std::endl;
				if(seen == seen_path.end() && this_name != next_name){
					std::pair<int,int> temp(startnode,*it);
					seen_path.insert(std::make_pair(temp, std::set<size_t>()));
					seen =  seen_path.find(std::make_pair(startnode,*it));
				}
				if(this_name != next_name){
					seen->second.insert(from_on_read);
					std::cout << "add to seen " << startnode << " to " << *it <<std::endl;
				}
	
				std::cout << "this adj! " << this_name <<std::endl;
				std::string refin;
				size_t ref_from = 0;
				size_t ref_to = 0;
				unsigned int ref_id;
			//	assert(this_name != next_name);
				if(this_name > 0){
					ref_id = rgraph.get_refid(refacc,this_name);
					if(this_name != next_name){
						ref_to = data.getSequence(ref_id).length()-1;
					}else{
						if(last_pos_on_ref >0){
							ref_to = last_pos_on_ref - 1;
						}else{
							REMAINDER = false;
							ref_to = 0;
						}
					}
				//	if(NO_REMAINDER == false || this_name != next_name){
					if(REMAINDER == true){
						refin = data.extract_seq_part(ref_id, ref_from, ref_to);
					}

				}else{
					int temp = -1*(this_name);
					ref_id = rgraph.get_refid(refacc,temp);
					if(this_name != next_name){
						ref_to = data.getSequence(ref_id).length()-1;
					}else{
						if(last_pos_on_ref >0){
							ref_to = last_pos_on_ref - 1;
						}else{
							REMAINDER = false;
						}
					}
				//	if(NO_REMAINDER == false || this_name != next_name){
					if(REMAINDER == true){
						std::string temp_seq_in = data.extract_seq_part(ref_id, ref_from, ref_to);
						std::string temp_seq_out;
						get_reverse_complement(temp_seq_in, temp_seq_out);
						refin = temp_seq_out;
					}
					size_t temp_ref_from = ref_from;
					ref_from = ref_to;
					ref_to = temp_ref_from;
				}
				std::string readin = data.extract_seq_part(read_id,from_on_read,to_on_read);
				std::string refout;
				std::string readout;
				pw_alignment nextal;
				if(this_name != next_name ||REMAINDER == true ){
					if(refin.length()>MAXGAP){
						std::cout << "reing length "<< refin.length() <<std::endl;
						continue;
					}
					assert(refin.length() <= MAXGAP);
					needleman<T> nl(data, model, readin, refin);
				//XXX	simpleNW nl( model, readin, refin);
					size_t type = 3;
					if(this_name == next_name) type = 1;
					nl.run_needleman(readacc,refacc,type,readout,refout);
					size_t readcounter = 0;
					for(size_t i = 0; i < readout.size();i++){
						if(readout.at(i)!= '-') readcounter++;
					}
					size_t this_to_on_read;
					if(readcounter>0 ){
						this_to_on_read = from_on_read + readcounter -1;
						assert(this_to_on_read <= to_on_read);
					}else{
						this_to_on_read = data.getSequence(read_id).length();
						if(from_on_read == to_on_read +1 ||( r2 == data.getSequence(read_id).length() ) ){//To deal with loops
							std::cout << "add to visited "<< *it <<std::endl;
							visited.insert(*it);
						}
					}
					pw_alignment nwal(refout,readout,ref_from,from_on_read,ref_to,this_to_on_read,ref_id,read_id);
					nextal = nwal;
					if(this_name == next_name){
						add_adjacencies(comp,nextal, neighbour_al , refacc, readacc);	
						//add the last al to adj
						std::cout << "add the last one! "<<std::endl;
					}
				}else{
					std::cout<<"shouldn't happen! "<<std::endl;
				//	exit(0);
					assert(this_name == next_name);
					assert(REMAINDER == false); 
					if(readin.length()==0){
						nextal = neighbour_al;
					}else{//Make an with only gap on ref(the neighbour node)
						std::string onlygap(readin.length(),'-');
						std::cout << "readin "<<readin << std::endl;
						std::cout<<"gap on ref graph "<<std::endl;
						pw_alignment onlygap_al(onlygap,readin,0,from_on_read,0,to_on_read,data.numSequences()+1,read_id);	
						nextal = onlygap_al;	
						add_adjacencies(comp,nextal, neighbour_al , refacc, readacc);		
					}
				}
				add_adjacencies(comp, this_al , nextal , refacc, readacc);//TODO I can do it locally in here, and later on which ever is chosen be added to the alignment graph		
				std::cout << "nextal "<<std::endl;
				nextal.print();
				std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator nwal_id = node_indices.at(comp).find(nextal);
				assert(nwal_id != node_indices.at(comp).end());
				std::map<size_t , double>::iterator wei = weight.find(nwal_id->second);
				if(wei == weight.end()){
					std::map<size_t,size_t>::iterator pre = cur_pre.find(nwal_id->second);
					assert(pre == cur_pre.end());
					pre = cur_pre.find(id->second);
//					assert(pre != cur_pre.end());
					double m1,m2;
					model.cost_function(nextal, m1, m2, refacc,readacc);
					double this_weight = m1 + pre_wei->second;
					weight.insert(std::make_pair(nwal_id->second,this_weight));
					cur_pre.insert(std::make_pair(nwal_id->second, id->second));
					al_distance.insert(std::make_pair(this_weight,nwal_id->second));
					std::cout << "add to seen indices "<< startnode << " " << *it << " "<< nwal_id->second << " " << this_weight << std::endl;
					seen_indices.insert(std::make_pair(std::pair<int, int>(startnode,*it) , std::pair<size_t , double>(nwal_id->second,this_weight)));
				}else{
				//	size_t WEIGHT = wei->second;
					std::multimap<double, size_t>::iterator test = al_distance.find(wei->second);
					if(test == al_distance.end()){
						al_distance.insert(std::make_pair(wei->second,nwal_id->second));						
					}else{ //TODO check if the second on is not the the same index and add it!

					}
					std::map<size_t,size_t>::iterator cur = cur_pre.find(nwal_id->second);
					if(cur == cur_pre.end()) std::cout << nwal_id->second << std::endl;
					assert(cur != cur_pre.end());
//					std::map<size_t,size_t>::iterator pre = cur_pre.find(id->second);
//					assert(pre != cur_pre.end());
					double m1,m2;
					model.cost_function(nextal, m1, m2, refacc,readacc);
					double this_weight = m1 + pre_wei->second;
					std::cout << "this weight "<< this_weight << " pre-weight "<< wei->second <<std::endl;
					if(this_weight < wei->second){
						size_t WEIGHT = wei->second;
						std::multimap<double, size_t>::iterator test = al_distance.find(wei->second);
						assert(test != al_distance.end());
						std::pair<std::multimap<double, size_t>::iterator, std::multimap<double, size_t>::iterator> findnode = al_distance.equal_range(wei->second);
						bool EXISTS = false;
						for(std::multimap<double, size_t>::iterator it2 = findnode.first; it2 != findnode.second; it2++){
							if(it2->second == nwal_id->second){
								al_distance.erase(it2);
								al_distance.insert(std::make_pair(this_weight,nwal_id->second));
								std::cout << "add to seen indices1 "<< startnode << " " << *it << " "<< nwal_id->second << " " << this_weight << std::endl;

	
								seen_indices.insert(std::make_pair(std::pair<int, int>(startnode,*it) , std::pair<size_t , double>(nwal_id->second,this_weight)));
								EXISTS = true;
								break;
							}
						}
						assert(EXISTS == true);
						cur->second = id->second;
						wei->second = this_weight;
					}else{
						std::cout << "weight is bigger! "<< this_weight<<std::endl;
					}
				}
			}
		}
}
		//Choosing the shortest path:
		std::cout << "al distance is "<<std::endl;
		for(std::multimap<double, size_t>::iterator test = al_distance.begin() ; test != al_distance.end(); test++){
			std::cout << "Weight "<< test->first << " index " << test->second <<std::endl;
		}
		std::cout << al_distance.begin()->first << " " << al_distance.begin()->second <<std::endl;
		assert(al_distance.size()>0);
		std::map<size_t, const pw_alignment>::iterator ID = indices_nodes.at(comp).find(al_distance.begin()->second);
		std::cout<< "the smallest is: " << al_distance.begin()->first<<std::endl;
		ID->second.print();
		assert(ID != indices_nodes.at(comp).end());
		this_al = ID->second;
		size_t ind = this_al.getreference1();
		if(ind != data.numSequences()+1){
			std::string tempstr = data.get_seq_name(ind);
			std::stringstream str;
			str<< tempstr;
			size_t temp;
			str>> temp;
			if(this_al.getbegin1()<this_al.getend1()){
				startnode = temp;
			}else if(this_al.getbegin1()> this_al.getend1()){
				startnode = -1 * temp;
			}else{
				assert(this_al.getbegin1()== this_al.getend1());
				unsigned int ref1 = this_al.getreference1();
				size_t b1 = this_al.getbegin1();
				size_t e1 = this_al.getend1() ; 
				std::cout << " from al " << this_al.get_al_ref1() << " from seq "<< data.extract_seq_part(ref1, b1 , e1) <<std::endl;
				std::string seq_from_al;
				for(size_t i =0; i < this_al.get_al_ref1().length();i++){
					if(this_al.get_al_ref1().at(i) != '-'){
						seq_from_al += this_al.get_al_ref1().at(i);
						break;
					}
				}
				if(seq_from_al == data.extract_seq_part(ref1, b1 , e1)){
					startnode = temp;
				}else{
					startnode = -1 * temp;
				}

			}
		}else{
			startnode = next_name;
		}
		al_distance.erase(al_distance.begin());
		pre_from_on_read == from_on_read;
//	}
	}
	//TODO The last al between read and and the remaining of the last node and the last node itself
	assert(startnode ==next_name);
	std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator id = node_indices.at(comp).find(this_al);
	assert(id != node_indices.at(comp).end());
	std::map<size_t , double>::iterator pre_wei = weight.find(id->second);
	assert(pre_wei != weight.end());
	size_t l2,r2;
	this_al.get_lr2(l2,r2);
	if(r2 != data.getSequence(read_id).length()){
			from_on_read = r2 + 1;
		}else{
			from_on_read = l2;
		}

/*	std::cout << "from on read "<< from_on_read << std::endl;
	std::cout << "to on read "<<to_on_read <<std::endl;
	assert(from_on_read <= to_on_read+1);
	size_t ref_from = 0;
	size_t ref_to = 0;
	pw_alignment remainder_al;
	if(last_pos_on_ref >0){
		std::cout <<"HERE!! " << last_pos_on_ref <<std::endl;
		ref_to = last_pos_on_ref - 1;	
		std::string refin;
		unsigned int ref_id;
		if(next_name>0){
			ref_id = rgraph.get_refid(refacc,next_name);
			refin = data.extract_seq_part(ref_id, ref_from, ref_to);			
		}else{
			int temp = -1 * next_name;
			ref_id = rgraph.get_refid(refacc,temp);
			std::string temp_seq_in = data.extract_seq_part(ref_id, ref_from, ref_to);
			std::string temp_seq_out;
			get_reverse_complement(temp_seq_in, temp_seq_out);
			refin = temp_seq_out;
		}
		if(from_on_read<=to_on_read){
			std::string readin = data.extract_seq_part(read_id,from_on_read,to_on_read);
			needleman<dynamic_mc_model> nl(data, model, readin, refin);
			size_t type = 1; //TODO check!
			std::string readout, refout;
			nl.run_needleman(readacc,refacc,type,readout,refout);
			pw_alignment nwal(refout,readout,ref_from,from_on_read,ref_to,to_on_read,ref_id,read_id);
			add_adjacencies(comp, this_al , nwal , refacc, readacc);			
			//weight!
			remainder_al = nwal;
		}else{//GAP only on read
			std::string onlygap(refin.length(),'-');
			size_t length = data.getSequence(read_id).length();
			pw_alignment onlygap_al(refin,onlygap,ref_from,from_on_read,ref_to,length,ref_id,read_id);
			add_adjacencies(comp, this_al , onlygap_al , refacc, readacc);			
			//weight!
			remainder_al = onlygap_al;
		}
		add_adjacencies(comp, remainder_al , neighbour_al , refacc, readacc);					
	}else{
		if(r2< to_on_read){ //GAP on ref only!
			std::string readin = data.extract_seq_part(read_id,from_on_read,to_on_read);
			std::string onlygap(readin.length(),'-');
			pw_alignment onlygap_al(onlygap,readin,0,from_on_read,0,to_on_read,data.numSequences()+1,read_id);	
			add_adjacencies(comp, this_al , onlygap_al , refacc, readacc);			
			add_adjacencies(comp, onlygap_al , neighbour_al , refacc, readacc);	

		}else{
			std::cout<< "THERE!! "<<std::endl;
			add_adjacencies(comp, this_al , neighbour_al , refacc, readacc);	
		}				
	}*/
	//	if(add_to_last == true){ //TODO Check the shortest path first and print it out, then if the last al is not the neighbour al add it to the path!
//
//	}
}
template<typename T>
void als_components<T>::finding_successors_on_subref(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors){
	//Information of the current alignment:
	size_t cur_l1,cur_r1,cur_l2,cur_r2;
	current_al.get_lr1(cur_l1,cur_r1);	
	current_al.get_lr2(cur_l2,cur_r2);

	unsigned int cur_ref1 = current_al.getreference1();
	unsigned int cur_ref2 = current_al.getreference2();
	size_t seqlength = data.get_seq_size(cur_ref1);//We use this seq from its r1+1 to length-1. It is part of the ref node that comes after current alignment.
	size_t readacc = data.accNumber(cur_ref2);
	size_t refacc = data.accNumber(cur_ref1);
	std::string current_ref_node_name = data.get_seq_name(current_al.getreference1());
	int current_name = std::stoi(current_ref_node_name);
	std::cout << "current node name "<< current_name << std::endl;
	std::string from_current_node;
	size_t from = cur_r1+1;
	size_t to = seqlength-1;
	from_current_node = data.extract_seq_part(cur_ref1,from,to);
	if(current_al.getbegin1() > current_al.getend1()){//If reverse
		std::string temp;
		get_reverse_complement(from_current_node,temp);
		from_current_node = temp;
		current_name = -1*current_name;
	}
	size_t current_remainder = to-from +1;
	std::cout << "current remainder is "<< to-from+1 << std::endl;
	for(std::multiset<pw_alignment, sort_pw_alignment_by_left>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
		//Information of a neighbor alignment:
		size_t nex_l1, nex_r1 ,nex_l2, nex_r2;
		const pw_alignment nex_al = *it;
		nex_al.get_lr1(nex_l1,nex_r1);
		nex_al.get_lr2(nex_l2,nex_r2);
		unsigned int nex_ref1 = nex_al.getreference1();

		std::string nex_seqname = data.get_seq_name(nex_al.getreference1());
		int nex_name = std::stoi(nex_seqname);
		if(nex_al.getbegin1()>nex_al.getend1()){
			nex_name = -1*nex_name;
		}
		std::cout << "neighbor's name is " << nex_name <<std::endl;
		nex_al.print();
		std::set<int>::iterator it1 = nodes.find(nex_name);
		if(it1 != nodes.end() && nex_name != current_name){//If they are neighbor but are on two separate node
			std::string from_next_node = "";
			from = 0;
			if(nex_l1 > 0){
				to = nex_l1-1;
				from_next_node = data.extract_seq_part(nex_ref1,from,to);//It could be that there are more than one node distance in between!!! they are all added.
				if(nex_al.getbegin1() > nex_al.getend1()){//If reverse
					std::cout<< "ref is reverse"<<std::endl;
					std::string temp;
					get_reverse_complement(from_next_node,temp);
					from_next_node = temp;
				}
			}
			if(from_next_node.length()>MAXGAP) continue;
			from = cur_r2 +1;
			to = nex_l2 -1;
			std::string on_read = data.extract_seq_part(cur_ref2,from,to);
			if(on_read.length() == 0 && from_next_node.length()== 0 && from_current_node.length()==0){
				add_adjacencies(comp, current_al, nex_al, refacc, readacc);
				continue;
			}
			std::vector<std::vector<size_t> >all_refs;
			std::vector<std::vector<int> >all_paths;
			std::vector<std::vector<std::string> >all_strings_from_ref_graph; //each string is content of a node on a path
			append_nodes(nodes,nex_name,current_name,nex_ref1, cur_ref1, refacc,all_refs,all_paths,all_strings_from_ref_graph);//add all the nodes to a string 
			assert(all_strings_from_ref_graph.size()==all_paths.size());
			assert(all_strings_from_ref_graph.size() == all_refs.size());
			assert(all_strings_from_ref_graph.size() > 0);
			for(size_t i = 0; i <all_strings_from_ref_graph.size();i++){
				std::string on_ref = from_current_node;
				for(size_t j = 0 ; j < all_strings_from_ref_graph.at(i).size();j++){//TODO if it worked, change it to string instead of a vector of string
					on_ref += all_strings_from_ref_graph.at(i).at(j);
				}
				on_ref += from_next_node;
				if(on_ref.length() > MAXGAP) continue;
				//Remove from the dynamic_adj in ref graph class:
				for(size_t j =0; j < all_paths.at(i).size()-1; j++){
					std::vector<int> this_edge;
					this_edge.push_back(all_paths.at(i).at(j));
					this_edge.push_back(all_paths.at(i).at(j+1));
					rgraph.delete_path(this_edge); 
				}
				if(from_current_node.length()==0){
					all_paths.at(i).erase(all_paths.at(i).begin());
					all_refs.at(i).erase(all_refs.at(i).begin());
				}
				if(from_next_node.length() == 0){
					all_paths.at(i).pop_back();
					all_refs.at(i).pop_back();
				}
				add_als(comp,refacc,readacc, all_paths.at(i),all_refs.at(i),on_ref, on_read , current_al, nex_al, from_current_node, from_next_node); // Use NW on the full length on_ref and on_read seq, cut it into pieces, make an al for each pieace, add to the graph
			}


		}
		if(it1 != nodes.end() && nex_name == current_name){//If they are neighbor and are on the same node on the ref graph, TODO should check if it is after or before!!

		}
	}
}
template<typename T>
void als_components<T>::add_als(size_t & comp, size_t & refacc, size_t & readacc, std::vector<int> & path, std::vector<size_t> & refs, std::string & onref , std::string & onread, const pw_alignment & cur_al , const pw_alignment & nex_al, std::string & from_current, std::string & from_next){
	size_t l2,r2;
	cur_al.get_lr2(l2,r2);
	size_t pos = r2 +1;
	std::cout << "pos "<<pos << std::endl;
	size_t cur_ref = cur_al.getreference1();
	size_t readid = cur_al.getreference2();
	size_t nex_ref = nex_al.getreference1();
	std::string read_out, ref_out;
	needleman<T> nl(data, model, onread, onref);
//XXX	simpleNW nl(model, onread, onref);
	size_t type = 1;
	nl.run_needleman(readacc,refacc,type,read_out,ref_out);
	std::vector<pw_alignment> als;
	size_t L2,R2;
	nex_al.get_lr2(L2,R2);
	size_t end_on_read_pos = L2-1;
	std::cout<< "END "<< L2-1<<std::endl;
	make_als(refacc, readacc, path, refs, ref_out, read_out, als, pos, cur_ref,nex_ref, from_current, from_next, readid, end_on_read_pos);
	add_to_graph(comp,als,cur_al,nex_al,refacc,readacc); 
}
template<typename T>
void als_components<T>::make_als(size_t & refacc, size_t & readacc, std::vector<int> & path, std::vector<size_t> & refs, std::string & onref , std::string & onread, std::vector<pw_alignment> & als, size_t & pos, size_t & cur_ref, size_t & nex_ref, std::string & from_current, std::string & from_next, size_t & readid , size_t & end_on_read_pos){
	size_t node_length = 0;
	size_t pos_on_read = pos;//Starting position on the read!
	std::cout<< "POS "<< pos_on_read <<std::endl;
	std::cout << "path size "<< path.size() <<std::endl;
	for(size_t i =0; i < path.size() ;i++){
		std::cout << "i "<< i << std::endl;
		size_t begin1 = 0;
		size_t begin2 = pos_on_read;
		node_length = data.getSequence(refs.at(i)).length();
		size_t end1 = node_length -1;
		if(i==0 &&refs.at(i)== cur_ref){
			node_length = from_current.length();
		}
		if(i == path.size()-1 && refs.at(i)== nex_ref){
			node_length = from_next.length();
		}
		std::string refout, readout;
		bool GAP = false;
		get_al_samples(onref,onread,node_length, refout, readout, pos_on_read , GAP);
		std::cout<< "onref " << onref << " on read "<< onread << " " << pos_on_read << std::endl;
		if(i == path.size()-1 && onref.length() != 0){//If there is only gap at the end
			assert(onref.length()== onread.length());
			refout.append(onref);
			readout.append(onread);
			size_t count = 0;
			for(size_t j =0; j < onread.length(); j++){
				if(onread.at(j)!='-'){
					count ++;
				}
			}
			pos_on_read = pos_on_read + count;
			if(count>0){
				GAP = false;
			}
		}
		if(i ==path.size()-1){
			std::cout << pos_on_read << " "<< end_on_read_pos<<std::endl;
			assert(pos_on_read-1 == end_on_read_pos);
		}

		if(GAP == false){
			size_t end2 = pos_on_read -1;
			if(path.at(i)>0){
				pw_alignment p(refout,readout,begin1,begin2,end1,end2,refs.at(i),readid);//Check if one sample is only gap!! I assume with this method it is not happing on the ref but only on the read!
				als.push_back(p);
			}else{
				pw_alignment p(refout,readout,end1,begin2,begin1,end2,refs.at(i),readid);
				als.push_back(p);
			}
		}else{
			size_t end2 = data.getSequence(readid).length();//Length of readid
			if(path.at(i)>0){
				pw_alignment p(refout,readout,begin1,begin2,end1,end2,refs.at(i),readid);
				als.push_back(p);
			}else{
				pw_alignment p(refout,readout,end1,begin2,begin1,end2,refs.at(i),readid);
				als.push_back(p);
			}
		}
	}

}
template<typename T>
void als_components<T>::get_al_samples(std::string & onref, std::string & onread, size_t & node_length , std::string & refout, std::string & readout, size_t & pos_on_read , bool & GAP){
	size_t counter = 0;
	size_t read_counter = 0;
	size_t pos = 0;
	std::string tempref;
	std::string tempread;
	for(size_t i = 0; i < onref.size(); i++){
		tempref += onref.at(i);
		tempread += onread.at(i);
		if(onref.at(i)!='-'){
			counter ++;
		}
		if(onread.at(i)!='-'){
			read_counter ++;
		}
		if(counter== node_length){
			pos = i;
			break;
		}
	}
	if(pos != onref.length()-1){
		std::cout << "pos "<< pos << " onref length "<< onref.length()<<std::endl;
		onref = onref.substr(pos+1);
		onread = onread.substr(pos+1);
	}else{
		onref.clear();
		onread.clear();
	}
	refout = tempref;
	readout = tempread;
	std::cout << "pos on read "<< pos_on_read << " counter "<< read_counter << std::endl;
	pos_on_read = pos_on_read + read_counter;
	if(read_counter == 0){
		GAP = true;
	}
}
template<typename T>
void als_components<T>::add_to_graph(size_t & comp, std::vector<pw_alignment> & als, const pw_alignment & cur_al , const pw_alignment & nex_al, size_t & refacc, size_t & readacc){
	for(size_t i =0; i < als.size()-1; i++){
		if(i == 0){
			add_adjacencies(comp,cur_al,als.at(0),refacc,readacc);
		}
		add_adjacencies(comp,als.at(i),als.at(i+1),refacc,readacc);
	}
	add_adjacencies(comp, als.back(), nex_al, refacc,readacc);
}
template<typename T>
std::multiset<pw_alignment,sort_pw_alignment_by_left> als_components<T>::find_successors(const pw_alignment & p, size_t & component){
	std::multiset<pw_alignment,sort_pw_alignment_by_left> successors;
	size_t l,r;
	p.get_lr2(l,r);
	std::map<unsigned int, std::vector<pw_alignment> > seen_refs;
	std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it = alignments.at(component).find(p);
	assert(it != alignments.at(component).end());
	for(std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it1 = it; it1 != alignments.at(component).end(); it1++){
		pw_alignment al = *it1;
		size_t l2,r2;
		al.get_lr2(l2,r2);
		if((l2>r) && (l2-r <=MAXGAP)){//Overlap less than MAXGAP,//XXX For the moment we dont consider the cases that there is an overlap less than MAXGAP
			std::cout<<"find one! "<<std::endl;
			al.print();
			double g1, g2;
			model.gain_function(al,g1,g2);
			std::map<unsigned int , std::vector<pw_alignment> >::iterator seen = seen_refs.find(al.getreference1());
			if(seen== seen_refs.end()){
			//	std::cout << "not seen yet"<<std::endl;
			//	successors.insert(al);
				std::vector<pw_alignment> temp;
				temp.push_back(al);
				seen_refs.insert(std::make_pair(al.getreference1(),temp));
			}else{ //Check for the overlap and pick the one with higher gain value
				std::vector<pw_alignment> seen_als = seen->second;
				for(size_t i = 0; i < seen_als.size(); i++){
					pw_alignment seen_al = seen_als.at(i);
				//	std::cout << "seen al"<<std::endl;
				//	seen_al.print();
					size_t sl2, sr2;
					seen_al.get_lr2(sl2,sr2);
					if(sr2 > l2 && sl2 < r2){
						double sg1, sg2;
						model.gain_function(seen_al,sg1, sg2);
					//	std::cout << "av gain "<< (sg1+sg2)/2 << " " << (g1+g2)/2 <<std::endl;
						if(((sg1 +sg2)/2) < ((g1+g2)/2)){
							seen->second.erase(seen->second.begin()+i);		
							seen->second.push_back(al);
						}
					}
				}
			}
		}
		if(l2 > r+MAXGAP){//since alignments are ordered by their left
			break;
		}
	}
//	std::cout << "add to successors: "<<std::endl;
	for(std::map<unsigned int, std::vector<pw_alignment> >::iterator it = seen_refs.begin() ; it != seen_refs.end() ; it++){
		for(size_t i = 0; i < it->second.size() ; i++){
		//	it->second.at(i).print();
			successors.insert(it->second.at(i));
		}
	}
//	std::cout << "successors size is "<< successors.size()<<std::endl;
	return successors;
}
template<typename T>
void als_components<T>::get_subgraph(const pw_alignment & p){//Return a sub graph that starts with reference1 of the current alignmnet and the nodes are closer than MAX_GAP to this first node.
	size_t l1,r1;
	p.get_lr1(l1,r1);
	p.print();
	size_t acc = data.accNumber(p.getreference1());
	std::string ref_accession = data.get_acc(acc);
	std::string seqname = data.get_seq_name(p.getreference1());
	std::stringstream temp;
	temp<< seqname;
	size_t name_temp;
	temp >> name_temp;
//	int name = std::stoi(seqname);
	int name = name_temp;
	if(p.getbegin1()>p.getend1()){
		name = -1*name;
	}
	std::cout << "name is " << name <<std::endl;
	clock_t t1;
	t1 = clock();
//	rgraph.deep_first_search(name, ref_accession,r1);
	t1 = clock() - t1;
  	printf ("It took me %d clicks (%f seconds).\n",t1,((float)t1)/CLOCKS_PER_SEC);

	clock_t t;
	t = clock();
	rgraph.bfs(name,ref_accession,r1);
	t = clock() - t;
	printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
//	clock_t t1;
//	t1 = clock();
//	rgraph.make_sub_graph(name,ref_accession,r1);
//	t1 = clock() - t1;
//  	printf ("It took me %d clicks (%f seconds).\n",t1,((float)t1)/CLOCKS_PER_SEC);
//	assert(t > t1);
/*	std::set<std::vector<int> > paths = rgraph.get_paths();
	std::cout<< "path size is "<< paths.size()<<std::endl;
	for(std::set<std::vector<int> >::iterator it = paths.begin(); it != paths.end(); it++){
		std::vector<int> this_p = *it;
		for(size_t j =0; j < this_p.size();j++){
			std::cout << this_p.at(j)<< " ";
		}
		std::cout<< " "<<std::endl;
	}*/
}
	/*XXX For the moment we keep overlapped cases out, they can be added later. In that case:
	on read:	________________|____|
					|______|______________
			 	     /	     /	
	node1fromref:   -------------|------|
					\     \
	node2fromref:			 |-----|-----------------	
	(if the length of these two pieces of ref nodes is less than MAXGAP)
	And using NeedlemanWunsch me create the centeral alignment. In fact we get three shorter als*/
template<typename T>
void als_components<T>::look_for_successors_on_paths(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors){//XXX Directions should be checked on ref 
	std::cout<<"look for neighbors on paths "<<std::endl;
	size_t from,to;
	size_t l1,r1,l2,r2;
	current_al.get_lr1(l1,r1);
	current_al.get_lr2(l2,r2);
	std::cout << l1 << " "<< l2 << " "<< r1 << " "<< r2 <<std::endl;
	unsigned int ref1 = current_al.getreference1();
	unsigned int ref2 = current_al.getreference2();
	size_t seqlength = data.get_seq_size(ref1);//We use this seq from its r1+1 to length-1. It is part of the ref node that comes after current alignment.
	size_t readacc = data.accNumber(ref2);
	size_t refacc = data.accNumber(ref1);
	std::string current_ref_node_name = data.get_seq_name(current_al.getreference1());
	int current_node_name = std::stoi(current_ref_node_name);
	std::cout << "current node name "<< current_node_name << std::endl;
	std::string from_current_node;
	from = r1+1;
	to = seqlength-1;
	from_current_node = data.extract_seq_part(ref1,from,to);
	if(current_al.getbegin1() > current_al.getend1()){//If reverse
		std::string temp;
		get_reverse_complement(from_current_node,temp);
		from_current_node = temp;
		current_node_name = -1*current_node_name;
	}
	size_t current_remainder = to-from +1;
	std::cout << "current remainder is "<< to-from+1 << std::endl;
	//should check if the first reference of the neighbors is one of the ints in the paths
//	for(std::multiset<pw_alignment, sort_pw_alignment_by_left>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
//	std::set<int> seen;
	for(std::multiset<pw_alignment, sort_pw_alignment_by_left>::iterator it = neighbors.begin(); it != neighbors.end(); it++){ //XXX changed it from reverse
	//	std::set<std::pair<int,int> > seen1;
		std::map<std::pair<int,int>, std::set<std::vector<int> > > seen1;

		alignments_on_each_ref_node.clear();
		size_t left1, right1 ,left2, right2;
		const pw_alignment pal = *it;
		pal.get_lr1(left1,right1);

		pal.get_lr2(left2,right2);
		std::cout << left1 << " "<< left2 << " "<< right1 << " "<< right2 <<std::endl;

		std::string seqname = data.get_seq_name(pal.getreference1());
		int name = std::stoi(seqname);
		if(pal.getbegin1()>pal.getend1()){
			name = -1*name;
		}
		std::cout << "neighbor's name is " << name <<std::endl;
		pal.print();
		if(pal.getreference1() == ref1){//If they share the same node on the reference //TODO also consider repeats as part of this condition!
			size_t onread_from = r2+1;
			size_t onread_to = left2-1;
			if(r1<left1 && left1-r1<=MAXGAP && name == current_node_name){//current happened on the node before its neighbor
			//if both are forwards or both are backwards
				from = r1+1;
				to = left1-1;
				if(name > 0){
					std::cout<<"when they share a node forward"<<std::endl;
					current_al.print();
					pal.print();
					if((r2 == left2-1) && (r1 != left1-1)){
						std::cout<<"only gap on read "<<std::endl;
						//only gap on read 
						onread_from = r2;
						onread_to = data.get_seq_size(ref2);
						make_al_on_a_node(comp,current_al,pal,true,ref1,from,to,ref2,onread_from,onread_to,refacc,readacc);
					}
					else if((r2 != left2-1) && (r1 == left1-1)){
						unsigned int temp = data.numSequences()+1;//only gap on ref
						from = 0;
						to = 0;
						make_al_on_a_node(comp,current_al,pal,true,temp,from, to,ref2,onread_from,onread_to,refacc,readacc);
					}
					else if((r2 != left2-1) && (r1 != left1-1)){
						make_al_on_a_node(comp,current_al,pal,true,ref1,from,to,ref2,onread_from,onread_to,refacc,readacc);
					}else{
						assert((r2 == left2-1) && (r1 == left1-1));
						add_adjacencies(comp,current_al,pal,refacc,readacc);

					}

				}else{//When both are backwards //TODO This is wrong!
					std::cout << "when they share a node backwards" <<std::endl;
					assert(name < 0);
					if((r2 = left2-1) && (r1 != left1-1)){
						//only gap on read
						onread_from = r2;
						onread_to = data.get_seq_size(ref2);
						make_al_on_a_node(comp,current_al,pal,false,ref1,from,to,ref2,onread_from,onread_to,refacc,readacc);
					}
					else if((r2 != left2-1) && (r1 == left1-1)){
						unsigned int temp = data.numSequences()+1;//only gap on ref
						from = 0;
						to = 0;
						make_al_on_a_node(comp, current_al,pal,false,temp,from,to,ref2,onread_from,onread_to,refacc,readacc);
					}
					else if((r2 != left2-1) && (r1 != left1-1)){
						make_al_on_a_node(comp,current_al,pal,false,ref1,from,to,ref2,onread_from,onread_to,refacc,readacc);
					}else{
						assert((r2 == left2-1) && (r1 == left1-1));						
						add_adjacencies(comp,current_al,pal,refacc,readacc);

					}
				}

			}else if(l1> right1 && l1-right1<=MAXGAP && name == current_node_name && name < 0){//TODO Do we even need that case???  what if one forwards and the other backwards? I think we might need to add this case instead
			//example: //XXX Nope! doesn't need to be added. It is wrong!
			//current: ref1 6 4
			//	   ref2 3 5
			//neighbor ref1 7 8
			//         ref2 7 8
			std::cout<< "neighbor before current "<<std::endl;
			//only if both ref1s are reverse

			}else{
				continue;
			}


		}else{//When they DONT share the same node on the ref graph
			std::set<int>::iterator it1 = nodes.find(name);
		//	std::set<int>::iterator s = seen.find(name);//XXX doesnt make sense if going forward!!
			if(it1 != nodes.end()){//if the neighbor is  aligned on a close node. NeedlemanWunsch is used to create an alignments for the gap part 
		//	if(it1 != nodes.end() && s == seen.end()){
				std::cout<<"here is one! "<<std::endl;
				unsigned int reference1 = pal.getreference1();
				std::string from_next_node;
				from = 0;
				if(left1 > 0){
					to = left1-1;
					from_next_node = data.extract_seq_part(reference1,from,to);//It could be that there are more than one node distance in between!!! they are all added.
					if(pal.getbegin1() > pal.getend1()){//If reverse
						std::cout<< "ref is reverse"<<std::endl;
						std::string temp;
						get_reverse_complement(from_next_node,temp);
						from_next_node = temp;
					}
				}else{//starts from the beginning of the last node on the path
					to = 0;
					from_next_node = "";
					std::cout<<"starts from the beginning of the last node on the path"<<std::endl;
				}
				if(from_next_node.length()>MAXGAP) continue;
				size_t last_remainder = left1;
				std::vector<std::vector<size_t> >all_refs;
				std::vector<std::vector<int> >all_paths;
				std::vector<std::vector<std::string> >all_strings_from_ref_graph; //each string is content of a node on a path
				append_nodes(nodes,name,current_node_name,reference1, ref1, refacc,all_refs,all_paths,all_strings_from_ref_graph);//add all the nodes to a string 
				assert(all_strings_from_ref_graph.size()==all_paths.size());
				assert(all_strings_from_ref_graph.size() == all_refs.size());
				assert(all_strings_from_ref_graph.size() > 0);
				std::cout << "size of strings "<< all_strings_from_ref_graph.size() << std::endl;
				for(size_t s = 0; s < all_strings_from_ref_graph.size(); s++){
					std::cout << all_strings_from_ref_graph.at(s).size()<<std::endl;
					std::cout<< "this path is "<<std::endl;
					for(size_t p = 0 ; p < all_paths.at(s).size(); p++){
						std::cout << all_paths.at(s).at(p) << " ";
						
					}
					std::cout << " "<<std::endl;
					bool SEEN = false;
					size_t counter = 0;
					from = r2+1;
					to = left2-1;//We dont need to check for left2==0 now since we already picked those which are bigger than r2 but if one day we look for overlapped one we should consider that too.
					std::string str2 = data.extract_seq_part(ref2,from,to);
					bool CONTINUE = false;
					std::vector<int> from_path;
					for(size_t p = 0 ; p < all_paths.at(s).size()-1; p++){//Exception is when size == 2
						from_path.push_back(all_paths.at(s).at(p));
						std::cout << "Count "<< p <<std::endl;
						std::pair<int,int> temp = std::make_pair(all_paths.at(s).at(p), all_paths.at(s).at(p+1));
						std::map<std::pair<int,int> , std::set<std::vector<int> > >::iterator this_edge = seen1.find(temp);
						bool READ_PATH; 
						if(this_edge != seen1.end()){
							//Find the path up to here! and if it doesnt exist we should add it and calculate the path again!
							std::set<std::vector<int> >::iterator pre_path = this_edge->second.find(from_path);
							if(pre_path == this_edge->second.end()){
								READ_PATH = true;
							}else{
								READ_PATH = false;
							}
						}
						if(this_edge == seen1.end()|| READ_PATH == true || CONTINUE == true){
							if(this_edge == seen1.end()){
								CONTINUE = true;
								std::set<std::vector<int> > this_set;
								this_set.insert(from_path);
								seen1.insert(std::make_pair(temp, this_set));
							}else{
								this_edge->second.insert(from_path);
							}
							counter ++;
						//	seen1.insert(std::make_pair(temp, from_path));
							std::vector<int> this_path;
							this_path.push_back(all_paths.at(s).at(p));
							this_path.push_back(all_paths.at(s).at(p+1));
							rgraph.delete_path(this_path); 
							std::string this_string;
							size_t onread_from = r2+1;
							size_t onref_from = r1+1;
							to = left2-1;
							pw_alignment nw_al;
							if(p == 0){
								this_string.append(from_current_node);
								std::string str2 = data.extract_seq_part(ref2,onread_from,to);
								if(this_string.length() != 0){//Gap only on ref is not allowed
									std::cout << "non empty!"<<std::endl;
									size_t type = 3;
									add_nodes_to_the_graph(type, str2, this_string,refacc,readacc,nw_al,ref2,ref1,onread_from,onref_from);
									nw_al.print();
									add_to_containers(comp,current_al,nw_al,refacc, readacc, current_node_name, current_node_name); //Add to adjacencies,indices(2), alignment_on_each_node
								}else{
									//Add the current al to the alignments_on_each_ref_node
									std::cout << "empty " << current_node_name <<std::endl;
									std::map<int , std::set<size_t> >::iterator it=alignments_on_each_ref_node.find(current_node_name);
									if(it == alignments_on_each_ref_node.end()){
										alignments_on_each_ref_node.insert(std::make_pair(current_node_name,std::set<size_t>()));
										it = alignments_on_each_ref_node.find(current_node_name);
									}
									std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator it1 = node_indices.at(comp).find(current_al);
									if(it1 == node_indices.at(comp).end()){
										node_indices.at(comp).insert(std::make_pair(current_al,index));
										indices_nodes.at(comp).insert(std::make_pair(index,current_al));
										index++;
										it1 = node_indices.at(comp).find(current_al);
									}
									it->second.insert(it1->second);
									std::cout<< "current al " << it1->second <<std::endl;
									current_al.print();

								}
							}
							if(p == all_paths.at(s).size()-2){
								this_string.clear();
								this_string.append(from_next_node);
								assert(all_refs.at(s).at(p+1) == reference1);
								unsigned int this_ref = all_refs.at(s).at(p+1);
								std::vector<pw_alignment> pre_als;
								std::vector<size_t> onread_froms;
								std::cout << "all_paths.at(s).at(p) " << all_paths.at(s).at(p) <<std::endl;
								get_previous_als(comp, all_paths.at(s).at(p),pre_als,onread_froms,to,refacc, readacc);
								for(size_t i = 0 ; i < pre_als.size() ; i++){
									std::cout<< "pre_al: "<<std::endl;
									std::cout<< "on read from "<< onread_froms.at(i)<<std::endl;
									str2 = data.extract_seq_part(ref2,onread_froms.at(i),to);
									onref_from = 0;
									if(str2.length() != 0 || this_string.length() != 0){
										std::cout << "str2 length "<< str2.length() <<std::endl;
										size_t type = 1; //TODO check!
										add_nodes_to_the_graph(type, str2, this_string,refacc,readacc,nw_al,ref2,this_ref,onread_froms.at(i),onref_from);
										add_to_containers(comp,pre_als.at(i),nw_al,refacc, readacc, all_paths.at(s).at(p),all_paths.at(s).at(p+1)); 
										size_t L2,R2;
										nw_al.get_lr2(L2,R2);
										nw_al.print();
										if(R2 != data.getSequence(ref2).length()){
											std::cout << "length " << data.getSequence(ref2).length()<< std::endl;
											assert(R2 == left2-1);
										}else{
											std::cout<< "GAP ONLY ON READ! L2: "<< L2 <<std::endl;
											assert(L2 == left2);
										}
										//Add new_al to pal
										add_to_containers(comp,nw_al,pal,refacc, readacc, all_paths.at(s).at(p+1),all_paths.at(s).at(p+1)); 
									}else{
										assert(str2.length() == 0 && this_string.length() == 0);
										add_to_containers(comp,pre_als.at(i),pal,refacc, readacc, all_paths.at(s).at(p),all_paths.at(s).at(p+1)); 
									}
								}		
								assert(all_paths.at(s).at(p+1) == name);
							}else{
						//	if(p != all_paths.at(s).size()-2 && p != 0 )
						//		assert(p != all_paths.at(s).size()-2 && p != 0 );
								this_string.clear();
								this_string.append(all_strings_from_ref_graph.at(s).at(p));
								assert(this_string.length() != 0);
								unsigned int this_ref = all_refs.at(s).at(p+1);
								std::vector<pw_alignment> pre_als;
								std::vector<size_t> onread_froms;

								get_previous_als(comp, all_paths.at(s).at(p),pre_als,onread_froms,to, refacc, readacc);
								std::cout << "pre als size "<< pre_als.size()<<std::endl;
								for(size_t i =0; i < pre_als.size();i++){
									std::cout << "this pre al" <<std::endl;
									pre_als.at(i).print();
									str2 = data.extract_seq_part(ref2,onread_froms.at(i),to);
									assert(str2.length() != 0 || this_string.length() != 0);
									onref_from = 0;
									size_t type = 3;
									add_nodes_to_the_graph(type , str2, this_string,refacc, readacc, nw_al, ref2, this_ref,onread_froms.at(i),onref_from);
									add_to_containers(comp,pre_als.at(i),nw_al,refacc, readacc, all_paths.at(s).at(p),all_paths.at(s).at(p+1));
									std::cout << all_paths.at(s).at(p) << " " << all_paths.at(s).at(p+1) <<std::endl;
								}

							}
						}
						/*else{
							SEEN = true;
							break;
						}*/
					}
				/*	if(SEEN == false){
					assert(counter == all_paths.at(s).size()-1);
					std::string str1 = from_current_node; //Remainder form the current node
					str1.append(all_strings_from_ref_graph.at(s));
					std::cout<< "str1 length "<< str1.length()<<std::endl;
				//	assert(str1.length()<= MAXGAP);
					std::cout << "from last node length "<< from_next_node.length()<<std::endl;
					str1.append(from_next_node);
					std::cout<< "str1 length at the end "<< str1.length()<<std::endl;
					from = r2+1;
					to = left2-1;//We dont need to check for left2==0 now since we already picked those which are bigger than r2 but if one day we look for overlapped one we should consider that too.
					std::string str2 = data.extract_seq_part(ref2,from,to);

					size_t type = 1; //Both ends are fixed.
					std::string read_out;
					std::string ref_out;

					if(str1.length() != 0 && str2.length() !=0 && str1.length()<=MAXGAP){
						std::cout<< "str2 " << str2.size() << " str1 "<< str1.size() <<std::endl;
						assert(str2.length()<=MAXGAP);
						std::cout << str1 << std::endl;
						std::cout << str2 << std::endl;
						needleman<dynamic_mc_model> nl(data, model, str2, str1);
						nl.run_needleman(readacc,refacc,type,read_out,ref_out);
						std::cout << "read out size "<<read_out.size() << " ref out size "<<ref_out.size()<<std::endl;
						std::cout<< "current " << current_remainder << " last "<< last_remainder<<std::endl;
						make_alignments(comp,from,to,current_al,p,read_out,ref_out, all_refs.at(s),all_paths.at(s),ref2,current_remainder,last_remainder,refacc,readacc);//makes als and save them in a container
					}else if(str1.length()==0 && str2.length()!=0){
						assert(str2.length()<=MAXGAP);
						std::string onlygap(str2.length(),'-');
						std::cout<<"gap on ref graph "<<std::endl;
						pw_alignment onlygap_al(onlygap,str2,0,from,0,to,data.numSequences()+1,ref2); //basically there is no ref for the gap, i set it to numseq+1.
						onlygap_al.print();
						std::cout <<"add adjs: "<<std::endl;
						add_adjacencies(comp, current_al,onlygap_al,refacc,readacc);
						add_adjacencies(comp,onlygap_al,p,refacc,readacc);
					}else if(str1.length()!=0&&str2.length()==0 && str1.length()<=MAXGAP){
						std::string onlygap(str1.length(),'-');
						std::cout<< str1.length()<<std::endl;
						from = r2;
						if(r2 == data.getSequence(ref2).length()){//It is unnecessary here i guess, since we always started with an ungapped alignment
							from = l2;
						}
						to = data.getSequence(ref2).length();
						std::cout<< "gap on read "<<std::endl;
						unsigned int virtual_ref = ref2;
						make_alignments(comp, from,to,current_al,p,onlygap,str1, all_refs.at(s),all_paths.at(s), virtual_ref,current_remainder,last_remainder,refacc,readacc);//XXX not so sure about this line, might need to write another make_al function to handle this situtation
					}else if(str1.length()==0 && str2.length()==0){
						std::cout<< "should be here!"<<std::endl;
						add_adjacencies(comp, current_al,p,refacc,readacc);

					}
					}*/
				}
			}
		}

	}

}
template<typename T>
void als_components<T>::get_previous_als(size_t & comp , int & index, std::vector<pw_alignment> & pre_als, std::vector<size_t> & onread_from,size_t & to, size_t & ref_acc, size_t & read_acc){//XXX Originally I was only checking the last alignment, but it was wrong. I need to check all of them! Becasue getting from different nodes to node index may create different als on it and all need to be taked into account.
	std::cout << "index "<< index <<std::endl;
	std::map<int , std::set<size_t> >::iterator it = alignments_on_each_ref_node.find(index);
	assert(it!=alignments_on_each_ref_node.end());
	std::set<size_t> al_indices = it->second;
	std::cout << "indices size "<< al_indices.size()<<std::endl;
	for(std::set<size_t>::iterator it2 = al_indices.begin() ; it2 != al_indices.end(); it2++){
		size_t this_index = *it2;
		std::cout << "this index "<< this_index <<std::endl;
		std::multimap<size_t, const pw_alignment>::iterator it1 =indices_nodes.at(comp).find(this_index);
		assert(it1 != indices_nodes.at(comp).end());
		pw_alignment pre_al;
		pre_al = it1->second;
		pre_al.print();
		size_t l1,r1,l2,r2;
		pre_al.get_lr1(l1,r1);
		pre_al.get_lr2(l2,r2);
		size_t position = r2 + 1;
		if(r2 == data.getSequence(pre_al.getreference2()).length()){
			std::cout<<"gap on read!"<<std::endl;
			position = l2;
		}
		if(pre_al.getreference1()< data.numSequences()){
			std::cout<< "ref1 "<< pre_al.getreference1()<< " ref2 "<<pre_al.getreference2() << "ref1 length "<< data.getSequence(pre_al.getreference1()).length() <<std::endl;
		}else{
			std::cout<< "ref1 "<< pre_al.getreference1()<< " ref2 "<<pre_al.getreference2() << std::endl;
		}
		if(pre_al.getreference1()< data.numSequences() && r1 == data.getSequence(pre_al.getreference1()).length()-1){
			pre_als.push_back(pre_al);
			onread_from.push_back(position);
		}
		if(pre_al.getreference1()< data.numSequences() && r1 < data.getSequence(pre_al.getreference1()).length()-1){
			//Make a new alignment till the end and add it to pre_als as well as the other containers
			std::cout << "on read to "<< to << std::endl;
			pw_alignment al;
			size_t type = 3;
			unsigned int read_id = pre_al.getreference2();
			unsigned int ref_id = pre_al.getreference1();
			std::string from_read = data.extract_seq_part(read_id, position,to);
			std::string from_ref;//Direction matters!
			size_t ref_from = r1+1;
			size_t ref_to = data.getSequence(pre_al.getreference1()).length()-1;
			if(index > 0){
				from_ref = data.extract_seq_part(ref_id,ref_from,ref_to);

			}else{
				from_ref = data.extract_reverse_seq_part(ref_id,ref_from,ref_to);
			}
			if(from_read.length()<=MAXGAP){
				add_nodes_to_the_graph(type,from_read,from_ref,ref_acc,read_acc,al,read_id,ref_id,position,ref_from);
				add_to_containers(comp, pre_al, al, read_acc, ref_acc, index, index);
				pre_als.push_back(al);
				al.get_lr2(l2,r2);
				position = r2 + 1;
				if(r2 == data.getSequence(al.getreference2()).length()){
					std::cout<<"gap on read!"<<std::endl;
					position = l2;
				}
				onread_from.push_back(position);
			}
		}
		if(pre_al.getreference1() == data.numSequences()+1){//TODO is it even happening? Seems I am not saving gap only on ref which means in that case I am losing some parts of the read or making another alignment on the same part more than once?!! Should think of this situtation !!
			std::cout<< "TODO!!"<<std::endl;
		}
	}
}
template<typename T>
void als_components<T>::add_to_containers(size_t & comp , const pw_alignment & p1 , const pw_alignment & p2 , size_t & ref_acc, size_t & read_acc, int & current_name, int & next_name){
	add_adjacencies(comp, p1 , p2, ref_acc, read_acc);
	std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator it1 = node_indices.at(comp).find(p1);
	assert(it1 != node_indices.at(comp).end());
	std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator it2 = node_indices.at(comp).find(p2);
	assert(it2 != node_indices.at(comp).end());
	std::map<int , std::set<size_t> >::iterator it=alignments_on_each_ref_node.find(current_name);
	if(it != alignments_on_each_ref_node.end()){
		std::set<size_t>::iterator it3 = it->second.find(it1->second);
		if(it3 == it->second.end()){
			it->second.insert(it1->second);
		}
	}else{
		alignments_on_each_ref_node.insert(std::make_pair(current_name,std::set<size_t>()));
		it = alignments_on_each_ref_node.find(current_name);
		it->second.insert(it1->second);
	}
//	if(it ==alignments_on_each_ref_node.end()){
//		alignments_on_each_ref_node.insert(std::make_pair(current_name,std::vector<size_t>()));
//		it = alignments_on_each_ref_node.find(current_name);
//	}
//	it->second.push_back(it1->second);
	it = alignments_on_each_ref_node.find(next_name);
	if(it ==alignments_on_each_ref_node.end()){
		alignments_on_each_ref_node.insert(std::make_pair(next_name,std::set<size_t>()));
		it = alignments_on_each_ref_node.find(next_name);
	}
	it->second.insert(it2->second);


}
template<typename T>
void als_components<T>::add_nodes_to_the_graph(size_t & type , std::string & seq_from_read, std::string & seq_from_ref, size_t & refacc , size_t & readacc, pw_alignment & nw_alignment, unsigned int & read_id , unsigned int & ref_id, size_t & onread_from, size_t & onref_from){//Always one node except on the first node
//	size_t type = 3; //Starts from begining of the read part, goes till end of the node from the ref graph. 
	std::string read_out;
	std::string ref_out;
	std::cout << "on read from " << onread_from << " "<< onref_from << std::endl;
	if(seq_from_read.length() != 0 && seq_from_ref.length() !=0){
		std::cout<< "str2 " << seq_from_read << " str1 "<< seq_from_ref.size() <<std::endl;
		assert(seq_from_read.length()<=MAXGAP);
		assert(seq_from_ref.length()<=MAXGAP);
		std::cout << seq_from_ref << std::endl;
		std::cout << seq_from_read << std::endl;
		needleman<T> nl(data, model, seq_from_read, seq_from_ref);
//XXX		simpleNW nl(model, seq_from_read, seq_from_ref);
		nl.run_needleman(readacc,refacc,type,read_out,ref_out);
		assert(read_out.length() == ref_out.length());
		size_t counter = 0;
		for(size_t i =0; i < read_out.size();i++){
			if(read_out.at(i)!= '-'){
				counter ++;
			}
		}
		size_t onread_to;
		if(counter != 0){
			onread_to = onread_from + counter-1;
		}else{
			onread_to = data.getSequence(read_id).length();
		}
		std::cout<< onread_from << " to: " << onread_to << " " << counter <<std::endl;
		counter = 0;
		for(size_t i =0; i < ref_out.size();i++){
			if(ref_out.at(i)!= '-'){
				counter ++;
			}
		}
		assert(counter != 0);
		size_t onref_to = onref_from + counter-1;
		pw_alignment p(ref_out,read_out,onref_from,onread_from,onref_to,onread_to, ref_id, read_id);
		nw_alignment = p;
	}else if(seq_from_ref.length()==0 && seq_from_read.length()!=0){//TODO can it only happen on the first or the last node on the path? how should i handle it?
		assert(seq_from_read.length()<=MAXGAP);
		std::string onlygap(seq_from_read.length(),'-');
		std::cout<<"gap on ref graph "<<std::endl;
		size_t onread_to = onread_from + seq_from_read.length()-1;
		pw_alignment onlygap_al(onlygap,seq_from_read,0,onread_from,0,onread_to,data.numSequences()+1,read_id); //basically there is no ref for the gap, i set it to numseq+1.
		onlygap_al.print();
		nw_alignment = onlygap_al;
	}
	else if(seq_from_ref.length()!=0&&seq_from_read.length()==0){
		assert(seq_from_ref.length()<=MAXGAP);
		std::string onlygap(seq_from_ref.length(),'-');
		std::cout<< seq_from_ref.length()<<std::endl;
		std::cout << "GAP ONLY ON READ"<<std::endl;
		size_t onref_to = onref_from + seq_from_ref.length()-1;
		size_t onread_to = data.getSequence(read_id).length();
		pw_alignment onlygap_al(seq_from_ref, onlygap, onref_from,onread_from,onref_to, onread_to,ref_id,read_id); //basically there is no ref for the gap, i set it to numseq+1.
		onlygap_al.print();
		nw_alignment = onlygap_al;

	}
	assert(seq_from_ref.length()!=0 || seq_from_read.length()!=0);
}
template<typename T>
void als_components<T>::make_al_on_a_node(size_t & comp, const pw_alignment& first_al , const pw_alignment& second_al, bool dir, unsigned int & ref1, size_t & left1, size_t & right1, unsigned int & ref2, size_t & left2, size_t & right2,size_t & refacc, size_t & readacc){
	std::string str1;
	size_t begin1=0;
	size_t end1=0;
	if(ref1 != data.numSequences()+1){
		str1 = data.extract_seq_part(ref1, left1,right1);
		begin1 = left1;
		end1 = right1;
		if(dir==false){//Make the reverse complement
			std::string temp;
			get_reverse_complement(str1, temp);
			begin1 = right1;
			end1 = left1;
			str1 = temp;
		}
	}
	std::string str2;
	assert(ref2 != data.numSequences()+2);//Just to get sure i fixed every parts in the code and removed seq+2
	if(right2 != data.getSequence(ref2).length()){
		str2 = data.extract_seq_part(ref2,left2,right2);
	}
	size_t type = 1; //Both ends are fixed.
	std::string read_out;
	std::string ref_out;
	if(str1.length() != 0 && str2.length() !=0){
		size_t readacc = data.accNumber(ref2);
		size_t refacc = data.accNumber(ref1);
		needleman<T> nl(data, model, str2, str1);
//XXX		simpleNW nl( model, str2, str1);
		nl.run_needleman(readacc,refacc,type,read_out,ref_out);
		pw_alignment p(ref_out,read_out,begin1,left2,end1,right2,ref1,ref2);
		std::cout << "1: "<<std::endl;
		p.print();
		add_adjacencies(comp,first_al,p,refacc,readacc);
		add_adjacencies(comp,p,second_al,refacc,readacc);
		first_al.print();
	}
	if(str1.length() != 0 && str2.length() ==0){//gap on read
		std::string temp(str1.length(),'-');
		str2 = temp;
		pw_alignment p(str1,str2,begin1,left2, end1,right2,ref1,ref2);
		std::cout << "2: "<<std::endl;
		p.print();
		add_adjacencies(comp, first_al,p,refacc,readacc);

		add_adjacencies(comp,p,second_al,refacc,readacc);

	}
	if(str1.length() == 0 && str2.length() !=0){//gap on ref
		std::string temp(str2.length(),'-');
		str1 = temp;
		pw_alignment p(str1,str2,begin1,left2, end1,right2,ref1,ref2);
		std::cout << "3: "<<std::endl;
		p.print();
		add_adjacencies(comp,first_al,p,refacc,readacc);
		add_adjacencies(comp,p,second_al,refacc,readacc);

	}
	assert(str1.length() != 0 || str2.length() != 0);

}
template<typename T>
void als_components<T>::get_paths(const std::set<int> & bfs_nodes, int & node_name , int & current_node_name, std::set<std::vector<int> > & paths){ //TODO do it with adjacencies also when they are read from the te stack add them to a seen container and never add the seen ones to the stack any more to avoid the loop!!
	std::set<int> seen;
	std::map<int, std::vector<std::vector<int> > > all_paths;
	std::vector<int> stack;
	stack.push_back(current_node_name);
	all_paths.insert(std::make_pair(current_node_name, std::vector<std::vector<int> >()));
	while(stack.size() != 0){
		int this_node;
		this_node = stack.at(0);
		stack.erase(stack.begin());
		seen.insert(this_node);
		std::set<int> adjs = rgraph.get_adjacencies(this_node);
		assert(adjs.size() != 0);
		std::cout << "thisnode "<< this_node<< std::endl;
		for(std::set<int>::iterator it = adjs.begin(); it != adjs.end() ; it++){
				std::map<int, std::vector<std::vector<int> > >::iterator node = all_paths.find(this_node);
				assert(node!=all_paths.end());
				std::set<int>::const_iterator bfsn = bfs_nodes.find(*it);
				std::set<int>::iterator see = seen.find(*it);
				if(*it != node_name && bfsn != bfs_nodes.end() && see == seen.end()){
					stack.push_back(*it);
				}
				if(bfsn != bfs_nodes.end()){
					for(size_t i =0 ; i < node->second.size(); i++){
						std::vector<int> p  = node->second.at(i);
						std::cout << "p size "<<p.size()<<std::endl;
						p.push_back(this_node);
						std::map<int, std::vector<std::vector<int> > >::iterator it1 = all_paths.find(*it);
						if(it1 == all_paths.end()){
							all_paths.insert(std::make_pair(*it, std::vector<std::vector<int> >()));
							it1 = all_paths.find(*it);
						}
						it1->second.push_back(p);
					}
					if(node->second.size() == 0){
						std::vector<int> p;
						p.push_back(this_node);
						std::map<int, std::vector<std::vector<int> > >::iterator it1 = all_paths.find(*it);
						assert(it1 == all_paths.end());
						all_paths.insert(std::make_pair(*it, std::vector<std::vector<int> >()));
						it1 = all_paths.find(*it);//TODO all_paths contains the same path several times and has to get fixed!!!
						it1->second.push_back(p);	
					}
				}
		}
	}
	std::map<int , std::vector<std::vector<int> > >::iterator it = all_paths.find(node_name);
	assert(it != all_paths.end());
	for(size_t i = 0; i < it->second.size();i++){
		paths.insert(it->second.at(i));
	}
}
template<typename T>
void als_components<T>::append_nodes(const std::set<int> & bfs_nodes ,int & name, int & current_node_name, unsigned int & next_ref_id, unsigned int & current_ref_id, size_t & refacc , std::vector<std::vector<size_t> >& all_refs, std::vector<std::vector<int> >& all_paths, std::vector<std::vector<std::string> >& all_strings){
	std::string ref_accession = data.get_acc(refacc);
	std::map< std::string, size_t> longname2seqidx;
	longname2seqidx = data.getLongname2seqidx();
	std::set<std::vector<int> > all_p = rgraph.get_paths();
	std::cout << "name !"<<name <<std::endl;
	std::set<std::vector<int> > PATH;
	std::cout << "all_p size: "<< all_p.size()<<std::endl; //TODO the same path is repeated several times! has to get fixed!
	for(std::set<std::vector<int> >::iterator it = all_p.begin(); it!= all_p.end() ; it++){
		std::vector<int> temp = *it;
		std::vector<int> this_p;
		for(size_t j = 0; j < temp.size(); j++){
			this_p.push_back(temp.at(j));
			if(temp.at(j)==name){
				PATH.insert(this_p);
				break;
			}
		}
	//	std::cout << this_p <<std::endl;
	//	PATH.insert(this_p);
	}
	std::cout << "PATH size "<< PATH.size()<<std::endl; //TODO check the paths and if they already exist dont add to all_path and all_string
	for(std::set<std::vector<int> >::iterator it =PATH.begin() ; it != PATH.end() ; it++){ //could be more than one path then more than one string !!
		std::vector<int> this_path = *it;
	//	this_path.push_back(name);
		std::vector<size_t> this_refs;
		this_refs.push_back(current_ref_id);
		std::string str1;
		std::vector<std::string> this_string;
		std::cout<< "this path size "<<this_path.size()<<std::endl;
		std::cout << this_path <<std::endl;
		if(this_path.size()>2){
			for(size_t j = 1; j < this_path.size()-1; j++){
				std::stringstream ss;
				int index = this_path.at(j);
				if(index<0){
					index = -1*index;
				}
				ss << index;
				std::string this_name = ss.str();
				std::string longname = ref_accession + ":"+this_name;
				std::cout<<longname<<std::endl;
				std::map<std::string, size_t>::iterator findseq = longname2seqidx.find(longname);
				assert(findseq != longname2seqidx.end());
				unsigned int ref = findseq->second;
				this_refs.push_back(ref);
				size_t from = 0;
				size_t to = data.getSequence(ref).length()-1;
				if(this_path.at(j)<0){
					std::string temp_in = data.extract_seq_part(ref,from,to);
					std::string temp_out;
					get_reverse_complement(temp_in,temp_out);
					str1.append(temp_out);
					this_string.push_back(temp_out);
					std::cout << "temp out size "<< temp_out<<std::endl;
				}else{
					std::string temp = data.extract_seq_part(ref,from, to);
					std::cout << "temp size "<< temp.size()<<std::endl;
					str1.append(data.extract_seq_part(ref, from, to));
					this_string.push_back(temp);
				}
			//	assert(str1.length() < MAXGAP);
			}
		}
		this_refs.push_back(next_ref_id);
		if(str1.length() >= MAXGAP){
			std::cout << "PATH LENGTH IS LONGER THAN MAXGAP! "<<std::endl;
		}
		if(str1.length()<MAXGAP){
			all_strings.push_back(this_string);
			all_refs.push_back(this_refs);
			all_paths.push_back(this_path);
		//	rgraph.delete_path(this_path); 
		}
	}
}
template<typename T>
void als_components<T>::make_alignments(size_t & comp, size_t & begin_on_read, size_t & end_on_read, const pw_alignment & first_al, const pw_alignment & last_al, std::string & read_out, std::string & ref_out, std::vector<size_t> &this_refs, std::vector<int> & this_path, unsigned int & ref2, size_t & first_length, size_t & last_length,size_t &refacc,size_t&readacc){//TODO Further improvement: maybe it makes sense to make another pw_alignment class with no ref number!
	assert(read_out.length() == ref_out.length());
	assert(this_refs.size()>=2);
	size_t counter = 0;
	size_t read_counter = 0;
	size_t ref_counter = 0;
	size_t current_pos = 0;
	std::string str1;
	std::string str2;
	size_t b1,e1,b2,e2;
	b2 = begin_on_read;
	e2 = begin_on_read;
	pw_alignment previous;
	//Making the first al
	std::cout << "first length " << first_length << std::endl;
//	if(this_refs.size()==2 && first_length>0 && last_length==0){
//		std::cout<< "all of it is located on the first ref"<<std::endl;
//	}
	if(first_length > 0){
		bool onlyfirst = false;
		if(this_refs.size()==2 && last_length==0){
			std::cout<< "all of it is located on the first ref"<<std::endl;
			onlyfirst = true;
		}
		compute_samples(onlyfirst,first_length,current_pos, read_counter, ref_counter , read_out, str2, ref_out, str1);
	//	for(size_t i = 0; i < ref_out.size();i++){
	//		if(counter < first_length){
	//			current_pos = i;
	//			str1+=ref_out.at(i);
	//			str2+=read_out.at(i);
	//			if(ref_out.at(i)!= '-'){
	//			//	std::cout<< counter << " "<<first_length<<std::endl;
	//				counter++;
	//			}
	//			if(read_out.at(i)!='-'){
	//			//	std::cout << read_counter <<std::endl;
	//				read_counter++;
	//			}
	//		}else{
	//			break;
	//		}
	//	}
	//	std::cout << "b2 " <<b2 << " "<< begin_on_read << " "<<read_counter-1 <<std::endl;
		if(read_counter >0){
			e2 = begin_on_read + read_counter-1;
		}else{
			assert(read_counter==0);//If only gap on the read
			size_t this_length = data.getSequence(ref2).length();
			std::cout<< "there is only gap on read " << b2 <<std::endl;
		//	b2 = this_length;
		//	b2 = e2; //we dont get only gap alignments one after each other here, cus it is just after an existing one
			e2 = this_length;
		}
		//If current forwards
		if(this_path.at(0)>0){
			b1 = data.getSequence(this_refs.at(0)).length()-first_length;
		//	e1= data.getSequence(this_refs.at(0)).length()-1;
			e1 = b1+ref_counter-1;
			assert(b1+ref_counter-1 == data.getSequence(this_refs.at(0)).length()-1);

		}else{//If current reverse, using the reverse complement!
			e1 = data.getSequence(this_refs.at(0)).length()-first_length;
		//	b1 = data.getSequence(this_refs.at(0)).length()-1;
			b1 = e1+ref_counter-1;
			assert(e1+ref_counter-1 == data.getSequence(this_refs.at(0)).length()-1);
		}
		std::cout<< "b2 "<< b2 << "e2 "<< e2 <<std::endl;
		pw_alignment p0(str1,str2,b1,b2,e1,e2,this_refs.at(0),ref2);
		p0.print();
		add_adjacencies(comp,first_al,p0,refacc,readacc);
		previous= p0;
		current_pos++;
		if(read_counter == 0){
			b2 = begin_on_read;
			e2 = begin_on_read;
		}else{
			b2 = e2+1;
		}
	}else{
		previous = first_al;
	}
	//Making all those are in between 
	for(size_t i = 1; i < this_refs.size()-1;i++){
		int edge = this_path.at(i);
		std::string sample1;
		std::string sample2;
		std::cout << "current pos "<< current_pos << std::endl;
		bool tillend = false;
		if(i == this_refs.size()-2 && last_length == 0){
			tillend=true;
		}
		compute_samples(tillend, data.getSequence(this_refs.at(i)).length(), current_pos, read_counter, ref_counter, read_out, sample2, ref_out, sample1);
		if(edge < 0){
		//	b1 = data.getSequence(this_refs.at(i)).length()-1;
			b1 = ref_counter-1;
			e1 = 0;
		}else{
			b1 = 0;
		//	e1 = data.getSequence(this_refs.at(i)).length()-1;//This could be wrong when we had gap at the end of an al
			e1 = ref_counter -1;
		}
		assert(data.getSequence(this_refs.at(i)).length()-1 == ref_counter-1);
		std::cout << "current pos after "<< current_pos << " read counter "<< read_counter << std::endl;
		if(read_counter >0){
			e2 = begin_on_read + read_counter-1;
			std::cout<< "b2 "<< b2 << " e2 "<<e2 <<std::endl;
			if(sample1.length() == 1){
				std::cout << "sample1 "<<sample1 <<std::endl;
				std::cout << "sample2 "<<sample2 <<std::endl;
			}
			pw_alignment p(sample1,sample2,b1,b2,e1,e2,this_refs.at(i),ref2);
			std::cout<<"p "<<std::endl;
			p.print();
			if(e2 != end_on_read){
				b2 = e2+1;
			}else{//TODO check it!
				std::cout<< "end of read part is reached! "<<std::endl;
				b2 = e2;
			}
			add_adjacencies(comp,previous,p,refacc,readacc);
			previous = p;
		}else{
			assert(read_counter==0);
			size_t this_length = data.getSequence(ref2).length();
			size_t left2,right2;
			previous.get_lr2(left2,right2);
			std::cout<<"l2 "<< left2 << " r2 "<<right2 <<std::endl;
			size_t coordinate = right2;
			if(right2 == this_length){
				coordinate = left2;
			}
			pw_alignment p(sample1,sample2,b1,coordinate,e1,this_length,this_refs.at(i),ref2);
			std::cout<<"p when all gap on read "<<std::endl;
			p.print();
			add_adjacencies(comp,previous,p,refacc,readacc);
			previous = p;
		}
		current_pos++;
	
	}
	//Making the last al
	if(last_length>0){// TODO check if it is possible read_counter be 0 at this part!
		std::cout<<"last length > 0"<<std::endl;
		str1.clear();
		str2.clear();
		counter =0;
		size_t r_counter = 0;
	//	if(begin_on_read != end_on_read){
	//		b2 = e2+1;
	//	}
		
		e2 = end_on_read;
	//	if(current_pos != 0){
	//		current_pos++;
	//	}
		for(size_t i = current_pos; i < ref_out.size();i++){
		//	if(counter < last_length){
				current_pos = i;
				str1+=ref_out.at(i);
				str2+=read_out.at(i);
				if(ref_out.at(i)!= '-'){
					counter++;
				}
				if(read_out.at(i)!='-'){
					read_counter++;
					r_counter++;
				}
		//	}else{
		//		break;
		//	}
		}
		assert(counter == last_length);
		if(this_path.back()>0){
			b1 = 0;
			e1= last_length-1;
		}else{
			b1 = last_length-1;
			e1 = 0;
		}
		if(r_counter== 0){
			std::cout << "last alignment has only gap on read! "<< b2 << " "<< begin_on_read << " " << end_on_read <<std::endl;
		//	size_t this_length = data.getSequence(ref2).length();
		//	size_t left2,right2;
		//	previous.get_lr2(left2,right2);
		//	std::cout<<"l2 "<< left2 << " r2 "<<right2 <<std::endl;
		//	size_t coordinate = right2;
		//	if(right2 == this_length){
		//		coordinate = left2;
		//	}
		//	std::cout<< "coordinate "<< coordinate <<std::endl;
		//	assert(b2 == end_on_read);
			e2 = data.getSequence(ref2).length();
		}
		pw_alignment p1(str1,str2,b1,b2,e1,e2,this_refs.back(),ref2);
		p1.print();
		add_adjacencies(comp,previous, p1,refacc,readacc);
		previous = p1;
		current_pos++;
	}else{//Dont need to do anything
		std::cout<<"last length = 0"<<std::endl;
	}
	add_adjacencies(comp,previous,last_al,refacc,readacc);
//	adjacencies.insert(std::make_pair(previous,last_al));
	std::cout << ref_out.length()<< " "<<current_pos<<std::endl;
	assert(ref_out.length()==current_pos);
	std::cout << begin_on_read << " "<<read_counter<< " "<< end_on_read <<std::endl;
	if(begin_on_read != end_on_read && read_counter != 0){
		assert(begin_on_read+read_counter-1 == end_on_read);
	}
}
template<typename T>
void als_components<T>::compute_samples(bool tillend, size_t node_length, size_t & current_pos, size_t & read_counter, size_t & ref_counter,std::string & read_in, std::string & read_out, std::string & ref_in, std::string & ref_out){
	ref_counter = 0;
//	if(tillend==false){
//		for(size_t i = current_pos; i< ref_in.length();i++){
//			if(counter < node_length){
//				current_pos = i;
//				ref_out+=ref_in.at(i);
//				read_out+=read_in.at(i);
//				if(ref_in.at(i)!= '-'){
//					counter++;
//				}
//				if(read_in.at(i)!='-'){
//					read_counter++;
//				}
//			}else{
//				break;
//			}
//		}
//	}else{
		for(size_t i = current_pos; i< ref_in.length();i++){
			current_pos = i;
			ref_out+=ref_in.at(i);
			read_out+=read_in.at(i);
			if(ref_in.at(i)!= '-'){
				ref_counter++;
			}
			if(read_in.at(i)!='-'){
				read_counter++;
			}
			if(ref_counter == node_length && tillend == false){
				break;
			}
		}
//	}
}
template<typename T>
void als_components<T>::get_reverse_complement(std::string & sequence_in , std::string & sequence_out){
	for(size_t i = sequence_in.size(); i >0; i--){
		sequence_out += dnastring::complement(sequence_in.at(i-1));
	}
}

template<typename T>
void als_components<T>::add_first_and_last_als(size_t & i, const pw_alignment & p, size_t & suc_size, size_t & first_left, size_t & last_right){//Create alignments that connect first and last alignments to the origin and the end.	
	size_t l2,r2;
	p.get_lr2(l2,r2);
	size_t length = data.getSequence(p.getreference2()).length();
	assert(p.getreference1()< data.numSequences() && p.getreference2() < data.numSequences() && l2 < length);
//	if(l2 != data.numSequences()+2){//Do not add gap only als to the begin nor end! (made no sense!)
		if(l2 < MAXGAP || (l2-(first_left-MAXGAP/2)<MAXGAP)){
			std::cout<< "add it to the beginning "<< std::endl;
			p.print();
			looking_for_first_al(i,p,first_left);	 
		}
	/*	if(p.getreference2()==1275){
			std::cout<< "length " << length <<std::endl;
			std::cout<< "r2 "<< r2 <<std::endl;
		}*/
		if((length - r2 < MAXGAP) || ((last_right-r2) + MAXGAP/2 <MAXGAP)){
			std::cout<< "add it to the end " << last_right << std::endl;
			p.print();
			if(suc_size == 0) get_subgraph(p);
			looking_for_last_al(i,p,last_right);
		}
//	}
}
template<typename T>
void als_components<T>::looking_for_first_al(size_t & comp, const pw_alignment & p, size_t & first_left){//We may create more than one al but keep the one with the highest gain from ref to read
	std::vector<std::vector<pw_alignment> >all_first_als;
//	size_t type = 3;
	size_t type = 2; 

	size_t from,to;
	size_t l1,l2,r1,r2;
	p.get_lr1(l1,r1);
	p.get_lr2(l2,r2);
	unsigned int ref2 = p.getreference2();
	unsigned int ref1 = p.getreference1();
	std::cout << "on node "<< ref1 <<std::endl;
	size_t readacc = data.accNumber(ref2);//this function is not used for gap only alignments so there is no problem using data.accNumber()
	size_t refacc = data.accNumber(ref1);
	std::string seq_from_ref;
/*	if(l1 != 0){
		from = 0;
		to = l1-1;
		if(p.getbegin1()<p.getend1()){
			seq_from_ref = data.extract_seq_part(ref1,from,to);
		}else{//Get the reverse if begin1 > end 1
			std::cout<<"looking for first als begin > end "<<std::endl;
			std::string temp_seq_in = data.extract_seq_part(ref1, from, to);
			std::string temp_seq_out;
			get_reverse_complement(temp_seq_in, temp_seq_out);
			seq_from_ref = temp_seq_out;
		}
	}*/
	size_t ref_from, ref_to;
	if(l1 != 0 && p.getbegin1() < p.getend1()){
		from = 0;
		to = l1-1;
		seq_from_ref = data.extract_seq_part(ref1,from,to);
		ref_from = 0;
		ref_to = l1-1;
	}
	if(r1 != data.getSequence(ref1).length()-1 && p.getbegin1()>p.getend1()){
		to = r1+1;
		from = data.getSequence(ref1).length()-1;
		std::cout<<"looking for first als begin > end "<<std::endl;
		std::string temp_seq_in = data.extract_seq_part(ref1, to, from);
		std::string temp_seq_out;
		get_reverse_complement(temp_seq_in, temp_seq_out);
		seq_from_ref = temp_seq_out;
		ref_from = data.getSequence(ref1).length()-1;
		ref_to = r1+1;
	}
	std::vector<std::vector<int> >refs;
	std::string seq;
	size_t new_al_begin;
	size_t new_al_end;

	if(l2 < MAXGAP){
		from = 0;
		if(l2 != 0){
			to = l2-1;
			seq = data.extract_seq_part(ref2, from, to);
			new_al_begin = 0;
			new_al_end = l2-1;
		}
	}else{
		std::cout << " l2 "<< l2 << " first left "<< first_left <<std::endl;
		assert(l2-(first_left-MAXGAP/2)<MAXGAP);
		from = first_left-MAXGAP/2;
		to = l2-1;
		seq = data.extract_seq_part(ref2, from, to);
		new_al_begin = from;
		new_al_end = l2-1;
	}
	if(seq.length() != 0){
		if(seq_from_ref.length()>=MAXGAP){
			//in this case we only use seq_from_ref and seq to make an alignemnt, pick MAXGAP baes of it and use NW to make an alignment
			std::string sub;
			if(p.getbegin1() < p.getend1()){
				sub = seq_from_ref.substr(seq_from_ref.length()-MAXGAP,MAXGAP);
				assert(sub.length() == MAXGAP);
				ref_from = l1-MAXGAP;
				ref_to = l1-1;
			}else{
				sub = seq_from_ref.substr(seq_from_ref.length()-MAXGAP,MAXGAP);
				assert(sub.length() == MAXGAP);

				ref_from = r1+MAXGAP;
				ref_to = r1+1;
			}
			seq_from_ref = sub;
			std::cout <<"ref seq" << seq_from_ref <<std::endl;
			std::string read_out;
			std::string ref_out;
			std::cout << "length of seq from ref "<< sub.length() << std::endl;
			needleman<T> nl(data, model, seq, seq_from_ref);
//XXX			simpleNW nl(model, seq, seq_from_ref);
			nl.run_needleman(readacc,refacc,type,read_out,ref_out);
			size_t count = 0;
			for(size_t i = 0; i < ref_out.length();i++){
				if(ref_out.at(i)!='-'){
					count++;
				}
			}
			if(p.getbegin1() < p.getend1() && count !=0){
				ref_from = l1 - count;
			}else if(p.getbegin1() > p.getend1() && count !=0 ){
				ref_from = r1+count;
			}else{
				assert(count == 0);
				ref1 = data.numSequences()+1;
				ref_from = 0;
				ref_to = 0;
			}
			pw_alignment p1(ref_out,read_out,ref_from,from,ref_to,to,ref1,ref2);
			add_adjacencies(comp,p1,p,refacc,readacc);
			add_edge_to_begin(comp,p1,refacc,readacc);

		}else{
			std::set<std::vector<int> > refs;//seq name it is not the ref id! //TODO what about loops?
//			size_t length = MAXGAP- seq_from_ref.length();
			size_t length = seq_from_ref.length();
			size_t length_on_read = seq.length();
			if(p.getbegin1()<p.getend1()){
				std::cout << "length is "<< length << "ref1 is "<< ref1 <<std::endl;
				refs=rgraph.get_predecessor(ref1,true,length,length_on_read);
				std::cout << refs.size() <<std::endl;
			}else{
				refs=rgraph.get_predecessor(ref1,false,length,length_on_read);
			}
			if(refs.size()!=0){
				make_first_nw_als(comp, p ,refs,seq_from_ref, ref2, from, to, readacc, refacc);
/*			//	std::reverse(seq.begin(),seq.end()); XXX Just commented
				std::cout<< "pre_refs: " << refs.size() <<std::endl;
			//	for(size_t i = 0 ; i < refs.size() ;i++)
				std::map<std::vector<int> ,size_t > seen_path;
				std::map<int, std::set<int> > seen_edge;
				for(std::set<std::vector<int> >::iterator it = refs.begin() ; it != refs.end() ; it++){
					std::string refin;
					std::vector<int> this_p = *it;
					std::vector<int> pre_path;
					std::vector<int> current_path;
					bool ALREADY_A_PATH = false;
					for(size_t j = this_p.size() ; j > 1; j--){//TODO Now it became so slow!! avoid making already existing alignments
						std::cout<< "this node: "<<this_p.at(j-1)<<std::endl;
						bool CONTINUE = false;	
						pre_path.push_back(this_p.at(j-1));
						std::map<int, std::set<int> >::iterator SEEN = seen_edge.find(this_p.at(j-1));
					//	if(SEEN != seen_edge.end()){
					//		std::map<std::vector<int>,size_t >::iterator SEEN_PATH = seen_path.find(pre_path);
					//		if(SEEN_PATH != seen_path.end()){
					//			ALREADY_A_PATH = true;
					//		//	current_pos_on_read = SEEN_PATH->second;
					//		}else{
					//			ALREADY_A_PATH = false;
					//		//	current_pos_on_read = 0;
					//		}
					//	}
					//	if(SEEN == seen_edge.end() || ALREADY_A_PATH == false||CONTINUE==true){ //CONTINUE = true; TODO
					//		CONTINUE = true;
							//Make an al between each member of refs and seq and pick the best
							current_path.push_back(this_p.at(j-1));
							from = 0;
							if(this_p.at(j-1)>0){
								unsigned int ref_id = rgraph.get_refid(refacc,this_p.at(j-1));//gets node name, retruns node id
								to = data.getSequence(ref_id).length()-1;
								std::cout<< data.extract_seq_part(ref_id, from, to)<<std::endl;
								refin += data.extract_seq_part(ref_id, from, to);
							}else{
								int temp = -1*(this_p.at(j-1));
								unsigned int ref_id = rgraph.get_refid(refacc,temp);
								to = data.getSequence(ref_id).length()-1;// Need to add  its reverse complement!
								std::string temp_seq_in = data.extract_seq_part(ref_id, from, to);
								std::string temp_seq_out;
								get_reverse_complement(temp_seq_in, temp_seq_out);
								refin += temp_seq_out;
							}
					//	}
					}
					std::cout << "seq from ref " << seq_from_ref <<std::endl;
					refin += seq_from_ref;//keep only the last MAXGAP bases 
					if(refin.length() > MAXGAP){
						refin.erase(refin.begin(),refin.end()-MAXGAP);
					}
					assert(refin.length()<=MAXGAP);
					std::cout<< "refin length "<< refin.length() << "seq length " << seq.length() <<std::endl;
				//	std::reverse(refin.begin(),refin.end()); Just commented
					std::string read_out, ref_out;
					needleman<dynamic_mc_model> nl(data, model, seq, refin);//TODO How much of seq is left?
					nl.run_needleman(readacc,refacc,type,read_out,ref_out);
					//Make all the als!
					size_t read_id = ref2;
					size_t current_node = ref1;
					std::vector<pw_alignment> first_als;
				//	make_first_als(this_p,ref_out,read_out,new_al_begin,new_al_end,seq_from_ref,l1,refacc,read_ref, current_node,first_als);
					if(seq_from_ref.length()==0) this_p.erase(this_p.begin());
					size_t remainder = seq_from_ref.length();
					make_first_als(this_p, ref_out , read_out, l2 , refacc, readacc, read_id, first_als,remainder); 
					all_first_als.push_back(first_als);
					std::cout<< "firsts are added! "<<std::endl;
				}
				//Choose the best set of als among those you created above (sum of their mod cost is the lowest)
				const std::vector<pw_alignment> best_als = find_the_best_als(all_first_als,refacc,readacc);
				std::cout<< "best als size "<< best_als.size() <<std::endl;
				add_adjacencies(comp, best_als.at(0), p,refacc,readacc);
				if(best_als.size()>=2){
					for(size_t i =0; i < best_als.size()-1;i++){
						add_adjacencies(comp,best_als.at(i+1),best_als.at(i),refacc,readacc);
					}
				}
				add_edge_to_begin(comp,best_als.at(best_als.size()-1),refacc,readacc); */
			}else{
				std::cout << "refs size is 0 "<<std::endl;
				//Gap + seqfromref
				if(seq_from_ref.length()!=0){
					std::reverse(seq.begin(),seq.end());
					std::string refin = seq_from_ref;
					std::reverse(refin.begin(),refin.end());
					std::cout<< "refin "<< refin << std::endl;
					std::string read_out, ref_out;
					needleman<T> nl(data, model, seq, refin);
			//XXX		simpleNW nl(model, seq, refin);
					nl.run_needleman(readacc,refacc,type,read_out,ref_out);//If gap only on the ref_out, change the ref to the virtual one! numseq+1
					std::reverse(read_out.begin(),read_out.end()); 
					std::reverse(ref_out.begin(),ref_out.end());
					std::cout << "ref out "<< ref_out << std::endl;
					bool gapOnly = true;
					for(size_t k = 0; k < ref_out.size() ; k++){
						if(ref_out.at(k)!='-'){
							gapOnly = false;
						}
					}
					if(gapOnly == false){
						pw_alignment p1(ref_out,read_out,ref_from,new_al_begin,ref_to,new_al_end,ref1,ref2);
						add_adjacencies(comp, p1, p,refacc,readacc);
						add_edge_to_begin(comp,p1,refacc,readacc);
					}else{
						unsigned int temp = data.numSequences()+1;//only gap on ref
						pw_alignment p1(ref_out,seq,0,new_al_begin,0,new_al_end,temp,ref2);
						add_adjacencies(comp, p1, p,refacc,readacc);
						add_edge_to_begin(comp,p1,refacc,readacc);
					}
				}else{//Only gap!
					unsigned int temp = data.numSequences()+1;//only gap on ref
					std::string str_ref(seq.length(),'-');
					pw_alignment p1(str_ref,seq,0,new_al_begin,0,new_al_end,temp,ref2);
					add_adjacencies(comp, p1, p,refacc,readacc);
					add_edge_to_begin(comp,p1,refacc,readacc);
				}
			}
		}
	}else{
		//There will no alignment before this node and it is directly connected to the start node
		//XXX HERE! 
		std::cout << "seq length 0 "<<std::endl;
		std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator it=node_indices.at(comp).find(p);
		p.print();
		if(it!=node_indices.at(comp).end()){
			std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(0);
			if(adj == adjacencies.at(comp).end()){
				adjacencies.at(comp).insert(std::make_pair(0, std::set<size_t>()));
				adj = adjacencies.at(comp).find(0);
			}
			adj->second.insert(it->second);
			//adjacencies.at(comp).insert(std::make_pair(0,it->second));
		}else{
			node_indices.at(comp).insert(std::make_pair(p,index));
			indices_nodes.at(comp).insert(std::make_pair(index,p));
			index++;
			it=node_indices.at(comp).find(p);
			std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(0);
			if(adj == adjacencies.at(comp).end()){
				adjacencies.at(comp).insert(std::make_pair(0, std::set<size_t>()));
				adj = adjacencies.at(comp).find(0);
			}
			adj->second.insert(it->second);
			adj = adjacencies.at(comp).find(it->second);
			assert(adj == adjacencies.at(comp).end());
			adjacencies.at(comp).insert(std::make_pair(it->second, std::set<size_t>()));


		//	adjacencies.at(comp).insert(std::make_pair(0,it->second));
		}

		double c1,c2,m1,m2;
		assert( p.getreference2() != data.numSequences()+2);
		assert(p.getreference1() != data.numSequences()+1); //XXX ??? cant it happen?? It could i guess!
		if(p.getreference1() != data.numSequences()+1){
		//	model.cost_function(p,c1,c2,m1,m2);
			model.cost_function(p, m1, m2, refacc, readacc);
		}else if(p.getreference1() == data.numSequences()+1){//m1 is the cost of insertion from ref accession to the read accession 
			m1 = 2.5 * p.alignment_length();//TODO change to insertion from ref accession to the read accession
		}
		weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(0,it->second),m1)); 
	}


}
template<typename T>
void als_components<T>::make_first_nw_als(size_t & comp, const pw_alignment & current_al, std::set<std::vector<int> > & refs, std::string & seq_from_ref, unsigned int & read_id, size_t & ini_read_from , size_t & read_to , size_t & read_acc, size_t & ref_acc){
	std::set<std::vector<int> > seen_path;
	std::map<int, std::set<int> > seen_node;
	alignments_on_each_ref_node.clear();
	std::string seq_from_read = data.extract_seq_part(read_id, ini_read_from,read_to);
	std::cout<< "on "<< read_id << " from "<< ini_read_from << " to "<< read_to <<std::endl;
	std::cout << seq_from_read <<std::endl;
	for(std::set<std::vector<int> >::iterator it = refs.begin() ; it != refs.end() ; it++){
		std::string refin;
		std::vector<int> this_p = *it;
		std::vector<int> pre_path;
		bool ALREADY_A_PATH = false;
		std::cout << "path size "<< this_p.size() << std::endl;
		for(size_t j = 0 ; j < this_p.size() ; j++){ //TODO break when no more base on read to go through
			std::cout<< "this node: "<<this_p.at(j)<<std::endl;
			bool CONTINUE = false;	
			if(j == 0){//TODO
				std::cout<< "seq from ref "<< seq_from_ref.length()<<std::endl;
				if(seq_from_ref.length()>0){
					assert(pre_path.size()==0);
					pre_path.push_back(this_p.at(j));
					std::set<std::vector<int> >::iterator a_path = seen_path.find(pre_path);
					if(a_path==seen_path.end()){
						pw_alignment this_pal;
						seen_path.insert(pre_path);
						std::string read_out, ref_out;
						needleman<T> nl(data, model, seq_from_read, seq_from_ref);
					//XXX	simpleNW nl(model, seq_from_read, seq_from_ref);
						size_t type = 5;
						nl.run_needleman(read_acc,ref_acc,type,read_out,ref_out);
						size_t read_counter = 0;
						for(size_t i = 0; i < read_out.size();i++){
							if(read_out.at(i)!= '-'){
								read_counter ++;
							}
						}
						size_t read_from = read_to- read_counter +1;
						if(read_counter == 0){
							read_to = data.getSequence(read_id).length();
						}
						size_t from = 0;
						size_t to =seq_from_ref.length() -1;
						unsigned int ref_id = rgraph.get_refid(ref_acc,this_p.at(j));//gets node name, retruns node id
						if(this_p.at(j)>0){ 
							std::cout << "is it here?? "<< this_p.at(j)<< std::endl;
							pw_alignment nw_al(ref_out,read_out,from, read_from ,to, read_to, ref_id, read_id);
							nw_al.print();
							this_pal = nw_al;
							add_to_containers(comp, nw_al, current_al, ref_acc, read_acc, this_p.at(j), this_p.at(j));
							if(read_from == ini_read_from){
								add_edge_to_begin(comp, nw_al, read_acc, ref_acc);
							}
						}else{
							from = data.getSequence(ref_id).length()-1;
							to = data.getSequence(ref_id).length()-1 - to; 
							std::cout << "from " << from << " to " << to <<std::endl;
							pw_alignment nw_al(ref_out,read_out,from, read_from ,to, read_to, ref_id, read_id);
							nw_al.print();
							this_pal = nw_al;
							add_to_containers(comp, nw_al, current_al , ref_acc, read_acc, this_p.at(j), this_p.at(j));
							if(read_from == ini_read_from){
								add_edge_to_begin(comp, nw_al, read_acc, ref_acc);
							}
						}
						if(this_p.size() ==1 && read_from != ini_read_from){
						//Gap only on ref or extending the previous one? TODO
							size_t l2,r2;
							this_pal.get_lr2(l2,r2);
							assert(l2 == read_from);
							std::string str_read = seq_from_read.substr(0, read_from);
							std::string str_ref(str_read.length(),'-');
							read_to = read_from -1;
							unsigned int gap_id = data.numSequences()+1;
							pw_alignment gap_only(str_ref, str_read, 0, ini_read_from, 0, read_to,gap_id,read_id);
							add_to_containers(comp, gap_only, this_pal , ref_acc, read_acc, this_p.at(j), this_p.at(j));
							add_edge_to_begin(comp, gap_only, read_acc, ref_acc);
						}
						//Add it to the current_al, create an al between the remaining parth and read, add it to seen path and pre_al ??
					}else{
						//DO nothing!
					}
				}else{
					assert(pre_path.size()==0);
					pre_path.push_back(this_p.at(j));
					std::set<std::vector<int> >::iterator a_path = seen_path.find(pre_path);
					if(a_path==seen_path.end()){
						seen_path.insert(pre_path);
						std::set<size_t> temp;
						std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator al_id = node_indices.at(comp).find(current_al);
						if(al_id == node_indices.at(comp).end()){
							node_indices.at(comp).insert(std::make_pair(current_al,index));
							indices_nodes.at(comp).insert(std::make_pair(index,current_al));
							index++;
							al_id = node_indices.at(comp).find(current_al);
						}
						assert(al_id != node_indices.at(comp).end());
						temp.insert(al_id->second);
						std::cout << "al id second "<< al_id->second << " " << this_p.at(j)<<std::endl;
						alignments_on_each_ref_node.insert(std::make_pair(this_p.at(j),temp));
						if(this_p.size() ==1){
							std::string str_ref(seq_from_read.length(),'-');
							unsigned int gap_id = data.numSequences()+1;
							pw_alignment gap_only(str_ref, seq_from_read, 0, ini_read_from, 0, read_to,gap_id,read_id);
							add_to_containers(comp, gap_only, current_al , ref_acc, read_acc, this_p.at(j), this_p.at(j));
							add_edge_to_begin(comp, gap_only, read_acc, ref_acc);
						}
					}else{
						//DO nothing!
					}
				}
			}else{
		//	if(j != 0 && j != this_p.size()-1){//j != 0 //Type is different from j = 0, in this case we wanna see all the ref, no matter if we get only gap on read!
				std::cout << pre_path <<std::endl;
				bool NOEDGE = false;
				std::map<int, std::set<int> >::iterator SEEN = seen_node.find(this_p.at(j));
				if(SEEN != seen_node.end()){
					std::set<int>::iterator edge = SEEN->second.find(this_p.at(j-1));
					if(edge!= SEEN->second.end()){
						std::set<std::vector<int> >::iterator SEEN_PATH = seen_path.find(pre_path);
						if(SEEN_PATH != seen_path.end()){
							ALREADY_A_PATH = true;
						}else{
							ALREADY_A_PATH = false;
						}
					}else{
						NOEDGE = true;
					}
				}
				std::cout << "HERE! "<<std::endl;
				if(SEEN == seen_node.end() || NOEDGE == true || ALREADY_A_PATH == false||CONTINUE==true){
					CONTINUE = true;
					if(SEEN == seen_node.end()){
						std::set<int> temp;
						temp.insert(this_p.at(j-1));
						seen_node.insert(std::make_pair(this_p.at(j),temp));
					}else if(NOEDGE == true){
						SEEN->second.insert(this_p.at(j-1));
					}
					if(ALREADY_A_PATH == false){
						pre_path.push_back(this_p.at(j));
						seen_path.insert(pre_path);					
					}
					std::map<int , std::set<size_t> >::iterator pre_al = alignments_on_each_ref_node.find(this_p.at(j-1));
					if(pre_al != alignments_on_each_ref_node.end()){
						size_t from = 0;
						size_t to;
						unsigned int ref_id;
						if(this_p.at(j)>0){
							ref_id = rgraph.get_refid(ref_acc,this_p.at(j));//gets node name, retruns node id
							to = data.getSequence(ref_id).length()-1;
							std::cout<< data.extract_seq_part(ref_id, from, to)<<std::endl;
							refin = data.extract_seq_part(ref_id, from, to);
						}else{
							int temp = -1*(this_p.at(j));
							ref_id = rgraph.get_refid(ref_acc,temp);
							to = data.getSequence(ref_id).length()-1;// Need to add  its reverse complement!
							std::string temp_seq_in = data.extract_seq_part(ref_id, from, to);
							std::string temp_seq_out;
							get_reverse_complement(temp_seq_in, temp_seq_out);
							refin = temp_seq_out;
						}
						if(refin.length() > MAXGAP){
							refin.erase(refin.begin(),refin.end()-MAXGAP);
							from = to - MAXGAP +1;
						}
						assert(refin.length()<=MAXGAP);
						std::cout<< "refin length "<< refin.length() <<std::endl;
						std::cout<< "pre size "<< pre_al->second.size() <<std::endl;
						for(std::set<size_t>::iterator pre_als = pre_al->second.begin(); pre_als != pre_al->second.end() ; pre_als ++){
							std::map<size_t, const pw_alignment>::iterator id =indices_nodes.at(comp).find(*pre_als);
							std::cout << "pre als: "<< *pre_als << std::endl;
							assert(id != indices_nodes.at(comp).end());
							const pw_alignment p = id->second;
							std::cout << "pre al "<<std::endl;
							p.print();
							size_t l2,r2;
							p.get_lr2(l2,r2);
							if(l2 == ini_read_from) continue;
							if(r2 == data.getSequence(read_id).length()){
								read_to = l2;
							}else{
								read_to = l2-1;
							}
							std::string readin = data.extract_seq_part(read_id , ini_read_from, read_to);
							std::cout << "read in "<< readin << std::endl;
							std::string read_out, ref_out;
							needleman<T> nl(data, model, readin, refin);
						//XXX	simpleNW nl(model, readin, refin);
							size_t type = 5;
							if(j == this_p.size()-1) type = 2;
							nl.run_needleman(read_acc,ref_acc,type,read_out,ref_out);
							size_t read_counter = 0;
							for(size_t i = 0; i < read_out.size();i++){
								if(read_out.at(i)!= '-'){
									read_counter ++;
								}
							}
							size_t read_from = read_to- read_counter +1;
							if(read_counter == 0){
								read_to = data.getSequence(read_id).length();
							}
							if(this_p.at(j)>0){ 
								pw_alignment nw_al(ref_out,read_out,from, read_from ,to, read_to, ref_id, read_id);
								nw_al.print();
								add_to_containers(comp, nw_al, p, ref_acc, read_acc, this_p.at(j), this_p.at(j-1));
								if(read_from == ini_read_from){
									std::cout << "first al is here! "<<std::endl;
									add_edge_to_begin(comp, nw_al, read_acc, ref_acc);
								}
							}else{
								pw_alignment nw_al(ref_out,read_out,to, read_from ,from, read_to, ref_id, read_id);
								nw_al.print();
								add_to_containers(comp, nw_al, p , ref_acc, read_acc, this_p.at(j), this_p.at(j-1));
								if(read_from == ini_read_from){
									std::cout << "first al is here1! "<<std::endl;
									add_edge_to_begin(comp, nw_al, read_acc, ref_acc);
								}
							}
						}
					}else{
						std::cout << "end of read is reached! "<<std::endl;
					}
					
				}
			}
			if(pre_path.back() != this_p.at(j)){
				pre_path.push_back(this_p.at(j));
			}

		}
	}
					
}

template<typename T>
void als_components<T>:: make_last_nw_als(size_t & comp, int& node_name, const pw_alignment & current_al, std::set<int> & nodes_on_paths ,std::string & seq_from_read, unsigned int & ref_id, unsigned int & read_id, size_t & read_from , size_t & read_to , size_t & read_acc, size_t & ref_acc){//TODO think of adding the posibility of handling loops
	std::cout << "nodes on path size "<< nodes_on_paths.size() <<std::endl;
	std::multimap<double, size_t> al_distance;
	std::map<size_t , double> weight;
	std::map<size_t,size_t> cur_pre;
	pw_alignment this_al = current_al;
	std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator id = node_indices.at(comp).find(this_al);
	if(id == node_indices.at(comp).end()){
		node_indices.at(comp).insert(std::make_pair(this_al , index));
		indices_nodes.at(comp).insert(std::make_pair(index, this_al));
		index++;
		id = node_indices.at(comp).find(this_al);
	}
	assert(id != node_indices.at(comp).end());
	weight.insert(std::make_pair(id->second, 0.0));
	std::string refin;
	int startnode = node_name;
	std::map<int , std::set<int> > adjacencies = rgraph.get_adjacencies();
	std::map<int, std::set<int> >::iterator adj = adjacencies.find(startnode);
	assert(adj != adjacencies.end());
	bool LASTAL = false;
	while(read_from != read_to+1){
		this_al.print();
		std::cout << "startnode "<< startnode <<std::endl;
		std::cout << read_from << " " << read_to <<std::endl;
		adj = adjacencies.find(startnode);
		if(adj == adjacencies.end()){
			unsigned int temp = data.numSequences()+1;//only gap on ref
			size_t length = read_to- read_from +1;
			std::string str_ref(length,'-');
			std::string readin = data.extract_seq_part(read_id,read_from,read_to);
			pw_alignment p1(str_ref,readin,0,read_from,0,read_to,temp,read_id);
			add_adjacencies(comp, this_al, p1,ref_acc,read_acc);
			add_edge_to_end(comp,p1);
			LASTAL = true;
			break;
		}
		id = node_indices.at(comp).find(this_al);
		assert(id != node_indices.at(comp).end());
		std::map<size_t , double>::iterator pre_wei = weight.find(id->second);
		assert(pre_wei != weight.end());
		std::set<int> this_set = adj->second;
		for(std::set<int>::iterator this_adj = this_set.begin() ; this_adj != this_set.end() ; this_adj++){
			int this_name = *this_adj;
			std::cout << "this_adj " << this_name << std::endl;
			std::set<int>::iterator onPaths = nodes_on_paths.find(this_name);
			if(onPaths == nodes_on_paths.end()){
				std::cout << "not on path! "<<std::endl;
				continue;
			}
			size_t ref_from = 0;
			size_t ref_to = 0;
			if(this_name > 0){
				ref_id = rgraph.get_refid(ref_acc,this_name);
				ref_to =  data.getSequence(ref_id).length()-1;
				refin = data.extract_seq_part(ref_id, ref_from, ref_to);
			}else{
				int temp = -1*(this_name);
				ref_id = rgraph.get_refid(ref_acc,temp);
				ref_to = data.getSequence(ref_id).length()-1;
				std::string temp_seq_in = data.extract_seq_part(ref_id, ref_from, ref_to);
				std::string temp_seq_out;
				get_reverse_complement(temp_seq_in, temp_seq_out);
				refin = temp_seq_out;
				size_t temp_ref_from = ref_from;
				ref_from = ref_to;
				ref_to = temp_ref_from;
			}
			std::string readin = data.extract_seq_part(read_id,read_from,read_to);
			std::string refout;
			std::string readout;
			pw_alignment nextal;
			size_t type = 3;
			if(refin.length() > MAXGAP ){
				//Take a sub string and finish it here!(type should be changed)
				refin = refin.substr(0,400);
				type = 4;//end of read should be reached
			}
			bool NOEDGE = true;
			std::map<int, std::set<int> >::iterator next_adjs = adjacencies.find(this_name);
			if(next_adjs == adjacencies.end()) NOEDGE = true;
			for(std::set<int>::iterator nex_adj = next_adjs->second.begin() ; nex_adj != next_adjs->second.end() ; nex_adj++){
				std::set<int>::iterator onPaths = nodes_on_paths.find(*nex_adj);
				if(onPaths != nodes_on_paths.end()){
					NOEDGE = false;
				}
			}
			if(NOEDGE == true){
				type = 4;
			}
			assert(refin.length() <= MAXGAP);
			needleman<T> nl(data, model, readin, refin);
			nl.run_needleman(read_acc,ref_acc,type,readout,refout);
			size_t readcounter = 0;
			for(size_t i = 0; i < readout.size();i++){
				if(readout.at(i)!= '-') readcounter++;
			}
			size_t this_to_on_read;
			if(readcounter>0 ){
				this_to_on_read = read_from + readcounter -1;
				std::cout << this_to_on_read << " " << readcounter << " " << read_from <<std::endl;
				assert(this_to_on_read <= read_to);
				if(this_to_on_read == read_to){ //XXX test
					type = 6;
					needleman<T> nl(data, model, readin, refin);
					refout.clear();
					readout.clear();
					nl.run_needleman(read_acc,ref_acc,type,readout,refout);
					size_t read_counter = 0;
					for(size_t i = 0; i < readout.size();i++){
						if(readout.at(i)!= '-') read_counter++;
					}
					std::cout << this_to_on_read << " " << read_counter << " " << read_from <<std::endl;
					assert(this_to_on_read == read_counter + read_from -1);
				}
			}else{
				this_to_on_read = data.getSequence(read_id).length();
			}
			pw_alignment nwal(refout,readout,ref_from,read_from,ref_to,this_to_on_read,ref_id,read_id);
			nextal = nwal;
			add_adjacencies(comp, this_al , nextal , ref_acc, read_acc);//TODO I can do it locally in here, and later on which ever is chosen be added to the alignment graph		
			std::cout << "nextal "<<std::endl;
			nextal.print();
			std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator nwal_id = node_indices.at(comp).find(nextal);
			assert(nwal_id != node_indices.at(comp).end());
			std::map<size_t , double>::iterator wei = weight.find(nwal_id->second);
			if(wei == weight.end()){
				std::map<size_t,size_t>::iterator pre = cur_pre.find(nwal_id->second);
				assert(pre == cur_pre.end());
				pre = cur_pre.find(id->second);
				double m1,m2;
				model.cost_function(nextal, m1, m2, ref_acc,read_acc);
				double this_weight = m1 + pre_wei->second;
				weight.insert(std::make_pair(nwal_id->second,this_weight));
				cur_pre.insert(std::make_pair(nwal_id->second, id->second));
				al_distance.insert(std::make_pair(this_weight,nwal_id->second));
			}else{
				std::multimap<double, size_t>::iterator test = al_distance.find(wei->second);
				if(test == al_distance.end()){
					al_distance.insert(std::make_pair(wei->second,nwal_id->second));						
				}else{ //TODO check if the second on is not the the same index and add it!

				}
				std::map<size_t,size_t>::iterator cur = cur_pre.find(nwal_id->second);
				if(cur == cur_pre.end()) std::cout << nwal_id->second << std::endl;
				assert(cur != cur_pre.end());
				double m1,m2;
				model.cost_function(nextal, m1, m2, ref_acc,read_acc);
				double this_weight = m1 + pre_wei->second;
				std::cout << "this weight "<< this_weight << " pre-weight "<< wei->second <<std::endl;
				if(this_weight < wei->second){
					size_t WEIGHT = wei->second;
					std::multimap<double, size_t>::iterator test = al_distance.find(wei->second);
					assert(test != al_distance.end());
					std::pair<std::multimap<double, size_t>::iterator, std::multimap<double, size_t>::iterator> findnode = al_distance.equal_range(wei->second);
					bool EXISTS = false;
					for(std::multimap<double, size_t>::iterator it2 = findnode.first; it2 != findnode.second; it2++){
						if(it2->second == nwal_id->second){
							al_distance.erase(it2);
							al_distance.insert(std::make_pair(this_weight,nwal_id->second));
							break;
						}
					}
					cur->second = id->second;
					wei->second = this_weight;
				}else{
					std::cout << "weight is bigger! "<< this_weight<<std::endl;
				}
			}		

		}	
		assert(al_distance.size()>0);
		std::map<size_t, const pw_alignment>::iterator ID = indices_nodes.at(comp).find(al_distance.begin()->second);
		std::cout<< "the smallest is: " << al_distance.begin()->first << " "<< al_distance.begin()->second <<std::endl;
		ID->second.print();
		assert(ID != indices_nodes.at(comp).end());
		this_al = ID->second;
		size_t l2,r2;
		this_al.get_lr2(l2,r2);
		size_t ind = this_al.getreference1();
		if(ind != data.numSequences()+1){
			std::string tempstr = data.get_seq_name(ind);
			std::stringstream str;
			str<< tempstr;
			size_t temp;
			str>> temp;
			if(this_al.getbegin1()<this_al.getend1()){
				startnode = temp;
			}else if(this_al.getbegin1()> this_al.getend1()){
				startnode = -1 * temp;
			}else{
				assert(this_al.getbegin1()== this_al.getend1());
				unsigned int ref1 = this_al.getreference1();
				size_t b1 = this_al.getbegin1();
				size_t e1 = this_al.getend1() ; 
				std::cout << " from al " << this_al.get_al_ref1() << " from seq "<< data.extract_seq_part(ref1, b1 , e1) <<std::endl;
				std::string seq_from_al;
				for(size_t i =0; i < this_al.get_al_ref1().length();i++){
					if(this_al.get_al_ref1().at(i) != '-'){
						seq_from_al += this_al.get_al_ref1().at(i);
						break;
					}
				}
				if(seq_from_al == data.extract_seq_part(ref1, b1 , e1)){
					startnode = temp;
				}else{
					startnode = -1 * temp;
				}

			}
			if(r2 != data.getSequence(this_al.getreference2()).length()){
				read_from = r2 +1;
			}else{
				//DO nothing
			}
		}else{
			//TODO
		}
		al_distance.erase(al_distance.begin());
	}
	if(LASTAL == false){
		std::cout << "add to the end : "<<std::endl;
		this_al.print();
		add_edge_to_end(comp,this_al);

	}
}

template<typename T>
double als_components<T>::get_cost(const pw_alignment & p, size_t & acc1, size_t & acc2){
	const std::map<std::string, std::vector<double>  > al_cost = model.get_al_cost(acc1, acc2);
}
template<typename T>
void als_components<T>::looking_for_last_al(size_t & comp, const pw_alignment & p, size_t & last_right){
	std::vector<std::vector<pw_alignment> >all_last_als;
	size_t type =4;
	size_t l1,l2,r1,r2;
	p.get_lr1(l1,r1);
	p.get_lr2(l2,r2);
	unsigned int ref2 = p.getreference2();
	unsigned int ref1 = p.getreference1();
	size_t readacc = data.accNumber(ref2);
	size_t refacc = data.accNumber(ref1);
	std::string seq_from_ref;
	size_t ref_from, ref_to;
	size_t length = data.getSequence(ref2).length();
	std::cout<< "length of seq "<< length<<std::endl;
/*	if(r1 < data.getSequence(ref1).length()-1){
		ref_from = r1+1;
		ref_to = data.getSequence(ref1).length()-1;
		seq_from_ref = data.extract_seq_part(ref1,ref_from,ref_to);
		if(p.getbegin1()>p.getend1()){//Get the reverse if begin1 > end 1
			std::cout<<"looking for last als begin > end "<<std::endl;

			std::string temp_seq_out;
			get_reverse_complement(seq_from_ref, temp_seq_out);
			seq_from_ref = temp_seq_out;
		}
	}*/
	if(l1 != 0 && p.getbegin1() > p.getend1()){
		ref_to = 0;
		ref_from = l1-1;
		std::string temp_seq_in = data.extract_seq_part(ref1,ref_to, ref_from);
		std::string temp_seq_out;
		get_reverse_complement(temp_seq_in, temp_seq_out);
		seq_from_ref = temp_seq_out;

	}
	if(r1 != data.getSequence(ref1).length()-1 && p.getbegin1()<p.getend1()){
		ref_from = r1+1;
		ref_to = data.getSequence(ref1).length()-1;
		seq_from_ref = data.extract_seq_part(ref1,ref_from,ref_to);
	}

	std::vector<std::vector<int> >refs;
	std::string seq;
	size_t new_al_begin;
	size_t new_al_end;
	if(length-r2<MAXGAP){//When it is close enough to the end of the read
		if(r2 != length-1){
			new_al_begin = r2+1;
			new_al_end = length-1;
			seq = data.extract_seq_part(ref2, new_al_begin, new_al_end);
		}
	}else{//When it is close enough to the end of a components of alignments on a read
		std::cout<<"if close to the end of component " << last_right << std::endl;
		assert( (last_right - r2) + MAXGAP/2  < MAXGAP);
		new_al_begin = r2+1;
		if(MAXGAP/2+last_right < length){
			new_al_end = MAXGAP/2 + last_right;
		}else{
			new_al_end = length-1;
		}	
		seq = data.extract_seq_part(ref2, new_al_begin, new_al_end);
	}
	
	if(seq.length() != 0){
		int current_node = std::stoi(data.get_seq_name(ref1)); 
		assert(current_node > 0);
		if(p.getbegin1() > p.getend1()){
			current_node = -1 * current_node;
		}
		std::cout << "seq length is "<< seq.length()<<std::endl;//Remained length from the read
		if(seq_from_ref.length()>=MAXGAP){
			//seq_from_ref.erase(seq_from_ref.begin()+MAXGAP,seq_from_ref.end());
			//assert(seq_from_ref.length()==MAXGAP);
			std::cout << "seq from ref length "<< seq_from_ref.length()<<std::endl;
			std::string sub = seq_from_ref.substr(0,MAXGAP);
			assert(sub.length() == MAXGAP);
			if(p.getbegin1() < p.getend1()){
				ref_from = r1+1;
				ref_to = r1+MAXGAP;
			}else{
				ref_from = l1-1;
				ref_to = l1-MAXGAP;
			}
			seq_from_ref = sub;
			std::string read_out;
			std::string ref_out;
			needleman<T> nl(data, model, seq, seq_from_ref);
		//XXX	simpleNW nl(model, seq, seq_from_ref);
			nl.run_needleman(readacc,refacc,type,read_out,ref_out);
			size_t count = 0;
			for(size_t i = 0; i < ref_out.length();i++){
				if(ref_out.at(i)!='-'){
					count++;
				}
			}
			if(p.getbegin1() < p.getend1()){
				ref_to = r1+ count;
			}else{
				ref_to = l1 - count;
			}

			pw_alignment p1(ref_out,read_out,ref_from,new_al_begin,ref_to,new_al_end,ref1,ref2);
			add_adjacencies(comp,p,p1,refacc,readacc);
			add_edge_to_end(comp,p1);

		}else{//Look for adjacent nodes
		/*	std::vector<std::vector<int> > refs;//The adjacent nodes from the current one on the ref graph(they are seq names not the ref ids!)
			size_t remainder = MAXGAP- seq_from_ref.length();
			std::cout << "length is "<< remainder << "ref1 is "<< ref1 <<std::endl;
			if(p.getbegin1()<p.getend1()){
				refs=rgraph.get_successor(ref1,true,remainder);
			}else{
				refs=rgraph.get_successor(ref1,false,remainder);
			}*/
		/*	std::set<std::vector<int> > paths = rgraph.get_paths();
		//	std::set<std::vector<int> >::iterator first = paths.begin();
		//	std::cout << "begin is "<< *first << std::endl;
		//	std::vector<int> first_path = *first;
			bool NOREF = false;*/
		//	assert(paths.size() != 0);
		//	if(paths.size()!= 0){
		//	std::set<int> nodes_on_paths = rgraph.get_subgraph_nodes();
			std::set<int> nodes_on_paths = rgraph.get_nodes();
			std::cout << "nodes size "<< nodes_on_paths.size() <<std::endl;
			if(nodes_on_paths.size() > 1){
				if(seq_from_ref.length() > 0){
					std::cout << "seq from node " << seq_from_ref.length() << std::endl;
					std::string read_out;
					std::string ref_out;
					needleman<T> nl(data, model, seq, seq_from_ref);
					type = 3;
					nl.run_needleman(readacc,refacc,type,read_out,ref_out);
					ref_from = r1 +1;
					ref_to = data.getSequence(p.getreference1()).length()-1;
					size_t count = 0;
					
					for(size_t i = 0; i < read_out.length();i++){
						if(read_out.at(i)!='-'){
							count++;
						}
					}
					if(p.getbegin1() > p.getend1()){
					//	ref_from = data.getSequence(p.getreference1()).length()-1;
					//	ref_to = r1 +1;
						ref_from = l1-1;
						ref_to = 0;
					}
					size_t al_end = new_al_begin + count -1 ;
					pw_alignment p1(ref_out,read_out,ref_from,new_al_begin,ref_to,al_end,ref1,ref2);
					add_adjacencies(comp,p,p1,refacc,readacc);
					new_al_begin = al_end +1;
					if(new_al_begin <= new_al_end){
						make_last_nw_als(comp, current_node ,p1 , nodes_on_paths ,seq,ref1,ref2,new_al_begin, new_al_end,readacc,refacc);
					}else{
						assert(new_al_begin = new_al_end +1);
						add_edge_to_end(comp,p1);
					}
				}else{
					std::cout << "seq from node " << seq_from_ref.length() << std::endl;					
					make_last_nw_als(comp, current_node ,p , nodes_on_paths ,seq,ref1,ref2,new_al_begin, new_al_end,readacc,refacc);
				}
			/*	std::set<std::vector<int> >::iterator first = paths.begin();
				std::cout << "begin is "<< *first << std::endl;
				for(std::set<std::vector<int> >::iterator it = paths.begin(); it != paths.end() ; it++){
		//	if(refs.size()!=0){
		//		std::cout<< "post nodes size "<< refs.size() <<std::endl;
		//		for(size_t i = 0 ; i < refs.size() ;i++){
		//			std::cout << refs.at(i)<<std::endl;
					std::string refin = seq_from_ref;
					std::vector<int> this_path = *it;
					std::cout << "first one " << this_path.at(0) << "remainder length of current node "<< seq_from_ref.length() <<std::endl;
					for(size_t j = 1; j < this_path.size(); j++){
						std::cout<< this_path.at(j)<<std::endl;		
						//Make an al between each member of refs and seq and pick the best
						size_t from = 0;
					//	if(refs.at(i).at(j)>0){
						if(this_path.at(j)>0){
							unsigned int ref_id = rgraph.get_refid(refacc,this_path.at(j));//gets node name, retruns node id
							size_t to = data.getSequence(ref_id).length()-1;
							refin += data.extract_seq_part(ref_id, from, to);
						}else{
							int temp = -1*this_path.at(j);
							unsigned int ref_id = rgraph.get_refid(refacc,temp);
							size_t to = data.getSequence(ref_id).length()-1;// Need to add  its reverse complement!
							std::string temp_seq_in = data.extract_seq_part(ref_id, from, to);
							std::string temp_seq_out;
							get_reverse_complement(temp_seq_in, temp_seq_out);
							refin += temp_seq_out;
						}
					}
					//keep only the first MAXGAP bases 
					if(refin.length() > MAXGAP){
						refin.erase(refin.begin()+MAXGAP,refin.end());
					}
					assert(refin.length()<=MAXGAP);
					std::cout<< "refin length "<< refin.length() << " seq length " << seq.length() <<std::endl;
					std::string read_out, ref_out;
					if(refin.length() != 0){
						needleman<T> nl(data, model, seq, refin);
						//XXX	simpleNW nl(model, seq, refin);
						nl.run_needleman(readacc,refacc,type,read_out,ref_out);
					}else{
						assert(paths.size() == 1);
						std::vector<int> first_path = *first;
						assert(first_path.size()==1);
						NOREF = true;
						add_edge_to_end(comp,p);
						break;
					}
					//Make all the als!
					size_t read_ref = ref2;
					size_t current_node = ref1;
					std::vector<pw_alignment> last_als;
					make_last_als(this_path,ref_out,read_out,new_al_begin,new_al_end,seq_from_ref,r1,refacc,read_ref, current_node,last_als); 
					all_last_als.push_back(last_als);
		//		}
				}*/
			/*	if(NOREF == false){
				//Choose the best set of als among those you created above (sum of their mod cost is the lowest)
				const std::vector<pw_alignment> best_als = find_the_best_als(all_last_als,refacc,readacc);
				add_adjacencies(comp,p,best_als.at(0),refacc,readacc);
				if(best_als.size()>2){
					std::cout <<"bigger than 2 " << best_als.size() <<std::endl;
					for(size_t i =0; i < best_als.size()-1;i++){
						add_adjacencies(comp,best_als.at(i),best_als.at(i+1),refacc,readacc);
					}
				}else if(best_als.size()==2){
					std::cout <<"equal 2 " << best_als.size() <<std::endl;
					add_adjacencies(comp,best_als.at(0),best_als.at(1),refacc,readacc);
				}
				add_edge_to_end(comp,best_als.at(best_als.size()-1));
				}*/
			}else{//Either only gap or the remainder from the ref node 
				if(seq_from_ref.length()!=0){
					std::string refin = seq_from_ref;
					std::string read_out, ref_out;
					needleman<T> nl(data, model, seq, refin);
				//XXX	simpleNW nl(model, seq, refin);
					nl.run_needleman(readacc,refacc,type,read_out,ref_out);
					size_t count = 0;
					for(size_t i = 0; i < ref_out.length();i++){
						if(ref_out.at(i)!='-'){
							count++;
						}
					}
					size_t ref_begin = r1+1;
					size_t ref_end = r1+count;
					if(p.getbegin1()<p.getend1()){ 
						if(count > 0){
							pw_alignment p1(ref_out,read_out,ref_begin,new_al_begin,ref_end,new_al_end,ref1,ref2);
							add_adjacencies(comp, p, p1,refacc,readacc);
							add_edge_to_end(comp,p1);
						}else{ //gap only
							unsigned int temp = data.numSequences()+1;//only gap on ref
							std::string str_ref(seq.length(),'-');
							pw_alignment p1(str_ref,seq,0,new_al_begin,0,new_al_end,temp,ref2);
							add_adjacencies(comp, p, p1,refacc,readacc);
							add_edge_to_end(comp,p1);
						}
					}else if(p.getbegin1()>p.getend1()){
						if(count > 0){
							ref_begin = l1 -1;
							ref_end = l1-count;
							assert(l1 >= count);
							pw_alignment p1(ref_out,read_out,ref_begin,new_al_begin,ref_end,new_al_end,ref1,ref2);
							add_adjacencies(comp, p, p1,refacc,readacc);
							add_edge_to_end(comp,p1);
						}else{
							unsigned int temp = data.numSequences()+1;//only gap on ref
							std::string str_ref(seq.length(),'-');
							pw_alignment p1(str_ref,seq,0,new_al_begin,0,new_al_end,temp,ref2);
							add_adjacencies(comp, p, p1,refacc,readacc);
							add_edge_to_end(comp,p1);
						}
					}else{
						assert(p.getbegin1()==p.getend1()); //TODO what if it is == ?
					}
				}else{//Only gap!
					unsigned int temp = data.numSequences()+1;//only gap on ref
					std::string str_ref(seq.length(),'-');
					pw_alignment p1(str_ref,seq,0,new_al_begin,0,new_al_end,temp,ref2);
					add_adjacencies(comp, p, p1,refacc,readacc);
					add_edge_to_end(comp,p1);
				}
			}

		}
	}else{
		//attach it directly to the end node
		assert(r2==length-1);
		std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator it=node_indices.at(comp).find(p);
		if(it!=node_indices.at(comp).end()){
			std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(it->second);
			if(adj == adjacencies.at(comp).end()){
				adjacencies.at(comp).insert(std::make_pair(it->second, std::set<size_t>()));
				adj = adjacencies.at(comp).find(it->second);
			}
			adj->second.insert(1);

		//	adjacencies.at(comp).insert(std::make_pair(it->second,1));
		}else{
			node_indices.at(comp).insert(std::make_pair(p,index));
			indices_nodes.at(comp).insert(std::make_pair(index,p));
			index++;
			it=node_indices.at(comp).find(p);
			std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(it->second);
			if(adj == adjacencies.at(comp).end()){
				adjacencies.at(comp).insert(std::make_pair(it->second, std::set<size_t>()));
				adj = adjacencies.at(comp).find(it->second);
			}
			adj->second.insert(1);

		//	adjacencies.at(comp).insert(std::make_pair(it->second,1));
		}
		weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(it->second,1),0));

	}
}
template<typename T>
void als_components<T>::make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout, size_t & left_on_current_node,size_t & refacc, size_t & readacc, size_t & read_id, std::vector<pw_alignment> & first_als, size_t & remainder_on_ref){
	size_t pos_on_ref = refout.size();
	size_t pos_on_read = left_on_current_node-1;
//	pw_alingnment pre_al = current_al;
	for(size_t i =0; i < nodes.size() ; i++){
		unsigned int ref_id;
		std::string onref;
		std::string onread;
		if(nodes.at(i)>0){
			ref_id = rgraph.get_refid(refacc,nodes.at(i));//gets node name, retruns node id
		}else{
			int temp = -1*nodes.at(i);
			ref_id = rgraph.get_refid(refacc,temp);
		}
		size_t node_length = data.getSequence(ref_id).length();
		std::cout << "on node "<< nodes.at(i) << "l "<< node_length << std::endl;
		if(i == 0 && remainder_on_ref > 0){
			//make the first al between the current node and read
			node_length = remainder_on_ref;
		}
		size_t counter =0;
		size_t read_counter = 0;
		for(int j = pos_on_ref-1; j >= 0; j--){
			std::cout << "j "<< j  << " "<< refout.at(j) <<std::endl;
			onref+=refout.at(j);
			onread+=readout.at(j);
			if(refout.at(j) != '-'){
				counter ++;
			}
			if(readout.at(j)!='-'){
				read_counter ++;
			}
			std::cout <<"c "<<counter << "l "<< node_length << std::endl;
			if(counter == node_length || j == 0){
				pos_on_ref = j;
				break;
			}
			if(counter == 0){//TODO gap only at the end of ref
				std::cout << "GAP ONLY AT THE BEGINNING!"<<std::endl;
			}
		}
		std::cout<<"HERE!!"<<std::endl;
		size_t from_on_read = 0;
		if(read_counter != 0){
			std::cout << pos_on_read << " " << read_counter<<std::endl;
			from_on_read = pos_on_read - read_counter+1;
			if(nodes.at(i)>0){
				pw_alignment p(onref,onread,0, from_on_read , node_length-1, pos_on_read,ref_id,read_id);
				first_als.push_back(p);

			}else{//If remainder >0 !! TODO
				pw_alignment p(onref,onread,node_length-1, from_on_read , 0, pos_on_read, ref_id,read_id);
				first_als.push_back(p);
			}
			if(from_on_read>0){
				pos_on_read = from_on_read - 1;
			}else{
				pos_on_read = from_on_read;
			}
			std::cout << "p on read" <<pos_on_read <<std::endl;
		}else{
			size_t to = data.getSequence(read_id).length();
			size_t ref_to = node_length -1;
			if(nodes.at(i)> 0){
				pw_alignment p(onref,onread,0, pos_on_read , ref_to, to, ref_id,read_id);
				first_als.push_back(p);
			}else{ //If remainder > 0!! TODO
				pw_alignment p(onref, onread ,ref_to, pos_on_read , 0, to, ref_id,read_id);
				first_als.push_back(p);
			}
		}	
		//add adjacencies p -> pre_p
	//	add_adjacencies(comp, p, pre_p);
		std::cout << "pos on ref " << pos_on_ref<<std::endl;
		if(pos_on_ref == 0) break;
	}
	std::cout << "pos on read "<< pos_on_read <<std::endl;
	assert(pos_on_read == 0);
	//add adjacencies 0 --> pre_al
//	add_edge_to_begin(comp, pre_p, refacc, readacc);
}
template<typename T>
void als_components<T>::make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & left_on_current_node,size_t & refacc, size_t & read_id, size_t & current_node, std::vector<pw_alignment> & first_als){//TODO !! It is wrong
	std::cout << "making first als "<<std::endl;
	std::string this_ref_part;
	std::string this_read_part;
	size_t ref_counter =0;
	size_t read_counter = 0;
	size_t position = 0;
	size_t current_position =0;
	size_t e2;
	e2 = last_end;
	size_t accu_ref_length = 0;
	if(seq_from_ref.length()>0){
		std::cout << "there is some left from the current node " << seq_from_ref.length() <<std::endl;
		for(size_t i = 0 ; i <refout.size() ;i++){
			if(refout.at(i)!='-'){
				ref_counter++;
			}
			if(readout.at(i)!='-'){
				read_counter++;
			}
			position++;
			this_ref_part += refout.at(i);
			this_read_part += readout.at(i);
			if(ref_counter == seq_from_ref.length()){
				break;
			}
		}
		//Make the al here: read readout and refout from size-position-1
		std::cout << "ref counter "<< ref_counter << " read counter "<<read_counter << std::endl;
		std::reverse(this_ref_part.begin(),this_ref_part.end());
		std::reverse(this_read_part.begin(),this_read_part.end());
		if(read_counter != 0){
			pw_alignment p(this_ref_part,this_read_part,seq_from_ref.length()-ref_counter,last_end-(read_counter-1), seq_from_ref.length()-1, last_end,current_node, read_id);
			std::cout << "al1" <<std::endl;
			p.print();
			assert(last_end-(read_counter-1)<=last_end);
			first_als.push_back(p);
			e2 = last_end-(read_counter-1)-1;
		}else{
			assert(read_counter == 0);
			size_t this_length = data.getSequence(read_id).length();
			pw_alignment p(this_ref_part,this_read_part,0,this_length, seq_from_ref.length()-1, this_length,current_node, read_id);
			std::cout << "al1" <<std::endl;
			p.print();
			first_als.push_back(p);
		}
		accu_ref_length += ref_counter;
		current_position = position;
		this_ref_part.clear();
		this_read_part.clear();
		read_counter = 0;
		ref_counter = 0;
	}

	for(size_t i = nodes.size()-1 ; i > 0 ;i--){
		if(refout.size() == current_position){
			break;
		}
		unsigned int ref_id;
		if(nodes.at(i)>0){
			ref_id = rgraph.get_refid(refacc,nodes.at(i));//gets node name, retruns node id
		}else{
			int temp = -1*nodes.at(i);
			ref_id = rgraph.get_refid(refacc,temp);
		}
		size_t current_node_length = data.getSequence(ref_id).length();
		std::cout << refout.size() << " " <<current_position <<std::endl;
		for(size_t j = current_position ; j <refout.size(); j++){
			if(refout.at(j)!='-'){
				ref_counter++;
			}
			if(readout.at(j)!='-'){
				read_counter++;
			}
			position++;
			this_ref_part += refout.at(j);
			this_read_part += readout.at(j);
			if(ref_counter == current_node_length){//When ref_counter is equal to the length of the current node
				break;
			}
		}
		assert(this_ref_part.length() != 0 || this_read_part.length() !=0);
		std::reverse(this_ref_part.begin(),this_ref_part.end());
		std::reverse(this_read_part.begin(),this_read_part.end());
		accu_ref_length += ref_counter;
		if(read_counter != 0){
			std::cout << current_node_length<< " " << ref_counter<<std::endl;
			pw_alignment p(this_ref_part,this_read_part,current_node_length-ref_counter,e2-(read_counter-1), current_node_length-1, e2,ref_id, read_id);
			std::cout << "read counter "<< read_counter << std::endl;
			std::cout << "al2" <<std::endl;
			p.print();
			assert(e2-(read_counter-1)<=e2);
			first_als.push_back(p);	
			e2 = e2-(read_counter-1)-1;
		}else{
			//only gap on read
			assert(read_counter == 0);//If this ref didnt map to any base on read.
			size_t this_length = data.getSequence(read_id).length();
			pw_alignment p(this_ref_part,this_read_part,current_node_length-ref_counter,e2, current_node_length-1, this_length,ref_id, read_id);//TODO Check if it is correct 
			first_als.push_back(p);
			std::cout << "al3" <<std::endl;
			p.print();
		}
		current_position = position;
		this_ref_part.clear();
		this_read_part.clear();
		read_counter = 0;
		ref_counter = 0;
	}
}
template<typename T>
void als_components<T>::make_last_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & right_on_current_node,size_t & refacc, size_t & read_id, size_t & current_node, std::vector<pw_alignment> & last_als){//TODO read counter == 0, and end position is reached condition!!
	std::cout << "making last als "<< refout.size() <<std::endl;
	std::string this_ref_part;
	std::string this_read_part;
	size_t ref_counter =0;
	size_t read_counter = 0;
//	size_t position = 0;
	size_t current_position =0;
	size_t b2;
	b2 = first_begin;
	size_t accu_ref_length = 0;
	size_t this_length = data.getSequence(read_id).length();
	if(seq_from_ref.length()>0){
		std::cout << "there is some left from the current node " << seq_from_ref.length() <<std::endl;
		bool onlyfirst = false;
		if(nodes.size()==0){
			std::cout<< "all of it is located on this node"<<std::endl;
			onlyfirst = true;
		}
		compute_samples(onlyfirst,seq_from_ref.length(),current_position, read_counter,ref_counter , readout, this_read_part, refout, this_ref_part);
	//	for(size_t i = 0 ; i <refout.size() ;i++){
	//		if(refout.at(i)!='-'){
	//			ref_counter++;
	//		}
	//		if(readout.at(i)!='-'){
	//			read_counter++;
	//		}
	//		position++;
	//		this_ref_part += refout.at(i);
	//		this_read_part += readout.at(i);
	//		if(ref_counter == seq_from_ref.length()){
	//			break;
	//		}
	//	}
		//Make the al here: read readout and refout from size-position-1
		std::cout<< right_on_current_node+ref_counter << " "<< data.getSequence(current_node).length() << " " <<ref_counter << " "<<read_counter<<std::endl;
	//	if(ref_counter == seq_from_ref.length()){//XXX It is correct only for + nodes not the - ones
	//		assert(right_on_current_node+ref_counter == data.getSequence(current_node).length()-1);
	//	}
	//	assert(right_on_current_node+ref_counter <= data.getSequence(current_node).length()-1);
		if(read_counter != 0){
			pw_alignment p(this_ref_part,this_read_part,right_on_current_node+1,first_begin,right_on_current_node+ref_counter, first_begin+(read_counter-1),current_node, read_id);
			last_als.push_back(p);
			std::cout << "al1" <<std::endl;
			p.print();
		}else{
			pw_alignment p(this_ref_part,this_read_part,right_on_current_node+1,first_begin-1,right_on_current_node+ref_counter,this_length, current_node, read_id);//TODO check if it is correct
			last_als.push_back(p);
			std::cout << "al1" <<std::endl;
			p.print();
		}
		accu_ref_length += ref_counter;
		current_position ++;
		this_ref_part.clear();
		this_read_part.clear();
		b2 = first_begin+read_counter;
		read_counter = 0;
		ref_counter = 0;
	}
	if(refout.size() > current_position){
	std::cout << "cur pos "<< current_position <<std::endl;
	std::cout<<"nodes size "<<nodes.size()<<std::endl;
	for(size_t i = 1 ; i < nodes.size() ;i++){
		std::cout<< "i "<< i <<std::endl;
		unsigned int ref_id;
		if(nodes.at(i)>0){
			ref_id = rgraph.get_refid(refacc,nodes.at(i));//gets node name, retruns node id
		}else{
			int temp = -1*nodes.at(i);
			ref_id = rgraph.get_refid(refacc,temp);
		}
		std::cout << "ref id " << ref_id <<std::endl;
		size_t current_node_length = data.getSequence(ref_id).length();
		bool tillend = false;
		if(i==nodes.size()-1){//When the ref sample ends with gaps/solution is not to break it when it is the last node in nodes
			tillend = true;
		}
		compute_samples(tillend, current_node_length, current_position, read_counter, ref_counter, readout, this_read_part, refout, this_ref_part);//TODO what if refcounter stays zero?
		std::cout << "size of ref out " <<refout.size() << " current pos " <<current_position <<std::endl;
	//	for(size_t j = current_position ; j <refout.size(); j++){
	//		if(refout.at(j)!='-'){
	//			ref_counter++;
	//		}
	//		if(readout.at(j)!='-'){
	//			read_counter++;
	//		}
	//		position++;
	//		this_ref_part += refout.at(j);
	//		this_read_part += readout.at(j);
	//		if(ref_counter == current_node_length){//When ref_counter is equal to the length of the current node
	//			break;
	//		}
	//	}
		accu_ref_length += ref_counter;
		std::cout<< "ref counter "<< ref_counter << std::endl;
		std::cout << current_node_length-1 << " "<< ref_counter-1 <<std::endl;
	//	if(accu_ref_length <=MAXGAP){
	//		assert(current_node_length-1==ref_counter-1);
	//	}
		if(read_counter != 0){//TODO what if ref_counter == 0???
		//	pw_alignment p(this_ref_part,this_read_part,0,b2, current_node_length-1, b2+read_counter-1 ,ref_id, read_id);
			pw_alignment p(this_ref_part,this_read_part,0,b2, ref_counter-1, b2+read_counter-1 ,ref_id, read_id);
			last_als.push_back(p);
			std::cout << "read counter "<< read_counter << std::endl;
			std::cout << "al2" <<std::endl;
			std::cout << "ref "<< this_ref_part << " read "<< this_read_part <<std::endl;
			p.print();
		}else{//If read_counter == 0
			pw_alignment p(this_ref_part,this_read_part,0,b2-1, ref_counter-1, this_length ,ref_id, read_id);//TODO check if it is correct
			last_als.push_back(p);
			std::cout << "read counter1 "<< read_counter << std::endl;
			std::cout << "al2" <<std::endl;
			p.print();
		}

	//	}else{
	//		pw_alignment p(this_ref_part,this_read_part,0,b2, ref_counter-1, b2+read_counter-1,ref_id, read_id);
	//		assert(b2 <= b2+read_counter-1);
	//		last_als.push_back(p);
	//		std::cout << "al3" <<std::endl;
	//		p.print();

	//	}
		current_position ++;
		this_ref_part.clear();
		this_read_part.clear();
		b2 = b2+read_counter;

		read_counter = 0;
		ref_counter = 0;
		std::cout<< "current position is "<< current_position << " refout length "<< refout.length() << std::endl;
		if(current_position == refout.length()){
			break;
		}

	}
}else{//TODO this shouldn't happen, check!
		assert(refout.length() == current_position);
		std::cout << "refout size "  << refout.length() << " current pos "<< current_position << std::endl;
		//Add it to end! TODO
	}

	std::cout<< b2 << " "<< last_end <<std::endl;
	assert(b2-1==last_end);
}
template<typename T>
const std::vector<pw_alignment> als_components<T>::find_the_best_als(std::vector<std::vector<pw_alignment> >& all_als,size_t & refacc, size_t & readacc)const{
	double min;
	size_t id = 0;
	for(size_t i = 0; i < all_als.size(); i++){
		std::cout<<"i "<< i << std::endl;
		double sum_of_cost = 0.0;
		for(size_t j =0; j < all_als.at(i).size();j++){
			pw_alignment p = all_als.at(i).at(j);
			p.print();
			double m1,m2;
			model.cost_function(p,m1,m2,refacc,readacc);//TODO both created al should be changed!
			sum_of_cost += m1;
		}
		if(i==0){
			min = sum_of_cost;
			id = i;
		}
		if(sum_of_cost< min){
			min = sum_of_cost;
			id = i;
		}
	}
	std::cout << "id of best al "<< id <<std::endl;
	return all_als.at(id);

}
template<typename T>
void als_components<T>::add_adjacencies(size_t & comp, const pw_alignment & p1 , const pw_alignment & p2, size_t & ref_acc, size_t & read_acc){
	p1.print();
	p2.print();
	size_t p1_index;
	std::map<const pw_alignment,size_t,compare_pw_alignment >::const_iterator it=node_indices.at(comp).find(p1);
	if(it == node_indices.at(comp).end()){
		std::cout << "here!"<<std::endl;
		node_indices.at(comp).insert(std::make_pair(p1,index));
		indices_nodes.at(comp).insert(std::make_pair(index,p1));
		index++;
		it=node_indices.at(comp).find(p1);
	}
	double c1,c2,m1,m2;
//	if(p2.getreference1() != data.numSequences()+1 && p2.getreference2() != data.numSequences()+2 && p2.getbegin2()!=p2.getend2() && p2.getbegin1() != p2.getbegin1()){
//		model.cost_function(p2,c1,c2,m1,m2);
//	}else if(p2.getreference1() == data.numSequences()+1&& p2.getbegin2()!=p2.getend2() && p2.getbegin1() != p2.getbegin1()){
//		model.cost_function(p2, m1, m2, ref_acc, read_acc);
//	}else if(p2.getreference2() == data.numSequences()+2&& p2.getbegin2()!=p2.getend2() && p2.getbegin1() != p2.getbegin1()){
//		model.cost_function(p2, m1, m2, ref_acc, read_acc);
//	}else{
//		assert(p2.getbegin2()==p2.getend2() || p2.getbegin1() == p2.getbegin1());
	model.cost_function(p2, m1, m2, ref_acc, read_acc);
	std::cout << "m1 "<< m1 << " m2 "<< m2 <<std::endl;
//		std::cout << "m1 for al of length one "<< m1 <<std::endl;
//	}
	p1_index = it->second;
	size_t seqsize = data.getSequence(p2.getreference2()).length();
	std::map<const pw_alignment,size_t,compare_pw_alignment>::const_iterator it1=node_indices.at(comp).find(p2);
	if((it1 == node_indices.at(comp).end())){//Just removed : ((p2.getbegin2()==p2.getend2())&&(p2.getend2()==seqsize))
		std::cout << "p1 index "<< p1_index << " p2 index "<<index<<std::endl;
		if(index == 241 || index == 247){
			std::cout<< "p2 is "<<std::endl;
			p2.print();
			std::cout << p2.get_al_ref1()<<std::endl;
			std::cout << p2.get_al_ref2() << std::endl;
		}
		node_indices.at(comp).insert(std::make_pair(p2,index));
		indices_nodes.at(comp).insert(std::make_pair(index,p2));
		size_t p2_index = index;
		index++;
		std::set<int> temp;
		temp.insert(p2_index);
		std::map<size_t , std::set<size_t> >::iterator it = adjacencies.at(comp).find(p1_index);
		if(it == adjacencies.at(comp).end()){
			adjacencies.at(comp).insert(std::make_pair(p1_index, std::set<size_t>()));
			it = adjacencies.at(comp).find(p1_index);
		}
		it->second.insert(p2_index);
		it = adjacencies.at(comp).find(p2_index);
		if(it == adjacencies.at(comp).end()){
			adjacencies.at(comp).insert(std::make_pair(p2_index, std::set<size_t>()));
		}
	//	adjacencies.at(comp).insert(std::make_pair(p1_index, temp));
		std::pair<size_t, size_t> this_pair;
		this_pair = std::make_pair(p1_index,p2_index);
		std::map<std::pair<size_t,size_t>,double>::iterator weight = weight_of_edges.at(comp).find(this_pair);
		if(weight == weight_of_edges.at(comp).end()){
			weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(p1_index, p2_index),m1));
		}
	}else{
		std::cout<<"else! "<<std::endl;
		size_t p2_index = it1->second;
		std::cout << p1_index << " to "<< p2_index<<std::endl;
		std::map<size_t , std::set<size_t> >::iterator it = adjacencies.at(comp).find(p1_index);
		if(it == adjacencies.at(comp).end()){
			adjacencies.at(comp).insert(std::make_pair(p1_index, std::set<size_t>()));
			it = adjacencies.at(comp).find(p1_index);
		}
		it->second.insert(p2_index);
		it = adjacencies.at(comp).find(p2_index);
		if(it == adjacencies.at(comp).end()){
			adjacencies.at(comp).insert(std::make_pair(p2_index, std::set<size_t>()));
		}

	//	weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(p1_index, p2_index),m1)); 
		std::pair<size_t, size_t> this_pair;
		this_pair = std::make_pair(p1_index,p2_index);
		std::map<std::pair<size_t,size_t>,double>::iterator weight = weight_of_edges.at(comp).find(this_pair);
		if(weight == weight_of_edges.at(comp).end()){//TODO
			weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(p1_index, p2_index),m1));
		}

	}
//	size_t p2_index = it1->second;
//	adjacencies.at(comp).insert(std::make_pair(p1_index, p2_index));
//	weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(p1_index, p2_index),m1));
	std::cout<<"done!"<<std::endl;
}
template<typename T>
void als_components<T>::add_edge_to_begin(size_t & comp , const pw_alignment & p, size_t & read_acc, size_t & ref_acc){
	std::map<const pw_alignment, size_t,compare_pw_alignment >::const_iterator it=node_indices.at(comp).find(p);
	assert(it != node_indices.at(comp).end());
	std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(0);
	if(adj == adjacencies.at(comp).end()){
		adjacencies.at(comp).insert(std::make_pair(0, std::set<size_t>()));
		adj = adjacencies.at(comp).find(0);
	}
	adj->second.insert(it->second);

//	adjacencies.at(comp).insert(std::make_pair(0, it->second));
	double c1,c2,m1,m2;
//	if(p.getreference1() != data.numSequences()+1 && p.getreference2() != data.numSequences()+2){
//		model.cost_function(p,c1,c2,m1,m2);
//	}else if(p.getreference1() == data.numSequences()+1){
//		model.cost_function(p, m1, m2, ref_acc, read_acc);
//	}else if(p.getreference2() == data.numSequences()+2){
		model.cost_function(p, m1, m2, ref_acc, read_acc);
//	}
	weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(0,it->second),m1));
}
template<typename T>
void als_components<T>::add_edge_to_end(size_t & comp , const pw_alignment & p){
	p.print();
	std::map<const pw_alignment ,size_t,compare_pw_alignment>::const_iterator it=node_indices.at(comp).find(p);
	if(it == node_indices.at(comp).end()){
		node_indices.at(comp).insert(std::make_pair(p,index));
		indices_nodes.at(comp).insert(std::make_pair(index,p));
		index++;
		it = node_indices.at(comp).find(p);
	}
	std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(it->second);
	if(adj == adjacencies.at(comp).end()){
		adjacencies.at(comp).insert(std::make_pair(it->second, std::set<size_t>()));
		adj = adjacencies.at(comp).find(it->second);
	}
	adj->second.insert(1);

//	adjacencies.at(comp).insert(std::make_pair(it->second,1));
	weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(it->second,1),0));
}
template<typename T>
const pw_alignment & als_components<T>::get_alignment(size_t & comp, size_t & nodeid)const{
	typename std::multimap<size_t , const pw_alignment>::const_iterator it = indices_nodes.at(comp).find(nodeid);
	assert(it != indices_nodes.at(comp).end());
	return it->second;
}
template<typename T>
void als_components<T>::add_expensive_edges(size_t & comp,size_t & refacc, size_t & readacc){
	for(std::map<size_t, const pw_alignment>::const_iterator it =indices_nodes.at(comp).begin(); it != indices_nodes.at(comp).end() ;it++){
		size_t node = it->first;
		std::cout << "this node "<< node <<std::endl;
		const pw_alignment p = it ->second;
		size_t left1,left2,right1,right2;
		p.get_lr1(left1,right1);
		p.get_lr2(left2,right2);
		p.print();
	//	std::cout<< data.getSequence(p.getreference2()).length()<<std::endl;
	/*	std::pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator> this_node = adjacencies.at(comp).equal_range(node);
		std::set<size_t> adjs;
		for(std::multimap<size_t,size_t>::iterator it1 = this_node.first ; it1 != this_node.second ;it1++){
			adjs.insert(it1->second);
			std::cout << it1->second << " ";
		}
		std::cout << " "<<std::endl;
		std::cout << "adjs size "<< adjs.size() <<std::endl;*/
		std::set<size_t> adjs;
		std::map<size_t , std::set<size_t> >::iterator this_node = adjacencies.at(comp).find(node);
		if(this_node == adjacencies.at(comp).end()){
			adjacencies.at(comp).insert(std::make_pair( node, std::set<size_t>()));
			this_node = adjacencies.at(comp).find(node);
		}
		adjs = this_node->second;
		for(std::map<size_t, const pw_alignment>::const_iterator it1 =indices_nodes.at(comp).begin(); it1 != indices_nodes.at(comp).end() ;it1++){
			const pw_alignment al = it1->second;
			size_t l1,l2, r1,r2;
			al.get_lr1(l1,r1);
			al.get_lr2(l2,r2);
			//if al has no overlap with p && is not in its adjs and started after the current one, then we add an edge between them
		//	if((left1 < r1 && right1 > l1) || (left2 < r2 && right2 > l2)){
		//		std::cout<<"there is overlap"<<std::endl;
		//	}
			std::set<size_t>::iterator adj = adjs.find(it1->first);
			std::cout << "HERE"<<std::endl;
		//	if((l1 > right1) && (left2> r2 || l2 > right2 ) && (adj == adjs.end())){
			if((p.getreference1()!= al.getreference1() || (l1 > right1|| r1 < left1)) && (l2 > right2) && (adj == adjs.end())){
				std::cout<< "expensive edge! "<<std::endl;
				al.print();
				double m1,m2;
				model.cost_function(al, m1, m2, refacc, readacc);
				std::cout << "m1 "<< m1 << " m2 "<< m2 <<std::endl;

			//	adjacencies.at(comp).insert(std::make_pair(node, it1->first));
				this_node->second.insert(it1->first);
				//These weight of these nodes are added to the weight map
				std::pair<size_t , size_t> this_edge(node,it1->first);
				double weight = (MAXGAP/2 * 2.5)+1000; 
				weight += (l2-right2-1)*2.5;
				weight += m1;
				weight_of_edges.at(comp).insert(std::make_pair(this_edge,weight));
			}
			
		}

	}
	
}
template<typename T>
void als_components<T>::add_to_maf(const pw_alignment & al, std::ofstream & output, bool & firstal){//TODO fix only gap alignments TODO thelengt of alignments are wrong!!
	size_t ref1,ref2;
	ref1 = al.getreference1();
	ref2 = al.getreference2();
	size_t l1,r1,l2,r2;
	al.get_lr2(l2,r2);
	al.get_lr1(l1,r1);
	size_t acc1,acc2;
	std::cout << "l1 "<< l1 << " r1 "<< r1 <<std::endl;
	acc2 = data.accNumber(ref2);
	std::string accname2 = data.get_acc(acc2);
	std::string seqname2 = data.get_seq_name(ref2);
	std::stringstream longname2;
	longname2 << accname2 << ':' << seqname2;
	assert(al.getbegin2()<= al.getend2());
//	std::cout<< "numseq + 1 "<< data.numSequences()+1 <<std::endl;
	assert(ref1 != data.numSequences()+1 || (l2 != data.getSequence(ref2).length() && r2 != data.getSequence(ref2).length()));
	if(ref1 != data.numSequences()+1){
	//	if(al.get_al_ref1().length()==0){
	//		al.print();
	//	}
		std::cout << "in add to maf "<< al.getbegin1() << " "<<al.getend1() << std::endl;
		assert(al.get_al_ref1().length()!=0 && al.get_al_ref2().length()!=0);
		output << "a " << std::endl; //TODO Maybe add a score by using the mod cost or gain of the al
		acc1 = data.accNumber(ref1);
		std::string accname1 = data.get_acc(acc1);
		std::string seqname1 = data.get_seq_name(ref1);
		std::stringstream longname1;
		longname1 << accname1 << ':' << seqname1;
		if(al.getbegin1()<al.getend1()){
			output << "s " << longname1.str() << " " << al.getbegin1() << " " << al.getend1()-al.getbegin1()+1 << " + " << data.getSequence(ref1).length() << " " << al.get_al_ref1() <<std::endl;
			previous_right1 = r1;
			previous_left1 = l1;
			previous_ref = ref1;
			forward = true;
		}else if(al.getbegin1() > al.getend1()){
			output << "s " << longname1.str() << " " << al.getend1() << " " << al.getbegin1()-al.getend1()+1 << " - " << data.getSequence(ref1).length() << " " << al.get_al_ref1() <<std::endl;
			previous_right1 = r1;
			previous_left1 = l1;
			previous_ref = ref1;
			forward = false;
		}else{
			assert(al.getbegin1()==al.getend1());
			std::cout<< "ref of length one " << al.get_al_ref1().size() <<std::endl;
			al.print();
		//	assert(al.get_al_ref1().size()==1);//It is possible that begin == end but the length is not zero. it happens when there are gaps around that single position. So this assertion can be wrong.
			char this_base;
			std::cout << al.get_al_ref1()<<std::endl;
			for(size_t i =0; i< al.get_al_ref1().size();i++){
				if(al.get_al_ref1().at(i)!='-'){
					this_base = al.get_al_ref1().at(i);
					break;
				}
			}
			std::cout <<"length " << data.getSequence(ref1).length() << "this base "<< this_base <<std::endl;
			if(this_base==data.getSequence(ref1).at(l1)){//To see if the sequence of length one is forward or backward
				output << "s " << longname1.str() << " " << al.getbegin1() << " " << al.getend1()-al.getbegin1()+1 << " + " << data.getSequence(ref1).length() << " " << al.get_al_ref1() <<std::endl;
				forward = true;
			}else{
			//	assert(al.get_al_ref1().at(0)== dnastring::complement(data.getSequence(ref1).at(l1)));
				std::cout << "dna: "<< data.getSequence(ref1).str()<<std::endl;
				assert(this_base == dnastring::complement(data.getSequence(ref1).at(l1)));

				output << "s " << longname1.str() << " " << al.getbegin1() << " " << al.getend1()-al.getbegin1()+1 << " - " << data.getSequence(ref1).length() << " " << al.get_al_ref1() <<std::endl;
				forward = false;
			}
			previous_right1 = r1;
			previous_left1 = l1;
			previous_ref = ref1;


			
		}
		if(al.getbegin2()< al.getend2() && r2 != data.getSequence(ref2).length()){
			output << "s " << longname2.str() << " " << al.getbegin2() << " " << al.getend2()-al.getbegin2()+1 << " + " << data.getSequence(ref2).length() << " " << al.get_al_ref2() <<std::endl;
			previous_right2 = r2;
		}else{
			if(r2 == data.getSequence(ref2).length()){//Alignments with only gap on the read
				std::cout<< "only gap on read " << std::endl;
				output << "s " << longname2.str() << " " << previous_right2+1 << " " << 0 << " + " << data.getSequence(ref2).length() << " " << al.get_al_ref2() <<std::endl;
			}else{
				assert(al.getbegin2()==al.getend2());
				output << "s " << longname2.str() << " " << al.getbegin2() << " " << al.getend2()-al.getbegin2()+1 << " + " << data.getSequence(ref2).length() << " " << al.get_al_ref2() <<std::endl;
				previous_right2 = r2;
			}
			
		}
	}else if(firstal == false){
		//Add to the previous alignment //TODO it is the case only when it is not the first al. In that case it has to be added to the next al. Also direction should be taken into account!! XXX It shouldnt be any gap only al in the maf maf!
		std::cout<<"only gap on ref"<<std::endl;
		acc1 = data.accNumber(previous_ref);
		std::string accname1 = data.get_acc(acc1);
		std::string seqname1 = data.get_seq_name(previous_ref);
		std::stringstream longname1;
		longname1 << accname1 << ':' << seqname1;
		if(forward == true){
			output << "s " << longname1.str() << " " << previous_right1 << " " << 0 << " + " << data.getSequence(previous_ref).length() << " " << al.get_al_ref1() <<std::endl;
		}else{
			output << "s " << longname1.str() << " " << previous_left1 << " " << 0 << " - " << data.getSequence(previous_ref).length() << " " << al.get_al_ref1() <<std::endl;

		}
		assert(r2 != data.getSequence(ref2).length());
	//	assert(al.getbegin2()< al.getend2());

		output << "s " << longname2.str() << " " << al.getbegin2() << " " << al.getend2()-al.getbegin2()+1 << " + " << data.getSequence(ref2).length() << " " << al.get_al_ref2() <<std::endl;
		previous_right2 = r2;	
	}else{
		assert(firstal == true);
		std::cout << "gap only at the beginning! "<<std::endl;
	}

}
template<typename T>
void als_components<T>::make_components(size_t & comp, std::vector<std::map<std::pair<size_t,size_t>,double> > & components, std::vector<std::set<size_t> > & comp_nodes, bool & FULLCOMP, size_t & no_comp){
	std::map<std::pair<size_t,size_t>, bool> visited;
	std::map<size_t, std::set<size_t> > temp_adj;
	for(std::map<size_t , std::set<size_t> >::iterator it =adjacencies.at(comp).begin() ; it != adjacencies.at(comp).end(); it++){
		for(std::set<size_t>::iterator this_adj = it->second.begin() ; this_adj != it->second.end() ;this_adj++){
			std::pair<size_t, size_t> this_pair(it->first, *this_adj);
	    		visited.insert(std::make_pair(this_pair,false));
		}
		temp_adj.insert(std::make_pair(it->first,it->second));
 	}
	while(temp_adj.size() != 0){
		std::map<std::pair<size_t,size_t>,double> component;
		std::set<size_t> c_nodes;
		std::map<size_t, std::set<size_t> >::iterator it = temp_adj.begin();
		std::list<size_t> queue;
		queue.push_back(it->first);
		assert(queue.size()==1);
		while(!queue.empty()){
			size_t this_node = queue.front();
			std::cout << "this node "<< this_node <<std::endl;
			queue.pop_front();
			it = temp_adj.find(this_node);
			if(it != temp_adj.end()){
				c_nodes.insert(this_node);
				for(std::set<size_t>::iterator this_adj = it->second.begin() ; this_adj != it->second.end() ;this_adj++){
					std::pair<size_t,size_t> this_pair(this_node, *this_adj);
					std::map<std::pair<size_t,size_t> , bool>::iterator vis = visited.find(this_pair);
					assert(vis != visited.end());
					if(vis->second == false){
						vis->second = true;
						queue.push_back(*this_adj);
						std::map<std::pair<size_t,size_t>,double>::iterator weight = weight_of_edges.at(comp).find(this_pair);
						assert(weight != weight_of_edges.at(comp).end());
						component.insert(std::make_pair(this_pair,weight->second));
						c_nodes.insert(*this_adj);
					}
				}
				temp_adj.erase(it);
			}		}
		components.push_back(component);
		std::set<size_t>::iterator it1 = c_nodes.find(0);
		std::set<size_t>::iterator it2 = c_nodes.find(1);
		if(it1 != c_nodes.end() && it2 != c_nodes.end()){
			no_comp = comp_nodes.size();	
			FULLCOMP = true;
		}else{
		//	FULLCOMP = false;
		}	
		comp_nodes.push_back(c_nodes);
	}
}


template<typename T>
void als_components<T>::filter_als(size_t & comp){
	size_t first_left, last_right;
	std::multiset<pw_alignment,sort_pw_alignment_by_right>::reverse_iterator revit = ordered_als_by_right.at(comp).rbegin(); //finds the last alignment of the component
	pw_alignment last_p = *revit;
	size_t lastr, lastl;
	last_p.get_lr2(lastl,lastr);
	last_right = lastr;
	size_t number = 0;
	for(std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it = alignments.at(comp).begin() ; it != alignments.at(comp).end(); it++){//For each alignment all of its successors are found
		pw_alignment p = *it;
		size_t l,r;
		p.get_lr2(l,r);
		if(it==alignments.at(comp).begin()){
			first_left = l;
		}
		if(l <= number*MAXGAP/2 + first_left && l >= (number-1)*MAXGAP/2 + first_left ){


		}else{


		}
	}


}
#endif
