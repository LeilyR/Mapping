#include "ref_graph.hpp"

void ref_graph::read_dot_file(std::string & refgraph, std::string & refacc){
	std::ifstream dotin(refgraph.c_str());
	if(!dotin){
		std::cerr << "Error : Cannot open " << refgraph.c_str() << std::endl;
		exit(1);
	}
	else{
		std::string line;
		while(getline(dotin,line)){
			if(line[0] == '/'){
				continue;
			}
			else{
				std::vector<std::string> parts;
				strsep(line, " ", parts);
				std::string from = parts.at(0);
				std::string to;
				if(parts.size()>2){
					to = parts.at(2);
					std::cout<<"from " << from << " to " << to << std::endl;
				}
				std::cout << "ref acc "<< refacc << std::endl;
				std::string dir1(from.end()-1,from.end());
				std::string name1(from.begin(),from.end()-1);
				std::string longname1 = refacc + ":"+name1;
				std::cout<<"dir1 "<< dir1 << " name1 "<< name1 << " "<< longname1<<std::endl;
				std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(longname1);
				std::cout << "seq id is "<< findseq1->second << " " << findseq1->first << std::endl;
				assert(findseq1!=longname2seqidx.end());
				if(findseq1==longname2seqidx.end()) {
					std::cerr << "Error: unknown sequence in dot File: " << longname1 << std::endl;
					exit(1);
				}
				size_t seqid1 = findseq1->second;
				if(parts.size() == 2){
					add_adjacencies(dir1,name1);
				}
				else{
					assert(to != "");
					std::string dir2(to.end()-1,to.end());
					std::string name2(to.begin(),to.end()-1);
					std::string longname2 = refacc + ":"+name2;
					std::cout<<"dir2 "<< dir2 << " name2 "<< name2 << " "<< longname2<<std::endl;

					add_adjacencies(dir1,dir2,name1,name2);
					add_adjacencies(dir2,name2);
					std::map<std::string, size_t>::iterator findseq2 = longname2seqidx.find(longname2);
					if(findseq2==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << longname2 << std::endl;
						exit(1);
					}
				}
			}
		}
	}
}

void ref_graph::add_adjacencies(std::string & dir1 , std::string & dir2 , std::string & name1 , std::string & name2){
//	std::cout << "n1 "<<name1 << " n2 "<<name2 <<std::endl;
	int node1 = std::stoi(name1);//TODO it shoud be changed to size_t!!!!
	int node2 = std::stoi(name2);
	assert(node1 > 0 && node2 > 0);
	if((dir1 == "+") && (dir2 == "+")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(node1);
		std::map<int, std::set<int> >::iterator dy_adj = dynamic_adjacencies.find(node1);

		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(node1, std::set<int>()));
			assert(dy_adj==dynamic_adjacencies.end());
			dynamic_adjacencies.insert(std::make_pair(node1, std::set<int>()));

			adj = adjacencies.find(node1);
			dy_adj = dynamic_adjacencies.find(node1);
		}
		adj->second.insert(node2);
		dy_adj->second.insert(node2);
		add_predecessor(node2,node1);
	}
	if((dir1 == "-") && (dir2 == "+")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(-1*node1);
		std::map<int, std::set<int> >::iterator dy_adj = dynamic_adjacencies.find(-1*node1);

		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(-1*node1, std::set<int>()));
			assert(dy_adj==dynamic_adjacencies.end());
			dynamic_adjacencies.insert(std::make_pair(-1*node1, std::set<int>()));

			adj = adjacencies.find(-1*node1);
			dy_adj = dynamic_adjacencies.find(-1*node1);

		}
		adj->second.insert(node2);
		dy_adj->second.insert(node2);

		int rev_node1 = -1*node1;
		add_predecessor(node2, rev_node1);
	}
	if((dir1 == "+") && (dir2 == "-")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(node1);
		std::map<int, std::set<int> >::iterator dy_adj = dynamic_adjacencies.find(node1);

		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(node1, std::set<int>()));
			assert(dy_adj==dynamic_adjacencies.end());
			dynamic_adjacencies.insert(std::make_pair(node1, std::set<int>()));

			adj = adjacencies.find(node1);
			dy_adj = dynamic_adjacencies.find(node1);

		}
		adj->second.insert(-1*node2);
		dy_adj->second.insert(-1*node2);
		int rev_node2 = -1*node2;

		add_predecessor(rev_node2, node1);

	}
	if((dir1 == "-") && (dir2 == "-")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(-1*node1);
		std::map<int, std::set<int> >::iterator dy_adj = dynamic_adjacencies.find(-1*node1);

		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(-1*node1, std::set<int>()));
			assert(dy_adj==dynamic_adjacencies.end());
			dynamic_adjacencies.insert(std::make_pair(-1*node1, std::set<int>()));

			adj = adjacencies.find(-1*node1);
			dy_adj = dynamic_adjacencies.find(-1*node1);

		}
		adj->second.insert(-1*node2);
		dy_adj->second.insert(-1*node2);

		int rev_node1 = -1*node1;
		int rev_node2 = -1*node2;
	//	std::cout << "rev: " << rev_node2 << " "<< rev_node1 <<std::endl;
		add_predecessor(rev_node2, rev_node1);

	}
}

void ref_graph::add_adjacencies(std::string & dir , std::string & name){
	int cent_id = std::stoi(name);
	assert(cent_id > 0);
	if(dir == "+"){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(cent_id);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(cent_id, std::set<int>()));
			dynamic_adjacencies.insert(std::make_pair(cent_id, std::set<int>()));

		}
	}else{
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(-1*cent_id);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair((-1*cent_id), std::set<int>()));
			dynamic_adjacencies.insert(std::make_pair((-1*cent_id), std::set<int>()));

		}
	}
}

void ref_graph::add_predecessor(int & this_node , int & pre_node){
	std::map<int, std::set<int> >::iterator pre = predecessors.find(this_node);
	if(pre == predecessors.end()){
		predecessors.insert(std::make_pair(this_node, std::set<int>()));
		pre = predecessors.find(this_node);
	}
	pre->second.insert(pre_node);
//	std::cout << "this node "<< this_node << " pre node " << pre_node<<std::endl;
	
}
const std::map<int , std::set<int> > & ref_graph::get_adjacencies()const{
	return adjacencies;
}


std::set<std::vector<int> > ref_graph::get_predecessor(unsigned int & ref_id, bool dir , size_t & left_on_ref, size_t & length_on_read){//int is ref number if negetive means reverse complement should be taken in to account
	std::cout << "ref name " << data.get_seq_name(ref_id) <<std::endl;
	int seqname = std::stoi(data.get_seq_name(ref_id));
	assert(seqname > 0);
//	size_t refacc = data.accNumber(ref_id);
	std::string acc_name = data.get_acc(data.accNumber(ref_id));
	if(dir == false){
		seqname = -1*seqname;
	}
	std::cout<<"seq name "<<seqname << " seq id "<< ref_id <<std::endl;
//	std::set<int>visited;
	std::set<std::vector<int> > pre_paths;
	std::map<int , std::set<int> >::const_iterator it = predecessors.find(seqname);
	if(it != predecessors.end()){
		std::cout<< "is not the name! "<<std::endl;
	//	size_t current_length = 0;
	//	look_for_predecessor(seqname,length,current_length,acc_name,visited,this_pre_nodes,all_pre_nodes);
		find_predecessors_bfs(seqname,acc_name ,left_on_ref, pre_paths, length_on_read);
	}
//	for(std::set<int>::iterator it = visited.begin(); it!= visited.end() ; it++){
//		std::cout<< " seen: "<< *it<<std::endl;
//	}
	return pre_paths;
}
void ref_graph::find_predecessors_bfs(int & startnode, std::string & refacc, size_t & left_on_ref, std::set<std::vector<int> > & pre_paths, size_t & length_on_read){
	std::set<int> seen;
	size_t max = MAXGAP;
	if(2*length_on_read < MAXGAP){
		max = 2*length_on_read;
	}
	std::map<int, std::set<int> >::const_iterator pre = predecessors.find(startnode);//It is possible that a node has no adjacent node
	if(pre != predecessors.end()){
		int begin_at = startnode;
		std::string seq_name = seqname(begin_at);
		size_t seqlength = seq_length(seq_name, refacc);
		std::cout << " "<<seqlength << " "<<left_on_ref<< std::endl;
		std::map<std::pair<int,int>, bool> visited;
		for(std::map<int, std::set<int> >::iterator it = predecessors.begin() ; it != predecessors.end() ;it++){
			for(std::set<int>::iterator this_pre = it->second.begin() ; this_pre != it->second.end(); this_pre ++){
      				visited.insert(std::make_pair(std::make_pair(it->first,*this_pre),false)); //Was( it->first ,false )
			}
 		}
		std::map<int, std::vector<std::vector<int> > > all_paths;
		all_paths.insert(std::make_pair(begin_at,std::vector<std::vector<int> >()));
		std::list<int> queue;
		std::map<int, std::vector<std::pair<std::vector<int> , size_t > > > length;
		std::vector<std::pair<std::vector<int>,size_t> > this_pair;
		std::vector<int> temp;
		this_pair.push_back(std::make_pair(temp, left_on_ref));
		length.insert(std::make_pair(begin_at, this_pair));
		size_t this_length = left_on_ref;
		if(left_on_ref >= max){
			std::vector<int> temp;
			temp.push_back(begin_at);
			pre_paths.insert(temp);
		//	pre_nodes_on_paths.insert(startnode);
			return ;
		}
		queue.push_back(begin_at);
		while(!queue.empty()){
			begin_at = queue.front();
			std::cout << "on node "<< begin_at<< std::endl;
			queue.pop_front();
			std::map<int, std::set<int> >::iterator it = predecessors.find(begin_at);
			if(it != predecessors.end()){
				std::set<int> this_predecessors = it->second;
				std::map<int, std::vector<std::vector<int> > >::iterator this_path = all_paths.find(begin_at);
				if(this_path !=all_paths.end()){
					std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator pre_len = length.find(begin_at);
					std::cout << pre_len->second << std::endl;
					assert(pre_len != length.end());
					bool EXPANDED = false;
					for(std::set<int>::iterator pred = this_predecessors.begin() ; pred != this_predecessors.end() ; pred++){
						std::cout << *pred << std::endl;
						std::map<std::pair<int,int>, bool>::iterator vis =visited.find(std::make_pair(begin_at,*pred));
						std::cout << "bool "<< int(vis->second)<<std::endl;
						if(vis->second == false){
							vis->second = true; 
							std::set<int>::iterator s = seen.find(*pred);
							if(s == seen.end()){
								queue.push_back(*pred);
							}
							seen.insert(*pred);
							int adjacent = *pred;
							std::string adj_name = seqname(adjacent);
							size_t adj_length = seq_length(adj_name, refacc);
							std::cout << adjacent << "adj length" << adj_length <<std::endl;
							std::vector<std::vector<int> >subpath = this_path->second;
							std::cout << subpath.size() <<std::endl;
							for(size_t i =0;i < subpath.size();i++){
								subpath.at(i).push_back(begin_at);
								std::cout<< "Here: "<<subpath.at(i)<<std::endl;
							}
							if(subpath.size()==0){
								std::vector<int> temp;
								temp.push_back(begin_at);
								subpath.push_back(temp);
							}
						//	std::cout << subpath <<std::endl;
							std::map<int, std::vector<std::vector<int> > >::iterator from_adj = all_paths.find(adjacent);
							if(from_adj == all_paths.end()){
								all_paths.insert(std::make_pair(adjacent,std::vector<std::vector<int> >()));
								from_adj = all_paths.find(adjacent);
							}
							assert(from_adj != all_paths.end());
							for(size_t i =0;i < subpath.size();i++){
								from_adj->second.push_back(subpath.at(i));
							}
							std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator len = length.find(adjacent);
							std::vector<std::pair<std::vector<int>,size_t> > this_pair;
							if(len == length.end()){
								length.insert(std::make_pair(adjacent, this_pair));
								len = length.find(adjacent);
							}
							assert(len !=length.end());
							this_pair = pre_len->second;
							for(size_t i = 0; i < this_pair.size() ;i++){
								std::pair<std::vector<int>, size_t> p = this_pair.at(i);
								p.first.push_back(begin_at);
								p.second+=adj_length;
								std::cout << p.second << std::endl;
								if(p.second< max){
									len->second.push_back(p);
								}
								if(p.second>=max){
									std::cout << "MAX reached! "<<std::endl;
									std::cout << p.first<<std::endl;
									std::vector<int> add_to_path = p.first;
								//	for(size_t n = 0; n < add_to_path.size() ; n++){
								//		pre_nodes_on_paths.insert(add_to_path.at(n));
								//	}
									add_to_path.push_back(adjacent);
								//	pre_nodes_on_paths.insert(adjacent);
									pre_paths.insert(add_to_path);
									bool SEEN = false;
									for(size_t j =0; j < from_adj->second.size() ;j++){
										if(from_adj->second.size()== 6){
											std::cout<< from_adj->second.at(j)<<std::endl;
										}
										if(from_adj->second.at(j)==p.first){
											from_adj->second.erase(from_adj->second.begin()+j);
											SEEN = true;
											break;
										}
									}
									if(SEEN!=true) std::cout << from_adj->second.size() <<std::endl;
									assert(SEEN == true);
								}
							}
							if(from_adj->second.size()==0){
								all_paths.erase(from_adj);
								//Remove it from the queue
							}
							EXPANDED = true;
						}
					}
					if(EXPANDED == true){
						all_paths.erase(this_path);
					}
				}
			}
		}
		std::cout << "HERE!!" << pre_paths.size() <<std::endl;
		for(std::set<vector<int> >::iterator it = pre_paths.begin() ; it != pre_paths.end() ; it++){
			std::cout << *it << std::endl;
			std::cout << " "<<std::endl;
		}
		for(std::map<int,std::vector<std::vector<int> > >::iterator it = all_paths.begin() ; it != all_paths.end() ; it++){
			std::vector<std::vector<int> >path = it->second;
			for(size_t i =0; i < path.size() ; i++){
				path.at(i).push_back(it->first);
				std::set<vector<int> >::iterator it1 = pre_paths.find(path.at(i));
				if(it1 == pre_paths.end()){
					pre_paths.insert(path.at(i));
				//	for(size_t j =0; j < path.at(i).size() ; j++){
				//		pre_nodes_on_paths.insert(path.at(i).at(j));
				//	}
				}
				std::cout << "this path "<< path.at(i)<<std::endl;
			}
		}
	}


}

const void ref_graph::look_for_predecessor(int & node , size_t & length , size_t & current_length, std::string & acc_name, std::set<int> & visited, std::vector<int> & this_pre_nodes, std::vector<std::vector<int> > & all_pre_nodes)const{
	std::map<int , std::set<int> >::const_iterator it = predecessors.find(node);//TODO Go bfs
	if(it != predecessors.end()){
		std::set<int> pre_nodes = it->second;
		if(it->second.size()==0){
			//Keep till here and current_legnth== 0;
			all_pre_nodes.push_back(this_pre_nodes);
			this_pre_nodes.clear();
			current_length = 0;
		}
		for(std::set<int>::iterator it1 = pre_nodes.begin() ; it1 != pre_nodes.end() ; it1++){
			int name = *it1;
			std::cout << "this pre: "<< name <<std::endl;
			std::set<int>::iterator seen= visited.find(name);
			if(seen==visited.end()){
				visited.insert(name);
				this_pre_nodes.push_back(name);
				int name1 = name;
				if(name < 0){
					name1 = -1*name;
				}
				std::stringstream str;
				str<<acc_name<<":"<<name1;
				std::cout<< str.str()<<std::endl;
				std::map<std::string, size_t>::const_iterator longname = longname2seqidx.find(str.str());
				assert(longname != longname2seqidx.end());
				size_t this_ref = longname->second;
				std::cout << "this length "<< data.getSequence(this_ref).length() <<std::endl;
				current_length += data.getSequence(this_ref).length();
				std::cout << "cur length "<< current_length << "length "<< length<<std::endl;
				if(current_length < length){
					look_for_predecessor(name,length,current_length,acc_name,visited,this_pre_nodes,all_pre_nodes);
				}else{
					//add to a container!
					std::cout << "pushed back! "<<std::endl;
					std::cout<< "this pre nodes size "<< this_pre_nodes.size() << std::endl;
					all_pre_nodes.push_back(this_pre_nodes);
					this_pre_nodes.clear();
					current_length = 0;
				}
			}
		}
	}else{
		//add to a container!
		std::cout << "pushed back1! "<<std::endl;
		std::cout<< "this pre nodes size "<< this_pre_nodes.size() << std::endl;
		all_pre_nodes.push_back(this_pre_nodes);
		this_pre_nodes.clear();
		current_length = 0;
	}
}
void ref_graph::read_gfa_for_adj(std::string & gfafile){
	std::ifstream in(gfafile.c_str());
	std::string line;
	while(getline(in,line)){
		if(line[0] == 'H'){//Skip the header
			continue;
		}else{
			assert(line[0] == 'S' || line[0]=='L' || line[0] == 'P');
		
			if(line[0]=='S'){
				std::vector<std::string> node;
				strsep(line, "\t" , node);
				std::string name = node.at(1);
				std::string content = node.at(2);
				nodes_content.insert(std::make_pair(std::stoi(name), content));
				assert(std::stoi(name) > 0);
			}
			
			if(line[0]== 'L'){
				std::vector<std::string> nodes;
				strsep(line, "\t" , nodes);
				std::string name1 , name2, dir1, dir2;
				name1 = nodes.at(1);
				name2 = nodes.at(3);
				dir1 = nodes.at(2);
				dir2 = nodes.at(4);
			//	std::cout << "L "<< line << std::endl;
			//	std::cout << "name1 "<< name1 << " name2 " << name2 << " dir1 " << dir1 << " dir2 " << dir2<<std::endl;
				add_adjacencies(dir1,dir2, name1, name2);
				add_adjacencies(dir2,name2);
				if(dir2 == "+"){
					dir2 = "-";
					add_adjacencies(dir2,name2);	
					if(dir1 == "+") dir1 = "-";
				//	else dir1 = "+";
					add_adjacencies(dir2,dir1, name2, name1);
					add_adjacencies(dir1,name1);	

				
				}else{
					dir2 = "+";
					add_adjacencies(dir2,name2);	
				}

			}
			if(line[0]=='P'){
				std::vector<std::string> nodes;
				strsep(line, "\t" , nodes);
				std::vector<std::string> node;
				strsep(nodes.at(2), "," , node);
				std::string seq = node.at(0).substr(0, node.at(0).length()-1);
				begin_nodes.insert(std::stoi(seq));
			}

		}
	}
}
void ref_graph::read_gfa_for_seq(std::string & gfafile , std::ofstream & fasta , std::ofstream & paths){
	std::ifstream in(gfafile.c_str());
	std::string line;
	while(getline(in,line)){
		if(line[0] == 'H'){//Skip the header
			continue;
		}else{
			assert(line[0] == 'S' || line[0]=='L' || line[0] == 'P');
		
			if(line[0]== 'S'){
				std::vector<std::string> node;
				strsep(line, "\t" , node);
				fasta<<">"<<node.at(1)<<std::endl;
				fasta<<node.at(2)<<std::endl;
			}
			if(line[0]=='P'){
				std::vector<std::string> node;
				strsep(line, "\t" , node);
				paths<<node.at(1)<<std::endl;
				paths<<node.at(2)<< std::endl;
			}
		}
	}


}
const std::vector<std::vector<int> > ref_graph::get_successor(unsigned int & ref_id, bool dir , size_t & length)const{//int is ref number if negetive means reverse complement should be taken in to account
	int seqname = std::stoi(data.get_seq_name(ref_id));
	std::string acc_name = data.get_acc(data.accNumber(ref_id));
	if(dir == false){
		seqname = -1*seqname;
	}
	std::cout<<"seq name "<<seqname << " seq id "<< ref_id <<std::endl;
	std::set<int>visited;
	std::vector<int> this_adjacencies;
	std::vector<std::vector<int> > all_adjacencies;
	std::map<int , std::set<int> >::const_iterator it = adjacencies.find(seqname);
	if(it != adjacencies.end() && it->second.size()!=0){
		std::cout<< "has " << it->second.size() <<" adjacent nodes! "<<std::endl;
		size_t current_length = 0;
		look_for_successor(seqname,length,current_length,acc_name,visited,this_adjacencies,all_adjacencies);
	}
	for(std::set<int>::iterator it = visited.begin(); it!= visited.end() ; it++){
		std::cout<< " seen: "<< *it<<std::endl;
	}
	std::cout << "all_adjacents size "<< all_adjacencies.size() << std::endl;
	return all_adjacencies;
}
const void ref_graph::look_for_successor(int & node, size_t & length, size_t & current_length, std::string & acc_name, std::set<int> & visited, std::vector<int> & this_adjacencies, std::vector<std::vector<int> > & all_adjacencies)const{
	std::map<int , std::set<int> >::const_iterator it = adjacencies.find(node);
	if(it != adjacencies.end()){
		std::set<int> post_nodes = it->second;
		if(it->second.size()==0){
			//Keep till here and current_legnth== 0;
			all_adjacencies.push_back(this_adjacencies);
			this_adjacencies.clear();
			current_length = 0;
		}
		for(std::set<int>::iterator it1 = post_nodes.begin() ; it1 != post_nodes.end() ; it1++){
			int name = *it1;
			std::cout << "this post: "<< name <<std::endl;
			std::set<int>::iterator seen= visited.find(name);
			if(seen==visited.end()){
				visited.insert(name);
				this_adjacencies.push_back(name);
				int name1 = name;
				if(name < 0){
					name1 = -1*name;
				}
				std::stringstream str;
				str<<acc_name<<":"<<name1;
				std::cout<< str.str()<<std::endl;
				std::map<std::string, size_t>::const_iterator longname = longname2seqidx.find(str.str());
				assert(longname != longname2seqidx.end());
				size_t this_ref = longname->second;
				std::cout << "this length "<< data.getSequence(this_ref).length() <<std::endl;
				current_length += data.getSequence(this_ref).length();
				std::cout << "cur length "<< current_length << "length "<< length<<std::endl;
				if(current_length < length){
					look_for_successor(name,length,current_length,acc_name,visited,this_adjacencies,all_adjacencies);
				}else{
					//add to a container!
					std::cout << "pushed back! "<<std::endl;
					std::cout<< "this pre nodes size "<< this_adjacencies.size() << std::endl;
					all_adjacencies.push_back(this_adjacencies);
					this_adjacencies.clear();
					current_length = 0;
				}
			}
		}
	}else{
		//add to a container!
		std::cout << "pushed back1! "<<std::endl;
		std::cout<< "this pre nodes size "<< this_adjacencies.size() << std::endl;
		all_adjacencies.push_back(this_adjacencies);
		this_adjacencies.clear();
		current_length = 0;
	}


}
const unsigned int ref_graph::get_refid(size_t & refacc, int & seqname)const{
	//Get accession name:
	std::string accname = data.get_acc(refacc);
	int name = seqname;
	//convert int to string:
	if( name < 0) name = -1*name;
	std::stringstream ss;
	ss << name;
	std::string longname = accname +":"+ss.str();
	std::cout << "long name "<< longname<<std::endl;
//	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
	std::map<std::string, size_t>::const_iterator findseq = longname2seqidx.find(longname);
	assert(findseq != longname2seqidx.end());
	unsigned int refid = findseq->second;
	return refid;
}
size_t ref_graph::seq_length(std::string & seq_name, std::string & accname){
	std::string longname = accname +":"+seq_name;
//	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
/*	std::cout << "long name "<<std::endl;
	for(std::map<std::string, size_t>::iterator tfindseq = longname2seqidx.begin(); tfindseq != longname2seqidx.end(); tfindseq++){
		std::cout << tfindseq->first << " " << tfindseq->second<<std::endl;
	}*/
	std::map<std::string, size_t>::iterator findseq = longname2seqidx.find(longname);
	assert(findseq != longname2seqidx.end());
	return data.get_seq_size(findseq->second);
}
std::string ref_graph::seqname(int & node){
	std::string temp;
	std::stringstream ss;
	if(node< 0){
		ss << (-1*node);
		temp = ss.str();
	}else{
		ss << node;
		temp = ss.str();
	}
	return temp;
}

void ref_graph::deep_first_search(int & startnode, std::string & refacc, size_t & right_on_ref){
	paths.clear();
	nodes_on_paths.clear();
	path_length.clear();
	parent_length.clear();
//	std::map<int , std::vector<std::pair<std::vector<int>,size_t> > > parent_length;
	std::map<int, std::set<int> >::iterator adj = adjacencies.find(startnode);//It is possible that a node has no adjacent node
	if(adj != adjacencies.end()){
		std::string seq_name = seqname(startnode);
		size_t seqlength = seq_length(seq_name, refacc);
		size_t remainder = seqlength-(right_on_ref+1);
		std::cout<< remainder << " "<<seqlength << " "<<right_on_ref<< std::endl;
		int accu_length = remainder;
		std::map<int, bool> visited;
		for(std::map<int, std::set<int> >::iterator it = adjacencies.begin() ; it != adjacencies.end() ;it++){
      			visited.insert(std::make_pair(it->first,false));
 		}
		//Recurssive dfs: 
		std::vector<int> apath;
		std::vector<int> temp;
		std::set<std::pair<std::vector<int>,size_t> > this_pair;
		this_pair.insert(std::make_pair(temp,remainder));
		
		parent_length.insert(std::make_pair(startnode, this_pair));
		if(remainder> MAXGAP){
			apath.push_back(startnode);
			nodes_on_paths.insert(startnode);
			paths.insert(apath);
		}else{
    			look_for_neighbors(startnode, visited, refacc);
		}
	}else{
		std::vector<int> apath;		
		apath.push_back(startnode);
		nodes_on_paths.insert(startnode);
		paths.insert(apath);
	}
	for(std::map<int , std::set<std::pair<std::vector<int>,size_t> > >::iterator it = parent_length.begin() ; it != parent_length.end() ; it++){
		std::cout<< "on node " << it->first << " " << it->second.size() << std::endl;
		for(std::set<std::pair<std::vector<int>,size_t> >::iterator it1 = it->second.begin() ; it1!=it->second.end() ; it1++){
			std::pair<std::vector<int>,size_t> this_pair = *it1;
			std::cout<< this_pair.first<<std::endl;
			std::cout<< "with length "<< this_pair.second <<std::endl;
		}
	}
}
void ref_graph::look_for_neighbors(int & node, std::map<int,bool> & visited , std::string & refacc){
	std::cout << "this node is "<< node << std::endl;
	nodes_on_paths.insert(node);
	std::map<int,bool>::iterator it = visited.find(node);
	assert(it != visited.end());
	std::string seq_name = seqname(node);
	size_t seqlength = seq_length(seq_name, refacc);
/*	std::cout << "its length is "<< seqlength<<std::endl;
	if(accu_length<=MAXGAP){
		apath.push_back(node);	//TODO wrong
		nodes_on_paths.insert(node);
	}
	std::cout << "apath size "<< apath.size() <<std::endl;
	for(size_t i =0; i < apath.size();i++){
		std::cout << apath.at(i)<< " ";
	}
	std::cout<<" "<<std::endl;
	//Look for the adjacent vertices: 
	std::map<int, std::set<int> >::iterator it1 = adjacencies.find(node);
	assert(it1 != adjacencies.end());
	std::cout<< "it has "<< it1->second.size() << " adjacents" <<std::endl;
	if(it1->second.size()==0){
		if(accu_length <= MAXGAP){
			paths.insert(apath);
			apath.pop_back();
			nodes_on_paths.erase(node);
		}
		accu_length -=seqlength;
	}
	std::cout << "accu length is "<< accu_length << std::endl; 
	std::cout<< "all the adjs " << it1->second.size() <<std::endl;
	for (std::set<int>::iterator adj = it1->second.begin(); adj != it1->second.end(); adj++){
		std::cout << *adj << " ";
	}
	std::cout<<std::endl;*/
	assert(it->second == false);
//	if(it->second == false){
		it->second = true;
		std::map<int, std::set<int> >::iterator it1 = adjacencies.find(node);
		assert(it1 != adjacencies.end());
		for(std::set<int>::iterator adj = it1->second.begin(); adj != it1->second.end(); adj++){
			int adjacent = *adj;
			std::string adj_name = seqname(adjacent);
			size_t adjlength = seq_length(adj_name, refacc);
			std::map<int , std::set<std::pair<std::vector<int>,size_t> > >::iterator len = parent_length.find(node);
			if(len != parent_length.end()){
				bool LESSTHANMAX = false;
				std::cout << "len second size "<< len->second.size() <<std::endl;
				for(std::set<std::pair<std::vector<int>,size_t> >::iterator len_it = len->second.begin() ; len_it != len->second.end() ;len_it++){
					std::pair<std::vector<int>,size_t> this_pair = *len_it;
					std::cout<<"this pair first  " << std::endl;
					std::cout << this_pair.first << std::endl;
					std::map<int,bool>::iterator it2 = visited.find(*adj);
					if(this_pair.second < MAXGAP){
						size_t length = this_pair.second + adjlength;
						std::cout << "len "<< length <<std::endl;
						std::map<int , std::set<std::pair<std::vector<int>,size_t> > >::iterator adj_len = parent_length.find(adjacent);
						std::vector<int> cur_path = this_pair.first;
						cur_path.push_back(node);
						std::pair<std::vector<int>, size_t> cur_pair(cur_path,length);
						if( adj_len == parent_length.end()){
							std::set<std::pair<std::vector<int>, size_t> > temp_set;
							parent_length.insert(std::make_pair(adjacent,temp_set));
							 adj_len = parent_length.find(adjacent);
						}
						assert(adj_len != parent_length.end());
						std::cout << "adjacent "<< adjacent << std::endl;
						std::cout<< "cur pair "<<cur_pair <<std::endl;
						adj_len->second.insert(cur_pair);
						if(it2->second == false){
						//	LESSTHANMAX = true;
							int current_node = *adj;
          						look_for_neighbors(current_node, visited,refacc);
						}
						//add it to parent length, keep the current ones
					}else{
						//Do nothing
					}
				}
	/*	accu_length += adjlength;
		
		if(adj == it1->second.begin()){
			parent_length = seqlength;
			std::cout << "parent length "<< parent_length <<std::endl;
		}
		std::cout<< "adj is "<< *adj << "adj length "<< adjlength << std::endl;
		std::set<int>::iterator adj1 =it1->second.end();
		std::map<int,bool>::iterator it2 = visited.find(*adj);
		assert(it2 != visited.end());
		if (it2->second == false && accu_length < MAXGAP){
			std::cout << "smaller"<<std::endl;*/
			//	if(LESSTHANMAX==true){
			//		int current_node = *adj;
          		//		look_for_neighbors(current_node, visited,refacc,parent_length);
			//	}
			}
/*		else if(it2->second == false && accu_length >= MAXGAP){ //XXX > only and the = belonged to the previous condition
			std::cout << "here! " <<std::endl;
			int current_node = *adj;
			it2->second = true; 
			apath.push_back(current_node);
			nodes_on_paths.insert(current_node);
			paths.insert(apath);
			std::cout<< "apath size: "<< apath.size() << " " << accu_length <<std::endl;
			for(size_t i =0; i < apath.size();i++){
				std::cout << apath.at(i)<< " ";
			}
			std::cout<<" "<<std::endl;
			accu_length -=adjlength;
			apath.pop_back();
			std::set<int>::iterator adj2 =it1->second.end();
			if(--adj2 == adj){
				std::cout<<"end of adjs"<<std::endl;
				std::cout << "accu before "<< accu_length <<std::endl;
				accu_length -= parent_length;
				std::cout <<"parent length "<<parent_length<<std::endl;
				std::cout<< "accu after "<< accu_length << std::endl;
			}

		}
		else if(it2->second == true){//TODO check this!
			std::cout << "if it is true! "<<std::endl;
			paths.insert(apath);//TODO Should change it, because if a node has more than one visited adjacencies the same path is added more than once!!
			accu_length -= adjlength;
		}
		else{
			std::cout<< "shouldn't happen!! "<<std::endl;
		}
	}
	apath.pop_back();
	std::cout<<"paths size is "<<paths.size()<<std::endl;*/
		}
//	}
}
void ref_graph::make_sub_graph(int & startnode, std::string & refacc, size_t & right_on_ref_node){
	sub_graph.clear();
	sub_graph_nodes.clear();
	pre_nodes_on_subgraph.clear();
	int first_node = startnode;
	std::map<int,std::set<int> >::iterator adj = adjacencies.find(startnode);
	if(adj != adjacencies.end()){
		std::string seq_name = seqname(startnode);
		size_t seqlength = seq_length(seq_name, refacc);
		size_t remainder = seqlength-(right_on_ref_node+1);
		if(remainder >= MAXGAP){
			sub_graph.insert(std::make_pair(startnode,std::set<int>()));
			sub_graph_nodes.insert(startnode);
			return;
		}
		std::map<int, bool> visited;
		std::cout << "adj size " << adjacencies.size() << " dy adj size "<< dynamic_adjacencies.size()<<std::endl;
		for(std::map<int, std::set<int> >::iterator it = adjacencies.begin() ; it != adjacencies.end() ;it++){
      			visited.insert(std::make_pair(it->first,false)); 
			for(std::set<int>::iterator this_adj = it->second.begin() ; this_adj != it->second.end(); this_adj ++){
				visited.insert(std::make_pair(*this_adj,false)); 
			}
 		}
		std::map<int, std::vector<std::pair<std::vector<int>,size_t> > > length;
		std::vector<std::pair<std::vector<int>,size_t> > this_pair;
		std::vector<int> temp;
		this_pair.push_back(std::make_pair(temp, remainder));
		length.insert(std::make_pair(startnode, this_pair));

		std::list<int> queue;
		queue.push_back(startnode);
		std::map<int,bool>::iterator vis=visited.find(startnode);
	//	std::cout << "start node "<< startnode <<std::endl;
		assert(vis != visited.end());
		vis->second = true;
		sub_graph_nodes.insert(startnode);
		while(!queue.empty()){
			startnode = queue.front();
		//	std::cout << "on node "<< startnode<< std::endl;
			queue.pop_front();
			adj = adjacencies.find(startnode);
			std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator pre_len = length.find(startnode);
			assert(pre_len != length.end());
			std::vector<std::pair<std::vector<int> , size_t > > parent_length = pre_len ->second;
			for(std::set<int>::iterator this_adj = adj->second.begin() ; this_adj != adj->second.end() ; this_adj++){
				vis=visited.find(*this_adj);
				assert(vis != visited.end());
				if(vis->second == false){
					vis->second = true;
					int adjacent_node = *this_adj;
					std::string adj_name = seqname(adjacent_node);
					size_t adj_length = seq_length(adj_name, refacc);
					std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator len = length.find(*this_adj);
					if(len == length.end()){
						std::vector<std::pair<std::vector<int>,size_t> > this_pair;
						std::vector<int> temp;
						this_pair.push_back(std::make_pair(temp, 0));
						length.insert(std::make_pair(*this_adj, this_pair));
						len = length.find(*this_adj);
					}
					for(size_t i = 0; i < parent_length.size() ; i++){
						std::pair<std::vector<int> , size_t > this_pair = parent_length.at(i);
						std::map<int, std::set<int> >::iterator it= sub_graph.find(startnode);
						if(it == sub_graph.end()){
							sub_graph.insert(std::make_pair(startnode,std::set<int>()));
							it = sub_graph.find(startnode);
						}
						it->second.insert(adjacent_node);
						sub_graph_nodes.insert(adjacent_node);

						std::map<int, std::set<int> >::iterator it1= pre_nodes_on_subgraph.find(adjacent_node);
						if(it1 == pre_nodes_on_subgraph.end()){
							pre_nodes_on_subgraph.insert(std::make_pair(adjacent_node,std::set<int>()));
							it1 = pre_nodes_on_subgraph.find(adjacent_node);
						}
						it1->second.insert(startnode);
						if(this_pair.second + adj_length >= MAXGAP){
							//Do nothing
							
						}else{
							queue.push_back(*this_adj);
							//Add it the length:
							std::vector<int> temp = this_pair.first;
							temp.push_back(startnode);
							len->second.push_back(std::make_pair(temp,this_pair.second+adj_length));
						}

					}
				//	queue.push_back(*this_adj);
				}else{
					int adjacent_node = *this_adj;
					std::map<int, std::set<int> >::iterator it= sub_graph.find(startnode);
					if(it == sub_graph.end()){
						sub_graph.insert(std::make_pair(startnode,std::set<int>()));
						it = sub_graph.find(startnode);
					}
					it->second.insert(adjacent_node);

					std::map<int, std::set<int> >::iterator it1= pre_nodes_on_subgraph.find(adjacent_node);
					if(it1 == pre_nodes_on_subgraph.end()){
						pre_nodes_on_subgraph.insert(std::make_pair(adjacent_node,std::set<int>()));
						it1 = pre_nodes_on_subgraph.find(adjacent_node);
					}
					it1->second.insert(startnode);

					std::set<int>::iterator it2 = sub_graph_nodes.find(adjacent_node);
				//	std::cout << "adj node " << adjacent_node << std::endl;
					assert(it2 != sub_graph_nodes.end());
				}

			}
		}
	}
/*	std::cout << " make_sub_graph for : " << first_node <<std::endl;
	for(std::map<int, std::set<int> >::iterator it = sub_graph.begin() ; it != sub_graph.end() ; it++){
		std::cout<< "from " << it->first << " to "<<std::endl;
		for(std::set<int>::iterator it1 = sub_graph_nodes.begin() ; it1 != sub_graph_nodes.end() ; it1++){
			std::cout << *it1<<std::endl;
		}
	}*/
}
void ref_graph::bfs(int & startnode, std::string & refacc, size_t & right_on_ref){//Need the a container to save the length till the current position
	paths.clear();
	nodes_on_paths.clear();
	path_length.clear();
	std::set<int> seen;
	seen.insert(startnode);
	std::map<int, std::set<int> >::iterator adj = dynamic_adjacencies.find(startnode);//It is possible that a node has no adjacent node 
	std::cout << "start node in bfs "<< startnode <<std::endl;
	if(adj != dynamic_adjacencies.end()){
		std::string seq_name = seqname(startnode);
		size_t seqlength = seq_length(seq_name, refacc);
		size_t remainder = seqlength-(right_on_ref+1);
	//	std::cout<< remainder << " "<<seqlength << " "<<right_on_ref<< std::endl;
		std::map<std::pair<int,int>, bool> visited;
		std::map<int, std::vector<std::vector<int> > > all_paths; //int is the last node on a path, and vector contains all the paths that end to this node
		std::vector<std::vector<int> > first_path;
		std::vector<int> temp;
		temp.push_back(0);
		first_path.push_back(temp);
		all_paths.insert(std::make_pair(startnode,first_path));
		std::list<int> queue;
	//	std::map<int,size_t> length;
		std::map<int, std::vector<std::pair<std::vector<int> , size_t > > > length;
		std::vector<std::pair<std::vector<int>,size_t> > this_pair;
		this_pair.push_back(std::make_pair(temp, remainder));
		length.insert(std::make_pair(startnode, this_pair));
		size_t this_length = remainder;
		if(remainder >= MAXGAP){
			std::vector<int> temp;
			temp.push_back(startnode);
			paths.insert(temp);
			nodes_on_paths.insert(startnode);
			return ;
		}
		for(std::map<int, std::set<int> >::iterator it = dynamic_adjacencies.begin() ; it != dynamic_adjacencies.end() ;it++){
			for(std::set<int>::iterator this_adj = it->second.begin() ; this_adj != it->second.end(); this_adj ++){
      				visited.insert(std::make_pair(std::make_pair(it->first,*this_adj),false)); //Was( it->first ,false )
			}
 		}

	//	length.insert(std::make_pair(startnode,remainder));
	//	std::map<int,bool>::iterator vis=visited.find(startnode);
	//	assert(vis != visited.end());
	//	vis->second = true;
		queue.push_back(startnode);
		while(!queue.empty()){
			startnode = queue.front();
			std::cout << "on node "<< startnode<< std::endl;
			queue.pop_front();
		//	std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator len = length.find(startnode); //TODO
		//	assert(len != length.end());
		//	size_t this_length = len->second;
		//	std::cout << "this length " <<this_length <<std::endl;
		//	if(this_length < MAXGAP)
		//	if(all_paths.size()!=0){//Its member are removed when they reach MAXGAP
			std::map<int, std::set<int> >::iterator it = dynamic_adjacencies.find(startnode);
			if(it == dynamic_adjacencies.end()) std::cout << "node "<< startnode << " doesnt exist " <<std::endl;
			assert(it != dynamic_adjacencies.end());
			std::set<int> this_adjs = it->second;
			std::map<int, std::vector<std::vector<int> > >::iterator this_path = all_paths.find(startnode);
			/*	if(this_path == all_paths.end()){
					for(std::map<int, std::vector<std::vector<int> > >::iterator test = all_paths.begin() ; test != all_paths.end() ; test++){
						std::cout << test->first<<std::endl;
						std::cout << test->second << std::endl;
					}
				}*/
			//	assert(this_path !=all_paths.end());
			if(this_path !=all_paths.end()){
				std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator pre_len = length.find(startnode);
			//	std::cout << pre_len->second << std::endl;
				assert(pre_len != length.end());//TODO!
				bool EXPANDED = false;
				for(std::set<int>::iterator adj = this_adjs.begin() ; adj != this_adjs.end() ; adj++){
					std::cout << "adj" << *adj << std::endl;
					std::map<std::pair<int,int>, bool>::iterator vis =visited.find(std::make_pair(startnode,*adj));
				//	std::cout << "bool "<< int(vis->second)<<std::endl;
					if(vis->second == false){
						vis->second = true; 
						std::set<int>::iterator s = seen.find(*adj); //XXX Just commented! See what happens! TODO cause a problem in case of loop :(
						if(s == seen.end()){
							std::cout << "add to queue " << *adj << std::endl;
							queue.push_back(*adj);
						}
						seen.insert(*adj);
						int adjacent = *adj;
						std::string adj_name = seqname(adjacent);
						size_t adj_length = seq_length(adj_name, refacc);
						std::cout << adjacent << "adj length" << adj_length <<std::endl;
					//	len = length.find(*adj);
					//	assert(len == length.end());
					//	length.insert(std::make_pair(*adj,this_length + adj_length));
						std::vector<std::vector<int> >subpath = this_path->second;
					//	std::cout << subpath.size() <<std::endl;
						for(size_t i =0;i < subpath.size();i++){
							subpath.at(i).push_back(startnode);
							std::cout<< "Here: "<<subpath.at(i)<<std::endl;
						}
						if(subpath.size()==0){
							std::vector<int> temp;
							temp.push_back(startnode);
							subpath.push_back(temp);
						}
					//	std::cout << subpath <<std::endl;
						std::map<int, std::vector<std::vector<int> > >::iterator from_adj = all_paths.find(adjacent);
						if(from_adj == all_paths.end()){
							all_paths.insert(std::make_pair(adjacent,std::vector<std::vector<int> >()));
							from_adj = all_paths.find(adjacent);
						}
						assert(from_adj != all_paths.end());
						for(size_t i =0;i < subpath.size();i++){
							from_adj->second.push_back(subpath.at(i));
						}
						std::cout << "from adj is " << from_adj->second.size()<<std::endl;
						std::cout << from_adj->second << std::endl;
						std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator len = length.find(adjacent);
						std::vector<std::pair<std::vector<int>,size_t> > this_pair;
						if(len == length.end()){
							length.insert(std::make_pair(adjacent, this_pair));
							len = length.find(adjacent);
						}
						assert(len !=length.end());
						this_pair = pre_len->second;
						for(size_t i = 0; i < this_pair.size() ;i++){
							std::pair<std::vector<int>, size_t> p = this_pair.at(i);
							p.first.push_back(startnode);
							p.second+=adj_length;
						//	std::cout << p.second << std::endl;
							if(p.second< MAXGAP){
								len->second.push_back(p);
							}
							if(p.second>=MAXGAP){
								std::cout << "MAX reached! "<<std::endl;
								std::cout << p.first<<std::endl;
								std::vector<int> add_to_path = p.first;
								for(size_t n = 0; n < add_to_path.size() ; n++){
									nodes_on_paths.insert(add_to_path.at(n));
								}
								add_to_path.push_back(adjacent);
								nodes_on_paths.insert(adjacent);
								add_to_path.erase(add_to_path.begin());
								paths.insert(add_to_path);
								bool SEEN = false;
								for(size_t j =0; j < from_adj->second.size() ;j++){
								//	if(from_adj->second.size()== 6){
										std::cout<< from_adj->second.at(j)<<std::endl;
								//	}
									if(from_adj->second.at(j)==p.first){
										from_adj->second.erase(from_adj->second.begin()+j);
										SEEN = true;
										break;
									}
								}
							//	if(SEEN!=true) std::cout << from_adj->second.size() <<std::endl;
								assert(SEEN == true);
							}
						}
						if(from_adj->second.size()==0){
							all_paths.erase(from_adj);
							//Remove it from the queue
						}
						EXPANDED = true;
					}else{//XXX Just added the 'else' hope it works , XXX It solved the issue for the graphs with no loop but not sure if it is not working for graphs with loop
						std::cout<< "it is seen! "<<std::endl;
						std::map<int, std::vector<std::vector<int> > >::iterator from_adj = all_paths.find(startnode);
						if(from_adj != all_paths.end()){
							std::cout << from_adj->second.size() <<std::endl;
							for(size_t j = 0; j < from_adj->second.size() ; j++){
								std::cout << from_adj->second.at(j)<<std::endl;
							}
							std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator len = length.find(startnode);
							assert(len != length.end());
							std::set<std::vector<int> > seenvectors;
							for(std::map<int,std::vector<std::vector<int> > >::iterator it = all_paths.begin() ; it != all_paths.end() ; it++){
								std::cout<< "it second size " << it->second.size() <<std::endl;
								for(size_t p = 0; p < it->second.size() ; p++){
									std::vector<int>::iterator position = std::find(it->second.at(p).begin(), it->second.at(p).end(), startnode);
									if(position != it->second.at(p).end()){
										//Get a sub vector from this position on 
										size_t index = std::distance(it->second.at(p).begin(), position);
										std::vector<int> v1(it->second.at(p).begin() + index, it->second.at(p).end());
										std::cout << "v1 "<< v1 <<std::endl;
										std::set<std::vector<int> >::iterator seen = seenvectors.find(v1);
										if(seen == seenvectors.end()){
											seenvectors.insert(v1);
											//Add it at the end of all from_adj->second members!
											std::cout << "from adj size is " << from_adj->second.size()<<std::endl;
											for(size_t j = 0; j < from_adj->second.size() ; j++){
												std::vector<int> concatenation = from_adj->second.at(j);
												concatenation.insert(concatenation.end(), v1.begin(), v1.end());
												std::cout << "concat: " << concatenation << std::endl;
												it->second.push_back(concatenation);
												//TODO check and then add the length!
												std::map<int, std::vector<std::pair<std::vector<int> , size_t > > >::iterator len1 = length.find(it->first);
												assert(len1 != length.end());
												size_t this_length = len1->second.at(0).second; //TODO
												len1->second.push_back(std::make_pair(concatenation,this_length));
											
											}
										}
									}
								}
							}
							all_paths.erase(from_adj);
						//	exit(1);
						}
					}
				}
				if(EXPANDED == true){
					all_paths.erase(this_path);
				}
			}
		//	}else{
		//		std::cout << "end of a path "<<std::endl;
		//	}
			std::cout << "queue " << queue << std::endl;
		}
	//	std::cout << "HERE!!" << paths.size() <<std::endl;
	//	for(std::set<vector<int> >::iterator it = paths.begin() ; it != paths.end() ; it++){
	//		std::cout << *it << std::endl;
	//		std::cout << " "<<std::endl;
	//	}
		for(std::map<int,std::vector<std::vector<int> > >::iterator it = all_paths.begin() ; it != all_paths.end() ; it++){
			std::vector<std::vector<int> >path = it->second;
			for(size_t i =0; i < path.size() ; i++){
				path.at(i).push_back(it->first);
				path.at(i).erase(path.at(i).begin());
				std::set<vector<int> >::iterator it1 = paths.find(path.at(i));
				if(it1 == paths.end()){
					paths.insert(path.at(i));
					for(size_t j =0; j < path.at(i).size() ; j++){
						nodes_on_paths.insert(path.at(i).at(j));
					}
				}else{
					std::vector<int> temp = *it1;
					assert(temp.at(0) != 0);
				}
			//	for(size_t j =0; j < path.at(i).size() ; j++){
			//		nodes_on_paths.insert(path.at(i).at(j));
			//	}
			//	paths.push_back(path.at(i));
			//	std::cout << "this path "<< path.at(i)<<std::endl;
			}
		//	nodes_on_paths.insert(it->first);		
		}
	}
	std::cout << "nodes on path size is " << nodes_on_paths.size()<<std::endl;

}
void ref_graph::delete_path(std::vector<int> & this_path){
	//DELETE THEM FROM DY_ADJ!
	std::cout << "delete "<< this_path <<std::endl;
	for(size_t i = 0; i < this_path.size()-1 ; i++){
		std::map<int , std::set<int> >::iterator it = dynamic_adjacencies.find(this_path.at(i));
		assert(it != dynamic_adjacencies.end());
		std::set<int>::iterator it1 = it->second.find(this_path.at(i+1));
		if(it1 != it->second.end()){
			std::cout << this_path.at(i) << " to "<< *it1 <<std::endl;
			it->second.erase(it1);

		}
	}

}
const std::set<std::vector<int> > ref_graph::get_paths()const{
/*	paths.clear();
	std::set<int>::iterator it = nodes_on_paths.find(startnode);
	assert(it != nodes_on_paths.end());
	std::vector<int> temp
	temp.push_back(startnode);
	paths.insert(temp);
	std::map<int, std::set<int> >::iterator it1 = adjacencies.find(startnode);
	for(std::set<int>::iterator adj = this_adjs.begin() ; adj != this_adjs.end() ; adj++){
		std::set<int>::iterator adj1 = nodes_on_paths.find(*adj);
		if(adj1 != nodes_on_paths.end()){
			temp.push_back(*adj);
			paths.insert(temp);
			temp.pop_back();
			qeue.push_back(*adj);
		}*/
	return paths;
}
const std::set<int> ref_graph::get_nodes()const{
	return nodes_on_paths;
}
const std::set<int> ref_graph::get_nodes_on_paths(const int & next)const{
	std::set<int> temp;
	temp.insert(next);
	std::map<int, std::set<int> >::const_iterator it = pre_nodes_on_subgraph.find(next);
	assert(it != pre_nodes_on_subgraph.end());
	int startnode = next;
	std::set<int> seen_nodes;
	std::list<int> queue;
	queue.push_back(startnode);
	while(!queue.empty()){
		startnode = queue.front();
		queue.pop_front();
		seen_nodes.insert(startnode);
		it = pre_nodes_on_subgraph.find(startnode);
		if(it != pre_nodes_on_subgraph.end()){
			std::set<int> pre_nodes = it->second;
			for(std::set<int>::iterator it1 = pre_nodes.begin() ; it1 != pre_nodes.end() ; it1++){
				std::set<int>::iterator seen = seen_nodes.find(*it1);
				if(seen == seen_nodes.end()){
					seen_nodes.insert(*it1);
					queue.push_back(*it1);
					temp.insert(*it1);
				}
			}
		}
	}
	return temp;

}
const std::set<int> ref_graph::get_subgraph_nodes()const{
	return sub_graph_nodes;
}


void Kgraph::make_Kgraph(){
	std::map<int, std::string > nodes = rgraph.get_nodes_content();
	std::map<int , bool> visited;
	for(std::map<int, std::string >::iterator it = nodes.begin() ; it != nodes.end(); it++){
		assert(it->first > 0);
		visited.insert(std::make_pair(it->first,false));
		visited.insert(std::make_pair(-1*it->first,false));
		
	}
///	for(std::map<int, std::string >::iterator it = nodes.begin() ; it != nodes.end(); it++){
		std::map<int, std::string >::iterator it = nodes.begin(); //TODO beginning of all paths!
	//	std::string remainder;
	//	std::string pre_pre_seq;
		std::map<int, bool>::iterator vis = visited.find(it->first);
		assert(vis != visited.end());
		std::map<std::pair<int,int>, std::vector<std::string> > remainder;
		std::string temp_seq;
		std::vector<std::string> temp;
		temp.push_back(temp_seq);
		remainder.insert(std::make_pair(std::make_pair(it->first,it->first),temp));
		if(vis->second == false){
			bfs(remainder , it->first , it->first, visited);
		}
	//	vis = visited.find(-1*it->first);
	//	if(vis->second == false){
	//		dfs(remainder,pre_pre_seq , -1*it->first, visited);
	//	}
///	}
}
void Kgraph::bfs(std::map<std::pair<int,int>, std::vector<std::string> > & remainder, int  this_node, int  pre_node, std::map<int,bool> & visited){
	//	std::map<int,bool>::iterator vis = visited.find(this_node);
	/*	for(std::map<std::string , std::set<std::string> >::iterator it = edges.begin() ; it != edges.end() ; it++){
			std::cout << "edges from "<< it->first << " to" <<std::endl;
			for(std::set<std::string>::iterator it1 = it->second.begin() ; it1 != it->second.end() ; it1++){
				std::cout << *it1<<std::endl;
			}
		}*/
		std::string pre_pre_seq;
		std::map<std::pair<int,int>, std::vector<std::string> >::iterator it = remainder.find(std::make_pair(this_node, pre_node));
		assert(it != remainder.end());
		assert(it->second.size() ==1);
		std::string leftover = it->second.at(0);
		remainder.erase(it); //XXX ??? why??
		std::list<std::pair<std::pair<int,int>,std::string> > queue;
		std::string temp_seq;
		std::pair<std::pair<int,int>,std::string> temp(std::make_pair(this_node,this_node),temp_seq);
		queue.push_back(temp);
		int start_node = this_node;
		while(!queue.empty()){
	//		vis->second = true;
			std::pair<std::pair<int,int>,std::string> start_pair = queue.front();
			std::string leftover = start_pair.second;
			start_node = start_pair.first.first;
		//	std::cout << "start_node " << start_node <<std::endl;
			queue.pop_front();
			std::string content;
			if(leftover.length()>0){
				content = leftover;
			}
			if(this_node > 0){
				content += rgraph.get_content(this_node);
			}else{ //TODO check if it is correct!
				int rev = -1*this_node;
				std::string temp = rgraph.get_content(rev);
				std::string tempseq;
				get_reverse_complement(temp, tempseq);
				content += tempseq;
			}
			if(content.length() >= 2*K){
				int pre_seq_value;
				std::string pre_seq;
			/*	if(remainder.length()== K){
					pre_seq = remainder;
					Hash hashVa(pre_seq);
					const int h = hashVa.get_hvalue();
					pre_seq_value = h;
				}else if(remainder.length() > 0 && remainder.length() < K){
					assert(pre_pre_seq.length()>0);
					pre_seq = pre_pre_seq;
					Hash hashVa(pre_seq);
					const int h = hashVa.get_hvalue();
					pre_seq_value = h;
				}*/
				for (size_t i = 0; i < content.length(); i += K){
    					std:: string seq = content.substr(i, K);
				//	Hash hashVa(seq);
				//	const int h = hashVa.get_hvalue();
					std::map<std::string , std::set<std::string> >::iterator it2 = edges.find(seq);
					if(it2 == edges.end()&&seq.length()==K){
						edges.insert(std::make_pair(seq , std::set<std::string>()));
					}
					if(pre_seq.length() != 0 && seq.length()==K){
						it2 = edges.find(pre_seq);					
						assert(it2 != edges.end());
						it2->second.insert(seq);
					//	std::cout << "add an edge from " << pre_seq << " to "<< seq <<std::endl;
					}
					pre_pre_seq = pre_seq;
					pre_seq = seq;
				//	pre_seq_value = h;
				}
				if(pre_seq.length()==K){
					leftover = pre_seq;
				}else{
					assert(pre_seq.length()<K);
					leftover = pre_pre_seq + pre_seq;
				}
			}else{
				leftover = content;	
			//	std::cout << "rem " << leftover << std::endl;	
			}
			std::set<int> adjs = rgraph.get_adjacencies(start_node);
		//	std::cout << "adj size " << adjs.size() <<std::endl;
			for(std::set<int>::iterator it = adjs.begin() ; it != adjs.end() ; it++){
		/*		std::map<std::pair<int,int>, std::vector<std::string> >::iterator it1 = remainder.find(std::make_pair(*it,this_node));
			//	std::cout << "this adj " << *it<<std::endl;
				assert(it1 == remainder.end());
				std::vector<std::string> temp;
				temp.push_back(leftover);
				remainder.insert(std::make_pair(std::make_pair(*it,this_node), temp));*/
				queue.push_back(std::make_pair(std::make_pair(*it,start_node),leftover));
			}
		//	for(std::set<int>::iterator it = adjs.begin() ; it != adjs.end() ; it++){
		//		dfs(remainder, *it, this_node , visited);
		//	}
		}
}

void ref_graph::merge_nodes(std::string & ref_acc , size_t & LENGTH){
	merge_nodes_bfs(ref_acc, LENGTH);
}

void ref_graph::merge_nodes_bfs(std::string & ref_acc , size_t & LENGTH){//XXX I am thinking of adding nodes together and creating a new graph! //TODO how to deal with loops here???!!
	std::map<int,std::set<int> > copyadj = adjacencies;
	std::map<int, std::set<int> >::iterator it = copyadj.begin(); //TODO later on needs to be replaced by begin of each path
	std::set<std::pair<int, int> > seen_edges;
	int name = copyadj.size();
	int startnode = 1;
	std::cout << "begin node "<< startnode << std::endl;
	std::map<int, std::set<int> >::iterator edge = new_edges.find(startnode);
	assert(edge == new_edges.end());
	new_edges.insert(std::make_pair(startnode,std::set<int>()));
	std::list<int> queue;
	std::map<int, std::set<std::pair<std::vector<int> , size_t > > > length;
	std::set<std::pair<std::vector<int> , size_t > > temp;
	temp.insert(std::make_pair(std::vector<int>(),0));
	length.insert(std::make_pair(startnode , temp));
	queue.push_back(startnode);
	while(!queue.empty()){
		startnode = queue.front();
		std::cout << "this startnode "<<startnode <<std::endl;
	//	if(startnode == 20) break;
		queue.pop_front();
		it = copyadj.find(startnode);
		assert(it != copyadj.end());
		if(it->second.size()!=0){
			std::map<int, std::set<std::pair<std::vector<int> , size_t > > >::iterator prelength = length.find(startnode);
			assert(prelength != length.end());
			/* raw thought:
			size_t adj_length = length of it->second.at(i);
			adj_length+ all prelength
			if < Length
				add
			else
				add to new_edge - naming start from adj.size+1
				Replace that part of path with its new name ?
			*/
			for(std::set<int>::iterator adj = it->second.begin(); adj != it->second.end(); adj++){
				int adjacent = *adj;
			//	std::cout << "adjacent "<< adjacent <<std::endl;
				std::map<int, std::set<std::pair<std::vector<int> , size_t > > >::iterator thislength = length.find(adjacent);
				if(thislength==length.end()){
					length.insert(std::make_pair(adjacent , std::set<std::pair<std::vector<int> , size_t > >()));
					thislength = length.find(adjacent);
				}
				std::string adj_name = seqname(adjacent);
				size_t adj_length = seq_length(adj_name, ref_acc);
			//	std::cout << "pre path size " << prelength->second.size() <<std::endl;
				for(std::set<std::pair<std::vector<int> , size_t > >::iterator this_set =prelength->second.begin(); this_set != prelength->second.end() ; this_set++){
					std::pair<std::vector<int> , size_t > thispair = *this_set;
					std::vector<int> nodes_on_path = thispair.first;
					size_t curlength = thispair.second;
					assert(curlength < LENGTH);
					curlength +=adj_length;
					std::vector<int> oldnodes = nodes_on_path;
					oldnodes.push_back(startnode);
					if(curlength >= LENGTH){
						int tempnode = oldnodes.at(0);
						oldnodes.push_back(adjacent);
						oldnodes.erase(oldnodes.begin());
						std::map<std::vector<int> , int >::iterator oldedge = old_edges.find(oldnodes);
						if(oldedge==old_edges.end()){
						//	std::cout << "name " << name << " contains " << oldnodes <<std::endl;
							old_edges.insert(std::make_pair(oldnodes,name)); //The new node and original nodes it contains
							std::set<std::pair<std::vector<int> , size_t > > temp;
							temp.insert(std::make_pair(std::vector<int>(),0));
							length.insert(std::make_pair(name , temp));
							new_edges.insert(std::make_pair(name, std::set<int>()));
							queue.push_back(name);
							name ++;
							oldedge = old_edges.find(oldnodes);
						}
						
						std::map<int,std::set<int> >::iterator adjofadj = copyadj.find(adjacent);
						assert(adjofadj != copyadj.end());
						copyadj.insert(std::make_pair(oldedge->second,adjofadj->second));
						edge = new_edges.find(tempnode);
						assert(edge != new_edges.end());
						edge->second.insert(oldedge->second);


					}else{
						std::set<std::pair<std::vector<int> , size_t > >::iterator check_path = thislength->second.find(std::make_pair(oldnodes,curlength));
						if(check_path ==thislength->second.end()){
							thislength->second.insert(std::make_pair(oldnodes,curlength));
						}
						std::set<std::pair<int, int> >::iterator seen = seen_edges.find(std::make_pair(startnode,adjacent));
						if(seen == seen_edges.end()){
							queue.push_back(adjacent);
							seen_edges.insert(std::make_pair(startnode,adjacent));
						}
					}
	
				}
			}

		}
	}
	std::cout << "new edges size is " << new_edges.size() << " old edges size is " << adjacencies.size() << std::endl;
	for(std::map<int, std::set<int> >::iterator it = new_edges.begin() ; it != new_edges.end() ; it++){
		std::cout << it->first << " to " << std::endl;
		for(std::set<int>::iterator it1 = it->second.begin() ;it1 != it->second.end(); it1++){
			std::cout << *it1 <<std::endl;
		}
	}
	for(std::map<std::vector<int>, int>::iterator it = old_edges.begin() ; it != old_edges.end() ; it++){
		std::cout << "old edges: " << it->first  << std::endl;
		std::cout << "index " << it->second <<std::endl;
	}

/*	std::cout << "copyadj "<< std::endl;
	for(std::map<int, std::set<int> >::iterator it = copyadj.begin() ; it != copyadj.end() ; it++){
		std::cout << it->first << " to " << std::endl;
		for(std::set<int>::iterator it1 = it->second.begin() ;it1 != it->second.end(); it1++){
			std::cout << *it1 <<std::endl;
		}
	}*/

	

}

void ref_graph::edge_contraction(){
	std::set<int> seen;
	for(std::set<int>::iterator it = begin_nodes.begin(); it != begin_nodes.end() ; it++){
		int begin_node = *it;
		find_d1_nodes(begin_node,seen);
	}
}
void ref_graph::find_d1_nodes(int & begin_node, std::set<int> & seen){
	int startnode = begin_node;
	std::list<int> queue;
	queue.push_back(startnode);
	while(!queue.empty()){
		startnode = queue.front();
		std::cout << "this startnode "<<startnode <<std::endl;
		queue.pop_front();
		std::map<int, std::set<int> >::iterator it = adjacencies.find(startnode);
		assert(it != adjacencies.end());
		if(it->second.size()==1){
			std::cout << "adj1"<<std::endl;
			std::set<int>::iterator adjs = it->second.begin();
			int adjacentnode = *adjs;
			it = adjacencies.find(adjacentnode);
			assert(it != adjacencies.end());
			std::map<int,std::set<int> >::iterator newedge = new_edges.find(startnode);
			assert(newedge == new_edges.end());
			new_edges.insert(make_pair(startnode,std::set<int>()));
			newedge = new_edges.find(startnode);
			if(it->second.size() != 0){
				newedge->second = it->second;
			}//content of startnode + content of adjacentnode
			std::map<int, std::string>::iterator content = nodes_content.find(startnode);	
			assert(content != nodes_content.end());
			std::string seq = content->second;
			content = nodes_content.find(adjacentnode);	
			assert(content != nodes_content.end());
			seq.append(content->second);
			std::cout << "add to content "<< startnode <<std::endl;
			new_nodes_content.insert(std::make_pair(startnode,seq));
			for(std::set<int>::iterator adjOfadj = it->second.begin() ; adjOfadj != it->second.end() ; adjOfadj++){
				std::set<int>::iterator SEEN = seen.find(*adjOfadj);
				if(SEEN == seen.end()){
					queue.push_back(*adjOfadj);
					seen.insert(*adjOfadj);
				}
			}
		}else{
			for(std::set<int>::iterator adjs = it->second.begin(); adjs != it->second.end() ; adjs++ ){
				std::set<int>::iterator SEEN = seen.find(*adjs);
				if(SEEN == seen.end()){
					queue.push_back(*adjs);
					seen.insert(*adjs);
				}
			}
				
			std::map<int,std::set<int> >::iterator newedge = new_edges.find(startnode);
			assert(newedge == new_edges.end());
			new_edges.insert(make_pair(startnode,std::set<int>()));
			newedge = new_edges.find(startnode);
			if(it->second.size() != 0){
				newedge->second = it->second;
			}
			std::map<int, std::string>::iterator content = nodes_content.find(startnode);	
			assert(content != nodes_content.end());
			std::string seq = content->second;
			std::cout << "add to content "<< startnode <<std::endl;
			new_nodes_content.insert(std::make_pair(startnode,seq));


		}
	}
	std::cout << "new edges size is " << new_edges.size() << " old edges size is " << adjacencies.size() << std::endl;
	for(std::map<int, std::set<int> >::iterator it = new_edges.begin() ; it != new_edges.end() ; it++){
		std::cout << it->first << " to " << std::endl;
		for(std::set<int>::iterator it1 = it->second.begin() ;it1 != it->second.end(); it1++){
			std::cout << *it1 <<std::endl;
		}
	}

}

void ref_graph::write_gfa(std::ofstream & gfaout){//TODO leave itlike this for the time being but later on add the path!
	gfaout << "H\tVN:Z:1.0" << std::endl; //TODO add the name os input sequences!
	for(std::map<int , std::set<int> >::iterator it = new_edges.begin();  it != new_edges.end(); it++){
		std::cout << "node is " << it->first << std::endl;
		gfaout << "S\t"<<it->first<<"\t";
		std::map<int, std::string>::iterator it1 = new_nodes_content.find(it->first);
		assert(it1 != new_nodes_content.end());
		gfaout<<it1->second<< std::endl;
		for(std::set<int>::iterator it2= it->second.begin() ; it2 != it->second.end() ; it2++){
			if(it->first> 0){
				gfaout << "L\t"<<it->first<<"\t+\t";			
			}else{
				//TODO is it even happening??
				gfaout << "L\t"<<std::abs(it->first)<<"\t-\t";							
			}
			if(*it2 > 0){
				gfaout <<*it2<<"\t+\t0M"<<std::endl;							
			}else{
				gfaout <<std::abs(*it2)<<"\t-\t0M"<<std::endl;
			}

		}


	}


}
