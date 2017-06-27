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
				std::cout << "seq id is "<< findseq1->second << std::endl;
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
	std::cout << "n1 "<<name1 << " n2 "<<name2 <<std::endl;
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
	std::cout << "this node "<< this_node << " pre node " << pre_node<<std::endl;
	
}
const std::map<int , std::set<int> > & ref_graph::get_adjacencies()const{
	return adjacencies;
}


std::set<std::vector<int> > ref_graph::get_predecessor(unsigned int & ref_id, bool dir , size_t & left_on_ref, size_t & length_on_read){//int is ref number if negetive means reverse complement should be taken in to account
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
		
			if(line[0]== 'L'){
				std::vector<std::string> nodes;
				strsep(line, "\t" , nodes);
				std::string name1 , name2, dir1, dir2;
				name1 = nodes.at(1);
				name2 = nodes.at(3);
				dir1 = nodes.at(2);
				dir2 = nodes.at(4);
				std::cout << "L "<< line << std::endl;
				std::cout << "name1 "<< name1 << " name2 " << name2 << " dir1 " << dir1 << " dir2 " << dir2<<std::endl;
				add_adjacencies(dir1,dir2, name1, name2);
				add_adjacencies(dir2,name2);
				if(dir2 == "+"){
					dir2 = "-";
					add_adjacencies(dir2,name2);					
				}else{
					dir2 = "+";
					add_adjacencies(dir2,name2);	
				}

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
	//convert int to string:
	std::stringstream ss;
	ss << seqname;
	std::string longname = accname +":"+ss.str();
	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
	std::map<std::string, size_t>::iterator findseq = longname2seqidx.find(longname);
	assert(findseq != longname2seqidx.end());
	unsigned int refid = findseq->second;
	return refid;
}
size_t ref_graph::seq_length(std::string & seq_name, std::string & accname){
	std::string longname = accname +":"+seq_name;
	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
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
	std::map<int, std::set<int> >::iterator adj = adjacencies.find(startnode);//XXX It is possible that a node has no adjacent node
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
		size_t parent_length = 0;

		if(accu_length> MAXGAP){
			apath.push_back(startnode);
			nodes_on_paths.insert(startnode);
			paths.insert(apath);
		}else{
    			look_for_neighbors(startnode, visited, refacc,accu_length,apath,parent_length);
		}
	}//TODO else?!
}
void ref_graph::look_for_neighbors(int & node, std::map<int,bool> & visited , std::string & refacc, int & accu_length, std::vector<int> & apath, size_t & parent_length){//TODO Maybe i need to think of setting the parent length in a better way though it seems it works fine for the time being
	std::cout << "this node is "<< node << std::endl;
	std::map<int,bool>::iterator it = visited.find(node);
	assert(it != visited.end());
	it->second = true;
	std::string seq_name = seqname(node);
	size_t seqlength = seq_length(seq_name, refacc);
	std::cout << "its length is "<< seqlength<<std::endl;
	if(accu_length<=MAXGAP){
		apath.push_back(node);	
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
	std::cout<<std::endl;
	for (std::set<int>::iterator adj = it1->second.begin(); adj != it1->second.end(); adj++){
		int adjacent = *adj;
		std::string adj_name = seqname(adjacent);
		size_t adjlength = seq_length(adj_name, refacc);
		accu_length += adjlength;
		if(adj == it1->second.begin()){
			parent_length = seqlength;
			std::cout << "parent length "<< parent_length <<std::endl;
		}
		std::cout<< "adj is "<< *adj << "adj length "<< adjlength << std::endl;
		std::set<int>::iterator adj1 =it1->second.end();
		std::map<int,bool>::iterator it2 = visited.find(*adj);
		assert(it2 != visited.end());
		if (it2->second == false && accu_length < MAXGAP){
			std::cout << "smaller"<<std::endl;
			int current_node = *adj;
          		look_for_neighbors(current_node, visited,refacc,accu_length,apath,parent_length);
		}
		else if(it2->second == false && accu_length >= MAXGAP){ //XXX > only and the = belonged to the previous condition
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
	std::cout<<"paths size is "<<paths.size()<<std::endl;
}
void ref_graph::bfs(int & startnode, std::string & refacc, size_t & right_on_ref){//Need the a container to save the length till the current position
	paths.clear();
	nodes_on_paths.clear();
	path_length.clear();
	std::set<int> seen;
	std::map<int, std::set<int> >::iterator adj = dynamic_adjacencies.find(startnode);//It is possible that a node has no adjacent node
	std::cout << "start node in bfs "<< startnode <<std::endl;
	if(adj != dynamic_adjacencies.end()){
		std::string seq_name = seqname(startnode);
		size_t seqlength = seq_length(seq_name, refacc);
		size_t remainder = seqlength-(right_on_ref+1);
	//	std::cout<< remainder << " "<<seqlength << " "<<right_on_ref<< std::endl;
		std::map<std::pair<int,int>, bool> visited;
		for(std::map<int, std::set<int> >::iterator it = dynamic_adjacencies.begin() ; it != dynamic_adjacencies.end() ;it++){
			for(std::set<int>::iterator this_adj = it->second.begin() ; this_adj != it->second.end(); this_adj ++){
      				visited.insert(std::make_pair(std::make_pair(it->first,*this_adj),false)); //Was( it->first ,false )
			}
 		}
		std::map<int, std::vector<std::vector<int> > > all_paths;
		all_paths.insert(std::make_pair(startnode,std::vector<std::vector<int> >()));
		std::list<int> queue;
	//	std::map<int,size_t> length;
		std::map<int, std::vector<std::pair<std::vector<int> , size_t > > > length;
		std::vector<std::pair<std::vector<int>,size_t> > this_pair;
		std::vector<int> temp;
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
						std::cout << "from adj is "<<std::endl;
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
								paths.insert(add_to_path);
								bool SEEN = false;
								for(size_t j =0; j < from_adj->second.size() ;j++){
								//	if(from_adj->second.size()== 6){
								//		std::cout<< from_adj->second.at(j)<<std::endl;
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
					}/*else{//XXX Just added the 'else' hope it works 
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
					}*/
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
				std::set<vector<int> >::iterator it1 = paths.find(path.at(i));
				if(it1 == paths.end()){
					paths.insert(path.at(i));
					for(size_t j =0; j < path.at(i).size() ; j++){
						nodes_on_paths.insert(path.at(i).at(j));
					}
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

