#include "map_read.hpp"
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
			assert(it != predecessors.end());
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
//	std::cout << "start node in bfs "<< startnode <<std::endl;
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
		//	std::cout << "on node "<< startnode<< std::endl;
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
				//	std::cout << *adj << std::endl;
					std::map<std::pair<int,int>, bool>::iterator vis =visited.find(std::make_pair(startnode,*adj));
				//	std::cout << "bool "<< int(vis->second)<<std::endl;
					if(vis->second == false){
						vis->second = true; 
						std::set<int>::iterator s = seen.find(*adj);
						if(s == seen.end()){
							queue.push_back(*adj);
						}
						seen.insert(*adj);
						int adjacent = *adj;
						std::string adj_name = seqname(adjacent);
						size_t adj_length = seq_length(adj_name, refacc);
					//	std::cout << adjacent << "adj length" << adj_length <<std::endl;
					//	len = length.find(*adj);
					//	assert(len == length.end());
					//	length.insert(std::make_pair(*adj,this_length + adj_length));
						std::vector<std::vector<int> >subpath = this_path->second;
					//	std::cout << subpath.size() <<std::endl;
						for(size_t i =0;i < subpath.size();i++){
							subpath.at(i).push_back(startnode);
						//	std::cout<< "Here: "<<subpath.at(i)<<std::endl;
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
							//	std::cout << "MAX reached! "<<std::endl;
							//	std::cout << p.first<<std::endl;
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
					}
				}
				if(EXPANDED == true){
					all_paths.erase(this_path);
				}
			}
		//	}else{
		//		std::cout << "end of a path "<<std::endl;
		//	}
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
void als_components::find_als_on_paths(std::ofstream & output,size_t & refacc, size_t & readacc){
//	add_to_start.resize(alignments.size());
//	add_to_end.resize(alignments.size());
	shortest_path.resize(alignments.size());//Shortest path is found for each component
	for(size_t i =0; i < alignments.size(); i++){//Over components of alignments between one specific read and the reference graph
		std::cout << "all als at comp "<< i << " are " << alignments.size() << std::endl;
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
				std::cout<< "nodes are "<<std::endl;
				for(std::set<int>::iterator it = nodes.begin(); it!=nodes.end();it++){
					std::cout<< *it<<std::endl;
				}
				if(nodes.size()!=0){
				//	look_for_successors_on_paths(i, p, nodes, successors); //TODO
				//	finding_successors_on_subref(i, p, nodes, successors);
					finding_successors_on_subref(i, p, nodes, successors,refacc,readacc);
				}
			}
			std::cout << "find first and last als"<<std::endl;
			add_first_and_last_als(i,p,first_left,lastr);//Make alignment to add the current al to the start and end node, use NW with the right type--> it can start at any cell on the last row and ends at any cell on the first row (type 3)
			std::cout << "adj size is "<< adjacencies.size()<<std::endl;
			for(std::map<size_t, std::set<size_t> >::iterator it = adjacencies.at(i).begin(); it != adjacencies.at(i).end(); it++){
				std::cout << it->first << " connects to: " << std::endl;
				for(std::set<size_t>::iterator s = it->second.begin() ; s != it->second.end() ; s++){
					std::cout << *s << " ";
				}
				std::cout << " "<<std::endl;

			}
		}
		add_expensive_edges(i,refacc,readacc);//i-->component //XXX Temporary commented out, needs to get fixed! TODO

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
		}
		dijkstra dij(data,model,weight_of_edges.at(i));
		dij.find_shortest_path(first_left,lastr,shortest_path.at(i));//Input: begin and end of the component, Output: shortest path
		std::cout<< "shortest path is: " << std::endl;
		double sum = 0.0;
		for(size_t j = 0; j < shortest_path.at(i).size(); j++){
			if(shortest_path.at(i).at(j) != 0 && shortest_path.at(i).at(j) !=1){
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
		for(size_t j = 0; j < shortest_path.at(i).size(); j++){
			std::cout << shortest_path.at(i).at(j) << " ";
			if(shortest_path.at(i).at(j) != 0 && shortest_path.at(i).at(j) !=1){
				std::multimap<size_t,const pw_alignment>::iterator it=indices_nodes.at(i).find(shortest_path.at(i).at(j));
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

//XXX Started thinking of creatig the best path rather than all the poosible ones to reduce the running time.
void als_components::finding_successors_on_subref(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors,size_t & refacc, size_t & readacc){
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
	std::set<std::vector<int> > all_p = rgraph.get_paths();
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
		//	std::set<std::vector<int> > all_p = rgraph.get_paths();
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
						//		std::cout << this_p << std::endl;
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
				needleman<dynamic_mc_model> nl(data, model, readin, from_current_node);
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
			std::cout << "PATH size "<< PATH.size()<<std::endl; 
			if(PATH.size()==1){
				std::set<std::vector<int> >::iterator test_it = PATH.begin();
				std::cout << *test_it <<std::endl;
			}
			if(nodes_on_paths.size() != 0){
				find_best_path_to_neighbour(comp,nex_name,current_name, cur_index , this_al,nex_al,last_pos_on_read, last_pos_on_ref,refacc ,readacc, nodes_on_paths,seen_path,seen_indices,weight, cur_pre);
			}
		}
		if(it1 != nodes.end() && nex_name == current_name){
			std::cout << "on the same node "<< std::endl; //TODO
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
void als_components::find_best_path_to_neighbour(size_t comp, int & next_name, int & current_name, size_t & cur_index, const pw_alignment & cur_al, const pw_alignment & neighbour_al , size_t & to_on_read, size_t & last_pos_on_ref , size_t & refacc, size_t & readacc, std::set<int> & nodes_on_paths, std::map<std::pair<int, int> , std::set<size_t> > & seen_path, std::map<std::pair<int, int> , std::pair<size_t, double> > & seen_indices, std::map<size_t , double> & weight, std::map<size_t , size_t>  & cur_pre){
	//TODO Do I need a visited container? YES!

//	std::map<size_t , size_t> cur_pre;
//	std::map<size_t , double> weight; //al index and it is accumulative weight from the beginning of the path
	weight.insert(std::make_pair(cur_index, 0));
	std::multimap<double, size_t> al_distance; //Double --> accumulative gain, size_t --> node index
	pw_alignment this_al = cur_al;
	unsigned int read_id = this_al.getreference2();
	int startnode = current_name;
	size_t from_on_read;
	size_t pre_from_on_read = data.getSequence(read_id).length();
	std::set<int> visited;
	//Create the first al from the remaining part of the current node!

	while(startnode!=next_name){
		std::cout << "start node "<< startnode <<std::endl;
		std::cout << "al distance size "<< al_distance.size()<<std::endl;
		std::set<int> adjs = rgraph.get_adjacencies(startnode); 
		assert(adjs.size() != 0); //TODO if erase the path then it should be removed!
		std::map<const pw_alignment, size_t, compare_pw_alignment>::iterator id = node_indices.at(comp).find(this_al);
		assert(id != node_indices.at(comp).end());
		std::map<size_t , double>::iterator pre_wei = weight.find(id->second);
		if(pre_wei == weight.end()){
			std::cout << id->second;
			this_al.print();
			std::cout << "-----"<<std::endl;
		}
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
			std::set<int>::iterator onPaths = nodes_on_paths.find(this_name);
			if(onPaths == nodes_on_paths.end()){
				std::cout << "not on path! "<<std::endl;
				continue;
			}
			bool DIFFERENT_POS = false;
			std::map<std::pair<int, int>, std::set<size_t> >::iterator seen = seen_path.find(std::make_pair(startnode,*it));
		 	if(seen != seen_path.end()){ //TODO doesn't work for the loops needs to be changed!! for which position on read? position should be checked
				std::set<size_t>::iterator pos = seen->second.find(from_on_read);
				if(pos != seen->second.end()){ 
					if(from_on_read == to_on_read +1 || (r2 == data.getSequence(read_id).length())){//To deal with loops
				//	if(from_on_read == to_on_read +1){//To deal with loops
						std::cout << "add to visited "<< *it <<std::endl;
						visited.insert(*it);
					}
					std::cout << "end of seen " << startnode << " to " << *it <<std::endl;
					std::map<std::pair<int, int> , std::pair<size_t, double> >::iterator seen_weight = seen_indices.find(seen->first);
					if(seen_weight == seen_indices.end()) continue;
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
					needleman<dynamic_mc_model> nl(data, model, readin, refin);
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
void als_components::finding_successors_on_subref(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors){
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
void als_components::add_als(size_t & comp, size_t & refacc, size_t & readacc, std::vector<int> & path, std::vector<size_t> & refs, std::string & onref , std::string & onread, const pw_alignment & cur_al , const pw_alignment & nex_al, std::string & from_current, std::string & from_next){
	size_t l2,r2;
	cur_al.get_lr2(l2,r2);
	size_t pos = r2 +1;
	std::cout << "pos "<<pos << std::endl;
	size_t cur_ref = cur_al.getreference1();
	size_t readid = cur_al.getreference2();
	size_t nex_ref = nex_al.getreference1();
	std::string read_out, ref_out;
	needleman<dynamic_mc_model> nl(data, model, onread, onref);
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
void als_components::make_als(size_t & refacc, size_t & readacc, std::vector<int> & path, std::vector<size_t> & refs, std::string & onref , std::string & onread, std::vector<pw_alignment> & als, size_t & pos, size_t & cur_ref, size_t & nex_ref, std::string & from_current, std::string & from_next, size_t & readid , size_t & end_on_read_pos){
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
void als_components::get_al_samples(std::string & onref, std::string & onread, size_t & node_length , std::string & refout, std::string & readout, size_t & pos_on_read , bool & GAP){
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
void als_components::add_to_graph(size_t & comp, std::vector<pw_alignment> & als, const pw_alignment & cur_al , const pw_alignment & nex_al, size_t & refacc, size_t & readacc){
	for(size_t i =0; i < als.size()-1; i++){
		if(i == 0){
			add_adjacencies(comp,cur_al,als.at(0),refacc,readacc);
		}
		add_adjacencies(comp,als.at(i),als.at(i+1),refacc,readacc);
	}
	add_adjacencies(comp, als.back(), nex_al, refacc,readacc);
}
std::multiset<pw_alignment,sort_pw_alignment_by_left> als_components::find_successors(const pw_alignment & p, size_t & component){
	std::multiset<pw_alignment,sort_pw_alignment_by_left> successors;
	size_t l,r;
	p.get_lr2(l,r);
	std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it = alignments.at(component).find(p);
	assert(it != alignments.at(component).end());
	for(std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it1 = it; it1 != alignments.at(component).end(); it1++){
		pw_alignment al = *it1;
		size_t l2,r2;
		al.get_lr2(l2,r2);
		if((l2>r) && (l2-r <=MAXGAP)){//Overlap less than MAXGAP,//XXX For the moment we dont consider the cases that there is an overlap less than MAXGAP
			std::cout<<"find one! "<<std::endl;
			al.print();
			successors.insert(al);
		}
		if(l2 > r+MAXGAP){//since alignments are ordered by their left
			break;
		}
	}
	return successors;
}
void als_components::get_subgraph(const pw_alignment & p){//Return a sub graph that starts with reference1 of the current alignmnet and the nodes are closer than MAX_GAP to this first node.
	size_t l1,r1;
	p.get_lr1(l1,r1);
	p.print();
	size_t acc = data.accNumber(p.getreference1());
	std::string ref_accession = data.get_acc(acc);
	std::string seqname = data.get_seq_name(p.getreference1());
	int name = std::stoi(seqname);
	if(p.getbegin1()>p.getend1()){
		name = -1*name;
	}
	std::cout << "name is " << name <<std::endl;
//	rgraph.deep_first_search(name, ref_accession,r1);
	rgraph.bfs(name,ref_accession,r1);
//	paths = rgraph.get_paths();
//	std::cout<< "path size is "<< paths.size()<<std::endl;
//	for(std::set<std::vector<int> >::iterator it = paths.begin(); it != paths.end(); it++){
//		std::vector<int> this_p = *it;
//		for(size_t j =0; j < this_p.size();j++){
//			std::cout << this_p.at(j)<< " ";
//		}
//		std::cout<< " "<<std::endl;
//	}
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

void als_components::look_for_successors_on_paths(size_t & comp, const pw_alignment & current_al ,std::set<int>& nodes, std::multiset<pw_alignment,sort_pw_alignment_by_left> & neighbors){//XXX Directions should be checked on ref 
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

				}else{//When both are backwards
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
			//example:
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
void als_components::get_previous_als(size_t & comp , int & index, std::vector<pw_alignment> & pre_als, std::vector<size_t> & onread_from,size_t & to, size_t & ref_acc, size_t & read_acc){//XXX Originally I was only checking the last alignment, but it was wrong. I need to check all of them! Becasue getting from different nodes to node index may create different als on it and all need to be taked into account.
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
void als_components::add_to_containers(size_t & comp , const pw_alignment & p1 , const pw_alignment & p2 , size_t & ref_acc, size_t & read_acc, int & current_name, int & next_name){
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
void als_components::add_nodes_to_the_graph(size_t & type , std::string & seq_from_read, std::string & seq_from_ref, size_t & refacc , size_t & readacc, pw_alignment & nw_alignment, unsigned int & read_id , unsigned int & ref_id, size_t & onread_from, size_t & onref_from){//Always one node except on the first node
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
		needleman<dynamic_mc_model> nl(data, model, seq_from_read, seq_from_ref);
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

void als_components::make_al_on_a_node(size_t & comp, const pw_alignment& first_al , const pw_alignment& second_al, bool dir, unsigned int & ref1, size_t & left1, size_t & right1, unsigned int & ref2, size_t & left2, size_t & right2,size_t & refacc, size_t & readacc){
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
		needleman<dynamic_mc_model> nl(data, model, str2, str1);
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
	if(str1.length() == 0 && str2.length() !=0){
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
void als_components::get_paths(const std::set<int> & bfs_nodes, int & node_name , int & current_node_name, std::set<std::vector<int> > & paths){ //TODO do it with adjacencies also when they are read from the te stack add them to a seen container and never add the seen ones to the stack any more to avoid the loop!!
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
void als_components::append_nodes(const std::set<int> & bfs_nodes ,int & name, int & current_node_name, unsigned int & next_ref_id, unsigned int & current_ref_id, size_t & refacc , std::vector<std::vector<size_t> >& all_refs, std::vector<std::vector<int> >& all_paths, std::vector<std::vector<std::string> >& all_strings){
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
void als_components::make_alignments(size_t & comp, size_t & begin_on_read, size_t & end_on_read, const pw_alignment & first_al, const pw_alignment & last_al, std::string & read_out, std::string & ref_out, std::vector<size_t> &this_refs, std::vector<int> & this_path, unsigned int & ref2, size_t & first_length, size_t & last_length,size_t &refacc,size_t&readacc){//TODO Further improvement: maybe it makes sense to make another pw_alignment class with no ref number!
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
void als_components::compute_samples(bool tillend, size_t node_length, size_t & current_pos, size_t & read_counter, size_t & ref_counter,std::string & read_in, std::string & read_out, std::string & ref_in, std::string & ref_out){
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
void als_components::get_reverse_complement(std::string & sequence_in , std::string & sequence_out){
	for(size_t i = sequence_in.size(); i >0; i--){
		sequence_out += dnastring::complement(sequence_in.at(i-1));
	}
}


void als_components::add_first_and_last_als(size_t & i, const pw_alignment & p,size_t & first_left, size_t & last_right){//Create alignments that connect first and last alignments to the origin and the end.	
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
		if(p.getreference2()==1275){
			std::cout<< "length " << length <<std::endl;
			std::cout<< "r2 "<< r2 <<std::endl;
		}
		if((length - r2 < MAXGAP) || ((last_right-r2) + MAXGAP/2 <MAXGAP)){
			std::cout<< "add it to the end " << last_right << std::endl;
			p.print();
			looking_for_last_al(i,p,last_right);		
		}
//	}
}
void als_components::looking_for_first_al(size_t & comp, const pw_alignment & p, size_t & first_left){//We may create more than one al but keep the one with the highest gain from ref to read
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
	if(l1 != 0 && p.getbegin1() < p.getend1()){
		from = 0;
		to = l1-1;
		seq_from_ref = data.extract_seq_part(ref1,from,to);
	}
	if(r1 != data.getSequence(ref1).length()-1 && p.getbegin1()>p.getend1()){
		to = r1+1;
		from = data.getSequence(ref1).length()-1;
		std::cout<<"looking for first als begin > end "<<std::endl;
		std::string temp_seq_in = data.extract_seq_part(ref1, to, from);
		std::string temp_seq_out;
		get_reverse_complement(temp_seq_in, temp_seq_out);
		seq_from_ref = temp_seq_out;
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
			size_t ref_from, ref_to;
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
			needleman<dynamic_mc_model> nl(data, model, seq, seq_from_ref);
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
			std::set<std::vector<int> > refs;//seq name it is not the ref id!
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
				//Gap + seqfromref
				if(seq_from_ref.length()!=0){
					std::reverse(seq.begin(),seq.end());
					std::string refin = seq_from_ref;
					std::reverse(refin.begin(),refin.end());
					std::cout<< "refin "<< refin << std::endl;
					std::string read_out, ref_out;
					needleman<dynamic_mc_model> nl(data, model, seq, refin);
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
						pw_alignment p1(ref_out,read_out,0,new_al_begin,seq_from_ref.length()-1,new_al_end,ref1,ref2);
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

void als_components::make_first_nw_als(size_t & comp, const pw_alignment & current_al, std::set<std::vector<int> > & refs, std::string & seq_from_ref, unsigned int & read_id, size_t & ini_read_from , size_t & read_to , size_t & read_acc, size_t & ref_acc){
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
						needleman<dynamic_mc_model> nl(data, model, seq_from_read, seq_from_ref);
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
							pw_alignment nw_al(ref_out,read_out,from, read_from ,to, read_to, ref_id, read_id);
							nw_al.print();
							this_pal = nw_al;
							add_to_containers(comp, nw_al, current_al, ref_acc, read_acc, this_p.at(j), this_p.at(j));
							if(read_from == ini_read_from){
								add_edge_to_begin(comp, nw_al, read_acc, ref_acc);
							}
						}else{
							pw_alignment nw_al(ref_out,read_out,to, read_from ,from, read_to, ref_id, read_id);
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
							needleman<dynamic_mc_model> nl(data, model, readin, refin);
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
double als_components::get_cost(const pw_alignment & p, size_t & acc1, size_t & acc2){
	const std::map<std::string, std::vector<double>  > al_cost = model.get_al_cost(acc1, acc2);
}
void als_components::looking_for_last_al(size_t & comp, const pw_alignment & p, size_t & last_right){
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
		std::cout << "seq length is "<< seq.length()<<std::endl;//Remained length from the read
		if(seq_from_ref.length()>=MAXGAP){
			//seq_from_ref.erase(seq_from_ref.begin()+MAXGAP,seq_from_ref.end());
			//assert(seq_from_ref.length()==MAXGAP);
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
			needleman<dynamic_mc_model> nl(data, model, seq, seq_from_ref);
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
			std::vector<std::vector<int> > refs;//The adjacent nodes from the current one on the ref graph(they are seq names not the ref ids!)
			size_t remainder = MAXGAP- seq_from_ref.length();
			std::cout << "length is "<< remainder << "ref1 is "<< ref1 <<std::endl;
			if(p.getbegin1()<p.getend1()){
				refs=rgraph.get_successor(ref1,true,remainder);
			}else{
				refs=rgraph.get_successor(ref1,false,remainder);//TODO
			}
			if(refs.size()!=0){
				std::cout<< "post nodes size "<< refs.size() <<std::endl;
				for(size_t i = 0 ; i < refs.size() ;i++){
					std::cout << refs.at(i)<<std::endl;
					std::string refin = seq_from_ref;
					for(size_t j = 0; j < refs.at(i).size(); j++){
						std::cout<< refs.at(i).at(j)<<std::endl;		
						//Make an al between each member of refs and seq and pick the best
						size_t from = 0;
						if(refs.at(i).at(j)>0){
							unsigned int ref_id = rgraph.get_refid(refacc,refs.at(i).at(j));//gets node name, retruns node id
							size_t to = data.getSequence(ref_id).length()-1;
							refin += data.extract_seq_part(ref_id, from, to);
						}else{
							int temp = -1*refs.at(i).at(j);
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
					needleman<dynamic_mc_model> nl(data, model, seq, refin);
					nl.run_needleman(readacc,refacc,type,read_out,ref_out);
					//Make all the als!
					size_t read_ref = ref2;
					size_t current_node = ref1;
					std::vector<pw_alignment> last_als;
					make_last_als(refs.at(i),ref_out,read_out,new_al_begin,new_al_end,seq_from_ref,r1,refacc,read_ref, current_node,last_als); 
					all_last_als.push_back(last_als);
				}
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
			}else{//Either only gap or the remainder from the ref node
				if(seq_from_ref.length()!=0){
					std::string refin = seq_from_ref;
					std::string read_out, ref_out;
					needleman<dynamic_mc_model> nl(data, model, seq, refin);
					nl.run_needleman(readacc,refacc,type,read_out,ref_out);
					size_t ref_begin = r1+1;
					size_t ref_end = data.getSequence(ref1).length()-1;
					if(p.getbegin1()<p.getend1()){ //TODO what is it is ==
						pw_alignment p1(ref_out,read_out,ref_begin,new_al_begin,ref_end,new_al_end,ref1,ref2);
						add_adjacencies(comp, p, p1,refacc,readacc);
						add_edge_to_end(comp,p1);
					}else{
						ref_begin = l1 -1;
						ref_end = 0;
						pw_alignment p1(ref_out,read_out,ref_begin,new_al_begin,ref_end,new_al_end,ref1,ref2);
						add_adjacencies(comp, p, p1,refacc,readacc);
						add_edge_to_end(comp,p1);
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
void als_components::make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout, size_t & left_on_current_node,size_t & refacc, size_t & readacc, size_t & read_id, std::vector<pw_alignment> & first_als, size_t & remainder_on_ref){
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
void als_components::make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & left_on_current_node,size_t & refacc, size_t & read_id, size_t & current_node, std::vector<pw_alignment> & first_als){//TODO !! It is wrong
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
void als_components::make_last_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & right_on_current_node,size_t & refacc, size_t & read_id, size_t & current_node, std::vector<pw_alignment> & last_als){//TODO read counter == 0, and end position is reached condition!!
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
	std::cout<<"nodes size "<<nodes.size()<<std::endl;
	for(size_t i = 0 ; i < nodes.size() ;i++){
		std::cout<< "i "<< i <<std::endl;
		unsigned int ref_id;
		if(nodes.at(i)>0){
			ref_id = rgraph.get_refid(refacc,nodes.at(i));//gets node name, retruns node id
		}else{
			int temp = -1*nodes.at(i);
			ref_id = rgraph.get_refid(refacc,temp);
		}
		size_t current_node_length = data.getSequence(ref_id).length();
		bool tillend = false;
		if(i==nodes.size()-1){//When the ref sample ends with gaps/solution is not to break it when it is the last node in nodes
			tillend = true;
		}
		compute_samples(tillend, current_node_length, current_position, read_counter, ref_counter, readout, this_read_part, refout, this_ref_part);//TODO what if refcounter stays zero?
		std::cout << refout.size() << " " <<current_position <<std::endl;
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
}
	std::cout<< b2 << " "<< last_end <<std::endl;
	assert(b2-1==last_end);
}

const std::vector<pw_alignment> als_components::find_the_best_als(std::vector<std::vector<pw_alignment> >& all_als,size_t & refacc, size_t & readacc)const{
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

void als_components::add_adjacencies(size_t & comp, const pw_alignment & p1 , const pw_alignment & p2, size_t & ref_acc, size_t & read_acc){
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

void als_components::add_edge_to_begin(size_t & comp , const pw_alignment & p, size_t & read_acc, size_t & ref_acc){
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
void als_components::add_edge_to_end(size_t & comp , const pw_alignment & p){
	p.print();
	std::map<const pw_alignment ,size_t,compare_pw_alignment>::const_iterator it=node_indices.at(comp).find(p);
	assert(it != node_indices.at(comp).end());
	std::map<size_t , std::set<size_t> >::iterator adj = adjacencies.at(comp).find(it->second);
	if(adj == adjacencies.at(comp).end()){
		adjacencies.at(comp).insert(std::make_pair(it->second, std::set<size_t>()));
		adj = adjacencies.at(comp).find(it->second);
	}
	adj->second.insert(1);

//	adjacencies.at(comp).insert(std::make_pair(it->second,1));
	weight_of_edges.at(comp).insert(std::make_pair(std::make_pair(it->second,1),0));
}
const pw_alignment & als_components::get_alignment(size_t & comp, size_t & nodeid)const{
	typename std::multimap<size_t , const pw_alignment>::const_iterator it = indices_nodes.at(comp).find(nodeid);
	assert(it != indices_nodes.at(comp).end());
	return it->second;
}
void als_components::add_expensive_edges(size_t & comp,size_t & refacc, size_t & readacc){
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
void als_components::add_to_maf(const pw_alignment & al, std::ofstream & output, bool & firstal){//TODO fix only gap alignments TODO thelengt of alignments are wrong!!
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
		//	if(al.get_al_ref1().size()!=1){
		//		std::cout<< al.get_al_ref1()<<std::endl;
		//	}
		//	assert(al.get_al_ref1().size()==1);//It is possible that begin == end but the length is not zero. it happens when there are gaps around that single position. So this assertion can be wrong.
			char this_base;
			std::cout << al.get_al_ref1()<<std::endl;
			for(size_t i =0; i< al.get_al_ref1().size();i++){
				if(al.get_al_ref1().at(i)!='-'){
					this_base = al.get_al_ref1().at(i);
					break;
				}
			}
			std::cout <<"length " << data.getSequence(ref1).length() <<std::endl;
			if(this_base==data.getSequence(ref1).at(l1)){//To see if the sequence of length one is forward or backward
				output << "s " << longname1.str() << " " << al.getbegin1() << " " << al.getend1()-al.getbegin1()+1 << " + " << data.getSequence(ref1).length() << " " << al.get_al_ref1() <<std::endl;
				forward = true;
			}else{
			//	assert(al.get_al_ref1().at(0)== dnastring::complement(data.getSequence(ref1).at(l1)));
				assert(this_base == dnastring::complement(data.getSequence(ref1).at(l1)));

				output << "s " << longname1.str() << " " << al.getbegin1() << " " << al.getend1()-al.getbegin1()+1 << " - " << data.getSequence(ref1).length() << " " << al.get_al_ref1() <<std::endl;
				forward = false;
			}
			previous_right1 = r1;
			previous_left1 = l1;
			previous_ref = ref1;


			
		}
		if(al.getbegin2()< al.getend2()){
			output << "s " << longname2.str() << " " << al.getbegin2() << " " << al.getend2()-al.getbegin2()+1 << " + " << data.getSequence(ref2).length() << " " << al.get_al_ref2() <<std::endl;
			previous_right2 = r2;
		}else{
			assert(al.getbegin2()==al.getend2());
			if(r2 == data.getSequence(ref2).length()){//Alignments with only gap on the read
				std::cout<< "only gap on read " << std::endl;
				output << "s " << longname2.str() << " " << previous_right2 << " " << 0 << " + " << data.getSequence(ref2).length() << " " << al.get_al_ref2() <<std::endl;//TODO check for the direction!
			}else{
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

void map_check::read_graph_maf(std::ifstream & graph_maf){
	std::cout<<"read ref_graph maf "<<std::endl;
	if(graph_maf) {
		std::string str;
		getline(graph_maf, str);
		std::cout << str.length()<<std::endl;
		while(str.at(0)=='#') { // skip headers
			getline(graph_maf, str);
		}
		assert(str.at(0)!='#');
		while(!graph_maf.eof()) {
			if(str.at(0)=='a'){
				//get the context and split it when '=' happened.
				std::vector<std::string> clusternumber;
				strsep(str, "=" , clusternumber);
				unsigned int id = atoi(clusternumber.at(1).c_str());
				std::cout<<"id "<< id<<" " <<std::endl;
				clusters.insert(std::make_pair(id,std::vector<std::string>()));
				getline(graph_maf, str);
				assert(str.at(0)=='s');
				bool IsCenter = true;
				std::string ref1,ref2;
				size_t begin1, begin2, end1, end2, left2, right2;
				size_t counter = 0;
				while(str.at(0)=='s'){
					counter ++;
					std::vector<std::string> parts;
					strsep(str, " ", parts);
					assert(parts.size()==7);
					std::stringstream sstr;
					if(IsCenter == true){
						ref1 = parts.at(6);
						ref2 = parts.at(6);
						IsCenter = false;
						unsigned int nodeLength = stoi(parts.at(3));
						unsigned int seq_length = stoi(parts.at(5));
						node_length.insert(std::make_pair(id, nodeLength));
						if(parts.at(4) == "+"){
							sstr<<parts.at(2);
							sstr >> begin1;
							end1 = begin1+nodeLength-1;
							begin2 = begin1;
							end2 = end1;
							left2 = begin2;
							right2 = end2;
						}else{//On the reverse strand the begin starts at end of the forward strand 
							sstr<<parts.at(2);
							sstr >> end1;
							left2 = seq_length-(end1+nodeLength-1)-1;
							right2 = left2 + nodeLength-1;
							begin1 = right2;
							end1 = left2;
							begin2 = begin1;
							end2 = end1;
						}
					}else{
						ref2 = parts.at(6);
						unsigned int length = stoi(parts.at(3));
						unsigned int seq_length = stoi(parts.at(5));
						if(parts.at(4) == "+"){
							sstr<<parts.at(2);
							sstr >> begin2;
							end2 = begin2+length-1;
							left2 = begin2;
							right2 = end2;
						}else{//On the reverse strand the begin starts at end of the forward strand 
							sstr<<parts.at(2);
							sstr >> end2;
							left2 = seq_length-(end2+length-1)-1;
							right2 = left2 + length-1;
							begin2 = right2;
							end2 = left2;
							std::cout << "b2 "<< begin2 << " e2 "<<end2 <<std::endl;

						}

					}
					std::vector<std::string>acc_parts;
					strsep(parts.at(1),":", acc_parts);
					std::stringstream curseq;
					if(id == 24 && acc_parts.at(0)=="nc000913.3"){
						std::ofstream out("maftest");
						for(size_t i = 0; i< ref2.length(); i++){
							if(ref2.at(i) != '-'){
								out<<ref2.at(i);
							}
						}
					}
					if(acc_parts.at(0)=="nc000913.3"){
						std::cout << "ref1 length " << ref1.length() << " ref2 length "<< ref2.length() <<std::endl;
						pw_alignment p(ref1, ref2, begin1, begin2, end1, end2, id, 0);
						al_from_graph.insert(std::make_pair(id,p));
		
					}
	
					//add them to the cluster container as name:clusternumber
					curseq<<acc_parts.at(0)<<":"<<id;
					std::cout <<curseq.str()<<std::endl;
					std::map<unsigned int, std::vector<std::string> >::iterator it = clusters.find(id);
					assert(it != clusters.end());
					it->second.push_back(curseq.str());
					//add them with their begin and end to the boundries container
				//	size_t start = atoi(parts.at(2).c_str());
				//	size_t length = atoi(parts.at(3).c_str());
				//	size_t end = length+start-1;
					size_t start = left2;
					size_t end = right2;
					boundries.insert(std::make_pair(curseq.str(),std::make_pair(start,end)));
					std::stringstream thisseq;
					thisseq<<acc_parts.at(0);
					test_boundries.insert(std::make_pair(thisseq.str(),std::make_pair(start,end)));
					std::cout<< start << " " << end << std::endl;
					getline(graph_maf, str);
					if(str.length()==0){
						break;
					}
				}
				if(counter > 2){
					std::cout << "big cluster "<< id << std::endl;

				}
			}else{
				std::cout<< "shouldnt happen! "<<std::endl;
				exit(1);
			}
		}
		graph_maf.close();
	}
}
void map_check::read_alignments(std::ifstream & alignments){
	std::string str;
	while(getline(alignments, str)) {
		if(str.at(0)!='#') { // skip headers
			if(str.at(0)=='a') { // next alignment 
				std::cout << str <<std::endl;
				std::string aline1;
				std::string aline2;
				getline(alignments, aline1);
				getline(alignments, aline2);
				getline(alignments, str); 

				std::vector<std::string> parts1;
				strsep(aline1, " ", parts1);
				assert(parts1.size()==7);
				std::vector<std::string>node_parts;
				strsep(parts1.at(1),":", node_parts);
				assert(node_parts.size()==2);

				std::vector<std::string> parts2;
				strsep(aline2, " ", parts2);
				assert(parts2.size()==7);
				std::vector<std::string>read_parts;
				strsep(parts2.at(1),":", read_parts);
				assert(read_parts.size()==2);

				std::stringstream temp;
				temp<<node_parts.at(1)<<":"<<read_parts.at(1);
				std::multimap<std::string , std::pair<unsigned int ,unsigned int> >::iterator it=nodes.find(temp.str());	
				if(it==nodes.end()){
					nodes.insert(std::make_pair(temp.str(),std::pair<unsigned int, unsigned int>()));
					it = nodes.find(temp.str());
				}
				unsigned int from = atoi(parts2.at(2).c_str());
				unsigned int length = atoi(parts2.at(3).c_str());
				unsigned int to = from+ length-1;
				it->second = std::make_pair(from,to);				
			}
		}
	}
	std::cout << "done! "	<<std::endl;
}
void map_check::check_output(std::ifstream & mapping_maf,const std::string & seqname){
	std::cout << "checking here! "<<std::endl;
	std::map<std::string,std::set<std::pair<size_t,size_t> > > sorted;
	sorted.insert(std::make_pair(seqname,std::set<std::pair<size_t,size_t> >()));
	for(std::multimap<std::string , std::pair<size_t, size_t> >::iterator it = test_boundries.begin(); it != test_boundries.end(); it++){
		if(it->first == seqname){
			std::map<std::string,std::set<std::pair<size_t,size_t> > >::iterator it1=sorted.find(seqname);
			assert(it1 != sorted.end());
			it1->second.insert(it->second);
		}
	}
	for(std::map<std::string,std::set<std::pair<size_t,size_t> > >::iterator it = sorted.begin() ; it != sorted.end() ; it++){
			std::cout<< it->first << " ";
			std::set<std::pair<size_t , size_t> > pairs = it->second;
			for(std::set<std::pair<size_t, size_t> >::iterator it1 = pairs.begin() ; it1 != pairs.end() ; it1++){
				std::pair<size_t,size_t>this_pair = *it1;
				std::cout<< this_pair.first << " "<<this_pair.second << std::endl;
			}
	}
	unsigned int refid;
	unsigned int readid = 0;
	size_t accu_length = 0;
	size_t count = 0;
	if(mapping_maf) {
		std::string str;
		while(getline(mapping_maf, str)) {
			if(str.at(0)!='#') { // skip headers
				if(str.at(0)=='a') { 
					count ++;
					std::string al1;
					std::string al2;
					getline(mapping_maf, al1);
					getline(mapping_maf, al2);
					std::vector<std::string> parts1;
				//	std::cout<< al1 <<std::endl;
					strsep(al1, " ", parts1);
					assert(parts1.size()==7);
					std::vector<std::string>ref_parts;
					strsep(parts1.at(1),":", ref_parts);
					refid = atoi(ref_parts.at(1).c_str());
					std::cout << "ref id "<< refid <<std::endl;
					std::vector<std::string> parts2;
					strsep(al2, " ", parts2);
					assert(parts2.size()==7);

					std::vector<std::string>read_parts;
					strsep(parts2.at(1),":", read_parts);
					std::string read_id = read_parts.at(1).c_str();
					std::cout << "read is "<< read_id << " its id is " << readid <<std::endl;

					size_t start = atoi(parts2.at(2).c_str());
					size_t length = atoi(parts2.at(3).c_str());
					accu_length += length;
					std::cout<<"length "<< length << " accu length "<< accu_length <<std::endl;

					size_t end = length+start-1;//TODO consider only gap als!
					std::cout << "from " << start << " to "<< end <<std::endl;
					size_t startonseq;
					size_t endonseq;
					if(readid>=1){
						startonseq = (readid*30030)+start;
						endonseq = length+startonseq-1;
					}else{
						startonseq = start;
						endonseq = end;
					}
					std::cout<< "startonseq "<< startonseq << " endonseq "<<endonseq<<std::endl;
					check_an_alignment(refid,seqname, startonseq ,endonseq);
					if(end == 30029){//XXX It is so specific and works in this case only since i set the lengths to the 30030, later on i should sum up the length of reads up here 
						readid++;
						accu_length = 0;
					}

				}

			}

		}
		mapping_maf.close();

	}
	std::cout <<"number of als "<< count<<std::endl;
	std::cout << "number of wrongly mapped bases: "<<error_counter << std::endl;


}
void map_check::check_an_alignment(unsigned int & ref1, const std::string & seqname, size_t & left2 , size_t & right2){
	std::map<unsigned int, std::vector<std::string> >::iterator it = clusters.find(ref1); //It reaches the end if it is a non aligned region 
	if(it != clusters.end()){
		std::stringstream str;
		str<< seqname<<":"<<ref1;
		std::cout << str.str()<<std::endl;
		std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator find_seq = boundries.find(str.str());
	//	assert(find_seq != boundries.end());
		std::pair<std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator, std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator> it1=boundries.equal_range(str.str());
		bool mapped = false;
		for(std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator intervals= it1.first ; intervals!=it1.second; intervals++){
			std::cout<<"here! " << intervals->first <<std::endl;
			std::pair<size_t,size_t>this_pair = intervals->second;
			std::cout<< this_pair.first << " "<< this_pair.second<< std::endl;
			if(left2>=this_pair.first && right2<=this_pair.second){
				std::cout<< "correctly mapped!"<<std::endl;
				std::cout<<"from " << left2 << " to "<< right2 <<std::endl;
				mapped =true;
				break;
			}
		}
		if(mapped == false){
			std::cout<< "mapping error!"<<std::endl;//TODO find this boundry in other clusters, if it is in another cluster that could be an issue :( slow test!
			error_counter += (right2 - left2);
			bool wronglyclustered = false;
			std::stringstream thisref;
			thisref<< ref1;
			for(std::multimap<std::string, std::pair<size_t, size_t> >::iterator it2 = boundries.begin() ; it2 != boundries.end() ; it2++){
				std::string seqlongname = it2->first;
				std::pair<size_t , size_t> this_interval = it2->second;
				std::vector<std::string> parts;
				strsep(seqlongname, ":", parts);
				if(parts.at(0)== seqname && parts.at(1)!=thisref.str() && this_interval.second<= 300000){
					std::cout<< "in cluster "<< parts.at(1)<< " " << this_interval.first << " "<<this_interval.second << std::endl;
				}
				if(parts.at(0)==seqname && left2>=this_interval.first && right2<=this_interval.second){
					std::cout << "clusterd to "<<parts.at(1)<<std::endl;
				}
			}
			
		}			
	}else{
		//If it might be mapped to a cluster of size one
		std::cout<< "mapped to a non aligned region!" << ref1 <<std::endl;
		std::map<int,std::pair<int, int> >::iterator it = non_aligned.find(ref1);
		if(it == non_aligned.end()){
			std::cout <<"mapping error"<<std::endl;
			error_counter += (right2 - left2);

		}
	//	assert(it != non_aligned.end());
		else{
			std::cout << "pos "<< left2 << " lower bound " << it->second.first << std::endl;
			//assert(left2 >= it->second.first);
			if(left2 < it->second.first){
				std::cout<< "mapping error" <<std::endl;
			}else{
				std::cout << "correctly mapped"<<std::endl;
			}
		}
	}

}
void map_check::read_txt_file_long_center(std::ifstream & txtfile, std::string & sequence_of_interest){//Read the text file that come along with the fasta output of the graph
	std::cout << "reading text file "<<std::endl;
	std::string str;
	getline(txtfile, str);
	while(!txtfile.eof()){
		std::vector<std::string> parts;
		strsep(str, ":", parts);
		if(parts.size()==4){
			std::vector<std::string>seq;
			strsep(parts.at(1),"|", seq);
			assert(seq.size()==4);
			if(seq.at(3)=="NC_000913.3"){
				std::vector<std::string>position;
				strsep(parts.at(2),"r", position);
				assert(position.size()==2);
				std::cout<<"insert"<< parts.at(0)<< " : " << position.at(0) << " l " << position.at(1)<<std::endl;
				assert(stoi(position.at(0))<=stoi(position.at(1)));
			//	int length = std::stoi(position.at(1).substr(1,position.at(1).length()-1));
				int length = std::stoi(position.at(1)) - std::stoi(position.at(0));
				non_aligned.insert(std::make_pair(std::stoi(parts.at(0)), std::make_pair(std::stoi(position.at(0)),length)));
			}

		}
		getline(txtfile, str);
	}
}
void map_check::read_txt_file(std::ifstream & txtfile, std::string & sequence_of_interest){//Read the text file that come along with the fasta output of the graph
	std::cout <<"read text file "<<std::endl;
//	sequence_of_interest = "NC_000913.3";
	std::string str;
	getline(txtfile, str);
	while(!txtfile.eof()){
		std::vector<std::string> parts;
		strsep(str, ":", parts);
		assert(parts.size()==4);
		std::vector<std::string>seq;
		strsep(parts.at(1),"|", seq);
		assert(seq.size()==4);
		if(seq.at(3)=="NC_000913.3"){
			std::cout << parts.at(2) << " " << parts.at(3)<<std::endl;
			int length = std::stoi(parts.at(3)); //TODO change them to size_t
			non_aligned.insert(std::make_pair(std::stoi(parts.at(0)), std::make_pair(std::stoi(parts.at(2)),length)));
		}
		getline(txtfile, str);
	}
}

const std::map<unsigned int , std::vector<std::string> > & map_check::get_clusters()const{
	return clusters;
}
const std::multimap<std::string , std::pair<size_t, size_t> > & map_check::get_intervals()const{
	return boundries;
}
const std::map<int,std::pair<int, int> > & map_check::get_nonaligned()const{
	return non_aligned;
}
const unsigned int map_check::get_node_length(unsigned int & node)const{
	std::map<unsigned int, unsigned int >::const_iterator it = node_length.find(node);
	assert(it != node_length.end());
	return it->second;

}
void test_sim_reads_mapping::read_sim_maffile(std::ifstream & sim_maffile, std::map<std::string , std::pair<size_t , size_t> > & onreads){//Read the maf file from simpb software
	size_t count = 0;
	std::string str;
	getline(sim_maffile,str);
	while(!sim_maffile.eof()){
		std::vector<std::string> parts1;
		std::vector<std::string> parts2;
		if(str == "a"){
			std::string al1;
			std::string al2;
			
			getline(sim_maffile,al1);//actual genome
			getline(sim_maffile,al2);//read
			count++;
			strsep(al1," ",parts1);
			assert(parts1.size()==15);
			size_t begin1 = size_t(std::stoi(parts1.at(10)));
			std::string str1 = parts1.at(14);
			std::cout << "length "<< str1.length()<<std::endl;
			std::string nogap;
			for(size_t i = 0; i < str1.size() ; i++){
				if(str1.at(i) != '-'){
					nogap += str1.at(i);
				}
			}
			std::stringstream temp;
			temp<< parts1.at(11);
			size_t length;
			temp >> length;
			size_t end1 = length + begin1-1; 
			reads_position_on_seq.insert(std::make_pair(count, std::make_pair(begin1,end1)));

			strsep(al2," ",parts2);
			assert(parts2.size()==7);
			size_t begin2 = size_t(std::stoi(parts2.at(2)));
			size_t end2 = size_t(std::stoi(parts2.at(3)))-1;
			std::string str2 = parts2.at(6);

			size_t read_begin, read_end , ref_begin, ref_end;
			read_begin = begin2;
			read_end = end2;
			if(parts2.at(4)== "-"){//TODO check the length and start from the begin of negative strand!!
				read_begin = end2;
				read_end = begin2;
			}
			ref_begin = begin1;
			ref_end = end1;
			if(parts1.at(12)=="-"){//TODO
				std::cout << "neg seq"<<std::endl;
				ref_begin = end1;
				ref_end = begin1;
			}
			pw_alignment p(str1,str2,ref_begin, read_begin, ref_end, read_end, 0, count);
			alignments.insert(std::make_pair(count,p));
			onreads.insert(std::make_pair(nogap, std::make_pair(ref_begin,ref_end)));
		}
		getline(sim_maffile,str);
	}

	for(std::map<size_t,std::pair<size_t,size_t> >::iterator it = reads_position_on_seq.begin() ; it!= reads_position_on_seq.end() ; it++){
		std::cout << "read " << it->first << " : "<< std::endl;
		std::cout << "from "<< it->second.first << " to "<< it->second.second <<std::endl;

	}
	assert(count == 30775);
}
void test_sim_reads_mapping::this_part_position_on_seq(bool & direction, size_t & read, size_t & start_on_read, size_t & end_on_read, size_t & start, size_t & end,size_t & read_length){
		std::cout << "begin on read "<< start_on_read << " end on read "<< end_on_read <<std::endl;
		std::map<size_t, const pw_alignment>::iterator it = alignments.find(read);// from simpb maf file
		assert(it != alignments.end());

		const pw_alignment & p = it->second;
		p.print();
		size_t l,r;
		p.get_lr2(l,r);
		read_length = r- l + 1;
		if((p.getbegin1()<p.getend1() && p.getbegin2() > p.getend2())||(p.getbegin2()<p.getend2() && p.getbegin1() > p.getend1())){
			std::cout << "change the dir"<<std::endl;
			direction = true;
		}
		std::string ref1 = p.get_al_ref1();
		std::string ref2 = p.get_al_ref2();
		assert(ref1.length() == ref2.length());
		int counter1 = -1;
		int counter2 = -1;
		int column = -1;
		size_t begin_col = 0;
		size_t end_col = 0;
		if(p.getbegin2()<p.getend2()){
			for(size_t i =0; i < ref2.length(); i++){//On read
				column ++;
				if(ref2.at(i)!= '-'){
					counter2 ++;
				}
				if(counter2 == start_on_read){
					begin_col = column ;
				}
				if(counter2 == end_on_read){
					end_col = column;
					break;
				}
			}
		}else{
			std::string reverse;
			for(size_t i = ref2.length()-1 ; i > 0; i--){
				reverse += dnastring::complement(ref2.at(i));
			}
			reverse += dnastring::complement(ref2.at(0));
			for(size_t i =0; i < reverse.length(); i++){//On read
			//	std::cout << reverse.at(i);
				column ++;
				if(reverse.at(i)!= '-'){
				//	std::cout << "no gap "<< reverse.at(i) <<std::endl;

					counter2 ++;
				}else{
				//	std::cout << "Gap" << counter1 << " " << column <<std::endl;
				}
				if(counter2 == start_on_read){
					begin_col = column ;
				}
				if(counter2 == end_on_read){
					end_col = column;
					break;
				}
			}
		//	std::cout << " "<<std::endl;

		}
		std::cout << "begin col "<< begin_col << " end col "<< end_col <<std::endl;
		column = -1;
		if(p.getbegin2()<p.getend2()){

			for(size_t i =0; i < ref1.length(); i++){
				column ++;
				if(ref1.at(i)!= '-'){
					counter1 ++;
				}
				if(column == begin_col){
					if(ref1.at(i)=='-'){ //TODO +1 could make more sense
						if(counter1 > 0){
							start = counter1-1;
						}else{	
						//	std::cout << counter1 <<std::endl;
							assert(counter1 ==-1);
							start = 0;
						}
					}else{
						start = counter1;
					}
				}
				if(column == end_col){
					if(ref1.at(i)=='-'){
						if(counter1 > 0){
						//	end = counter1-1;
							end = counter1; //XXX Just changed it!
						}else{	
							assert(counter1 ==-1);
							end = 0;
						}
					}else{
						end = counter1;
					}
					break;
				}
			}
		}else{
			for(size_t i =0; i < ref1.length(); i++){
				column ++;
				if(ref1.at(i)!= '-'){
					counter1 ++;
				}
				if(column == ref1.length()-1 - end_col){
					std::cout << "pos "<< ref1.length() << " " << end_col << std::endl;
					if(ref1.at(i)=='-'){
						std::cout << "is on gap " << column <<std::endl;
						if(counter1 >= 0){
							start = counter1+1;
							std::cout << counter1 <<std::endl;
						}else{	
						//	std::cout << counter1 <<std::endl;
							assert(counter1 ==-1);
							start = 0;
						}
					}else{
						start = counter1;
					}
				}
				if(column == ref1.length()-1-begin_col){
					if(ref1.at(i)=='-'){
						std::cout << "the end is on gap " << column  << " " << counter1 <<std::endl;

						if(counter1 > 0){
						//	end = counter1-1;
							end = counter1; //XXX Just changed it
						}else{	
						//	assert(counter1 ==-1);
							end = 0;
						}
					}else{
						end = counter1;
					}
					break;
				}
			}
		}
		std::cout << "start "<<start << " end "<< end << std::endl;
}

void test_sim_reads_mapping::check_output(std::ifstream & mapping_maf){ //Read the mapping maf file
	if(mapping_maf) {
		std::string str;
		while(getline(mapping_maf, str)) {
			if(str.at(0)!='#') { //Skip headers
				if(str.at(0)=='a') { 
					std::string al1;
					std::string al2;
					getline(mapping_maf, al1);
					getline(mapping_maf, al2);
	
					//al2 is from read
					std::vector<std::string> parts2;
					strsep(al2," ",parts2);
					assert(parts2.size()==7);
					std::vector<std::string> readname_parts;
					strsep(parts2.at(1),":",readname_parts);
					assert(readname_parts.size()==2);
					size_t read_id = size_t(stoi(readname_parts.at(1)));
					std::cout << "read id "<< read_id <<std::endl;
					size_t begin_on_ref2 = size_t(stoi(parts2.at(2)));
					std::string ref2 = parts2.at(6);	
					size_t length = size_t(stoi(parts2.at(3)));
					size_t end_on_ref2 = length + begin_on_ref2 -1;
					std::map<size_t , std::pair<size_t , size_t> >::iterator it = reads_position_on_seq.find(read_id);
					assert(it != reads_position_on_seq.end());
					std::cout <<"position from the simpb maf file "<< it->second.first << " "<<it->second.second <<std::endl;
					size_t start_on_seq, end_on_seq;
					bool direction = false;
					size_t read_length;
					this_part_position_on_seq(direction, read_id, begin_on_ref2, end_on_ref2, start_on_seq, end_on_seq,read_length);
				//	if(direction == false){
						start_on_seq += it->second.first ;
						end_on_seq += it->second.first;
				//	}else{
				//		std::cout << "here! "<< std::endl;
				//		start_on_seq = it->second.second - start_on_seq ;
				//		end_on_seq = it->second.second - end_on_seq;

				//	}
					std::cout << "position on seq: " << start_on_seq << " "<< end_on_seq <<std::endl; 
					assert(end_on_seq <=it->second.second && start_on_seq <= it->second.second);
				//	if(read_id == 38) exit(0);
					//al1 is from the reference graph
					std::vector<std::string> parts1;
					strsep(al1," ",parts1);
					assert(parts1.size()==7);
					std::vector<std::string> refname_parts;
					strsep(parts1.at(1),":",refname_parts);
					assert(refname_parts.size()==2);
					unsigned int ref_id = stoi(refname_parts.at(1));
					size_t begin_on_ref1 = size_t(std::stoi(parts1.at(2)));
					std::string ref1 = parts1.at(6);	
					size_t end_on_ref1 = begin_on_ref1 + size_t(std::stoi(parts1.at(3))) -1; 
					std::cout<< begin_on_ref1 << " "<<end_on_ref1<<std::endl; //Begin and end on the center itself. We might need to check the corresponding member!
					std::map<unsigned int , std::vector<std::string> > clusters = mp_check.get_clusters();
					std::cout << "centers size "<< clusters.size()  << " id "<< ref_id <<std::endl;
					std::map<unsigned int , std::vector<std::string> >::iterator center = clusters.find(ref_id);	
					if(center != clusters.end()){
						std::stringstream curseq;
						curseq<<"nc000913.3"<<':'<<ref_id; //XXX too specific
						std::cout << curseq.str()<<std::endl;
						bool IsAMember = false;
						std::cout << center->second.size() <<std::endl;
						for(size_t i = 0 ; i < center->second.size(); i++){
							std::cout << center->second.at(i) <<std::endl;
							if(center->second.at(i)==curseq.str()){
								IsAMember = true;
								break;
							}
						}
						if(IsAMember == false){
							std::cout<<"mapping error!"<<std::endl;
						}else{ //Should check and see if begin and end from the read are in these boundries. //TODO how to check if they were not shifted?
	
							std::multimap<std::string , std::pair<size_t, size_t> > boundries = mp_check.get_intervals();
							std::multimap<std::string , std::pair<size_t, size_t> >::iterator member = boundries.find(curseq.str());
							assert(member != boundries.end()); 
							std::pair<std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator, std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator> it1=boundries.equal_range(curseq.str());
							bool mapped = false;
							for(std::multimap<std::string , std::pair<size_t, size_t> >::const_iterator intervals= it1.first ; intervals!=it1.second; intervals++){
								std::pair<size_t,size_t>this_pair = intervals->second;
								std::cout << "from "<< this_pair.first << " to "<< this_pair.second <<std::endl;
								size_t start = start_on_seq;
								size_t end = end_on_seq;
							//	if(start > end){
							//		start = end_on_seq;
							//		end = start_on_seq;
							//	}
								if(start>=this_pair.first && end<= this_pair.second){
									std::cout << "correct!"<<std::endl;
									mapped = true;
									int pos = -1;
									find_position_on_member(ref_id, begin_on_ref1,end_on_ref1,this_pair.first , this_pair.second, pos);

									std::cout << start << " - " <<this_pair.first << "  == " << pos <<std::endl;;
									if(this_pair.first + pos != start){
										if(begin_on_ref2==0 && end_on_ref2<400){
											std::cout << "begin_shifted"<<std::endl;
										}
										else if(begin_on_ref2 >= read_length -400 && end_on_ref2 == read_length-1){
											std::cout<<"end_shifted"<<std::endl;
										}
										else{
											int shift = start - this_pair.first - pos;
											std::cout << "Shifted!" << shift <<std::endl;  //TODO check if they are begin or end of a component!
 
										}
									}
									break;
								}
							}
							if(mapped == false){
								std::cout << "mapping error!"<<std::endl;
							}
						}
					}else{//TODO
						std::cout << "mapped to a non-aligned region! "<<std::endl;
						std::map<int,std::pair<int, int> > non_aligned = mp_check.get_nonaligned();
						std::map<int,std::pair<int, int> >::iterator it = non_aligned.find(ref_id);
						if(it == non_aligned.end()){
							std::cout <<"mapping error on nonaligned"<<std::endl; //Means it mapped to the non aligned region from another genome sequence.
						}else if(start_on_seq >= it->second.first && end_on_seq < it->second.first + it->second.second){//Check the borders! //TODO check it more precise!
							std::cout << "correct! on non aligned"<<std::endl;
							size_t this_begin = begin_on_ref1 + it->second.first;
							size_t this_end = end_on_ref1 + it->second.first;
							if(start_on_seq == this_begin && end_on_seq == this_end){
								std::cout << "precise"<<std::endl;
							}else{
								std::cout << "Shifted! on non aligned "<<std::endl;
								std::cout << start_on_seq << "!= " <<  this_begin << " or " << end_on_seq <<" != "<< this_end<<std::endl;
							}
						}else{
							std::cout << " non aligned Shifted!" <<std::endl; //TODO! check
						}

					}
					

				}
			}
		}

	}

}

void test_sim_reads_mapping::find_position_on_member(unsigned int & center, size_t & begin_on_center , size_t & end_on_center, size_t & begin_on_mem , size_t & end_on_mem, int & position){
	std::multimap<size_t , const pw_alignment> al_from_graph = mp_check.get_als();
	std::multimap<size_t , const pw_alignment>::iterator it = al_from_graph.find(center);
	assert(it != al_from_graph.end());
	std::pair<std::multimap<size_t , const pw_alignment >::const_iterator, std::multimap<size_t , const pw_alignment >::const_iterator> it1=al_from_graph.equal_range(center); //From output maf file of ref-graph
	for(std::multimap<size_t , const pw_alignment >::const_iterator als= it1.first ; als!=it1.second; als++){
		pw_alignment p = als->second;
		size_t l2,r2;
		p.get_lr2(l2,r2);

		p.print();
		//Check the coordinats of the second reference:
		if(l2 >= begin_on_mem && r2 << end_on_mem && p.getreference2()==0){
			std::string ref1 = p.get_al_ref1();
			std::string ref2 = p.get_al_ref2();
			size_t column = p.alignment_length();
			int counter = -1;
			int counter2 = -1;
			for(size_t i =0; i < ref1.length(); i++){
				if(ref1.at(i)!='-'){
					counter ++;
				}else{
					std::cout << "a gap on node "<<std::endl;
				}
				if(counter == begin_on_center){
					column = i;
					break;
				}
			}
			std::cout << "till "<< column <<std::endl;
			assert(column < p.alignment_length());
			for(size_t i = 0; i <= column; i++){
				if(ref2.at(i)!='-'){
					counter2 ++;
				}
				else{
					std::cout << "a gap on ref after pos " << counter2 << " at column "<< i <<std::endl;
				}
			}
			position = counter2;
			std::cout << "pos on ref " << position <<std::endl;
			assert(position>=0);
			break;
		}
	}
}
void test_reveal::read_gfa(std::ifstream & gfa){
	std::string str;
	size_t from = 0;
	size_t to = 0;
	while(getline(gfa,str)){
		if(str.at(0)=='S'){
			std::vector<std::string> parts;
			strsep(str,"\t",parts);
			assert(parts.size()==7);
			std::stringstream ststream(parts.at(1));
			size_t id;
			ststream >> id;
			std::vector<std::string> ref_parts;
			strsep(parts.at(4),":",ref_parts);
			to = from + parts.at(2).length()-1;
			if(ref_parts.at(2).at(0)=='0'){
				ref_graph_nodes.insert(std::make_pair(id,std::make_pair(from,to)));
				from = to + 1;
			}
		}

	}
	for(std::map<size_t , std::pair<size_t,size_t> >::iterator it = ref_graph_nodes.begin() ; it != ref_graph_nodes.end() ; it++){
		std::cout<< "on "<< it->first << " from "<< it->second.first << " to " << it->second.second << std::endl;
	}

}
void test_reveal::read_the_result(std::ifstream & mapping_maf){ //Read the mapping maf file
	if(mapping_maf) {
		std::string str;
		while(getline(mapping_maf, str)) {
			if(str.at(0)!='#') { //Skip headers
				if(str.at(0)=='a') { 
					std::string al1;
					std::string al2;
					getline(mapping_maf, al1);
					getline(mapping_maf, al2);
	
					//al1 is from ref
					std::vector<std::string> parts1;
					strsep(al1," ",parts1);
					assert(parts1.size()==7); 
				/*	if(parts1.size() != 7){
						std::cout << "not 7 "<<std::endl;
						std::cout<< parts1 << std::endl;
						continue;
					}*/
					std::vector<std::string> refname_parts;
					strsep(parts1.at(1),":",refname_parts);
					assert(refname_parts.size()==2);
					std::stringstream sstream(refname_parts.at(1));
					size_t ref_id;
					sstream >> ref_id;
				//	std::string name_and_dir = refname_parts.at(1);
				//	name_and_dir += parts1.at(4);
				//	nodes.push_back(name_and_dir);

					//al1 is from read
					std::vector<std::string> parts2;
					strsep(al2," ",parts2);
					assert(parts2.size()==7); 
					std::vector<std::string> readname_parts;
					strsep(parts2.at(1),":",readname_parts);
					std::string sub = readname_parts.at(1).substr(4);
					std::stringstream sstream3(sub);
					size_t read_id;
					sstream3 >> read_id;
					std::cout << "read id " << read_id << std::endl;
					std::stringstream sstream1(parts2.at(2));
					size_t read_pos;
					sstream1 >> read_pos;
				//	std::cout << "read pos "<<read_pos<<std::endl;
					std::stringstream sstream2(parts2.at(3));
					size_t al_length;
					sstream2 >> al_length;
					size_t from = read_id*30030 + read_pos;
					size_t to = read_id*30030 + read_pos + al_length-1;
					from_mapping_output.insert(std::make_pair(ref_id, std::make_pair(from,to)));
					std::cout<< "on "<< ref_id << " from "<< from << " to " << to << std::endl;

				}
			}
		}
	}
//	std::cout<< "-------------------------------"<<std::endl;
//	for(std::map<size_t , std::pair<size_t,size_t> >::iterator it = from_mapping_output.begin() ; it != from_mapping_output.end() ; it++){
//		std::cout<< "on "<< it->first << " from "<< it->second.first << " to " << it->second.second << std::endl;
//	}

}
void test_reveal::compare_with_reveal(){
	size_t total_error_base = 0;
	int error_begin = 0;
	int error_end = 0;
	bool ERROR = false;
	for(std::multimap<size_t , std::pair<size_t,size_t> >::iterator it = from_mapping_output.begin() ; it != from_mapping_output.end() ; it++){
		std::cout<< "from mapping on "<< it->first << " from "<< it->second.first << " to " << it->second.second << std::endl;
		std::map<size_t , std::pair<size_t,size_t> >::iterator it1 = ref_graph_nodes.find(it->first);
		if(it1 != ref_graph_nodes.end()){
			std::cout<< "from gfa file "<< it1->first << " from "<< it1->second.first << " to " << it1->second.second << std::endl;

		}else{
			std::cout << "ERROR! on ref"<<std::endl;
			ERROR = true;
			error_begin = it->second.second - it->second.first + 1;
			error_end = 0;
		}
		if(it->second.first>= it1->second.first && it->second.first <= it1->second.second){
			
		}else{
			std::cout << "ERROR!"<<std::endl;
			if(ERROR == false){
				ERROR = true;
				if(it->second.first < it1->second.first){
					error_begin = it->second.first - it1->second.first;
					error_end = 0;
					error_begin = std::abs(error_begin);
				}
				else{
					assert(it->second.first > it1->second.second);
					error_end = it->second.second-it1->second.second;
					error_begin = 0;
					error_end = std::abs(error_end);
				}
				std::cout << "b "<< error_begin << " e "<< error_end <<std::endl;
			}

		}
		if(it->second.second>= it1->second.first && it->second.second <= it1->second.second){

		}else{//TODO
			std::cout << "ERROR! at end"<<std::endl;
			if(ERROR == false){
				ERROR = true;
				error_begin = it->second.second - it->second.first + 1;
				error_end = 0;
			}

		}
		if(ERROR == true){
		//	total_error_base += (error_end - error_begin);
			std::cout << "b1 "<< error_begin << " e1 "<< error_end <<std::endl;			
			total_error_base += (error_end + error_begin);
			ERROR = false;
		}
	}
	std::cout << "total error base "<< total_error_base <<std::endl;


}
void test_reveal::compare_with_path(std::ifstream & path){//XXX This is a quick and dirty test!
	//path is reveal output path.
	std::string str;
	getline(path, str);
	getline(path,str);
	assert(str.at(0) != 'n');
	std::vector<std::string> path_parts;
	strsep(str,",",path_parts);
	std::cout << "path size "<< path_parts.size()<<std::endl;
	for(size_t i =0; i < path_parts.size() ; i++){
		std::cout << path_parts.at(i)<<std::endl;
	}
	std::cout << "---------------"<<std::endl;
	for(size_t i =0; i < path_parts.size() ; i++){
		if(i<nodes.size()){
			if(path_parts.at(i)== nodes.at(i)){
				std::cout << "correct"<<  path_parts.at(i) << " "<< nodes.at(i) <<std::endl;
			}else{
				std::cout << "wrong! " << path_parts.at(i) << " "<< nodes.at(i) <<std::endl;
			}
		}else{
			std::cout<< "ERROR: size of two container are not the same"<<std::endl;
			break;
		}

	}


}
