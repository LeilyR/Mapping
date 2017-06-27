#ifndef TEST_MAPPING_CPP
#define TEST_MAPPING_CPP

#include "test_mapping.hpp"

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

void test_sim_reads_mapping::check_output(std::ifstream & mapping_maf, std::map<size_t, std::pair<size_t,size_t> > & nodes_on_ref_graph){ //Read the mapping maf file
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
					std::cout << "position of this part of read on seq: " << start_on_seq << " "<< end_on_seq <<std::endl; 
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
					std::cout<<"on "<< ref_id << " : "<< begin_on_ref1 << " "<<end_on_ref1<<std::endl; //Begin and end on the center itself. We might need to check the corresponding member!
				//	compare_with_my_graph(ref_id, start_on_seq , end_on_seq, begin_on_ref1,end_on_ref1, begin_on_ref2, end_on_ref2, read_length);
					size_t temp_ref_id = ref_id;
					compare_with_reveal_graph(nodes_on_ref_graph, temp_ref_id);	
				}
			}
		}

	}

}
void test_sim_reads_mapping::compare_with_my_graph(unsigned int & ref_id, size_t & start_on_seq, size_t & end_on_seq,size_t & begin_on_ref1, size_t & end_on_ref1, size_t & begin_on_ref2, size_t & end_on_ref2, size_t & read_length){
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
						std::cout << start << " - " <<this_pair.first << "  == " << pos <<std::endl;
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
void test_sim_reads_mapping::compare_with_reveal_graph(std::map<size_t, std::pair<size_t,size_t> > & nodes_on_ref_graph, size_t & ref_id){
	std::map<size_t, std::pair<size_t,size_t> >::iterator it = nodes_on_ref_graph.find(ref_id);
	if(it == nodes_on_ref_graph.end()){
		std::cout<< "mapping error! "<<std::endl;
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
	size_t error = 0;
	for(std::multimap<size_t , std::pair<size_t,size_t> >::iterator it = from_mapping_output.begin() ; it != from_mapping_output.end() ; it++){
		std::cout<< "from mapping on "<< it->first << " from "<< it->second.first << " to " << it->second.second << std::endl;
		std::map<size_t , std::pair<size_t,size_t> >::iterator it1 = ref_graph_nodes.find(it->first);
		if(it1 != ref_graph_nodes.end()){
			std::cout<< "from gfa file "<< it1->first << " from "<< it1->second.first << " to " << it1->second.second << std::endl;
			if((it->second.first>= it1->second.first && it->second.first <= it1->second.second) && (it->second.second>= it1->second.first && it->second.second <= it1->second.second)){

			}else if(it->second.second < it1->second.first || it->second.first > it1->second.second){
				std::cout << "ERROR!"<<std::endl;
				error = it->second.second - it1->second.first +1;
				total_error_base += error;
			}else if(it->second.first < it1->second.first && it->second.second <= it1->second.second){
				std::cout << "ERROR!"<<std::endl;
				error = it1->second.first - it->second.first;
				total_error_base += error;
			}else if(it->second.first >= it1->second.first && it->second.second > it1->second.second){
				std::cout << "ERROR!"<<std::endl;
				error = it->second.second - it1->second.second + 1;
				total_error_base += error;
			}else{
				assert(it->second.first < it1->second.first && it->second.second > it1->second.second);
				std::cout << "ERROR!"<<std::endl;
				error = it1->second.first - it->second.first + 1;
				error += it->second.second - it1->second.second + 1;
				total_error_base += error;
			}

		}else{
			std::cout << "ERROR! on ref"<<std::endl;
			error = it->second.second - it->second.first + 1;
			total_error_base += error;
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

#endif
