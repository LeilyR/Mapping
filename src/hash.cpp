#ifndef HASH_CPP
#define HASH_CPP



#include "hash.hpp"

void readmap::find_best_match(){
	size_t n = read.length()/(mismatches+1);//There is at least  a perfect fragment
	for(size_t i =0; i < read.length() ;i += n){
		std:: string seq = read.substr(i,n);
		calculate_first_kmers(seq);
	}

}
void readmap::calculate_first_kmers(std::string & seq){
	assert(seq.length() >= 2*K-1);
	std::map<std::string , std::set<std::string> > edges = kgraph.get_edges();
	for(size_t i =0; i < K ; i++){
		std::string subseq = seq.substr(i, K);
		std::map<std::string , std::set<std::string> >::iterator it = edges.find(seubseq);
		if(it != edges.end){
			//TODO  a perfect fragment is found! Make alignment
			make_fitting_alignment(subseq,position);
		}
	}


}
