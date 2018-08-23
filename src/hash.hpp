#ifndef HASH_HPP
#define HASH_HPP

#include "ref_graph.hpp"
#include <cassert>

class Hash{

	public:
		Hash(std::string key){
			std::hash<std::string> hash_value;
			value = hash_value(key);
		}
		const int get_hvalue()const{
			return value;
		}
	private:
		int value;

};


class readmap{
	public:
		readmap(const Kgraph & kg, std::string & read , size_t & error_rate):kgraph(kg){
			this->read = read;
			this->error_rate = error_rate;
			mismatches = read.length()*error_rate/100;
		}
		void find_best_match();
		void calculate_first_kmers(std::string);
	private:
		const Kgraph & kgraph;
		std::string read;
		size_t error_rate;
		size_t mismatches;


};
#endif
