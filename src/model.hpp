#ifndef MODEL_HPP
#define MODEL_HPP

#include "pw_alignment.hpp"
#include"data.hpp"
#include "dynamic_mc.hpp"


//Interval tree:
#include "alignment_index.hpp"

#include <map>
#include <vector>
#include <cassert>
#include <omp.h>

class compute_overlapped_reads_avl {
	public:
	compute_overlapped_reads_avl(const std::vector<pw_alignment> & als_in, size_t num_sequences, size_t num_threads):alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
	//	RIGHT = 0;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
		for(size_t  i =0; i < als_in.size();i++){
			const pw_alignment * al = &als_in.at(i);
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::vector<const pw_alignment*> sorted_als(sorter.size());
		size_t num = 0;
		for(typename std::multimap<size_t, const pw_alignment*>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
		//	al->print();
		//	std::cout << "add one"<<std::endl;
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	}
	~compute_overlapped_reads_avl(){
		delete alind;
	}
	//This two are used for mapping reads against a graph
	void compute_on_second_ref(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & ccs){
		std::cout << "compute CC on " << alignments.size() << std::endl;
		std::set <const pw_alignment*, compare_pointer_pw_alignment > seen;
		std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
		for(typename std::set<const pw_alignment *, compare_pointer_pw_alignment >::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		//	std::cout << "compute_cc" <<std::endl;
			const pw_alignment * al = *it;
		//	al->print();
			typename std::set< const pw_alignment*, compare_pointer_pw_alignment >::iterator seenal = seen.find(al);
	//		std::cout << " seen " << seen.size() << std::endl;
			if(seenal == seen.end()) {
			//	std::cout << " getcc" << std::endl;
				std::set< const pw_alignment*, compare_pointer_pw_alignment > cc;
				get_cc_on_second_ref(*al, cc, seen);
			//	std::cout << "FOUND CC size " << cc.size() << std::endl;
			//	for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::const_iterator it = cc.begin(); it != cc.end();it++){
			//		const pw_alignment * test = *it;
			//		test->print();
			//	}
				sorter.insert(std::make_pair(cc.size(), cc));
			}	
		}
		for(typename std::multimap<size_t , std::set<const pw_alignment *, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
			ccs.push_back(it->second);
		}

	}
	void get_cc_on_second_ref(const pw_alignment & al , std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment *, compare_pointer_pw_alignment> & seen) {
		std::vector<size_t> left(2);
		std::vector<size_t> right(2);
		al.get_lr1(left.at(0), right.at(0));
		al.get_lr2(left.at(1), right.at(1));
		std::vector<size_t>reference(2);
		reference.at(0) = al.getreference1();
		reference.at(1) = al.getreference2();
		cc_step_on_second_ref(&al, reference.at(1), left.at(1), right.at(1), cc, seen);
	}
	void cc_step_on_second_ref(const pw_alignment * al,size_t current, size_t left, size_t right,std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc,std::set <const pw_alignment *,compare_pointer_pw_alignment>  & seen){
		seen.insert(al);
		cc.insert(al);
		std::vector<const pw_alignment *> result;
	
		alind->search_overlap(current, left, right,result);
		for(size_t i =0; i < result.size();i++){
			const pw_alignment * this_result = result.at(i);
			size_t ref = this_result->getreference2();
			size_t l,r;
			this_result->get_lr2(l,r);
			std::set <const pw_alignment* ,compare_pointer_pw_alignment>::iterator it=seen.find(this_result);
			if(it == seen.end()){
				cc_step_on_second_ref(this_result, ref, l, r, cc, seen);
			}
		}
	}

	private:
	std::set<const pw_alignment*, compare_pointer_pw_alignment> alignments; 
	alignment_index  * alind;
	size_t num_threads;

};

#include "intervals.cpp"


#endif



