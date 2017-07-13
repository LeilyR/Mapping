#include <iostream>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <ostream>
#include <vector>

#include "pw_alignment.hpp"
#include "data.hpp"
#include "dynamic_mc.hpp"
#include "model.hpp"
#include "alignment_index.hpp"
#include "intervals.hpp"

//Aligning reads:
#include "simple_NW_algo.hpp"
#include "needleman_wunsch.hpp"
#include "dijkstra.hpp"
#include "ref_graph.hpp"
#include "map_read.hpp"
#include "test_mapping.hpp"



#define VERSION "0.0.1" 
#define NUMVERISON 0
#define MIN_READ_VERSION 0

#define TIMING 1


#define TEST 0



/**
General TODO

Look at diversity within clusters. Can this diversity be used to differentiate between orthologous and paralogous groups,
should they be treated differently in mapping? Some more stringent approach to place repeat instances? Discuss this with Felix again?

Maybe another round of AP-clustering after all reads are aligned? Maybe include not only reads, but also all origingal sequences in that cluster
maybe do multiple alignments?

**/

void usage() {
	std::cerr << "Mapping program" << std::endl;
	std::cerr << "Usage:" << std::endl;
}

int do_fasta_prepare(int argc, char * argv[]) {

	if(argc < 4) {
		usage();
		std::cerr << "Program fasta_prepare" << std::endl; 
		std::cerr << "Combine fasta files from different species/accessions into a common file" << std::endl;
		std::cerr << "Parameters: " << std::endl;
		std::cerr << "output file -- fasta file" << std::endl;
		std::cerr << "<list of fasta files> -- input file, one file per species/accession, no colon in file name"<< std::endl;
		return 1;
	}


	std::ofstream outf(argv[2]);
	if(outf) {
		std::vector<std::ifstream*> inputs;
		std::vector<std::string> names;
		for(int i=3; i<argc; ++i) {
			std::string istr(argv[i]);
			std::vector<std::string> slashparts;
			strsep(istr, "/", slashparts);
			std::string afterslash = slashparts.at(slashparts.size()-1);
			std::vector<std::string> periodparts;
			strsep(afterslash, ".", periodparts);
			std::string name = periodparts.at(0);
			if(periodparts.size()>1) {
				std::stringstream sstr;
				for(size_t j=0; j<periodparts.size()-1; j++) {
					if(j!=0) sstr << ".";
					sstr << periodparts.at(j);
				}
				name = sstr.str();
			}
			for(size_t j=0; j<name.length(); ++j) {
				if(':' == name.at(j)) {
					std::cerr << "Error: " << argv[i] << " file name contains a colon" << std::endl;
					exit(1);
				}
			}
			std::ifstream * iin = new std::ifstream(argv[i]);
			if(!(*iin)) {
				std::cerr << "Error: cannot read: " << argv[i] << std::endl;
				exit(1);
			}
			inputs.push_back(iin);
			names.push_back(name);
		}

		for(size_t i=0; i<inputs.size(); ++i) {
			std::string str;
			size_t j = 0;
			while(getline((*(inputs.at(i))), str)) {
				if(str.length()>0) { // ignore empty lines
				if(str.at(0)=='>') {
					std::string nname = str.substr(1);
					std::stringstream nhead;
				//	nhead << ">" << i << ":" << j;
					nhead << ">" << names.at(i) << ":" << nname;
					outf << nhead.str() << std::endl;
				//	outf << nhead.str();
					j = j +1;

						
						
/*
	for(size_t col = 0; col < alignment_length(); col++) {
		char c1;
		char c2;
		alignment_col(col, c1, c2);
		std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
	if(col>50) break;	
	}
*/	

				} else {
				//	outf << str;
					outf << str << std::endl;
				}
				}
			}
			inputs.at(i)->close();
			delete inputs.at(i);
		}
		outf.close();
	} else {
		std::cerr << "Error: cannot write: " << argv[2] << std::endl;
		exit(1);
	}

	outf.close();

	return 0;
}

/*
    src -- The name of one of the source sequences for the alignment. For sequences that are resident in a browser assembly, the form 'database.chromosome' allows automatic creation of links to other assemblies. Non-browser sequences are typically reference by the species name alone.
    start -- The start of the aligning region in the source sequence. This is a zero-based number. If the strand field is '-' then this is the start relative to the reverse-complemented source sequence.
    size -- The size of the aligning region in the source sequence. This number is equal to the number of non-dash characters in the alignment text field below.
    strand -- Either '+' or '-'. If '-', then the alignment is to the reverse-complemented source.
    srcSize -- The size of the entire source sequence, not just the parts involved in the alignment.
    text -- The nucleotides (or amino acids) in the alignment and any insertions (dashes) as well.
*/

void write_maf_record(std::ostream & out, const std::string & src, size_t start, size_t size, char strand, size_t srcSize, const std::string & alignment_part) {
	out << "s " << src << " " << start << " " << size << " " << strand << " " << srcSize << " " << alignment_part << std::endl;
}


/*
	write half a pairwise alignment to a maf file 
	(this write the s-line only)

*/
void write_maf_record(std::ostream & out, const all_data & data, const pw_alignment & al, size_t reference) {
	assert(reference==0 || reference==1);
	// get the requested half alignment
	size_t print_seq;
	size_t print_start;
	size_t print_end;
	std::string print_al;
	// TODO print size not end
	if(reference==0) {
		print_seq = al.getreference1();
		print_start = al.getbegin1();
		print_end = al.getend1();
		print_al = al.get_al_ref1();
	} else {
		print_seq = al.getreference2();
		print_start = al.getbegin2();
		print_end = al.getend2();
		print_al = al.get_al_ref2();
	}
	// translate reference sequence index number to corresponding long name (accession:contig)
	size_t acc = data.accNumber(print_seq);
	std::string accname = data.get_acc(acc);
	std::string seqname = data.get_seq_name(print_seq);
	std::stringstream write_longname;
	write_longname << accname << ':' << seqname;
	char strand = '+';
	size_t size = data.get_seq_size(print_seq);
	size_t al_on_ref_size = print_end - print_start + 1;
	if(print_end < print_start) {
		al_on_ref_size = print_start - print_end + 1;
		strand = '-';
		print_start = size - print_start -1;


		
	} 
	write_maf_record(out, write_longname.str(), print_start, al_on_ref_size, strand, size, print_al);
	

}
std::string get_reverse(const std::string & seq){
	std::stringstream temp;
	std::vector<std::string> center_parts;
	strsep(seq, ":" , center_parts);
	unsigned int dir = std::atoi(center_parts.at(0).c_str());
	unsigned int ref = std::atoi(center_parts.at(1).c_str());
	unsigned int left = std::atoi(center_parts.at(2).c_str());
	if(dir == 0){
		temp<< 1 <<":"<<ref<<":"<<left;
	}else{
		temp << 0 << ":"<<ref<<":"<<left;
	}
	return temp.str();
}

int make_fasta_file(int argc, char* argv[]){
	if(argc< 4){
		std::cerr<< "Program: make_fasta"<<std::endl;
		std::cerr<< "Parameters: "<<std::endl;
		std::cerr << "* input: txt file"<<std::endl;
		std::cerr << "* output: fasta file"<<std::endl;
		return 1;
	}
	std::string textfile(argv[2]);
	std::string fastafile(argv[3]);
	all_data data;
	std::ifstream in(textfile.c_str());
	std::string str;
	std::stringstream curseq;
	while(getline(in, str)){
		curseq<<str;
	}
	std::vector<std::string> sequence;
	sequence.push_back(curseq.str());
	std::ofstream outs(fastafile.c_str());
	data.make_fasta(outs,sequence);
	return 0;
}
int make_long_reads(int argc, char * argv[]){
	if(argc < 3){
		std::cerr<< "Program: sim_perfect_reads"<<std::endl;
		std::cerr<< "Parameters: "<<std::endl;
		std::cerr << "* input: genome fasta file "<<std::endl;
		std::cerr << "* output: reads fasta file " << std::endl;
		return 1;
	}
	std::string inputfastafile(argv[2]);
	std::string outputfastafile(argv[3]);

	std::ifstream fastain(inputfastafile.c_str());
	std::map<std::string , std::string> contigs; // first string is the name, the second one is the sequence content
	if(fastain) {
		std::string str;
		std::stringstream curseq;
		std::string name;
		std::string pre_name;
		std::cout<<"check point "<<std::endl;
		while(getline(fastain, str)){
			if(str.at(0)=='>'){ //Separate with | keep the first two parts as its name
				std::vector<std::string> parts;
				strsep(str,"|",parts);
				name = parts.at(1);
				if(curseq.str().size()!=0){
					contigs.insert(std::make_pair(pre_name,curseq.str()));
					curseq.str("");
				}else{
					std::cerr << "Warning: A read of length zero is skipped " <<std::endl;
				}
			}
			else {
				curseq << str;
				pre_name = name;
			}
		}
		// store last reas
		if(curseq.str().size()!=0){
			contigs.insert(std::make_pair(pre_name,curseq.str()));
		}else{
			std::cerr << "Warning: A read of length zero is skipped " <<std::endl;
		}
		curseq.str("");
		fastain.close();
	} else {
		std::cerr << "Error: cannot read " << std::endl;
		exit(1);
	}
	std::multimap<std::string, std::string> reads; //name, content
	for(std::map<std::string,std::string>::iterator it = contigs.begin() ; it != contigs.end() ; it++){
		//Break the sequence into shorter pieces
		std::string contig = it->second;
		if(contig.length() > 15000){
			std::cout << "long contig"<< it->first << std::endl;
			for(size_t i = 0; i < contig.length()/15000 ; i++){
				std::string this_read = contig.substr(i*15000, (i+1)*15000-1);
				std::stringstream name;
				name << it->first << "|" << i*15000<<"|"<<(i+1)*15000-1;
				reads.insert(std::make_pair(name.str(),this_read));
			}
			std::string this_read = contig.substr((contig.length()/15000)*15000, contig.length()-1);
			std::stringstream name;
			name<< it->first << "|"<< (contig.length()/15000)*15000 << "|" << contig.length()-1;
			reads.insert(std::make_pair(name.str(), this_read));

		}else{
			std::stringstream name;
			name<<it->first<< "|" << 0 << "|" << contig.length()-1;
			reads.insert(std::make_pair(name.str(),contig));
		}
	}
	size_t count = 0;
	std::ofstream fastaout(outputfastafile.c_str());

	for(std::multimap<std::string, std::string>::iterator it = reads.begin() ; it != reads.end() ; it++){
		size_t rest = 0;
		fastaout<<">"<<count<< "|" << it->first <<std::endl;
		count++;
		if(it->second.length() > 70){
			for(size_t i = 0; i < it->second.length()-70; i++){
				for(size_t j = i; j < i+70;j++){
					fastaout<<it->second.at(j);
				}
				i +=69;
				rest = i+1;
				fastaout<<std::endl;
			}
			size_t length = 0;
			for(size_t i = rest ; i < it->second.length();i++){
				fastaout<<it->second.at(i);
				length++;
			}
			fastaout<<std::endl;
			assert(length <= 70);
		}else{
			fastaout<<it->second<<std::endl;
		}
	}


}
int do_read_data(int argc, char * argv[]) {
	if(argc < 4) {
		usage();
		std::cerr << "Program: read_data" << std::endl;
		std::cerr << "Parameters:" << std::endl;
		std::cerr << "* fasta file from fasta_prepare" << std::endl;
		std::cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		return 1;

	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);

//	std::string samfile(argv[3]);

// Read all data
	all_data data;
	typedef overlap overlap_type;
	data.read_fasta_maf(fastafile, maffile);

//	data.read_accknown_fasta_sam(fastafile,samfile);
//	std::vector<pw_alignment> alignments = data.getAlignments();
	size_t seq_length = 0;
	for(size_t i =0; i < data.numSequences();i++){
		seq_length += data.getSequence(i).length();
	}
	std::cout << "total base number is " << seq_length <<std::endl;
	return 0;
}

int gfa_to_fasta(int argc, char * argv[]){
	if(argc < 4){
		std::cerr << "Program: gfa_to_fasta"<<std::endl;
		std::cerr << "Parameters:"<<std::endl;
		std::cerr<< "Input: GFA file"<<std::endl;
		std::cerr <<"Output: fasta file"<<std::endl;
		std::cerr<<"Output: text file"<<std::endl;
		return 1;

	}
	std::string gfafile(argv[2]);
	std::string fastafile(argv[3]);
	std::string textfile(argv[4]);
	std::ofstream fastaout(fastafile.c_str());
	std::ofstream txtout(textfile.c_str());
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
				fastaout<<">"<<node.at(1)<<std::endl;
				fastaout<<node.at(2)<<std::endl;
			}
			if(line[0]=='L'){
				
			}
			if(line[0]=='P'){
				
				std::vector<std::string> node;
				strsep(line, "\t" , node);
				txtout<<node.at(1)<<std::endl;
				std::vector<std::string> parts;
				strsep(node.at(2), "M" , parts);
				txtout<<parts.at(0)<< std::endl;
				break; // XXX Cus in my test case i only need the first path
			}
		}
	}
}

int do_read_map(int argc, char * argv[]){
	/*
	1. Gets read fasta, graph fasta, alignment maf file and read+graph fasta as an input
	2. For each read makes alignment graph and detect the best path on the graph.
	3. Note that all the fasta files are output of fasta_prepare
	4. Also since LAST is used for creating pairwise alignment, I already know that reads are always on the second reference. This should be checked in case of using any other aligner.
	*/

	typedef dynamic_mc_model use_model;
//	typedef overlap_interval_tree<located_alignment> overlap_type;
	if(argc < 9){ 
		std::cerr << "Program: map_read" << std::endl;
		std::cerr << "Parameters:" << std::endl;
		std::cerr << "* Ref(graph) gfa file ('nogfa' if not using gfa)" << std::endl;
		std::cerr << "* Ref(graph) dot file ('nodot' if not using dot)" << std::endl;
		std::cerr << "* Read+Ref fasta file" << std::endl; //after running fasta_prepare on the original ref and read file(before fasta_prepare)
		std::cerr << "* alignment maf file "<<std::endl; //I used LAST (make alignemnt between Ref.fasta and Read.fasta)
		std::cerr << "* output maf file " << std::endl; 
		std::cerr << "genomes and reads fasta "<<std::endl;
		std::cerr << "genomes and reads maf file "<< std::endl;
		std::cerr << "* num_threads(optional, defult is 1)" <<std::endl;
		return 1;
	}
	size_t num_threads = 1;
	if(argc ==11){
		num_threads = std::atoi(argv[9]);
	}
	std::string refgfa(argv[2]);
	std::string refdotfile(argv[3]);
	std::string allfasta(argv[4]);
	std::string alignmentmaf(argv[5]);
	std::string outputmaf(argv[6]);
	std::ofstream output(outputmaf.c_str());
	output << "##maf version=1 " << std::endl;
	output << "# Result of mapping over a graph " << std::endl;
	output << "#" << std::endl;
	std::string genomes_and_reads_fasta(argv[7]);
	std::string genomes_and_reads_maf(argv[8]);
	
/*	all_data data1; //Just to test if model returns more reliable cost!
	data1.read_fasta_maf_forward_read_only(genomes_and_reads_fasta, genomes_and_reads_maf);
	std::cout<< "data 1 numacc "<<data1.numAcc()<<std::endl;
	wrapper wrap;
	use_model m(data1,wrap, num_threads); //TODO take the wrapper out from the model!
	m.train();
//	model m(data1);
//	m.train();*/
	//--------------------
	all_data data;
//TODO read a fasta file containing genomes have been used in making graph and reads, alignments between full genomes and reads. Consider all the genomes as on accession and all the reads as an other accession
	data.read_fasta_maf_forward_read_only(allfasta,alignmentmaf);//All the alignments are read in a way the 'read' reference (second ref) of it is forward.
	wrapper wrap;
	use_model m(data,wrap, num_threads); //TODO take the wrapper out from the model!
//	model m(data);
	m.train();

	//TODO anthoer data object, read allfasta and alignmentmaf in here. use this object from this point on

	std::cout << "numseq: "<<data.numSequences()<<std::endl;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_with_pos_gain; 
	std::string ref_accession;
	const pw_alignment & al = data.getAlignment(0);
	assert(al.getbegin2() < al.getend2());
	size_t ref_acc = data.accNumber(al.getreference1());
	std::cout << "ref accession is " << ref_acc <<std::endl;
	ref_accession = data.get_acc(ref_acc);
#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < data.numAlignments();i++){
		std::cout << "al " << i << std::endl;
		const pw_alignment & al = data.getAlignment(i);
		assert(al.getbegin2() < al.getend2());
		size_t read_acc = data.accNumber(al.getreference2());
	//	if(i == 0){
	//		size_t acc = data.accNumber(al.getreference1());
	//		std::cout << "ref accession is " << acc <<std::endl;
	//		ref_accession = data.get_acc(acc);
	//	}
		double g1 ,g2;
		std::cout << "ref1 " <<al.getreference1() << " numAcc "<< data.numAcc() << std::endl;
	
	//	m.gain_function(al, ref_acc, read_acc, g1,g2);
		double av_gain = (g1+g2)/2 ;
		if(av_gain > 0){
#pragma omp critical(al_insert)
{
			al_with_pos_gain.insert(&al);
}
		}
	}

//	std::map< std::string, size_t> longname2seqidx;
//	longname2seqidx = data.getLongname2seqidx();
//	std::cout<<longname2seqidx.size()<<std::endl;
	ref_graph rgraph(data);
	if(refdotfile != "nodot"){
		std::cout<< "read dot file "<<std::endl;
		rgraph.read_dot_file(refdotfile,ref_accession);//Fill in the vertix and edge container with the reference graph nodes and adjacencies
	}
	if(refgfa != "nogfa"){
		std::cout<< "read gfa file "<<std::endl;
		rgraph.read_gfa_for_adj(refgfa);
	}
	std::cout << "after gfa "<<std::endl;
	int test_node = 31888;
	std::set<int> temp_adj = rgraph.get_adjacencies(test_node);
	for(std::set<int>::iterator it = temp_adj.begin() ; it != temp_adj.end() ; it++){
		std::cout << *it <<std::endl;
	}
	size_t acc;
	size_t num = 0;
	std::vector<std::vector<pw_alignment> > all_als;
	for(size_t i =0; i< data.numAlignments();i++){//Read all the als and fill in all_als vector with als of each read
		const pw_alignment & p = data.getAlignment(i);
	//	p.print();
//	for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = al_with_pos_gain.begin() ; it != al_with_pos_gain.end() ;it++){
//		const pw_alignment * p = *it;
		size_t ref2 = p.getreference2(); 
		if(i == 0){
			size_t acc2 = data.accNumber(ref2);
			all_als.resize(data.getAcc(acc2).size());
			std::cout<< "all_als size is "<<all_als.size()<<std::endl;
		}
		acc = data.accNumber(p.getreference1());
		if(acc == 1){
			all_als.at(ref2).push_back(p);	
		}else{
			assert(acc == 0);
			size_t num = data.get_seq_number_of_acc(acc);
			std::cout<< "num "<< num<< "ref2 "<< ref2 <<std::endl;
			all_als.at(ref2-num).push_back(p);
		}
	}
	for(size_t i =0; i < all_als.size();i++){
		std::cout <<"read i  "<< i <<std::endl;
		std::vector<pw_alignment> als = all_als.at(i);
		std::cout<<als.size()<<std::endl;
	}
//Goes through all the reads and for each of them find the best path using Dijkestra graph, the best path should be save in the output file
	std::map<int, std::set<int> > adjacencies = rgraph.get_adjacencies(); //Reference graph
//	deep_first df(data, adjacencies);//Deep first search algorithm on reference graph nodes to finde sub graph of length MAXGAP
//	for(size_t i =0; i < all_als.size();i++){
//	size_t i = 227;
	for(size_t i =0; i < all_als.size();i++){
		std::cout <<"read i  "<< i  <<std::endl;
		std::vector<pw_alignment> alignments = all_als.at(i);
	/*	std::vector<pw_alignment> alignments;
		std::map<double , pw_alignment> ordered_als;
		for(size_t j =0; j < all_als.at(i).size(); j++){
			const pw_alignment al = all_als.at(i).at(j);
			double g1 ,g2;	
			m.gain_function(al,g1,g2);
			double av_gain = (g1+g2)/2 ;
			ordered_als.insert(std::make_pair(av_gain, al));
		}
		size_t counter = 0;
		for(std::map<double, pw_alignment>::reverse_iterator it = ordered_als.rbegin() ; it != ordered_als.rend(); it++){
			if(counter <= 200){
				alignments.push_back(it->second);
			}else break;
			counter ++;
		}*/
		dnastring current_read = data.getSequence(i);
		size_t readacc = data.accNumber(i);
		std::cout << "als size "<< alignments.size() <<std::endl;
		if(alignments.size() != 0){
			compute_overlapped_reads_avl connected_comps(alignments, data.numSequences(), num_threads);
			std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > no_fraction_ccs;
			connected_comps.compute_on_second_ref(no_fraction_ccs); //We are interested in overlaps only on second reference of an al
			std::cout << "number of components is " << no_fraction_ccs.size() <<std::endl;
			std::cout<<no_fraction_ccs.at(0).size()<<std::endl;
			als_components<use_model> ccs(data, rgraph, m, no_fraction_ccs);//check the distance between the last right and the first left of the next component
			ccs.merg_components();
			std::cout<< "components are merged! "<<std::endl;
		//	std::cout<< "sequences are 2: "<<data.get_seq_name(2)<<" 3: "<<data.get_seq_name(3)<<" 8: "<<data.get_seq_name(8)<< " 17: "<< data.get_seq_name(17) << " 18: " <<data.get_seq_name(18)<<std::endl;
		//	std::cout << data.get_seq_name(193)<< " " << data.get_seq_name(303) <<std::endl;
			ccs.find_als_on_paths(output,acc,readacc);//It also calls dijkstra algorithm inside
		}
//		exit(0); //XXX Temp
	}
	std::cout << "done! "<<std::endl;




	return 0;
}
int check_mapping_result(int argc , char * argv[]){//Use it to test the mapping result over the reads i made out of one certain accession among the existing one on the graph
	if(argc<7){
		std::cerr<<"Program: test_mapping" <<std::endl;
		std::cerr << "Parameters:" << std::endl;
		std::cerr<<"maf file from reference graph"<<std::endl;//graph.maf
		std::cerr<<"maf file from mapping"<<std::endl;
		std::cerr<<"accession of interest"<<std::endl; //Name of the sequence that I made reads out of it
		std::cerr<<"maf file alignment between the reference and reads"<<std::endl;
		std::cerr<< "* text file contains center's info"<<std::endl;
		return 1;
	}
	std::string graphmaffile(argv[2]);
	std::string mapmaffile(argv[3]);
	std::string sequence(argv[4]);
	std::string alignments(argv[5]);
	std::string centerids(argv[6]);

	std::ifstream als(alignments.c_str());
	std::ifstream graphmaf(graphmaffile.c_str());
	std::ifstream mappingoutput(mapmaffile.c_str());
	std::ifstream centeridtxt(centerids.c_str());
	map_check mc;
//	mc.read_txt_file_long_center(centeridtxt,sequence);

	mc.read_txt_file(centeridtxt,sequence);
	mc.read_graph_maf(graphmaf);//XXX This is a bit too specific since I use it for the maffile i produce and know that each accession includes only one full length genome rather than couple of contiges. 
//	exit(1);
	mc.read_alignments(als); //alignments between reference and reads
	std::cout<< "als were read"<<std::endl;
	mc.check_output(mappingoutput,sequence);
//	mc.print_nodes();
//	mc.print_graph();

	return 0;
}
int check_sim_mapping_result(int argc, char * argv[]){
	if(argc<6){
		std::cerr<< "Program: test_sim_mapping "<<std::endl;
		std::cerr<< "Parameters: "<<std::endl;
		std::cerr<< " maffile from simulation program "<<std::endl;
		std::cerr<< " maffile from mapping reads against the graph "<<std::endl;
		std::cerr<< "maffile contains referecne graph "<< std::endl;
		std::cerr<< "textfile contains center's info"<<std::endl;
		std::cerr<< "gfafile contains ref graph " << std::endl;
		return 1;
	}
	std::string simmaf(argv[2]);
	std::ifstream maffile(simmaf.c_str());
	std::string mapmaf(argv[3]);
	std::ifstream map_maffile(mapmaf.c_str());
	std::string graph_maf(argv[4]);
	std::ifstream graphmaf(graph_maf.c_str());
	std::string centersinfo(argv[5]);
	std::ifstream centers(centersinfo.c_str());
	std::string refgfa(argv[6]);
	std::ifstream gfain(refgfa.c_str());
	test_reveal test;	
	test.read_gfa(gfain); //Read the reference 
	std::map<size_t, std::pair<size_t,size_t> > nodes_on_ref_graph;
	nodes_on_ref_graph = test.get_ref_graph_nodes();
/*	std::ifstream in("fastatest");//Only center of the cluster 24, read from fasta file containing both nodes and reads
	char c;
	std::string fasta;
	in.get(c);
	fasta +=c;
	while(!in.eof()){
		in.get(c);
		fasta +=c;
	}
	std::cout<< fasta.length()<<std::endl;
	map_check mc;
	mc.read_graph_maf(graphmaf);

	std::ifstream in1("maftest");//includes the member of cluster 24 that contains the sequence of interest from graph.maf
	char c1;
	std::string maf;
	in1.get(c1);
	maf +=c1;
	while(!in1.eof()){
		in1.get(c1);
		maf +=c1;
	}
	std::cout<< maf.length()<<std::endl;
	std::cout <<"node 24 from 35659 to 36799 "<<std::endl;
	for(size_t i = 35659; i<= 36799;i++){
		std::cout << maf.at(i);
	}
	std::cout << " "<<std::endl;

	std::ifstream in2("testnc");//Full nc sequence
	char c2;

	std::string seq;
	in2.get(c2);
	seq +=c2;
	while(!in2.eof()){
		in2.get(c2);
		seq +=c2;
	}
	std::cout<< seq.length()<<std::endl;
	std::cout << "seq at 2332525 and 2332526 "<< seq.at(2332525) << seq.at(2332526) <<std::endl;
	std::cout<< "at 4189250 "<< seq.at(4189250) << seq.at(4189251) <<std::endl;
	std::cout<< "at 4153591 "<< seq.at(4153591) << seq.at(4153592) <<std::endl;
	std::cout<< "at 4153591+35659 "<< seq.at(4153591+35659) << seq.at(4153592+35660) <<std::endl;
	std::cout<< "at 4153591+35658 "<< seq.at(4153591+35658) << seq.at(4153592+35659) <<std::endl;
	std::cout << "node at 35658 " << maf.at(35658)<< maf.at(35659)<<std::endl;
	for(size_t i =0; i < maf.length() ; i++){
		if(maf.at(i) != seq.at(i+4153591)){
			std::cout << "warning at "<< i <<std::endl;
			std::cout << maf.at(i) << " "<< seq.at(i+4153591)<<std::endl;//TODO it has to be fixed in creating graph!!
		}
	//	assert(maf.at(i) == seq.at(i+4153591));
	}
	std::ifstream in12("cluster12");
	std::string center12;
	in12.get(c);
	center12 +=c;
	while(!in12.eof()){
		in12.get(c);
		center12 +=c;
	}
	std::cout << "center12 " << center12.at(13672)<<center12.at(13673)<<center12.at(13674) <<std::endl;*/
//	exit(0);
	std::map<std::string , std::pair<size_t , size_t> > onreads;
	map_check mc;
//	mc.read_graph_maf(graphmaf);
	std::string sequence; //TODO
//	mc.read_txt_file(centers,sequence);
	test_sim_reads_mapping test_map(mc);
	test_map.read_sim_maffile(maffile, onreads);//read tha maf file from simpb 
//compare onreads with the original seq from testnc:
	for(std::map<std::string , std::pair<size_t,size_t> >::iterator it = onreads.begin() ; it != onreads.end() ; it++){//onreads contains informations from sim.maf
		std::string this_part = it->first;
		size_t from = it->second.first;
		size_t to = it->second.second;
		std::cout << "f " << from << " t "<< to <<std::endl;
	/*	if(from < to){
			for(size_t j = from ; j <=to ; j++){
				assert(seq.at(j)==this_part.at(j-from));
			}
		}else{//compare it with reverse complement of the seq
			std::string reverse;
			for(size_t i = from ; i > to ; i--){
				reverse+=dnastring::complement(seq.at(i));
			}
			reverse+=dnastring::complement(seq.at(to));

			for(size_t j = to ; j <=from ; j++){
				assert(reverse.at(j)==this_part.at(j));
			}

		}*/
	}
	std::cout << "position and content were checked" << onreads.size() <<std::endl;
	test_map.check_output(map_maffile,nodes_on_ref_graph);

}
int do_test_reveal(int argc , char* argv[]){
	if(argc< 3){
		std::cerr<< "Program: test_reveal"<<std::endl;
		std::cerr<< "*input: mapping maf file"<<std::endl;
		std::cerr << "*input: reveal GFA" <<std::endl;
		return 1;
	}

	std::string maffile(argv[2]);
	std::string gfafile(argv[3]);
	std::ifstream mafin(maffile.c_str());
	std::ifstream gfain(gfafile.c_str());
	test_reveal test;
	
	test.read_gfa(gfain);
	std::cout << "gfa was read! "<<std::endl;
	test.read_the_result(mafin);
	test.compare_with_reveal();

}
int check_reads(int argc, char* argv[]){
	if(argc < 4){
		std::cerr << "Program: check_read" << std::endl;
		std::cerr << "* fasta file from fasta_prepare" << std::endl;
		std::cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		std::cerr << "read fasta"<<std::endl;
		return 1;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string readfasta(argv[4]);
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	data.read_fasta(readfasta);
	std::ofstream choose_read;
	choose_read.open ("thisRead.txt");
	std::string first_read = data.get_read(0);
	std::cout<<"read lenght "<< first_read.length()<<std::endl;
	std::string this_read = data. get_read(6);
	size_t sequence_counter = 6;
	choose_read.write(reinterpret_cast<char*>(&sequence_counter),sizeof(uint32_t));
	for(size_t c = 0; c< this_read.length(); c++){
		choose_read<< this_read.at(c);
	}
	choose_read.close();
	std::ifstream thisread("thisRead.txt");
	size_t thisref = 0;
	size_t from = 180180;
	size_t to = 210209;

	data.compare_seq_with_read(thisread, thisref ,from, to);
}
int check_graph_out_put(int argc , char* argv[]){
	if(argc < 6){
		usage();
		std::cerr << "Program: check_graph" << std::endl;
		std::cerr << "* fasta file from fasta_prepare" << std::endl;
		std::cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		std::cerr << "* fasta file contains graphs nodes "<<std::endl;
		std::cerr << "* text file contain nodes info. " << std::endl;
		std::cerr << "* maf file contains graph output. "<< std::endl;
		return 1;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string nodesfasta(argv[4]);
	std::string textfile(argv[5]);
	std::string graphmaffile(argv[6]);
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	std::ifstream txtfile(textfile.c_str());
	std::map<size_t, std::vector<std::pair<unsigned int , size_t> > > nodes;
	std::map<std::pair<unsigned int, size_t> , size_t > rights;
	std::string str;
	getline(txtfile, str);
	while(!txtfile.eof()){
		std::vector<std::string> parts;
		strsep(str, ":", parts);
		for(size_t i = 1 ; i < parts.size()-1;i++){
			unsigned int seqid = 6;
			if(parts.at(i)=="gi|556503834|ref|NC_000913.3|"){
				seqid = 0;
			}
			if(parts.at(i)=="gi|15829254|ref|NC_002695.1|"){
				seqid = 1;
			}
			if(parts.at(i)=="gi|170079663|ref|NC_010473.1|"){
				seqid = 2;
			}
			if(parts.at(i)=="gi|254160123|ref|NC_012967.1|"){
				seqid = 3;
			}
			if(parts.at(i) =="gi|26111730|gb|AE014075.1|"){
				seqid = 4 ;
			}
			assert(seqid < 5);
			std::string position = parts.at(i+1);
			std::vector<std::string> possep;
			strsep(position,"r",possep);
			assert(possep.size()==2);
			size_t left = size_t(std::stoi(possep.at(0)));
			size_t right = size_t(std::stoi(possep.at(1)));
			std::map<size_t, std::vector<std::pair<unsigned int , size_t> > >::iterator it = nodes.find(size_t(std::stoi(parts.at(0))));
			if(it ==nodes.end()){
				nodes.insert(std::make_pair(size_t(std::stoi(parts.at(0))), std::vector<std::pair<unsigned int , size_t> >()));
				 it = nodes.find(size_t(std::stoi(parts.at(0))));
			}
			it->second.push_back(std::make_pair(seqid,left));
			rights.insert(std::make_pair(std::make_pair(seqid,left),right));
			i++;
		}
		getline(txtfile, str);
	}
	std::map<size_t , std::string> longcenters;
	std::ifstream node(nodesfasta.c_str());
	getline(node,str);
	while(!node.eof()){//Read each node check if the content is equal to what is on the actual sequence
		assert(str.at(0)=='>');
		std::cout << str<<std::endl;
		std::string temp;
		temp = str.substr(1, str.length()-1);
		std::stringstream sstream(temp);
		size_t id;
		sstream >> id;
		std::cout << "node "<< id <<std::endl;
		std::string content;//from graph
		getline(node,str);
		while(str.at(0) != '>' && !node.eof()){
			std::cout << str.length() << " "<< content.length() << std::endl;
			content.append(str);
			if(id == 99 && str.length() ==37){
				std::cout << "HERE!!"<<std::endl;
				break;
			}
			getline(node,str);
		}
		longcenters.insert(std::make_pair(id, content));
		std::string sequence;
		std::map<size_t, std::vector<std::pair<unsigned int , size_t> > >::iterator it = nodes.find(id);
		assert(it != nodes.end());
		for(size_t j =0; j < it->second.size() ; j++){
			size_t left = it->second.at(j).second;
			unsigned int refid = it->second.at(j).first;
			std::map<std::pair<unsigned int, size_t> , size_t >::iterator it1 = rights.find(it->second.at(j));
			assert(it1 != rights.end());
			size_t right = it1->second;
			std::string this_seq;
			if(left< right || (left==right)){
				this_seq = data.extract_seq_part(refid, left, right);
			}else{
				assert(right< left);
				this_seq = data. extract_reverse_seq_part(refid,right,left);
			}
			sequence.append(this_seq);

		}
		for(size_t i = 0; i < sequence.size() ;i++){
		//	std::cout << "at "<< i << " on seq " << sequence.at(i) << " on the node "<< content.at(i) <<std::endl;
			assert(sequence.at(i)==content.at(i));
		}
		if(id == 99 && str.length() ==37){
			break;
		}

	
	}
	std::cout<<"compare full length sequences: "<<std::endl;
	std::ifstream nc("testnc");
	char c;
	std::string nc_str;
	nc.get(c);
	while(!nc.eof()){
		nc_str += c;
		nc.get(c);
	}
	unsigned int refId = 0;
	size_t start = 0;
	size_t end = data.getSequence(0).length()-1;
	std::string nc_seq = data.extract_seq_part(refId,start,end);
	assert(nc_seq.length() == nc_str.length());
	for(size_t i =0 ; i < nc_str.size();i++){
		assert(nc_str.at(i)== nc_seq.at(i));
	}
	std::cout<< "comparison is done!"<<std::endl;
	std::ifstream graphmaf(graphmaffile.c_str());
	getline(graphmaf, str);
	while(str.at(0)=='#'){
		getline(graphmaf,str);
	}
	assert(str.at(0)=='a');
	while(!graphmaf.eof()){//Only read the first line of each al and check the content!
		std::string str1;
		if(str.at(0)=='a'){
			std::vector<std::string> name_split;
			strsep(str," ",name_split);
			assert(name_split.size()==2);
			std::string temp = name_split.at(1).substr(8,name_split.at(1).length()-1);
			std::stringstream sstr;
			sstr<< temp;
			size_t clusterid;
			sstr >> clusterid;
			getline(graphmaf,str1);
			std::vector<std::string> split;
			strsep(str1," ",split);
			assert(split.size()==7);
			std::string centerseq = split.at(6);
			std::string nogap;
			for(size_t i = 0; i < centerseq.size();i++){
				if(centerseq.at(i)!='-'){
					nogap +=centerseq.at(i);
				}
			}
			if(split.at(4)=="-"){
				//GET the reverse_compl
				std::string reverse;
				for(size_t i = nogap.size(); i > 0; i --){
					reverse += dnastring::complement(nogap.at(i-1));
				}
				nogap = reverse;
			}
			std::map<size_t , std::string>::iterator it = longcenters.find(clusterid);
			assert(it != longcenters.end());
			std::cout << "id "<< clusterid << " " << nogap.length() << " "<< it->second.length() << std::endl;
			assert(nogap.length() == it->second.length());
			for(size_t i =0; i < nogap.length();i++){
				if(nogap.at(i)!=it->second.at(i)){
					std::cout<<"at "<< i << " : " << nogap.at(i) << " " << it->second.at(i) << std::endl;
				}
				assert(nogap.at(i)==it->second.at(i));
			}
		}
		getline(graphmaf,str);
	}
}
int needleman_wunsch(int argc, char * argv[]){
	if(argc < 4){
		usage();
		std::cerr << "Program: needleman_test" << std::endl;
		std::cerr << "* fasta file from fasta_prepare" << std::endl;
		std::cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;

		
		return 1;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);

	size_t num_threads = 1;
//	std::string ref = "ACTGCTGATTC";
//	std::string read = "GACTGGATTC";
	std::string ref = "GAAGCTGCTGAAAAAGCCAA";
	std::string read = "GAAGCTGCTGAAAA";

	all_data data;
	data.read_fasta_maf(fastafile, maffile);

	wrapper wrap;
//	std::ofstream outs("testoutput");
	dynamic_mc_model m(data,wrap, num_threads);
//	model m(data);
	m.train();
//	m.write_parameters(outs);

	needleman<dynamic_mc_model> n(data, m, read, ref);
//	simpleNW n(m,read,ref);
	size_t readid = 1;
	size_t refid = 0;
	n.compute_matrix(readid, refid);
	n.print_score_matrix();
	n.print_path();
	size_t type=3;
	std::string readstr;
	std::string refstr;
	n.find_the_best_path(type, readstr, refstr);



	std::map<int, std::set<int> > adjacencies;
//	adjacencies.insert(make_pair(1,2));
//	adjacencies.insert(make_pair(1,3));
//	adjacencies.insert(make_pair(3,4));
//	adjacencies.insert(make_pair(3,5));
//	adjacencies.insert(make_pair(4,5));
//	adjacencies.insert(make_pair(4,6));
//	adjacencies.insert(make_pair(4,7));




//	deep_first df(data, adjacencies);

	return 0;
}

#if !TEST
	
int main(int argc, char * argv[]) {

	if(argc < 2) {
		usage();
		exit(1);
	}
	std::string program(argv[1]);

	if(0==program.compare("fasta_prepare")) {
		return do_fasta_prepare(argc, argv);
	}
	else if(0==program.compare("make_fasta")){
		return make_fasta_file(argc,argv);
	}
	else if (0== program.compare("read_data")){
		return do_read_data(argc,argv);
	}
	else if(0==program.compare("gfa_to_fasta")){
		return gfa_to_fasta(argc,argv);
	}
	else if(0==program.compare("map_read")){
		return do_read_map(argc,argv);
	}
	else if(0==program.compare("needleman_test")){
		return needleman_wunsch(argc,argv);
	}
	else if(0==program.compare("test_mapping")){
		return check_mapping_result(argc,argv);
	}
	else if(0 == program.compare("test_reveal")){
		return do_test_reveal(argc,argv);
	}
	else if(0 == program.compare("test_sim_mapping")){
		return check_sim_mapping_result(argc,argv);
	}
	else if(0== program.compare("sim_perfect_reads")){
		return make_long_reads(argc,argv);
	}
	else {
		usage();
	}

	
	return 1;
}



#endif
#include "needleman_wunsch.cpp"
#include "dynamic_mc.cpp"
#include "pw_alignment.cpp"
#include "intervals.cpp"
#include "alignment_index.cpp"
#include "map_read.cpp"




















