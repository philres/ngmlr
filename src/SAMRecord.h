/*
 * SAMrecord.h
 *
 *  Created on: Apr 4, 2014
 *      Author: fritz
 */

#ifndef SAMRECORD_H_
#define SAMRECORD_H_
#include <string.h>
#include <stdint.h>
#include "api/BamMultiReader.h"
//#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;
class SAMRecord {

private:

	int mapped_flag;
	string chr;
	int mapping_pos;
	int pair_number;

	vector<CigarOp> cigar;
	string read_name;
	uint32_t mq;
	string seq;
	string qualities;
	string tags;

	vector<CigarOp> convert_cigar(string CIGAR) {
		int num = -1;
		vector<CigarOp> tmp;
		for (size_t i = 0; i < CIGAR.size(); i++) {
			if (num == -1) {
				num = atoi(&CIGAR[i]);
			} else if (num != -1 && atoi(&CIGAR[i]) == 0 && CIGAR[i] != '0') {
				CigarOp op;
				op.Length = num;
				op.Type = CIGAR[i];
				tmp.push_back(op);
				num = -1;
			}
		}
		return tmp;
	}

public:
	SAMRecord() {
		mapping_pos = -1;
		chr = "NOT VALID";
		read_name = "NOT SET";
		mq = -1;
		mapped_flag = 4; //read is not mapped
	}
	~SAMRecord() {

	}
// ############# Setter: ################
	void set_mapped_flag(int flag) {
		this->mapped_flag = flag;
	}
	void set_chr(string chr) {
		this->chr = chr;
	}
	void set_mapping_pos(int mapping_pos) {
		this->mapping_pos = mapping_pos;
	}
	void set_CIGAR(vector<CigarOp> CIGAR) {
		this->cigar = CIGAR;
	}
	void set_CIGAR(string CIGAR) {
		this->cigar = convert_cigar(CIGAR);
	}
	void set_read_name(string name) {
		this->read_name = name;
	}
	void set_mapping_quality(uint32_t mq) {
		this->mq = mq;
	}
	void set_sequence(string sequence) {
		this->seq = sequence;
	}
	void set_qualities(string qualities) {
		this->qualities = qualities;
	}
	void set_tags(string tags) {
		this->tags = tags;
	}

	void set_pair_number(int pair_number) {
		this->pair_number = pair_number;
	}

// ############# Getter: ################
	int get_mapped_flag() {
		return this->mapped_flag;
	}
	string get_chr() {
		return this->chr;
	}
	int get_mapping_pos() {
		return this->mapping_pos;
	}
	string get_CIGAR_str() {

		ostringstream ss;
		for (size_t i = 0; i < this->cigar.size(); i++) {
			ss << this->cigar[i].Length;
			ss << this->cigar[i].Type;
		}
		return string(ss.str());
	}
	vector<CigarOp> get_CIGAR() {
		return this->cigar;
	}
	string get_read_name() {
		return this->read_name;
	}
	uint32_t get_mapping_quality() {
		return this->mq;
	}
	string get_sequence() {
		return this->seq;
	}
	string get_qualities() {
		return this->qualities;
	}
	string get_tags() {
		return this->tags;
	}

	int get_pair_number() {
		return this->pair_number;
	}

// ############# Query Flag: ################

	bool is_mapped() {
		return !(this->get_mapped_flag() & 0x4);
	}
	bool is_proper_pair() {
		return !(this->get_mapped_flag() & 0x2);
	}
	bool is_plus_strand() {
		return (this->get_mapped_flag() & 0x10);
	}

	bool is_mate_mapped() {
		return (this->get_mapped_flag() & 0x8);
	}
	bool is_mate_plus_strand() {
		return (this->get_mapped_flag() & 0x20);
	}
	bool is_paired() {
		return (this->get_mapped_flag() & 0x1);
	}
	bool is_primaryAlignment() {
		return (this->get_mapped_flag() & 0x100);
	}
	bool is_second_mate() {
		return (this->get_mapped_flag() & 0x40);
	}
};

#endif /* SAMRECORD_H_ */
