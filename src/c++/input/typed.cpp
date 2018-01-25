#include <algorithm>
#include <iterator>
#include <vector>
#include <stdint.h>
#include "global/common.h"
#include "input/typed.h"
using namespace std;

namespace ONETOOL {

void make_typed_string(vector<char> &out, const string &_in, bool typed)
{
	vector<char> tmp_vector;
	out.resize(0);

	if (_in == "." || _in == " " || _in == "") {
		if (typed == false)
			return;

		int8_t tmp = (int8_t)0;
		tmp = tmp << 4;
		tmp = tmp | (int8_t)7;
		out.push_back(tmp);

		return;
	}

	if (typed == true) {
		if (_in.length() >= 15) {
			int8_t tmp = (int8_t)15;
			tmp = tmp << 4;
			tmp = tmp | (int8_t)7;
			out.push_back(tmp);

			make_typed_int(tmp_vector, (int)_in.length(), typed);
			out.insert(out.end(), tmp_vector.begin(), tmp_vector.end());
		} else {
			int8_t tmp = (int8_t)_in.length();
			tmp = tmp << 4;
			tmp = tmp | (int8_t)7;
			out.push_back(tmp);
		}
	}
	out.reserve(out.size()+_in.size());
	copy(_in.begin(), _in.end(), back_inserter(out));
}

void make_typed_int(vector<char> &out, const int &_in, bool typed)
{
	vector<char> tmp_char;
	out.resize(0);

	int type;
	int8_t size_type = (int8_t)1;
	if (_in < 127 && _in >-127)
		type = 1;
	else if (_in < 32767 && _in>-32767)
		type = 2;
	else
		type = 3;

	make_int(tmp_char, _in, type);

	if (typed == true) {
		size_type = size_type << 4;
		size_type = size_type | type;
		out.push_back(size_type);
	}
	out.insert(out.end(), tmp_char.begin(), tmp_char.end());
}

void make_typed_string_vector(vector<char> &out, const vStr &_in, int number)
{
	vector<char> tmp_char;
	int max_val = 0;
	int8_t size_type;
	out.resize(0);

	if (number == -1) {
		for (unsigned int ui=0; ui<_in.size(); ui++) {
			if ((int)_in[ui].size() > max_val)
				max_val = (int)_in[ui].size();
		}
	} else
		max_val = number;

	if (max_val < 15) {
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)7;

		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)7;

		out.push_back(size_type);

		make_typed_int(tmp_char, max_val, true);
		out.insert(out.end(), tmp_char.begin(), tmp_char.end());
	}

	for (unsigned int ui=0 ; ui<_in.size() ; ui++)
		for (unsigned int uj=0 ; (int)uj<max_val ; uj++)
			if (_in[ui] == ".")
				out.push_back('\0');
			else if (uj<_in[ui].size())
				out.push_back(_in[ui][uj]);
			else
				out.push_back('\0');
}

void make_typed_GT_vector(vector<char> &out, vStr &_in)
{
	vector<char> tmp_vector;
	int8_t size_type;
	int max_ploidy = 0;
	out.resize(0);

//FIXME	max_ploidy = *max_element(ploidy.begin(), ploidy.end());

	if (max_ploidy < 15) {
		size_type = (int8_t)max_ploidy;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back(size_type);

		make_typed_int(tmp_vector, max_ploidy, true);
		out.insert(out.end(), tmp_vector.begin(), tmp_vector.end());
		tmp_vector.resize(0);
	}

	for (unsigned int ui=0; ui<_in.size(); ui++) {
//FIXME		encode_genotype(tmp_vector, _in[ui], max_ploidy);
		out.insert(out.end(), tmp_vector.begin(), tmp_vector.end());
		tmp_vector.resize(0);
	}
	out.insert(out.end(), tmp_vector.begin(), tmp_vector.end());
}

void encode_genotype(vector<char> &out, string &_in, int exp_size)
{
	int8_t tmp_int = 0;
	int8_t phased = 0;
	out.resize(exp_size);
	int idx = 0;

	for (unsigned int ui=0; ui<_in.length(); ui++) {
		if (_in[ui] =='|')
			phased = 1;
		else if (_in[ui] == '/')
			phased = 0;
		else {
			if (_in[ui] != '.') {
//FIXME				tmp_int = header::str2int(_in.substr(ui, 1));
				tmp_int++;
				tmp_int = tmp_int << 1;
				tmp_int = tmp_int | phased;
			} else
				tmp_int = (int8_t)0x80;

			out[idx] = (int8_t)tmp_int;
			idx++;
		}
	}
	while (idx<exp_size) {
		out[idx] = (int8_t)0x81;
		idx++;
	}
}

void make_typed_int_vector(vector<char> &out, const string &_in, int number)
{
	vector<char> tmp_char;
	vector<int> tmp_ints;
	vStr split_string;
	int converted = 0, type;
	int8_t size_type;
	unsigned int max = 0;
	unsigned int max_val = 0;
	out.resize(0);

	if (_in == " " || _in == "." || _in == "") {
		size_type = (int8_t)0;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back(size_type);
		return;
	}

//FIXME	header::tokenize(_in, ',', split_string);
	if (number == -1) {
		if (split_string.size() > max_val)
			max_val = (wsUint)split_string.size();
	} else
		max_val = number;

	for (unsigned int ui=0; ui<max_val; ui++) {
		if (ui<split_string.size()) {
//FIXME			converted = header::str2int(split_string[ui], 0x80000000);

			if ((abs(converted) >(int)max) && (converted != (int)0x80000000))
				max = abs(converted);
		} else
			converted = 0x80000001;

		tmp_ints.push_back(converted);
	}

	if (max < 127)
		type = 1;
	else if (max < 32767)
		type = 2;
	else
		type = 3;

	if (max_val < 15) {
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back(size_type);

		make_typed_int(tmp_char, max_val, true);
		out.insert(out.end(), tmp_char.begin(), tmp_char.begin());
	}

	for (unsigned int ui=0; ui<tmp_ints.size(); ui++) {
		make_int(tmp_char, tmp_ints[ui], type);
		out.insert(out.end(), tmp_char.begin(), tmp_char.end());
	}
}

void make_typed_int_vector(vector<char> &out, const vStr &_in, int number)
{
	vector<char> tmp_char;
	vector<int> tmp_ints;
	vStr split_string;
	int converted = 0, type;
	int8_t size_type;
	unsigned int max = 0;
	unsigned int max_val = 0;
	out.resize(0);

	if (number == -1) {
		unsigned int tmp_int = 0;
		for (unsigned int ui=0; ui<_in.size(); ui++) {
			tmp_int = (wsUint)count(_in[ui].begin(), _in[ui].end(), ',');

			if (tmp_int > max_val)
				max_val = tmp_int;
		}
		max_val++;
	} else
		max_val = number;

	for (unsigned int ui=0; ui<_in.size(); ui++) {
//FIXME		header::tokenize(_in[ui], ',', split_string);
		for (unsigned int uj=0; uj<max_val; uj++) {
			if (uj<split_string.size()) {
//FIXME				converted = header::str2int(split_string[uj], 0x80000000);

				if ((abs(converted) >(int)max) && (converted != (int)0x80000000))
					max = abs(converted);
			} else
				converted = 0x80000001;

			tmp_ints.push_back(converted);
		}
	}

	if (max < 127)
		type = 1;
	else if (max < 32767)
		type = 2;
	else
		type = 3;
	if (max_val < 15) {
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back(size_type);

		make_typed_int(tmp_char, max_val, true);
		out.insert(out.end(), tmp_char.begin(), tmp_char.begin());
	}

	for (unsigned int ui=0; ui<tmp_ints.size(); ui++) {
		make_int(tmp_char, tmp_ints[ui], type);
		out.insert(out.end(), tmp_char.begin(), tmp_char.end());
	}
}

void make_typed_int_vector(vector<char> &out, const vector<int> &_in)
{
	vector<char> tmp_char;
	int type;
	int8_t size_type;
	unsigned int max = 0;
	out.resize(0);

	for (unsigned int ui=0; ui<_in.size(); ui++) {
		if ((abs(_in[ui]) >(int)max) && ((int8_t)_in[ui] != (int8_t)0x80))
			max = abs(_in[ui]);
	}

	if (max < 127)
		type = 1;
	else if (max < 32767)
		type = 2;
	else
		type = 3;

	if (_in.size() < 15) {
		size_type = (int8_t)_in.size();
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back(size_type);

		make_typed_int(tmp_char, (int)_in.size(), true);
		out.insert(out.end(), tmp_char.begin(), tmp_char.begin());
	}

	for (unsigned int ui=0; ui<_in.size(); ui++) {
		make_int(tmp_char, _in[ui], type);
		out.insert(out.end(), tmp_char.begin(), tmp_char.end());
	}
}

void make_int(vector<char> &out, const int &_in, int type)
{
	out.resize(0);
	if (type == 1) {
		int8_t tmp_int;
		if (_in == (int)0x80000000 || _in >= 128)
			tmp_int = (int8_t)0x80;
		else if (_in == (int)0x80000001)
			tmp_int = (int8_t)0x81;
		else
			tmp_int = (int8_t)_in;
		out.push_back((int8_t)tmp_int);
	} else if (type == 2) {
		int16_t tmp_int;

		if (_in == (int)0x80000000 || _in >= 32768)
			tmp_int = (int16_t)0x8000;
		else if (_in == (int)0x80000001)
			tmp_int = (int8_t)0x8001;
		else
			tmp_int = (int16_t)_in;

		int8_t split;
		for (unsigned int ui=0; ui<2; ui++) {
			split = tmp_int & (int16_t)0x00FF;//0000000011111111
			out.push_back(split);
			tmp_int = tmp_int >> 8;
		}
	} else {
		int32_t tmp_int;
		tmp_int = (int32_t)_in;

		int8_t split;
		for (unsigned int ui=0; ui<4; ui++) {
			split = tmp_int & (int32_t)0x0000FF;
			out.push_back((int8_t)split);
			tmp_int = tmp_int >> 8;
		}
	}
}

void make_typed_float_vector(vector<char> &out, const string &_in, int number)
{
	vStr split_string;
	int8_t size_type;
	int max_val = 0;
	out.resize(0);

	if (_in == " " || _in == "." || _in == "") {
		size_type = (int8_t)0;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back(size_type);
		return;
	}

//FIXME	header::tokenize(_in, ',', split_string);
	if (number == -1)
		max_val = (int)split_string.size();
	else
		max_val = number;

	if (max_val < 15) {
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back(size_type);

		vector<char> size_vector;
		make_typed_int(size_vector, max_val, true);
		out.insert(out.end(), size_vector.begin(), size_vector.end());
	}

	float value;
	char missing[4] ={ char(0x01), char(0x00), char(0x80), char(0x7F) };
	char end[4] ={ char(0x02), char(0x00), char(0x80), char(0x7F) };

	for (unsigned int ui=0; (int)ui<max_val; ui++) {
		if (ui < split_string.size()) {
//FIXME			value = (float)header::str2double(split_string[ui], 0x7F800001);
		} else
			value = (float)0x7F800002;

		char *p = (char *)&value;

		for (unsigned int uj=0; uj<sizeof(value); uj++)
			if (value == (float)0x7F800001)
				out.push_back(missing[uj]);
			else if (value == (float)0x7F800002)
				out.push_back(end[uj]);
			else
				out.push_back(p[uj]);
	}
}

void make_typed_float_vector(vector<char> &out, const vStr &_in, int number)
{
	vStr split_string;
	int8_t size_type;
	unsigned int max_val = 0;
	out.resize(0);

	if (number == -1) {
		unsigned int tmp_int = 0;
		for (unsigned int ui=0; ui<_in.size(); ui++) {
			tmp_int = (wsUint)count(_in[ui].begin(), _in[ui].end(), ',');
			if (tmp_int > max_val)
				max_val = tmp_int;
		}
		max_val++;
	} else
		max_val = number;

	if (max_val < 15) {
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back(size_type);
	} else {
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back(size_type);

		vector<char> size_vector;
		make_typed_int(size_vector, max_val, true);
		out.insert(out.end(), size_vector.begin(), size_vector.end());
	}

	float value;
	char missing[4] ={ char(0x01), char(0x00), char(0x80), char(0x7F) };
	char end[4] ={ char(0x02), char(0x00), char(0x80), char(0x7F) };
	for (unsigned int ui=0; ui<_in.size(); ui++) {
//FIXME		header::tokenize(_in[ui], ',', split_string);
		for (unsigned int uj=0; uj<max_val; uj++) {
			if (uj < split_string.size()) {
//FIXME				value = (float)header::str2double(split_string[uj], 0x7F800001);
			} else
				value = (float)0x7F800002;

			char *p = (char *)&value;

			for (unsigned int uk=0; uk<sizeof(value); uk++) {
				if (value == float(0x7F800001))
					out.push_back(missing[uk]);
				else if (value == float(0x7F800002))
					out.push_back(end[uk]);
				else
					out.push_back(p[uk]);
			}
		}
	}
}

} // End namespace ONETOOL
