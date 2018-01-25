#pragma once
#ifndef __WISARD_TYPED_H__
#define __WISARD_TYPED_H__

namespace ONETOOL {

	using namespace std;

void make_typed_string(vector<char> &out, const string &_in, bool typed); //
void make_typed_int(vector<char> &out, const int &_in, bool typed); //
void make_int(vector<char> &out, const int &_in, int type);
void make_typed_int_vector(vector<char> &out, const vStr &_in, int number = -1);
void make_typed_int_vector(vector<char> &out, const string &_in, int number = -1);
void make_typed_int_vector(vector<char> &out, const vector<int> &_in);
void make_typed_float_vector(vector<char> &out, const string &_in, int number = -1);
void make_typed_float_vector(vector<char> &out, const vStr &_in, int number = -1);
void make_typed_string_vector(vector<char> &out, const vStr &_in, int number = -1);
void make_typed_GT_vector(vector<char> &out, vStr &_in);
void make_type_size(vector<char> &out, const unsigned int &type, const unsigned int &size);

} // End namespace ONETOOL

#endif
