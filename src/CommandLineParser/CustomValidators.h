/*
 * CustomValidators.h
 *
 *  Created on: Jan 11, 2014
 *      Author: sebastian
 */
/**
 * Here we overload the validate function in order to accept only valid
 * command line input parameters.
 *
 * For a detailed description on overloading the default validator see:
 * http://www.boost.org/doc/libs/1_55_0/doc/html/program_options/howto.html#idp163429032
 */
#ifndef CUSTOMVALIDATORS_H_
#define CUSTOMVALIDATORS_H_

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <string>
using namespace boost;
using namespace boost::program_options;
using namespace std;

struct input_PDB_filename {
	input_PDB_filename(string const & val) :
			filename(val) {	}
	string filename;
	friend ostream& operator <<(ostream& s, const input_PDB_filename& idpbf) {
		s << idpbf.filename;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * PDB input filenames.
 *
 * The function takes four parameters. The first is the storage for
 * the value, and in this case is either empty or contains an instance
 * of the magic_number class. The second is the list of strings found
 * in the next occurrence of the option. The remaining two parameters
 * are needed to workaround the lack of partial template specialization
 * and partial function template ordering on some compilers.
 *
 * The function first checks that we don't try to assign to the same
 * option twice. Then it checks that only a single string was passed
 * in. Next the string is verified. If that test is passed, the parsed
 * value is stored into the v variable.
 */

void validate(any& v, vector<string> const& values,
		input_PDB_filename* /* target_type */, int) {
	static regex r("[^[.NUL.]]+(?:\\.pdb|\\.PDB|\\.xyzr|\\XYZR|\\.xyzrn|\\XYZRN)");
	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);
	// Do regex match
	if (regex_match(s, r))
		v = any(input_PDB_filename(s));
	 else
		throw validation_error(validation_error::invalid_option_value);
};
struct filename {
	filename(string const & outname) :
		fname(outname) {	}
	string fname;
	friend ostream& operator <<(ostream& s, const filename& out) {
		s << out.fname;
		return s;
	}
};
void validate(any& v, vector<string> const& values,
		filename* /* target_type */, int) {
	static regex r("[^[.NUL.]]+");
	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);
	// Do regex match
	if (regex_match(s, r))
		v = any(filename(s));
	 else
		throw validation_error(validation_error::invalid_option_value);
};

struct probe_radius {
public:
	probe_radius(float p) :
			p(p) { }
	float p;
	friend ostream& operator <<(ostream& s, const probe_radius& pr) {
		s << pr.p;
		return s;
	}
};
/**
 * Here we overload the validate function in order to accept only valid
 * probe radius parameters.
 */

void validate(any& v, vector<string> const& values,
		probe_radius* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0)
			v = any(probe_radius(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct resolution_param {
public:
	resolution_param(float r) :
			r(r) { }
	float r;
	friend ostream& operator <<(ostream& s, const resolution_param& rp) {
		s << rp.r;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		resolution_param* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0)
			v = any(resolution_param(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct out_type {
public:
	out_type(int st) :
			st(st) { }
	int st;
	friend ostream& operator <<(ostream& s, const out_type& type) {
		s << type.st;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		out_type* /* target_type */, int) {
	static regex r("[123456]");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		v = any(out_type(lexical_cast<int>(s)));
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct surf_type {
public:
	surf_type(int st) :
			st(st) { }
	int st;
	friend ostream& operator <<(ostream& s, const surf_type& type) {
		s << type.st;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		surf_type* /* target_type */, int) {
	static regex r("[123]");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		v = any(surf_type(lexical_cast<int>(s)));
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct num_iterations {
public:
	num_iterations(int iter) :
			iter(iter) { }
	int iter;
	friend ostream& operator <<(ostream& s, const num_iterations& type) {
		s << type.iter;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		num_iterations* /* target_type */, int) {
	static regex r("[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);
	if (regex_match(s, r)) {
		v = any(num_iterations(lexical_cast<int>(s)));
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};


struct shrink_param {
public:
	shrink_param(float lambda) :
		lambda(lambda) { }
	float lambda;
	friend ostream& operator <<(ostream& s, const shrink_param& rp) {
		s << rp.lambda;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		shrink_param* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0)
			v = any(shrink_param(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct inflate_param {
public:
	inflate_param(float mu) :
		mu(mu) { }
	float mu;
	friend ostream& operator <<(ostream& s, const inflate_param& rp) {
		s << rp.mu;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		inflate_param* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp < 0)
			v = any(inflate_param(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

#endif /* CUSTOMVALIDATORS_H_ */
