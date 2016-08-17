//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
//
// Copyright (C) 2016  biospi Laboratory, University of Bristol, UK
//
// This file is part of seaMass.
//
// seaMass is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef _SEAMASSRESTORATION_BASISFUNCTIONS_HPP_
#define _SEAMASSRESTORATION_BASISFUNCTIONS_HPP_


#include "core.hpp"
#include "seaMass.hpp"


struct CoeffsMetadata
{
	CoeffsMetadata(ii d);
	~CoeffsMetadata();

	li size() const;
	void print(std::ostream& out) const;
	void operator=(const CoeffsMetadata& cm);

	ii d;        // dimension of the output coefficients
	std::vector<ii> l; // coefficient dyadic level for each dimension
	std::vector<ii> o; // coefficient offset
	std::vector<ii> n; // number of coefficients per set
};


////////////////////////////////////////////////////////////////////////////////


class Basis
{
protected:
    ii index;       // index of this basis in the serialised tree
    Basis* parent;  // parent node
    ii child_count; // how many children synthesise to this node
    bool transient; // if transient, coefficients not part of fitting
    
    CoeffsMetadata cm;
  
    double volume; // statistic of the last error() call
    double discrep; // statistic of the last error() call
    double erro; // statistic of the last error() call
    double maxerr; // statistic of the last error() call

public:
	Basis(std::vector<Basis*>& bases,
          ii dimensions,
          Basis* parent = NULL, bool
          transient = false);
    virtual ~Basis() {}
    
	virtual void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true) = 0;
	virtual void analysis(std::vector<fp>& es, const std::vector<fp>& fs) = 0;
	virtual void l2norm(std::vector<fp>& es, const std::vector<fp>& fs) = 0;
	virtual void error(std::vector<fp>& fs, const std::vector<fp>& gs);
    
    ii get_index() { return index; }
    Basis* get_parent() { return parent; }
    bool is_transient() { return transient; }
    
    const CoeffsMetadata& get_cm() { return cm; }
    
    // some statistics of the last error() call
    double get_volume() { return volume; }
    double get_discrep() { return discrep; }
    double get_error() { return erro; }
    double get_maxerror() { return maxerr; }
};


////////////////////////////////////////////////////////////////////////////////


class BasisResampleMZ : public Basis
{
protected:
	const std::vector<li>& is;
    
    // CSR sparse A basis matrices and their transposes
	std::vector<ii> nnz;
	std::vector<ii> m;
	std::vector< std::vector<fp> > a;
	std::vector< std::vector<ii> > ia, ja;
    //vector< vector<fp> > at;
    //vector< vector<ii> > iat, jat;

	double mz_min, mz_max;
    
public:
	BasisResampleMZ(std::vector<Basis*>& bases,
		const std::vector< std::vector<double> >& mzs,
		const std::vector<fp>& gs,
		const std::vector<li>& is,
		const std::vector<ii>& js,
        ii rc,
        ii order = 3,
         bool transient = false);
    
    ~BasisResampleMZ();
    
	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true);
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs);
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs);

	double get_min() const { return mz_min; }
	double get_max() const { return mz_max; }
};


////////////////////////////////////////////////////////////////////////////////


class BasisResampleRT : public Basis
{
protected:
	ii nnz; ii m; std::vector<fp> a; std::vector<ii> ia, ja; // CSR sparse A basis matrix
	std::vector<fp> at; std::vector<ii> iat, jat;

	double rt_min, rt_max;
    
public:
	BasisResampleRT(std::vector<Basis*>& bases,
                    Basis* parent,
					const std::vector<double>& rts,
					const std::vector<ii>& js,
					const std::vector<double>& exposures,
                    ii rc,
                    ii order = 3,
                    bool transient = false);
    
    ~BasisResampleRT();
    
	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true);
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs);
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs);

	double get_min() const { return rt_min; }
	double get_max() const { return rt_max; }
};


////////////////////////////////////////////////////////////////////////////////


class BasisDyadicScale : public Basis
{
protected:
	ii nnz; ii m, n; std::vector<fp> a; std::vector<ii> ia, ja; // CSR sparse A basis matrix
	std::vector<fp> at; std::vector<ii> iat, jat;
    ii dim;
    
public:
	BasisDyadicScale(std::vector<Basis*>& bases,
                     Basis* parent,
                     ii dim,
                     ii order = 3,
                     bool transient = false);
    
    ~BasisDyadicScale();
    
	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true);
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs);
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs);
};

#endif // _SEAMASSRESTORATION_BASISFUNCTIONS_HPP_

