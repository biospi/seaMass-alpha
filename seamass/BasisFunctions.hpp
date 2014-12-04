//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Bioinformatics Laboratory, University of Manchester, UK
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
	Basis(vector<Basis*>& bases,
          ii dimensions,
          Basis* parent = NULL, bool
          transient = false);
    virtual ~Basis() {}
    
    virtual void synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum = true) = 0;
    virtual void analysis(vector<fp>& es, const vector<fp>& fs) = 0;
    virtual void l2norm(vector<fp>& es, const vector<fp>& fs) = 0;
    virtual void error(vector<fp>& fs, const vector<fp>& gs);
    
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
    const vector<li>& is;
    
    // CSR sparse A basis matrices and their transposes
    vector<ii> nnz;
    vector<ii> m;
    vector< vector<fp> > a;
    vector< vector<ii> > ia, ja;
    //vector< vector<fp> > at;
    //vector< vector<ii> > iat, jat;

	double mz_min, mz_max;
    
public:
	BasisResampleMZ(vector<Basis*>& bases,
                    const vector< vector<double> >& mzs,
                    const vector<fp>& gs,
                    const vector<li>& is,
                    const vector<ii>& js,
                    ii rc,
                    ii order = 3,
                    bool transient = false);
    
    ~BasisResampleMZ();
    
    void synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum = true);
    void analysis(vector<fp>& es, const vector<fp>& fs);
    void l2norm(vector<fp>& es, const vector<fp>& fs);

	double get_min() const { return mz_min; }
	double get_max() const { return mz_max; }
};


////////////////////////////////////////////////////////////////////////////////


class BasisResampleRT : public Basis
{
protected:
    ii nnz; ii m; vector<fp> a; vector<ii> ia, ja; // CSR sparse A basis matrix
    vector<fp> at; vector<ii> iat, jat;

	double rt_min, rt_max;
    
public:
	BasisResampleRT(vector<Basis*>& bases,
                    Basis* parent,
                    const vector<double>& rts,
                    const vector<ii>& js,
                    ii rc,
                    ii order = 3,
                    bool transient = false);
    
    ~BasisResampleRT();
    
    void synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum = true);
    void analysis(vector<fp>& es, const vector<fp>& fs);
    void l2norm(vector<fp>& es, const vector<fp>& fs);

	double get_min() const { return rt_min; }
	double get_max() const { return rt_max; }
};


////////////////////////////////////////////////////////////////////////////////


class BasisDyadicScale : public Basis
{
protected:
    ii nnz; ii m, n; vector<fp> a; vector<ii> ia, ja; // CSR sparse A basis matrix
    vector<fp> at; vector<ii> iat, jat;
    ii dim;
    
public:
	BasisDyadicScale(vector<Basis*>& bases,
                     Basis* parent,
                     ii dim,
                     ii order = 3,
                     bool transient = false);
    
    ~BasisDyadicScale();
    
    void synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum = true);
    void analysis(vector<fp>& es, const vector<fp>& fs);
    void l2norm(vector<fp>& es, const vector<fp>& fs);
};

#endif // _SEAMASSRESTORATION_BASISFUNCTIONS_HPP_

