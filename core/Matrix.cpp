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


#include "Matrix.hpp"

#include "MatrixSparse.hpp"

#include <iomanip>


using namespace std;


Matrix::Matrix()
	: m_(0), n_(0), vs_(0), isOwned_(false)
{
}


/*Matrix::Matrix(const Matrix& a)
{
	init(a.m, a.n);
	memcpy(vs_, a.vs_, sizeof(fp) * a.size());
}*/


Matrix::~Matrix()
{
	free();
}


void Matrix::init(li m, ii n)
{
	free();

	m_ = m;
	n_ = n;
	s_ = n;
	vs_ = new fp[m_ * n_];
	isOwned_ = true;
}


void Matrix::init(li m, ii n, fp v)
{
	init(m, n);

	for (ii i = 0; i < m_ * n_; i++)
	{
		vs_[i] = v;
	}
}


void Matrix::init(li m, ii n, ii s, fp* vs)
{
	free();

	m_ = m;
	n_ = n;
	s_ = s;
	vs_ = vs;
}


void Matrix::init(const Matrix& a, li i, ii j, li m, ii n, ii s)
{
	free();
	
	m_ = m;
	n_ = n;
	s_ = (s == 0) ? a.s_ : s;
	vs_ = &a.vs_[i * s_ + j];
	isOwned_ = false;
}


void Matrix::free()
{
	m_ = 0;
	n_ = 0;

	if (isOwned_)
	{
		delete[] vs_;
		isOwned_ = false;
	}
}


li Matrix::nnz() const
{
	li nnz = 0;
	for (li i = 0; i < (li)m_ * n_; i++)
	{
		if (vs_[i] != 0.0) nnz++;
	}
	return nnz;
}


li Matrix:: size() const
{
	return (li)m_ * n_;
}


li Matrix::mem() const
{
	return sizeof(*this) + (li)m_ * n_ * sizeof(fp);
}


li Matrix::m() const
{
	return m_;
}


ii Matrix::n() const
{
	return n_;
}


bool Matrix::operator!() const
{ 
	return vs_ == 0; 
}


void Matrix::mul(const MatrixSparse& a, const Matrix& x, bool accumulate, bool transposeA)
{
	if (!*this)
	{
		// initialise output matrix if not already done
		init(transposeA ? a.n_ : a.m_, x.n_);
		accumulate = false;
	}

	static fp alpha = 1.0;
	fp beta = (fp)(accumulate ? 1.0 : 0.0);
	mkl_scsrmm(transposeA ? "T" : "N", &a.m_, &x.n_, &a.n_, &alpha, "G**C", a.vs_, a.js_, a.is_, &a.is_[1], x.vs_, &x.n_, &beta, vs_, &x.n_);

	cout << "  Y" << *this << (accumulate ? " += A" : " = A") << (transposeA ? "t" : "") << a << " . X" << x << endl;
}


void Matrix::elementwiseMul(const Matrix& a, const Matrix& b)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

	// split into tasty chunks for openmp and when using 32bit MKL_INT
	static ii chunk_size = 0x10000;
	ii chunks = (ii)size() / chunk_size + 1;
	ii last_chunk = size() % chunk_size;
	if (last_chunk == 0)
	{
		chunks--;
		last_chunk = chunk_size;
	}

	#pragma omp parallel for
	for (ii k = 0; k < chunks; k++)
	{
		vsMul((k == chunks - 1) ? last_chunk : chunk_size, &a.vs_[k * chunk_size], &b.vs_[k * chunk_size], &vs_[k * chunk_size]);
	}

	cout << "  Y" << *this << " = A" << a << " ./ B" << b << endl;
}


void Matrix::elementwiseMul(fp scale, const Matrix& a)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

	#pragma omp parallel for
	for (ii i = 0; i < a.size(); i++)
	{
		vs_[i] = scale * a.vs_[i];
	}

	cout << "  Y" << *this << " = " << scale << " * " << a << endl;
}


void Matrix::elementwiseDiv(const Matrix& n, const Matrix& d)
{
	// initialise output matrix if not already done
	if (!*this) init(d.m_, d.n_);

	// split into tasty chunks for openmp and when using 32bit MKL_INT
	static ii chunk_size = 0x10000;
	ii chunks = (ii) size() / chunk_size + 1;
	ii last_chunk = size() % chunk_size;
	if (last_chunk == 0)
	{
		chunks--;
		last_chunk = chunk_size;
	}

	#pragma omp parallel for
	for (ii k = 0; k < chunks; k++)
	{
		vsDiv((k == chunks - 1) ? last_chunk : chunk_size, &n.vs_[k * chunk_size], &d.vs_[k * chunk_size], &vs_[k * chunk_size]);
	}

	cout << "  Y" << *this << " = N" << n << " ./ D" << d << endl;
}


void Matrix::elementwiseSqrt(const Matrix& a)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

	// split into tasty chunks for openmp and when 32bit using MKL_INT
	static ii chunk_size = 0x10000;
	ii chunks = (ii) size() / chunk_size + 1;
	ii last_chunk = size() % chunk_size;
	if (last_chunk == 0)
	{
		chunks--;
		last_chunk = chunk_size;
	}

	#pragma omp parallel for
	for (ii k = 0; k < chunks; k++)
	{
		vsSqrt((k == chunks - 1) ? last_chunk : chunk_size, &a.vs_[k * chunk_size], &vs_[k * chunk_size]);
	}

	cout << "  Y" << *this << " = sqrt(A" << a << ")" << endl;
}


ostream& operator<<(ostream& os, const Matrix& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
		os << "[" << a.m() << "," << a.n() << "]:" << a.nnz() << "/" << a.size() << ":" << defaultfloat << setprecision(3) << 100.0 * a.nnz() / (double) a.size() << "%";
	}

	return  os;
}

