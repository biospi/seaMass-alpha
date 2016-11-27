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
#include <iostream>
#include <cmath>

using namespace std;


Matrix::Matrix()
	: m_(0), n_(0), vs_(0), isOwned_(false)
{
}


Matrix::~Matrix()
{
	free();
}


void Matrix::init(li m, ii n)
{
	free();

	m_ = m;
	n_ = n;
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


void Matrix::init(li m, ii n,fp* vs)
{
	free();

	m_ = m;
	n_ = n;
	vs_ = vs;
}


void Matrix::init(const Matrix& a, li i, ii j, li m, ii n)
{
	free();

	m_ = m;
	n_ = n;
	vs_ = &a.vs_[i * n_ + j];
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


void Matrix::mul(const MatrixSparse& a, const Matrix& x, bool accumulate, bool transposeA, bool transposeXY)
{
#ifndef NDEBUG
	cout << "  " << (transposeXY ? "(A" : "A") << (transposeA ? "t" : "") << a << " . X" << (transposeXY ? "t" : "") << x << (transposeXY ? ")t" : "") << (accumulate ? " =+ " : " = ") << flush;
#endif

	if (!*this)
	{
		// initialise output matrix if not already done
		li m = transposeA ? a.n_ : a.m_;
		li n = transposeXY ? x.m_ : x.n_;
		init(transposeXY ? n : m, transposeXY ? m : n);
		accumulate = false;
	}

	static fp alpha = 1.0;
	fp beta = (fp)(accumulate ? 1.0 : 0.0);
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	mkl_sparse_s_mm(
		transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE,
		alpha,
		a.mat_,
		descr,
		transposeXY ? SPARSE_LAYOUT_COLUMN_MAJOR : SPARSE_LAYOUT_ROW_MAJOR,
		x.vs_,
		transposeXY ? x.m_ : x.n_,
		x.n_,
		beta,
		vs_,
		transposeXY ? (transposeA ? a.n_ : a.m_) : x.n_
	);

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::copy(const Matrix& a)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  A" << a << " = " << flush;
#endif

	#pragma omp parallel for
	for (ii i = 0; i < a.size(); i++)
	{
		vs_[i] = a.vs_[i];
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwiseAdd(const Matrix& a, fp beta)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  A" << a << " + " << beta << " = " << flush;
#endif

	#pragma omp parallel for
	for (ii i = 0; i < a.size(); i++)
	{
		vs_[i] = a.vs_[i] + beta;
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwiseMul(const Matrix& a, fp beta)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  A" << a << " * " << beta << " = " << flush;
#endif

	#pragma omp parallel for
	for (ii i = 0; i < a.size(); i++)
	{
		vs_[i] = a.vs_[i] * beta;
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwiseMul(const Matrix& a, const Matrix& b)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  A" << a << " ./ B" << b << " = " << flush;
#endif

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

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwiseDiv(const Matrix& n, const Matrix& d)
{
	// initialise output matrix if not already done
	if (!*this) init(d.m_, d.n_);

#ifndef NDEBUG
	cout << "  N" << n << " ./ D" << d << " = " << flush;
#endif

	#pragma omp parallel for
	for (li i = 0; i < d.size(); i++)
	{
		if (d.vs_[i] != 0.0)
		{
			vs_[i] = n.vs_[i] / d.vs_[i];
		}
		else
		{
			vs_[i] = 0.0;
		}
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwiseSqrt(const Matrix& a)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  sqrt(A" << a << ") = " << flush;
#endif

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

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwisePow(const Matrix& a, fp power)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  (A" << a << ")^" << power << " = " << flush;
#endif

	// split into tasty chunks for openmp and when 32bit using MKL_INT
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
		vsPowx((k == chunks - 1) ? last_chunk : chunk_size, &a.vs_[k * chunk_size], power, &vs_[k * chunk_size]);
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::elementwiseLn(const Matrix& a)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  ln(A" << a << ") = " << flush;
#endif

	#pragma omp parallel for
	for (li i = 0; i < a.size(); i++)
	{
		if (a.vs_[i] != 0.0)
		{
			vs_[i] = log(a.vs_[i]);
		}
		else
		{
			vs_[i] = 0.0;
		}
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void Matrix::prune(const Matrix& a, fp threshold)
{
	// initialise output matrix if not already done
	if (!*this) init(a.m_, a.n_);

#ifndef NDEBUG
	cout << "  (A" << a << " >= ";
	cout.unsetf(ios::floatfield); 
	cout << setprecision(8) << threshold << ") ? A : 0.0 = " << flush;
#endif

	#pragma omp parallel for
	for (li i = 0; i < size(); i++)
	{
		vs_[i] = (a.vs_[i] >= threshold) ? a.vs_[i] : (fp)0.0;
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


double Matrix::sum() const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum(Y" << *this << ") = " << flush;
#endif

	#pragma omp parallel for reduction(+:sum)
	for (li i = 0; i < size(); i++)
	{
		sum += vs_[i];
	}

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


double Matrix::sumSqrs() const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum((Y" << *this << ")^2) = " << flush;
#endif

	#pragma omp parallel for reduction(+:sum)
	for (li i = 0; i < size(); i++)
	{
		sum += vs_[i] * vs_[i];
	}

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


double Matrix::sumSqrDiffs(const Matrix& a) const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum((Y" << *this << " - A" << a << ")^2) = " << flush;
#endif

	#pragma omp parallel for reduction(+:sum)
	for (li i = 0; i < size(); i++)
	{
		sum += (vs_[i] - a.vs_[i]) * (vs_[i] - a.vs_[i]);
	}

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


const fp* Matrix::getVs() const
{
	return vs_;
}


ostream& operator<<(ostream& os, const Matrix& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
		os << "[" << a.m() << "," << a.n() << "]:" << a.nnz() << "/" << a.size() << ":";
		os.unsetf(ios::floatfield);
		os << setprecision(3) << 100.0 * a.nnz() / (double) a.size() << "%";
	}

	return  os;
}
