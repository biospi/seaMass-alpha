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


#include "seaMass.hpp"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <limits>

using namespace std;


void
seaMass::
notice()
{
    cout << endl;
    cout << "seaMass - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
    cout << endl;
}


seaMass::
seaMass(Input& input, const std::vector<ii>& scales)
{
	init(input, scales);
}


seaMass::
seaMass(Input& input, const Output& seed)
{
	vector<ii> scales(seed.scales.size());
	init(input, scales);

	// todo
	/*cout << seed.weights.size() << endl;
	for (ii i = 0; i < seed.weights.size(); i++)
	{
		//cout << seed.scales[0][i] << "," << seed.scales[1][i] << ":" << seed.offsets[0][i] << "," << seed.offsets[1][i] << ":" << seed.weights[i] << endl;
	}*/
}


void
seaMass::
init(Input& input, const std::vector<ii>& scales)
{
	// for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
	// this is conservative, 4 times might be ok, but 2 times isn't enough
	//merge_bins(mzs, intensities, 0.125 / rc_mz);

	////////////////////////////////////////////////////////////////////////////////////
	// INIT BASIS FUNCTIONS

	// Create our tree of bases
	ii order = 3; // B-spline order

	if (input.spectrum_index.size() <= 2)
	{
		dimensions = 1;

		BasisResampleMZ* bResampleMZ = new BasisResampleMZ(bases, input.bin_counts, input.spectrum_index, input.bin_edges, scales[0], order);
		while (bases.back()->get_cm().n[0] > order + 1)
		{
			new BasisDyadicScale(bases, bases.back(), 0, order);
		}
	}
	else
	{
		dimensions = 2;

		BasisResampleMZ* bResampleMZ = new BasisResampleMZ(bases, input.bin_counts, input.spectrum_index, input.bin_edges, scales[0], order, true);
		BasisResampleRT* bResampleRT = new BasisResampleRT(bases, bResampleMZ, input.start_times, input.finish_times, input.exposures, scales[1], order);
		Basis* last = bResampleRT;
		while (last->get_cm().n[1] > order + 1)
		{
			last = new BasisDyadicScale(bases, last, 1, order);
			while (bases.back()->get_cm().n[0] > order + 1)
			{
				new BasisDyadicScale(bases, bases.back(), 0, order);
			}
		}
	}

	optimiser = new OptimiserASRL(bases, input.bin_counts, 2);
}


bool
seaMass::
step(double shrinkage, double tolerance)
{
	////////////////////////////////////////////////////////////////////////////////////
	// OPTIMISATION
	double threshold = 0.0;// 000000001; // L0 threshold
	double grad = numeric_limits<double>::max();

	// l1
	li nc = 0;
	for (ii j = 0; j < (ii)bases.size(); j++)
		if (!bases[j]->is_transient())
			nc += bases[j]->get_cm().size();
	cout << " L1 nc=" << nc << " shrinkage=" << shrinkage << ":" << fixed << setprecision(2) << shrinkage << " tolerance=" << tolerance << ":" << setprecision(6) << tolerance << endl;

	for (ii i = 0; grad > tolerance; i++)
	{
		grad = optimiser->step(i, shrinkage);
		optimiser->threshold(threshold);

		li nnz = 0;
		for (ii j = 0; j < (ii)bases.size(); j++)
			if (!bases[j]->is_transient())
				for (ii i = 0; i < (ii)bases[j]->get_cm().size(); i++)
					if (optimiser->get_cs()[j][i] > 0.0)
					{
						nnz++;
					}

		cout << "  f: " << fixed << setw(5) << i;
		cout << "  err: " << setw(8) << setprecision(5) << bases.front()->get_error();
		cout << "  max: " << setw(8) << setprecision(1) << bases.front()->get_maxerror();
		cout << "  dis: " << setw(8) << setprecision(5) << bases.front()->get_discrep();
		cout << "  vol: " << setw(8) << setprecision(5) << bases.front()->get_volume();
		cout << "  nnz: " << setw(8) << setprecision(5) << nnz;
		cout << "  grad: " << setiosflags(ios::fixed) << setprecision(6) << setw(8) << grad << endl;
	}

	// l0
	cout << " L0 threshold=" << fixed << setprecision(2) << threshold << " tolerance=" << tolerance << ":" << setprecision(6) << tolerance << endl;

	grad = numeric_limits<double>::max();
	for (ii i = 0; i < 20; i++)
	{
		grad = optimiser->step(i, 0.0);

		li nnz = 0;
		for (ii j = 0; j < (ii)bases.size(); j++)
			if (!bases[j]->is_transient())
				for (ii i = 0; i < (ii)bases[j]->get_cm().size(); i++)
					if (optimiser->get_cs()[j][i] > 0.0)
					{
						nnz++;
					}

		cout << "  f: " << fixed << setw(5) << i;
		cout << "  err: " << setw(8) << setprecision(5) << bases.front()->get_error();
		cout << "  max: " << setw(8) << setprecision(1) << bases.front()->get_maxerror();
		cout << "  dis: " << setw(8) << setprecision(5) << bases.front()->get_discrep();
		cout << "  vol: " << setw(8) << setprecision(5) << bases.front()->get_volume();
		cout << "  nnz: " << setw(8) << setprecision(5) << nnz;
		cout << "  grad: " << setiosflags(ios::fixed) << setprecision(6) << setw(8) << grad << endl;
	}

	return false;
}


void
seaMass::
get_output(Output& output) const
{
	output.scales.resize(dimensions);
	output.offsets.resize(dimensions);
	output.baseline_size.resize(dimensions);
	output.baseline_scale.resize(dimensions);
	output.baseline_offset.resize(dimensions);

	for (ii i = 0; i < dimensions; i++)
	{
		output.baseline_size[i] = bases[dimensions - 1]->get_cm().n[i];
		output.baseline_scale[i] = bases[dimensions - 1]->get_cm().l[i];
		output.baseline_offset[i] = bases[dimensions - 1]->get_cm().o[i];
	}

	for (ii j = 0; j < (ii)bases.size(); j++)
	{
		if (!bases[j]->is_transient())
		{
			for (ii i = 0; i < (ii)bases[j]->get_cm().size(); i++)
			{
				double c = optimiser->get_cs()[j][i];
				if (c > 0.0)
				{
					output.weights.push_back(c);

					ii index = i;
					for (ii d = 0; d < dimensions; d++)
					{
						output.scales[d].push_back(bases[j]->get_cm().l[d]);
						output.offsets[d].push_back(bases[j]->get_cm().o[d] + index % bases[j]->get_cm().n[d]);
						index /= bases[j]->get_cm().n[d];
					}
				}
			}
		}
	}
}


void
seaMass::
get_output_bin_counts(std::vector<fp>& bin_counts) const
{
	optimiser->synthesis(bin_counts);
}


void
seaMass::
get_output_control_points(ControlPoints& control_points) const
{
	optimiser->synthesis(control_points.coeffs, dimensions - 1);

	for (ii i = 0; i < dimensions; i++)
	{
		control_points.size[i] = bases[dimensions - 1]->get_cm().n[i];
		control_points.scale[i] = bases[dimensions - 1]->get_cm().l[i];
		control_points.offset[i] = bases[dimensions - 1]->get_cm().o[i];
	}
}
