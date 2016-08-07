//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Laboratory for Medical Bioinformatics, University of Manchester, UK
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


#include "seamass.hpp"
#include "BasisFunctions.hpp"
#include "OptimiserASRL.hpp"


#include <iostream>


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
seaMass(std::vector<double>& bin_counts,
        std::vector<double>& bin_locations,
        std::vector<li>& spectrum_index,
        std::vector<double>& start_times,
        std::vector<double>& finish_times,
        std::vector<double>& exposures,
        std::vector<char>& resolutions, double shrinkage, double tolerance)
{
    // difference between carbon12 and carbon13
    double rc_mz = pow(2.0, (double) mz_res) * 60 / 1.0033548378;

    // for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
    // this is conservative, 4 times might be ok, but 2 times isn't enough
    //merge_bins(mzs, intensities, 0.125 / rc_mz);

    // Convert intensities into format used in algorithm for gs
    vector<fp> gs; vector<li> is; vector<ii> js;
    create_gs(gs, is, js, bin_counts);
	for (ii j = 0; j < (ii)bin_counts.size(); j++) vector<double>().swap(bin_counts[j]);
    
     // Create our tree of bases
    vector<Basis*> bases;
    ii order = 3; // B-spline order
    
    // Construct BasisResampleMZ root node
    cout << endl << "Spectrometry rc_mz=" << rc_mz << ":" << mz_res << endl;
	BasisResampleMZ* bResampleMZ = new BasisResampleMZ(bases, bin_locations, gs, is, js, rc_mz, order, true);
	double mz_min = bResampleMZ->get_min();
	double mz_max = bResampleMZ->get_max();

	/*ostringstream mzh5;
	mzh5 << "/" << config_id << "/" << rc0_mz << "/";
	if (debug)
	{
		h5out->write_cdata(mzh5.str(), mzs, "SpectrumMZ");
	}*/

	for (ii j = 0; j < (ii)bin_locations.size(); j++) vector<double>().swap(bin_locations[j]);
    //while (bases.back()->get_cm().n[0] > order + 1)
    //{
    //    new BasisDyadicScale(bases, bases.back(), 0, order, true);
    //}
    ii n_core_bases = bases.size();

    double duration = 0.0;

	double start = omp_get_wtime();
    double rc_st = pow(2.0, (double) st_res);
    cout << endl << "Scan time rc_st=" << rc_st << ":" << st_res << endl;
        
    ////////////////////////////////////////////////////////////////////////////////////
    // INIT OTHER BASIS FUNCTIONS
        
    for (ii i = n_core_bases; i < (ii) bases.size(); i++) delete bases[i];
    bases.resize(n_core_bases);
        
    // BasisResampleRT
    BasisResampleRT* bResampleRT = new BasisResampleRT(bases, bResampleMZ, scan_times, js, exposures, rc_st, order);
	double rt_min = bResampleRT->get_min();
	double rt_max = bResampleRT->get_max();
	ii cs_basis = bResampleRT->get_index();
    while (bases.back()->get_cm().n[0] > order + 1)
    {
        new BasisDyadicScale(bases, bases.back(), 0, order);
    }
        
    // The dyadic scale bases
    Basis* last = bResampleRT;
    while (last->get_cm().n[1] > order + 1)
    {
        last = new BasisDyadicScale(bases, last, 1, order);
        while (bases.back()->get_cm().n[0] > order + 1)
        {
            new BasisDyadicScale(bases, bases.back(), 0, order);
        }
    }
        
    ////////////////////////////////////////////////////////////////////////////////////
    // OPTIMISATION
	double thres = 0.0;// 000000001; // L0 threshold
        
    //cout << "a1" << endl;
    OptimiserASRL* optimiser = new OptimiserASRL(bases, gs, 2);
    //cout << "a2" << endl;

    double shr = pow(2.0, (double) shrinkage);
    double tol = pow(2.0, (double) tolerance);
    double grad = DBL_MAX;
     
    // l1
    li nc = 0;
    for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
        nc += bases[j]->get_cm().size();
    cout << " L1 nc=" << nc << " shrinkage=" << shrinkage << ":" << fixed << setprecision(2) << shr << " tolerance=" << tolerance << ":" << setprecision(6) << tol << endl;

    //cout << "a3" << endl;
    /*if (debug >= 2)
    {
        ostringstream oss;
        oss << "/" << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/init";
 		vector<fp> fs;
		optimiser->synthesis(fs, h5out, oss.str(), cs_basis);
	}*/

    for (ii i = 0; grad > tol; i++)
    {
        grad = optimiser->step(i, shr);
		optimiser->threshold(thres);

        li nnz = 0;
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
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
                
        /*if (debug >= 2)
        {
            ostringstream oss;
            oss << "/" << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/L1/" << setfill('0') << setw(8) << i << "/";
 			vector<fp> fs;
			optimiser->synthesis(fs, h5out, oss.str(), cs_basis);
		}*/

	}
            
    // l0
    cout << " L0 threshold=" << fixed << setprecision(2) << thres << " tolerance=" << tolerance << ":" << setprecision(6) << tol << endl;
            
    grad = DBL_MAX;
    for (ii i = 0; i < 20; i++)
    {
        grad = optimiser->step(i, 0.0);
                
        li nnz = 0;
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
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
                
        /*if (debug >= 2)
        {
            ostringstream oss;
            oss << "/"  << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/L0/" << setfill('0') << setw(8) << i;
			vector<fp> fs;
			optimiser->synthesis(fs, h5out, oss.str(), cs_basis);
		}*/
    }
	cout << "Duration: " << (omp_get_wtime() - start) << "seconds" << endl;

    //////////////////////////////////////////////////////////////////////////////////
    // OUTPUT
			
	/*if (debug)
	{
		// write fc, fcs, fs seamass output
		ostringstream oss;
		oss  << "/" << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/";
		vector<fp> fs;
		optimiser->synthesis(fs, h5out, oss.str(), cs_basis);

		// write gs original spectrum intensities
		h5out->write_cdata(mzh5.str(), gs, "SpectrumCount");
		h5out->write_cdata(mzh5.str(), is, "SpectrumCountIndex");
		h5out->write_cdata(mzh5.str(), js, "SpectrumIndex");
		h5out->write_cdata(mzh5.str(), rts, "SpectrumScanTimes");
		h5out->write_cdata(mzh5.str(), exposures, "SpectrumExposure");
	}*/

    // write smv viz r-tree
	/*start = omp_get_wtime();
	ostringstream oss;
	oss << config_id << "_" << rc0_mz << "_" << rcr << "_" << shr << "_" << tol;
	vizout->write_cs(oss.str(), bases, n_core_bases, optimiser->get_cs(), mz_min, mz_max, rt_min, rt_max, optimiser->compute_norm_max_counts(n_core_bases));
 	cout << "Duration: " << (omp_get_wtime() - start)/60.0 << "mins" << endl;*/

    //ostringstream oss2; oss2 << id << "_" << config_id << "_" << rc0_mz << "_" << rcr << "_" << shr << "_" << tol << ".error.csv";
    //optimiser->calc_error(oss2.str());           
    delete optimiser;
	//delete vizout;
    
    // output mzML
    /*for (ii j = 0; j < js.size(); j++)
    {
        intensities[js[j]].resize(is[j+1]-is[j]);
        for (ii i = 0; i < is[j+1]-is[j]; i++)
        {
            intensities[js[j]][i] = gs[is[j]+i];
        }
    }*/
    
    ////////////////////////////////////////////////////////////////////////////////////
    for (ii j = 0; j < bases.size(); j++) delete bases[j];
	
	//if (debug) delete h5out; 
	//omp_set_num_threads(_threads);
}

bool
seaMass::
iteration()
{
	return false;
}


void create_gs(vector<fp>& gs,
	vector<li>& is,
	vector<ii>& js,
	const vector< vector<double> >& intensities)
{
	ii nj = 0;
	for (ii j = 0; j < intensities.size(); j++)
		if (intensities[j].size() > 0) nj++;

	js.resize(nj);
	is.resize(nj + 1);
	is.front() = 0;
	for (ii i = 0, j = 0; j < intensities.size(); j++)
		if (intensities[j].size() > 0)
		{
			js[i] = j;
			is[i + 1] = is[i] + intensities[j].size();
			i++;
		}

	gs.resize(is.back());
	#pragma omp parallel for
	for (ii j = 0; j < nj; j++)
		for (ii i = 0; i < intensities[js[j]].size(); i++)
		{
			gs[is[j] + i] = intensities[js[j]][i];
		}

	cout << "Raw gs=[" << js.size() << "/" << intensities.size() << "]:" << gs.size() << " mem=" << fixed << setprecision(2) << (is.size()*sizeof(li) + js.size()*sizeof(ii) + gs.size()*sizeof(fp)) / (1024.0*1024.0) << "Mb";
	cout << endl;
}