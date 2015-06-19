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

#include <fstream>
#include <sstream>

#include "seamass.hpp"
#include "BasisFunctions.hpp"
#include "OptimiserASRL.hpp"
#include "SMOWriter.hpp"
#include "SMVWriter.hpp"

namespace seamass
{

using namespace std;

void notice()
{
    cout << endl;
    cout << "seaMass - Copyright (C) 2015 - biospi Laboratory, EEE, University of Liverpool, UK" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
    cout << endl;
}

void process(const std::string& id,
             const std::string& config_id,
			 int instrument_type,
             vector<double>& rts,
             vector< vector<double> >& mzs,
             vector< vector<double> >& intensities,
             int rc0_mz, int rc1_mz,
             int rc0_rt, int rc1_rt,
             int shrinkage0, int shrinkage1,
             int tolerance0, int tolerance1,
		     int threads, int debug)
{
	int _threads = omp_get_num_threads();
	omp_set_num_threads(threads);

	cout << "Config threads=" << threads << " debug=" << debug << endl;
    cout << "Input id=" << id << " config_id=" << config_id << " instrument_type=" << instrument_type << endl;
    cout << endl;

	// create SMO if debug >= 1
	ostringstream oss; oss << id << ".smo";
	SMOWriter* h5out = 0;
	if (debug) h5out = new SMOWriter(oss.str()); 

	// create SMV directory
	ostringstream oss2; oss2 << id << ".smv";
	SMVWriter* vizout = new SMVWriter(oss2.str());

    ////////////////////////////////////////////////////////////////////////////////////
    // INIT RAW DATA AND MZ BASIS
    
    // difference between carbon12 and carbon13
    double rc_mz = pow(2.0, (double) rc0_mz) * 60 / 1.0033548378;

    // Ensure the raw data is in binned format and compute exposures
	vector<fp> exposures;
    bin_mzs_intensities(mzs, intensities, instrument_type, exposures);
    
    // for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
    // this is conservative, 4 times might be ok, but 2 times isn't enough
    merge_bins(mzs, intensities, 0.125 / rc_mz);

    // Convert intensities into format used in algorithm for gs
    vector<fp> gs; vector<li> is; vector<ii> js;
    create_gs(gs, is, js, intensities);
    for (ii j = 0; j < (ii) intensities.size(); j++) vector<double>().swap(intensities[j]);
    
     // Create our tree of bases
    vector<Basis*> bases;
    ii order = 3; // B-spline order
    
    // Construct BasisResampleMZ root node
    cout << endl << "Spectrometry rc_mz=" << rc0_mz << ":" << rc_mz << endl;
    BasisResampleMZ* bResampleMZ = new BasisResampleMZ(bases, mzs, gs, is, js, rc0_mz, order, true);
	double mz_min = bResampleMZ->get_min();
	double mz_max = bResampleMZ->get_max();

	ostringstream mzh5;
	mzh5 << "/" << config_id << "/" << rc0_mz << "/";
	h5out->write_cdata(mzh5.str(), mzs,"SpectrumMZ");

    for (ii j = 0; j < (ii) mzs.size(); j++) vector<double>().swap(mzs[j]);
    //while (bases.back()->get_cm().n[0] > order + 1)
    //{
    //    new BasisDyadicScale(bases, bases.back(), 0, order, true);
    //}
    ii n_core_bases = bases.size();

    double duration = 0.0;
    for (ii rcr = rc0_rt; rcr <= rc1_rt; rcr++)
    {
		double start = omp_get_wtime();
        double rc_rt = pow(2.0, (double) rcr);
        cout << endl << "Chromatography rc_rt=" << rcr << ":" << rc_rt << endl;
        
        ////////////////////////////////////////////////////////////////////////////////////
        // INIT OTHER BASIS FUNCTIONS
        
        for (ii i = n_core_bases; i < (ii) bases.size(); i++) delete bases[i];
        bases.resize(n_core_bases);
        
        // BasisResampleRT
        vector<ii> scale_bases(1, n_core_bases);
        BasisResampleRT* bResampleRT = new BasisResampleRT(bases, bResampleMZ, rts, js, rcr, order);
		double rt_min = bResampleRT->get_min();
		double rt_max = bResampleRT->get_max();
        while (bases.back()->get_cm().n[0] > order + 1)
        {
            new BasisDyadicScale(bases, bases.back(), 0, order);
        }
        
        // The dyadic scale bases
        Basis* last = bResampleRT;
        while (last->get_cm().n[1] > order + 1)
        {
            scale_bases.push_back(bases.size());
            last = new BasisDyadicScale(bases, last, 1, order);
            while (bases.back()->get_cm().n[0] > order + 1)
            {
                new BasisDyadicScale(bases, bases.back(), 0, order);
            }
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // OPTIMISATION
        double thres = 1.0; // L0 threshold
        
        for (ii shr = shrinkage0; shr <= shrinkage1; shr++)
        for (ii tol = tolerance1; tol >= tolerance0; tol--)
        {
            //cout << "a1" << endl;
            OptimiserASRL* optimiser = new OptimiserASRL(bases, gs, 2);
            //cout << "a2" << endl;

            double shrinkage = pow(2.0, (double) shr);
            double tolerance = pow(2.0, (double) tol);
            double grad = DBL_MAX;
     
            // l1
            li nc = 0;
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
                nc += bases[j]->get_cm().size();
            cout << " L1 nc=" << nc << " shrinkage=" << shr << ":" << fixed << setprecision(2) << shrinkage << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;

            //cout << "a3" << endl;

            for (ii i = 0; grad > tolerance; i++)
            {
                grad = optimiser->step(i, shrinkage);
                
                li nnz = 0;
                for (ii j = 0; j < (ii) bases.size(); j++)
                if (!bases[j]->is_transient())
                for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
                if (optimiser->get_cs()[j][i] > thres)
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
                
                if (debug >= 2)
                {
                    ostringstream oss;
                    oss << "/" << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/L1/" << setfill('0') << setw(8) << i;
                    optimiser->write_h5(*h5out, oss.str(), scale_bases, is, js, exposures);
                }
            }
            
            // l0
            cout << " L0 threshold=" << fixed << setprecision(2) << thres << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;
            
            optimiser->threshold(thres);
            grad = DBL_MAX;
            for (ii i = 0; grad > tolerance; i++)
            {
                grad = optimiser->step(i, 0.0);
                
                li nnz = 0;
                for (ii j = 0; j < (ii) bases.size(); j++)
                if (!bases[j]->is_transient())
                for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
                if (optimiser->get_cs()[j][i] > thres)
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
                
                if (debug >= 2)
                {
                    ostringstream oss;
                    oss << "/"  << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/L0/" << setfill('0') << setw(8) << i;
                    optimiser->write_h5(*h5out, oss.str(), scale_bases, is, js, exposures);
                }
            }
			cout << "Duration: " << (omp_get_wtime() - start)/60.0 << "mins" << endl;

            //////////////////////////////////////////////////////////////////////////////////
            // OUTPUT
			
			// write smo
			if (debug)
			{
				ostringstream oss;
				oss  << "/" << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/";
				optimiser->write_h5(*h5out, oss.str(), scale_bases, is, js, exposures);
			}

            // write smv viz r-tree
			start = omp_get_wtime();
			ostringstream oss;
			oss << config_id << "_" << rc0_mz << "_" << rcr << "_" << shr << "_" << tol;
			vizout->write_cs(oss.str(), bases, n_core_bases, optimiser->get_cs(), mz_min, mz_max, rt_min, rt_max, optimiser->compute_norm_max_counts(n_core_bases));
 			cout << "Duration: " << (omp_get_wtime() - start)/60.0 << "mins" << endl;

            //ostringstream oss2; oss2 << id << "_" << config_id << "_" << rc0_mz << "_" << rcr << "_" << shr << "_" << tol << ".error.csv";
            //optimiser->calc_error(oss2.str());           
            delete optimiser;
        }
        
    }
	delete vizout;
    
    // output mzML
    for (ii j = 0; j < js.size(); j++)
    {
        intensities[js[j]].resize(is[j+1]-is[j]);
        for (ii i = 0; i < is[j+1]-is[j]; i++)
        {
            intensities[js[j]][i] = gs[is[j]+i];
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    for (ii j = 0; j < bases.size(); j++) delete bases[j];
	
	if (debug) delete h5out; 
	omp_set_num_threads(_threads);
}

}


