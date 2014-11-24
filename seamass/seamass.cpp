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

namespace seamass
{


using namespace std;


void notice()
{
    cout << endl;
    cout << "seaMass - Copyright (C) 2013 - CADET Bioinformatics Laboratory, University of Manchester, UK" << endl;
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

	ostringstream oss; oss << id << ".smo";
	HDF5File* h5out = 0;
	if (debug) h5out = new HDF5File(oss.str()); 

    ////////////////////////////////////////////////////////////////////////////////////
    // INIT RAW DATA AND MZ BASIS
    
    // difference between carbon12 and carbon13
    double rc_mz = pow(2.0, (double) rc0_mz) * 60 / 1.0033548378;

    // Ensure the raw data is in binned format
    bin_mzs_intensities(mzs, intensities, instrument_type);
    
    // for speed only, merge bins if rc_mz is set higher than twice bin width
    merge_bins(mzs, intensities, 0.5 / rc_mz);

    // Convert intensities into format used in algorithm for gs
    vector<fp> gs; vector<li> is; vector<ii> js;
    create_gs(gs, is, js, intensities);
    for (ii j = 0; j < (ii) intensities.size(); j++) vector<double>().swap(intensities[j]);
    
     // Create our tree of bases
    vector<Basis*> bases;
    ii order = 3; // B-spline order
    
    // Construct BasisResampleMZ root node
    cout << endl << "Spectrometry rc_mz=" << rc0_mz << ":" << rc_mz << endl;
    BasisResampleMZ* bResampleMZ = new BasisResampleMZ(bases, mzs, gs, is, js, rc_mz, order);
    for (ii j = 0; j < (ii) mzs.size(); j++) vector<double>().swap(mzs[j]);
    while (bases.back()->get_cm().n[0] > order + 1)
    {
        new BasisDyadicScale(bases, bases.back(), 0, order);
    }
    ii n_core_bases = bases.size();

    double duration = 0.0;
    for (ii rcr = rc0_rt; rcr <= rc1_rt; rcr++)
    {
        double rc_rt = pow(2.0, (double) rcr);
        cout << endl << "Chromatography rc_rt=" << rcr << ":" << rc_rt << endl;
        
        ////////////////////////////////////////////////////////////////////////////////////
        // INIT OTHER BASIS FUNCTIONS
        
        for (ii i = n_core_bases; i < (ii) bases.size(); i++) delete bases[i];
        bases.resize(n_core_bases);
        
        // BasisResampleRT
        vector<ii> scale_bases(1, n_core_bases);
        BasisResampleRT* bResampleRT = new BasisResampleRT(bases, bResampleMZ, rts, js, rc_rt, order);
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
            OptimiserASRL* optimiser = new OptimiserASRL(bases, gs, 2);
            
            double shrinkage = pow(2.0, (double) shr);
            double tolerance = pow(2.0, (double) tol);
            double grad = DBL_MAX;
     
            // l1
            ii nc = 0;
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
                nc += bases[j]->get_cm().size();
            cout << " L1 nc=" << nc << " shrinkage=" << shr << ":" << fixed << setprecision(2) << shrinkage << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;

            for (ii i = 0; grad > tolerance; i++)
            {
                double start = omp_get_wtime();
                grad = optimiser->step(i, shrinkage);
                duration += omp_get_wtime() - start;
                
                ii nnz = 0;
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
                    oss << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/L1/" << setfill('0') << setw(8) << i;
                    optimiser->write_h5(*h5out, oss.str(), scale_bases, is, js);
                }
            }
            
            // l0
            cout << " L0 threshold=" << fixed << setprecision(2) << thres << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;
            
            optimiser->threshold(thres);
            grad = DBL_MAX;
            for (ii i = 0; grad > tolerance; i++)
            {
                grad = optimiser->step(i, 0.0);
                
                ii nnz = 0;
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
                    oss << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/_debug/L0/" << setfill('0') << setw(8) << i;
                    optimiser->write_h5(*h5out, oss.str(), scale_bases, is, js);
                }
             }
            
            //////////////////////////////////////////////////////////////////////////////////
            // OUTPUT
			if (debug)
			{
				ostringstream oss;
				oss << config_id << "/" << rc0_mz << "/" << rcr << "/" << shr << "/" << tol << "/";
				optimiser->write_h5(*h5out, oss.str(), scale_bases, is, js);
			}

            // output coeffs csv
            ostringstream oss3; oss3 << id << "_" << config_id << "_" << rc0_mz << "_" << rcr << "_" << shr << "_" << tol << ".csv";
            cout << "Writing " << oss3.str() << endl;
            ofstream ofs(oss3.str().c_str());
            for (ii j = n_core_bases; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            for (ii y = 0; y < bases[j]->get_cm().n[1]; ++y)
            for (ii x = 0; x < bases[j]->get_cm().n[0]; ++x)
            if (optimiser->get_cs()[j][x + y*bases[j]->get_cm().n[0]] > 0.0)
            {
                ofs << bases[j]->get_cm().l[0] << ","
                    << bases[j]->get_cm().l[1] << ","
                    << bases[j]->get_cm().o[0] + x << ","
                    << bases[j]->get_cm().o[1] + y << ","
                    << optimiser->get_cs()[j][x + y*bases[j]->get_cm().n[0]] << endl;
            }
            ofs.close();
            
            ostringstream oss2; oss2 << id << "_" << config_id << "_" << rc0_mz << "_" << rcr << "_" << shr << "_" << tol << ".error.csv";
            //optimiser->calc_error(oss2.str());
            
            delete optimiser;
        }
        
    }
    
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
   
    ii minutes = (int) (duration/60);
    cout << "Finished, Duration: " << minutes << "m" << setprecision(3) << (duration - 60*minutes) << "s" << endl;
	
	if (debug) delete h5out; 
	omp_set_num_threads(_threads);
}

}


