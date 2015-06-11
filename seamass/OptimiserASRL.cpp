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


#include "OptimiserASRL.hpp"
#include <fstream>
#include <algorithm>


OptimiserASRL::
OptimiserASRL(const vector<Basis*>& _bases,
              vector<fp>& _gs,
              ii _accell) :
    bases(_bases),
    gs(_gs),
    accell(_accell),
    iteration(0)
{
    // pre-compute weights
    wcs.resize(bases.size());
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        wcs[j].resize(bases[j]->get_cm().size());
        if (j == 0)
        {
            vector<fp> ts(gs.size(), 1.0);
            bases[j]->analysis(wcs[0], ts);
        }
        else
        {
            bases[j]->analysis(wcs[j], wcs[bases[j]->get_parent()->get_index()]);
        }
    }
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient()) vector<fp>().swap(wcs[j]);
    }
    
    // pre-compute l2norm
    l2.resize(bases.size());
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        l2[j].resize(bases[j]->get_cm().size());
        if (j == 0)
        {
            vector<fp> ts(gs.size(), 1.0);
            bases[j]->l2norm(l2[0], ts);
        }
        else
        {
            bases[j]->l2norm(l2[j], l2[bases[j]->get_parent()->get_index()]);
        }
    }
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient())
        {
            vector<fp>().swap(l2[j]);
        }
        else
        {   // sqrt
            vsSqrt(bases[j]->get_cm().size(), &(l2[j][0]), &(l2[j][0]));
        }
    }
    
    // starting cs from ng
    cs.resize(bases.size());
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        cs[j].resize(bases[j]->get_cm().size());
        if (j == 0)
        {
            bases[j]->analysis(cs[0], gs);
        }
        else
        {
            bases[j]->analysis(cs[j], cs[bases[j]->get_parent()->get_index()]);
        }
    }
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient()) vector<fp>().swap(cs[j]);
    }
    // any bases that are too small must be deleted to avoid numerical problems
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (!bases[j]->is_transient())
        for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
        if (wcs[j][i] < 0.001) cs[j][i] = 0.0;
    }
    
    // temporaries required for acceleration
    if (accell >= 1)
    {
        c0s.resize(bases.size());
        u0s.resize(bases.size());
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        {
            c0s[j].resize(bases[j]->get_cm().size());
            u0s[j].resize(bases[j]->get_cm().size());
        }
    }
    if (accell >= 2)
    {
        q0s.resize(bases.size());
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        {
            q0s[j].resize(bases[j]->get_cm().size());
        }
    }
    
    // how much memory are we using?
    ii size = 0;
    for (ii j = 0; j < cs.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < wcs.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < l2.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < c0s.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < u0s.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < q0s.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    cout << endl << "Vector Extrapolated Sparse Richardson-Lucy mem=" << fixed << setprecision(2) << (size*sizeof(fp) + sizeof(this))/(1024.0*1024.0) << "Mb" << endl;
}


OptimiserASRL::
~OptimiserASRL()
{
}


void
OptimiserASRL::
threshold(double thresh)
{
    for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
    #pragma omp parallel for
    for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
    {
        if (cs[j][i] < thresh) cs[j][i] = 0.0;
    }
    
    double volume = bases.front()->get_volume();
    
    for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
    #pragma omp parallel for
    for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
    {
        cs[j][i] /= volume;
    }
}


double
OptimiserASRL::
step(ii iteration, double shrinkage)
{
    double synthesis_start = omp_get_wtime();
    // synthesis except root
    vector< vector<fp> > ts(bases.size());
    for (ii j = (ii) bases.size() - 1; j > 0; j--)
    {
        ii pj = bases[j]->get_parent()->get_index();
        if (ts[pj].size() == 0)
        if (bases[pj]->is_transient())
        {
            ts[pj].resize(bases[pj]->get_cm().size(), 0.0);
        }
        else
        {
            ts[pj].resize(bases[pj]->get_cm().size());
            for (ii i = 0; i < (ii) bases[pj]->get_cm().size(); i++) ts[pj][i] = cs[pj][i];
        }

        if (ts[j].size() == 0)
        {
            bases[j]->synthesis(ts[pj], cs[j]);
        }
        else
        {
            bases[j]->synthesis(ts[pj], ts[j]);
            
        }
        vector<fp>().swap(ts[j]);
    }
    //synthesis root
    vector<fp> fs(gs.size());
    if (ts.front().size() == 0)
    {
        bases.front()->synthesis(fs, cs.front());
    }
    else
    {
        bases.front()->synthesis(fs, ts.front());
    }
    vector<fp>().swap(ts.front());
    // timing
    double synthesis_duration = omp_get_wtime() - synthesis_start;
    static double culm_synthesis_duration = 0.0;
    culm_synthesis_duration += synthesis_duration;
    
    double error_start = omp_get_wtime();
    // Poisson likelihood error
    bases.front()->error(fs, gs);
    // timing
    double error_duration = omp_get_wtime() - error_start;
    static double culm_error_duration = 0.0;
    culm_error_duration += error_duration;
    
    double analysis_start = omp_get_wtime();
    // analysis root
    vector< vector<fp> > es(bases.size());
    es.front().resize(bases.front()->get_cm().size());
    bases.front()->analysis(es.front(), fs);
    vector<fp>().swap(fs);
    // analysis except root
    for (ii j = 1; j < (ii) bases.size(); j++)
    {
        es[j].resize(bases[j]->get_cm().size());
        bases[j]->analysis(es[j], es[bases[j]->get_parent()->get_index()]);
        // use child count here to delete transient es sooner
    }
    // timing
    double analysis_duration = omp_get_wtime() - analysis_start;
    static double culm_analysis_duration = 0.0;
    culm_analysis_duration += analysis_duration;

    double shrinkage_start = omp_get_wtime();
    // shrinkage
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient())
        {
            vector<fp>().swap(es[j]);
        }
        else
        {
            #pragma omp parallel for
            for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
            if (es[j][i] >= FLT_MIN && cs[j][i] >= FLT_MIN)
            {
                es[j][i] *= cs[j][i] / (shrinkage * l2[j][i] + wcs[j][i]);
            }
            else
            {
                es[j][i] = 0.0;
            }
        }
    }
    double shrinkage_duration = omp_get_wtime() - shrinkage_start;
    static double culm_shrinkage_duration = 0.0;
    culm_shrinkage_duration += shrinkage_duration;
    
    double accel_start = omp_get_wtime();
    // unaccelerated update
    double sum = 0.0;
    double sumd = 0.0;
    fp a = 0.0;
    if (accell == 0)
    {
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        #pragma omp parallel for reduction(+:sum,sumd)
        for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
        {
            sum += cs[j][i] * cs[j][i];
            sumd += (es[j][i] - cs[j][i])*(es[j][i] - cs[j][i]);

            cs[j][i] = es[j][i];
        }        
    }
    else // accelerated update    
    {
        // init/update u0s and compute acceleration factor a
        if (iteration == 0)
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            #pragma omp parallel for
            for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
            if (cs[j][i] > 0.0)
            {
                u0s[j][i] = es[j][i]/cs[j][i];
            }
            a = 0.0;
        }
        else
        {
            double sum_u0u = 0.0;
            double sum_u0u0 = 0.0;
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            #pragma omp parallel for reduction(+:sum_u0u,sum_u0u0)
            for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
            if (cs[j][i] > 0.0)
            {
                double old_u0 = u0s[j][i];
                u0s[j][i] = es[j][i]/cs[j][i];
                
                old_u0 = old_u0 > 0.0 ? c0s[j][i]*log(old_u0) : 0.0;
                sum_u0u += old_u0 * (u0s[j][i] > 0.0 ? es[j][i]*log(u0s[j][i]) : 0.0);
                sum_u0u0 += old_u0 * old_u0;
            }
            a = sqrt(sum_u0u/sum_u0u0);
            a = a > 0.0f ? a : 0.0f;
            a = a < 1.0f ? a : 1.0f;
        }    
        
        if (iteration == 0) // unaccelerated this time, but save es as c0s
        {            
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
                if (cs[j][i] > 0.0)
                {
                    sum += cs[j][i] * cs[j][i];
                    sumd += (es[j][i] - cs[j][i])*(es[j][i] - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = es[j][i];
                    // for next itteration
                    c0s[j][i] = es[j][i];
                }
                vector<fp>().swap(es[j]);
            }
        }
        else if (accell == 1) // linear vector extrapolation
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
                if (cs[j][i] > 0.0)
                {
                    fp c1 = es[j][i] * powf(es[j][i]/c0s[j][i], a);
                    
                    sum += cs[j][i] * cs[j][i];
                    sumd += (c1 - cs[j][i])*(c1 - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = c1;
                    
                    // for next itteration
                    c0s[j][i] = es[j][i];
                }
                vector<fp>().swap(es[j]);
            }
        }
        else if (iteration == 1) // linear vector extrapolation this time, but save the qs
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
                if (cs[j][i] > 0.0)
                {
                    fp q = es[j][i]/c0s[j][i];
                    fp c1 = es[j][i] * powf(es[j][i]/c0s[j][i], a);
                    
                    sum += cs[j][i] * cs[j][i];
                    sumd += (c1 - cs[j][i])*(c1 - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = c1;
                    
                    // for next itteration
                    c0s[j][i] = es[j][i];
                    q0s[j][i] = q;
                }
                vector<fp>().swap(es[j]);
            }
        }
        else // quadratic vector extrapolation
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
                if (cs[j][i] > 0.0)
                {
                    fp q = es[j][i]/c0s[j][i];
                    fp c1 = es[j][i] * powf(q, a) * powf(q/q0s[j][i], 0.5f*a*a);
                    
                    sum += cs[j][i] * cs[j][i];
                    sumd += (c1 - cs[j][i])*(c1 - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = c1;
                    
                    // for next itteration
                    c0s[j][i] = es[j][i];
                    q0s[j][i] = q;
                }
                vector<fp>().swap(es[j]);
            }
        }
    }
    double accel_duration = omp_get_wtime() - accel_start;
    static double culm_accel_duration = 0.0;
    culm_accel_duration += accel_duration;
    
    //cout << "Iteration   Durations: synthesis=" << setprecision(2) << synthesis_duration << " error=" << error_duration << " analysis=" << analysis_duration << " shrinkage=" << shrinkage_duration << " accel=" << accel_duration << " all=" << synthesis_duration+error_duration+analysis_duration+shrinkage_duration+accel_duration << endl;
    //cout << "Culminative Durations: synthesis=" << culm_synthesis_duration << " error=" << culm_error_duration << " analysis=" << culm_analysis_duration << " shrinkage=" << culm_shrinkage_duration << " accel=" << culm_accel_duration<< " all=" << culm_synthesis_duration+culm_error_duration+culm_analysis_duration+culm_shrinkage_duration +culm_accel_duration << endl;
    
    return sqrt(sumd)/sqrt(sum);
}


fp
OptimiserASRL::
compute_norm_max_counts(ii n_core_bases)
{
    // synthesis except root
	fp max_counts = 0.0;
    vector< vector<fp> > ts(bases.size());
    for (ii j = (ii) bases.size() - 1; j >= n_core_bases; j--)
    {
        ii pj = bases[j]->get_parent()->get_index();
        if (ts[pj].size() == 0)
        if (bases[pj]->is_transient())
        {
            ts[pj].resize(bases[pj]->get_cm().size(), 0.0);
        }
        else
        {
            ts[pj].resize(bases[pj]->get_cm().size());
            for (ii i = 0; i < (ii) bases[pj]->get_cm().size(); i++) ts[pj][i] = cs[pj][i];
        }

        if (ts[j].size() == 0)
        {
            bases[j]->synthesis(ts[pj], cs[j]);
        }
        else
        {
            bases[j]->synthesis(ts[pj], ts[j]);
            
        }
        vector<fp>().swap(ts[j]);

		if (j == n_core_bases)
		{
			for (ii i = 0; i < ts[pj].size(); ++i)
			{
				max_counts = max_counts > ts[pj][i] ? max_counts : ts[pj][i];
			}
			// normalise
			max_counts *=  pow(2.0, bases[j]->get_cm().l[0]) * pow(2.0, bases[j]->get_cm().l[1]);
			break;
		}
    }
	return max_counts;
}


void
OptimiserASRL::
write_h5(const SMOWriter& file, const string& datafilename, const vector<ii>& scale_bases,
         const vector<li>& is, const vector<ii>& js, const vector<fp>& gains)
{
    vector<double> tic(js.size(),0);
    vector<double> rtic(js.size(),0);

    //ii k = scale_bases.size() - 1;
    ii k = 0;
    vector< vector<fp> > ts(bases.size());
    for (ii j = (ii) bases.size() - 1; j >= 0; j--)
    {
        ii pj = 0;
        if (j > 0)
        {
            pj = bases[j]->get_parent()->get_index();
            if (ts[pj].size() == 0)
            if (bases[pj]->is_transient())
            {
                ts[pj].resize(bases[pj]->get_cm().size(), 0.0);
            }
            else
            {
                ts[pj].resize(bases[pj]->get_cm().size());
                for (ii i = 0; i < (ii) bases[pj]->get_cm().size(); i++)
                    ts[pj][i] = cs[pj][i];
            }
        }

        vector<fp>* p = (ts[j].size() == 0) ? &cs[j] : &ts[j];

        if (j == scale_bases[k])
        {
            // write 2D B-spline coefficients
            ostringstream oss;
            oss << datafilename << "/cs";
            file.write_cs(oss.str(), bases[j]->get_cm(), *p);
            if (k > 0) k--;
        }

        if (j > 0)
        {
            bases[j]->synthesis(ts[pj], *p);
        }
        else
        {
            // write residual 1D B-spline coefficients
            ostringstream oss;
            oss << datafilename << "fcs";
            file.write_cs(oss.str(), bases[0]->get_cm(), *p);

            // write fs
            vector<fp> fs(gs.size(),0);
            if (ts.front().size() == 0)
            {
                bases.front()->synthesis(fs, cs.front());
            }
            else
            {
                bases.front()->synthesis(fs, ts.front());
            }
            // write fs spectrum intensities
            ostringstream oss2;
            oss2 << datafilename;
            file.write_cdata(oss2.str(), fs,"fSpectrumCount");

            // write gs spectrum intensities
            ostringstream oss3;
            oss3 << datafilename.substr(0, datafilename.find("/",datafilename.find("/",1)+1)+1);
            file.write_cdata(oss3.str(), gs,"gSpectrumCount");

            // write mz scan index for spectrum intensities
            file.write_cdata(oss3.str(), is,"SpectrumIndex");

			// write mz scan index for spectrum intensities
			file.write_cdata(oss3.str(), gains, "SpectrumGain");

            for (ii j = 0; j < js.size(); j++)
            for (ii i = is[j]; i < is[j+1]; i++)
            {
                tic[j] += fs[i];
            }
            //for (ii j = 0; j < js.size(); j++)
            //    cout << tic[j] << endl;

            // detect and write centroids (todo)
            /*ostringstream oss; oss << id << "_d.h5";

            CoeffsMetadata dcm = bases[j]->get_cm();
            dcm.n[0]--;
            vector<fp> ds(dcm.size());

            cout << bases[j]->get_cm().n[0] << "," << bases[j]->get_cm().n[1] << endl;
            cout << dcm.n[0] << "," << dcm.n[1] << endl;

            for (ii y = 0; y < dcm.n[1]; y++)
            for (ii x = 0; x < dcm.n[0]; x++)
            {
                ds[x+y*dcm.n[0]] = (*p)[x+1+y*bases[j]->get_cm().n[0]] - (*p)[x+y*bases[j]->get_cm().n[0]];
            }
            write_h5(oss.str(), "Image", dcm, ds);*/
        }

        vector<fp>().swap(ts[j]);
    }

    // residuals
    /*
    for (ii j = scale_bases.front() - 1; j >= 0; j--)
    {
        ii pj = 0;
        if (j > 0)
        {
            pj = bases[j]->get_parent()->get_index();
            if (ts[pj].size() == 0)
            if (bases[pj]->is_transient() || pj > scale_bases.front())
            {
                ts[pj].resize(bases[pj]->get_cm().size(), 0.0);
            }
            else
            {
                ts[pj].resize(bases[pj]->get_cm().size());
                for (ii i = 0; i < bases[pj]->get_cm().size(); i++)
                    ts[pj][i] = cs[pj][i];
            }
        }

        vector<fp>* p = (ts[j].size() == 0) ? &cs[j] : &ts[j];

        if (j > 0)
        {
            bases[j]->synthesis(ts[pj], *p);
        }
        else
        {
            // write residual 1D B-spline coefficients
            ostringstream oss;
            oss << datafilename << "rcs";
            file.write_cs(oss.str(), bases[0]->get_cm(), *p);

            vector<fp> rs(gs.size(),0);
            if (ts.front().size() == 0)
            {
                bases.front()->synthesis(rs, cs.front());
            }
            else
            {
                bases.front()->synthesis(rs, ts.front());
            }
            ostringstream oss2;
            oss2 << datafilename << "rs";
            file.write_fs(oss2.str(), rs, is, js);

            for (ii j = 0; j < js.size(); j++)
            for (ii i = is[j]; i < is[j+1]; i++)
            {
                rtic[j] += rs[i];
            }
            //for (ii j = 0; j < js.size(); j++)
            //    cout << tic[j] << ":" << rtic[j] << endl;

       }

        vector<fp>().swap(ts[j]);
    }*/
}


void
OptimiserASRL::
calc_error(const std::string& id)
{
    vector<fp> tosort;
    for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
        tosort.insert(tosort.end(),cs[j].begin(),cs[j].end());
    std::sort(tosort.begin(), tosort.end());

    cout << "Writing " << id << endl;
    ofstream ofs(id.c_str());

    for (ii i = tosort.size() % 1000; i < (ii) tosort.size(); i += 1000)
    {
        double thresh = tosort[i];
        if (thresh == 0.0) continue;
        
        ii n = 0;
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        //#pragma omp parallel for reduction(+:sum_u0u,sum_u0u0)
        for (ii i = 0; i < (ii) bases[j]->get_cm().size(); i++)
        if (cs[j][i] < thresh)
        {
            cs[j][i] = 0.0;
        }
        else
        {
            n++;
        }
    
        vector< vector<fp> > ts(bases.size());
        for (ii j = (ii) bases.size() - 1; j >= 0; j--)
        {
            ii pj = 0;
            if (j > 0)
            {
                pj = bases[j]->get_parent()->get_index();
                if (ts[pj].size() == 0)
                if (bases[pj]->is_transient())
                {
                    ts[pj].resize(bases[pj]->get_cm().size(), 0.0);
                }
                else
                {
                    ts[pj].resize(bases[pj]->get_cm().size());
                    for (ii i = 0; i < bases[pj]->get_cm().size(); i++)
                        ts[pj][i] = cs[pj][i];
                }
            }
            
            vector<fp> * p = (ts[j].size() == 0) ? &cs[j] : &ts[j];
            if (j > 0) bases[j]->synthesis(ts[pj], *p);
            
            else if (j==0)
            {
                vector<fp> fs(gs.size());
                if (ts.front().size() == 0)
                {
                    bases.front()->synthesis(fs, cs.front());
                }
                else
                {
                    bases.front()->synthesis(fs, ts.front());
                }
                vector<fp>().swap(ts.front());
                
                ii size_d = 0;
                double err = 0.0;
                //#pragma omp parallel for simd reduction(+:err,size_d)
                for (ii i = 0; i < (ii) gs.size(); i++)
                {
                    if (fs[i] > 0.0 && gs[i] >= 0.0)
                    {
                        err += fabs(gs[i]-fs[i]);
                        size_d++;
                    }
                }
                err = err / size_d;
                ofs << thresh << "," << n << "," << err << endl;
            }
            
            vector<fp>().swap(ts[j]);
        }
    }
    ofs.close();

}




