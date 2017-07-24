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


#include <limits>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
#include "../core/DatasetSeamass.hpp"
#include "../peak/SMData.hpp"
#include "../peak/MathOperator.hpp"
#include "../peak/BsplineData.hpp"
#include "../peak/PeakOperator.hpp"
#include "../peak/PeakData.hpp"
#include "../peak/PeakManager.hpp"
#include "seamass-restore.hpp"

using namespace std;
using namespace kernel;
namespace po = boost::program_options;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
    try
#endif
    {
        string filePathIn;
        vector<char> scales(2);
        int debugLevel;
        bool peakWidth;
        double threshold;
        ii samplingRate;
        bool reconstruct;
        bool centroid;
        bool deconvolve;

        // *******************************************************************

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Restore an mzMLb/smv file to mzMLb/smb with optional reconstruction or centroiding.\n"
            "\n"
            "seamass-restore [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input file in mzMLv or smv format produced by 'seamass'.")
            ("sampling_rate,s",po::value<int>(&samplingRate)->default_value(6),
             "Number of data points to generate per unit b-spline. Default=6")
            ("deconvolve,v", po::bool_switch(&deconvolve)->default_value(false),
             "Output deconvolved peaks.")
            ("reconstruct,r", po::bool_switch(&reconstruct)->default_value(false),
             "Reconstruct using original m/z locations (ignores sample_rate, deconvolve, centroid).")
            ("centroid,c", po::bool_switch(&centroid)->default_value(false),
             "Output centroided data.")
            ("centroid_threshold,t", po::value<double>(&threshold)->default_value(1000.0),
             "Minimum ion count density to keep centroid. Default is 1000.")
            ("centroid_widths,w",po::bool_switch(&peakWidth)->default_value(false),
             "Output 'SMW' file containing centroid information including location,"
             "width.")
            ("debug,d", po::value<int>(&debugLevel)->default_value(0),
             "Debug level.")
        ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "seaMass-restore : Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;
        initKernel(debugLevel);

        Subject::setDebugLevel(debugLevel);
        Observer* observer = 0;
        if (debugLevel % 10 >= 1)
            Subject::registerObserver(observer = new Observer());

        ObserverMatrix* observerMatrix = 0;
        ObserverMatrixSparse* observerMatrixSparse = 0;
        if (debugLevel / 10 >= 1)
        {
            SubjectMatrix::registerObserver(observerMatrix = new ObserverMatrix());
            SubjectMatrixSparse::registerObserver(observerMatrixSparse = new ObserverMatrixSparse());
        }

        if (vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        string fileStemOut = boost::filesystem::path(filePathIn).stem().string();
        Dataset* dataset = FileFactory::createFileObj(filePathIn, fileStemOut, Dataset::WriteType::Input);
        if (!dataset)
            throw runtime_error("ERROR: Input file is missing or incorrect");

        Seamass::Input input;
        Seamass::Output output;
        string id;

        while (dataset->read(input, output, id))
        {
            if (debugLevel % 10 == 0)
                cout << "Processing " << id << endl;

            // load back into Seamass
            Seamass seamassCore(input, output);

            if (reconstruct)
            {
                seamassCore.getOutputBinCounts(input.counts);
                input.type = Seamass::Input::Type::Binned;
            }
            else if (centroid)
            {
                Seamass::ControlPoints contpts;
                seamassCore.getOutputControlPoints1d(contpts, deconvolve);

                // Now preform 1D Centroid
                cout << "Performing 1D centoiding of scans"<<endl;
                VecMat<double> mzPeak;
                VecMat<float> pkPeak;
                vector<size_t> mzpkVecSize;

                uli dims[2];
                double mzRes;
                double rtRes;
                vector<ii> offset=contpts.offset;
                mzRes=double(contpts.scale[0]);
                if (contpts.scale.size() > 1)
                    rtRes=double(contpts.scale[1]);
                else
                    rtRes=0;

                VecMat<float> rawCoeff(contpts.extent[1],contpts.extent[0],contpts.coeffs);
                rawCoeff.getDims(dims);

                SMData2D<OpUnit> A(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
                SMData2D<OpNablaH> dhA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);
                SMData2D<OpNabla2H> d2hA(dims,&offset[0],mzRes,rtRes,rawCoeff.v);

                for (size_t i = 0; i < A.rt.size(); ++i)
                {
                    A.rt[i] = input.startTimes[i];
                    dhA.rt[i] = input.startTimes[i];
                    d2hA.rt[i] = input.startTimes[i];
                }

                BsplineData<> bsData(A,dhA,d2hA);

                PeakManager<PeakData,BsplineData,Centroid2D> centriodPeak(bsData,threshold);
                centriodPeak.execute();
                centriodPeak.peak->getPeakMat(mzPeak, pkPeak, contpts.extent[1], mzpkVecSize);

                input.type = Seamass::Input::Type::Centroided;
                vector<double>().swap(input.locations);
                vector<fp>().swap(input.counts);
                vector<li>().swap(input.countsIndex);
                input.countsIndex.push_back(0);

                uli peakDims[2];
                mzPeak.getDims(peakDims);
                for (ii i = 0; i < mzpkVecSize.size(); i++)
                {
                    if (mzpkVecSize[i] > 0)
                    {
                        li idxOffset=li(i*peakDims[1]);
                        input.locations.insert(input.locations.end(), mzPeak.v.begin()+idxOffset,
                                               mzPeak.v.begin()+idxOffset+mzpkVecSize[i]);
                        input.counts.insert(input.counts.end(), pkPeak.v.begin()+idxOffset,
                                            pkPeak.v.begin()+idxOffset+mzpkVecSize[i]);
                    }
                    input.countsIndex.push_back(input.counts.size());
                }
                if (peakWidth == true)
                {
                    centriodPeak.peak->writePeakWidth(fileStemOut,NC_DOUBLE);
                }
            }
            else
            {
                cout << "Performing high resolution output of seaMass." << endl;

                Seamass::ControlPoints contpts;
                seamassCore.getOutputControlPoints1d(contpts, deconvolve);

                vector<double>().swap(input.locations);
                vector<fp>().swap(input.counts);
                vector<li>().swap(input.countsIndex);
                input.countsIndex.push_back(0);

                int csCol=contpts.extent[0];
                int csRow=contpts.extent[1];

                ii m = samplingRate-1;
                ii k = 4;
                ii n = ii(csCol) - k + 1;

                float Mval[16] = {0.1666667, 0.6666667, 0.1666667, 0.0, -0.5, 0.0,
                                  0.5, 0.0, 0.5, -1.0, 0.5, 0.0, -0.1666667, 0.5,
                                  -0.5, 0.1666667};

                float *M = alcMat(M, k, k);
                float *T = alcMat(T, m, k);
                float *TM = alcMat(TM,m,k);

                float dt = 1 / (float(m));
                vector<float> t(m);

                for (int i = 0; i < t.size(); ++i)
                    t[i] = float(i) * dt;

                for (int i = 0; i < k * k; ++i)
                    M[i] = Mval[i];

                int tidx = 0;
                for (int i = 0; i < (m * k); i = i + k, ++tidx)
                {
                    T[i] = 1.0;
                    T[i + 1] = t[tidx];
                    T[i + 2] = t[tidx] * t[tidx];
                    T[i + 3] = t[tidx] * t[tidx] * t[tidx];
                }

                matDmul(T,M,TM,m,k,k);

                input.counts.reserve(m*n*csRow);
                vector<double>().swap(input.locations);

                input.type = Seamass::Input::Type::Sampled;

                //ii cptr=0;
                for(int idx = 0; idx < csRow; ++idx)
                {
                    /* M = 1/6*[1,4,1,0;-3,0,3,0;3,-6,3,0;-1,3,-3,1]
                     * M =1/6 * | 1     4     1     0 |
                     *          |-3     0     3     0 |
                     *          | 3    -6     3     0 |
                     *          |-1     3    -3     1 |
                     *
                     * P=T*M*C
                     * TM=T*M
                     * P=TM*C
                     */

                    float *C = alcMat(C, k, n);
                    float *P = alcMat(P, m, n);
                    vector<double> mz;

                    ii crow,ccol,csize;
                    crow = k;
                    ccol = csCol - k + 1;
                    csize= n*k;

                    float **Cidx;
                    vector<float*> matidx;
                    matidx.reserve(crow);

                    for (ii i = 0; i < csize; i++) {
                        C[i] = 0;
                    }

                    for (int i=0; i < crow; ++i)
                    {
                        matidx.push_back(&C[i*ccol]);
                    }
                    Cidx=matidx.data();


                    ii cptr=idx*csCol;
                    for (ii j = 0; j < ccol; ++j)
                    {
                        for(ii i = 0; i < crow; ++i)
                        {
                            Cidx[i][j]=contpts.coeffs[cptr];
                            ++cptr;
                        }
                        cptr=cptr-crow+1;
                    }

                    matDmul(TM,C,P,m,k,n);

                    genMZAxis(mz,contpts,m*n,samplingRate-1);

                    input.locations.insert(input.locations.end(), mz.begin(), mz.end());
                    for (int j = 0; j < n; ++j)
                    {
                        for (int i = 0; i < m; ++i)
                        {
                            input.counts.push_back(P[j+i*n]);
                        }
                    }
                    input.countsIndex.push_back(input.locations.size());

                    delMat(C);
                    delMat(P);
                }
                delMat(M);
                delMat(T);
                delMat(TM);

            }

            dataset->write(input, id);

            if (debugLevel % 10 == 0)
                cout << endl;
        }

        delete dataset;
        cout << endl;
    }
#ifdef NDEBUG
    catch(exception& e)
    {
        cerr << e.what() << endl;
        cout << endl;
        return 1;
    }
#endif
    return 0;
}
