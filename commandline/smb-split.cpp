//
// Original author: Ranjeet Bhamber <ranjeet.bhamber <a.t> bristol.ac.uk>
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


#include "../io/FileNetcdf.hpp"
#include <kernel.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
using namespace std;
using namespace kernel;
namespace po = boost::program_options;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
    try
#endif
    {
        string filePath;
        int sections;
        int tiles;

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Splits a smb file into many smaller smb file for HPC processing and use for combining many smaller smb files back into origonal smb\n"
            "\n"
            "smb-split [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<string>(&filePath),
             "Input file in smb format. Use mzmlb2smb to convert from mzMLb.")
            ("tiles,s", po::value<int>(&sections)->default_value(6),
             "The number of sections the RT will be split into. This will generate (2*N-1) Tiles inorder to cover the RT range.")
        ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "smb-split : Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;

        if(vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        //boost::filesystem::path smbPathStem = boost::filesystem::path(filePath).stem();
        //cout << smbPathStem.string()<<endl;
        boost::filesystem::path fileIn(filePath);
        cout << fileIn.string()<<endl;
        cout << fileIn.filename().string()<<endl;
        cout << fileIn.stem().string()<<endl;
        cout << fileIn.extension().string()<<endl;




        boost::filesystem::create_directories(fileIn.stem());
        boost::filesystem::path pathOut;
        pathOut = fileIn.stem() / fileIn.stem();

        cout<<"file path out:"<<pathOut.string()<<endl;

        vector<li> csIndex;
        vector<li> mzIndex;
        vector<double> rtTimes;

        //FileNetcdf fileSmbIn(fileIn.filename().string());
        FileNetcdf fileSmbIn(fileIn.string());
        fileSmbIn.read_VecNC("startTimes",rtTimes);
        fileSmbIn.read_VecNC("countsIndex",csIndex);

        for(li i=0; i < csIndex.size(); ++i)
        {
            mzIndex.push_back(csIndex[i]+i);
        }


        size_t N = rtTimes.size();
        size_t dN = size_t(ceil(double(N)/double(sections)));
        //size_t dN = (N/tiles)+1;

        tiles = 2*sections-1;

        vector<double> binLocations;
        vector<fp> counts;
        vector<li> countsIndex;
        vector<fp>  exposures;
        vector<double> finishTimes;
        vector<double> startTimes;

        bool overlap = false;
        //for(int i = 0 ; i < sections ; ++i)
        int idx = 0;
        for(int i = 0 ; i < tiles; ++i)
        {
            //string fileOut = fileIn.stem().string()+"."+to_string(i)+".smb";
            string fileOut = pathOut.string()+"."+to_string(i)+".smb";
            cout<<"Generating file: "<<fileOut<<endl;
            FileNetcdf fileSmbOut(fileOut,NC_NETCDF4);

            binLocations.clear();
            counts.clear();
            countsIndex.clear();
            exposures.clear();
            finishTimes.clear();
            startTimes.clear();

            size_t beginIdx, endIdx;
            size_t csBeginIdx, csLen;
            size_t mzBeginIdx, mzLen;
            size_t len;
            size_t lenP1;

            if (overlap)
            {
                beginIdx = dN * idx + dN/2;
                //if (i == sections - 1)
                //{
                //    endIdx = N;
                //    len = N - beginIdx;
                //}
                //else
                //{
                endIdx = beginIdx + dN;
                len = dN;
                //}
                lenP1 = len + 1;
                ++idx;
                overlap = false;
            }
            else
            {
                beginIdx = dN*idx;
                if ( i == tiles - 1)
                {
                    endIdx   = N;
                    len = N-beginIdx;
                }
                else
                {
                    //endIdx   = dN*idx+dN;
                    endIdx   = beginIdx + dN;
                    len = dN;
                }
                lenP1=len+1;
                overlap = true;
            }

            csBeginIdx = size_t(csIndex[beginIdx]);
            csLen      = size_t(csIndex[endIdx]) - csBeginIdx;
            mzBeginIdx = size_t(mzIndex[beginIdx]);
            mzLen      = size_t(mzIndex[endIdx]) - mzBeginIdx;

            fileSmbIn.read_HypVecNC("startTimes",startTimes,&beginIdx,&len);
            fileSmbIn.read_HypVecNC("finishTimes",finishTimes,&beginIdx,&len);
            fileSmbIn.read_HypVecNC("exposures",exposures,&beginIdx,&len);
            fileSmbIn.read_HypVecNC("countsIndex",countsIndex,&beginIdx,&lenP1);
            fileSmbIn.read_HypVecNC("counts",counts,&csBeginIdx,&csLen);
            fileSmbIn.read_HypVecNC("binLocations",binLocations,&mzBeginIdx,&mzLen);

            fileSmbOut.write_VecNC("binLocations",binLocations,NC_DOUBLE);
            fileSmbOut.write_VecNC("counts",counts,NC_FLOAT);
            fileSmbOut.write_VecNC("countsIndex",countsIndex,NC_INT64);
            fileSmbOut.write_VecNC("exposures",exposures,NC_FLOAT);
            fileSmbOut.write_VecNC("finishTimes",finishTimes,NC_DOUBLE);
            fileSmbOut.write_VecNC("startTimes",startTimes,NC_DOUBLE);

            fileSmbOut.close();
        }

    }
#ifdef NDEBUG
    catch(exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
#endif

    return 0;
}
