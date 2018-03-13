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
#include <boost/program_options.hpp>
#include "smb-split.hpp"
#include <algorithm>

using namespace std;
using namespace kernel;
namespace po = boost::program_options;


int main(int argc, const char * const * argv)
{
#ifdef NDEBUG
    try
#endif
    {
        string fileDirPath;
        int sections;
        bool restore;

        po::options_description general(
            "Usage\n"
            "-----\n"
            "Splits a smb file into many smaller smb file for HPC processing and use for combining many smaller smb files back into origonal smb\n"
            "\n"
            "smb-split [OPTIONS...] <file>\n"
        );

        general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<string>(&fileDirPath),
             "Input file in smb format or directory containing smb files. Use mzmlb2smb to convert from mzMLb.")
            ("tiles,s", po::value<int>(&sections)->default_value(6),
             "The number of sections the RT will be split into. This will generate (2*N-1) Tiles in order to cover the RT range.")
            ("restore,r", po::bool_switch(&restore)->default_value(false),
             "Restore split smb files into single file. The --file, -f, switch must be set and be a valid directory containing smb files."
             "Convert smv files to smb using seamass-restore.")
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
        boost::filesystem::path fileIn(fileDirPath);
        //cout << fileIn.string()<<endl;
        //cout << fileIn.filename().string()<<endl;
        //cout << fileIn.stem().string()<<endl;
        //cout << fileIn.extension().string()<<endl;

        if (restore == true && !boost::filesystem::is_directory(fileDirPath))
        {
            cout<<"Error: Invalid directory given for smb restore."<<endl;
            exit(1);
        }

        if (boost::filesystem::is_regular_file(fileDirPath))
        {
            cout<<"processing File: "<<fileIn.filename()<<endl;
            createSmbTiles(sections, fileIn);
        }
        else if (restore == false && boost::filesystem::is_directory(fileDirPath))
        {
            for (boost::filesystem::directory_entry &file : boost::filesystem::directory_iterator(fileDirPath))
            {
                cout<<"Checking directory for smb"<<endl;
                if (file.path().extension().string() == ".smb")
                {
                    cout<<"processing File: "<<file.path()<<endl;
                    createSmbTiles(sections,file.path());
                }
                else
                {
                    cout<<"Not a smb file: "<<file.path()<<endl;
                }
            }
        }
        else if (restore == true && boost::filesystem::is_directory(fileDirPath))
        {
            FileNetcdf fileSmbTile;
            vector<SmbTile> smbTileList;
            cout<<"Restoring Split smb file into single file"<<endl;
            for (boost::filesystem::directory_entry &file : boost::filesystem::directory_iterator(fileDirPath))
            {
                cout<<"Checking directory for smb"<<endl;
                if (file.path().extension().string() == ".smb")
                {
                    vector<double> rtTimes;
                    vector<li> csIdx;
                    cout<<"processing File: "<<file.path().filename()<<endl;
                    fileSmbTile.open(file.path().string());
                    fileSmbTile.read_VecNC("startTimes",rtTimes);
                    fileSmbTile.read_VecNC("countsIndex",csIdx);

                    //SmbTile newData = {file.path(), rtTimes, csIdx, int(smbTileList.size()),0,0,0};
                    SmbTile newData(file.path(), rtTimes,csIdx);
                    smbTileList.push_back(newData);

                    smbTileList.back().id=getSmbIndex<int>(smbTileList.back().fileName.stem().string());
                    fileSmbTile.close();
                }
                else
                {
                    cout<<"Not a smb file: "<<file.path()<<endl;
                }
            }

            sort(smbTileList.begin(),smbTileList.end(),
                 [](const SmbTile &x, const SmbTile &y)->bool
                 {
                     return x.id < y.id;
                 }
            );

            for (SmbTile i: smbTileList)
            {
                cout<<"Name of tile smbFile: "<<i.fileName.filename().string();
                cout<<"    rt Begin: " <<i.startTime.front()<<"\trt End: "<<i.startTime.back()<<"    Size: "<<i.startTime.size()<<endl;
            }


            // Find out what mode data is in, re-sampled or centroid mode
            bool centroidMode;
            fileSmbTile.open(smbTileList.front().fileName.string());
            fileSmbTile.search_Group("binLocations");
            fileSmbTile.search_Group("centroidLocations");
            vector<InfoGrpVar> dataSetInfo = fileSmbTile.get_Info();
            fileSmbTile.close();

            if (dataSetInfo.size() == 1)
            {
                if (dataSetInfo.front().varName == "centroidLocations")
                {
                    cout<<"Processing Centroid mode data..."<<endl;
                    centroidMode = true;
                }
                else if (dataSetInfo.front().varName == "binLocations")
                {
                    cout<<"Processing Profile re-sampled mode data..."<<endl;
                    centroidMode = false;
                }
                else
                {
                    cout<<"Error: Unsupported Input mode data format..."<<endl;
                    exit(1);
                }
            }
            else
            {
                cout<<"Error: Invalid Input data files."<<endl;
                cout<<"       Input Data files need to be in either profile or centroid mode format."<<endl;
                exit(1);
            }

            li drt = smbTileList.front().startTime.size()/4;

            li offsetTile=0;

            smbTileList.front().begIdx=0;
            smbTileList.front().endIdx=drt*3;
            smbTileList.front().len=smbTileList.front().endIdx - smbTileList.front().begIdx+1;
            smbTileList.front().offset=smbTileList.front().countsIndex.front();
            offsetTile=smbTileList.front().offset;

            smbTileList.front().csBeg=smbTileList.front().countsIndex[smbTileList.front().begIdx];
            smbTileList.front().csEnd=smbTileList.front().countsIndex[smbTileList.front().endIdx+1];
            smbTileList.front().csLen=smbTileList.front().csEnd - smbTileList.front().csBeg;
            smbTileList.front().csIdxLen=smbTileList.front().len;
            if (centroidMode == true)
            {
                smbTileList.front().mzBeg=smbTileList.front().countsIndex[smbTileList.front().begIdx];
                smbTileList.front().mzEnd=smbTileList.front().countsIndex[smbTileList.front().endIdx+1];
                smbTileList.front().mzLen=smbTileList.front().csEnd - smbTileList.front().csBeg;
            }
            else
            {
                smbTileList.front().mzBeg=smbTileList.front().csBeg + smbTileList.front().begIdx;
                smbTileList.front().mzEnd=smbTileList.front().csEnd + smbTileList.front().endIdx;
                smbTileList.front().mzLen=smbTileList.front().mzEnd - smbTileList.front().mzBeg+1;
            }


            for (li i = 1; i < smbTileList.size() - 1; ++i)
            {
                ptrdiff_t posIdx = distance(smbTileList[i].startTime.begin(),
                                            find(smbTileList[i].startTime.begin(),
                                                 smbTileList[i].startTime.end(),
                                                 smbTileList[i-1].startTime[smbTileList[i-1].endIdx]));

                ptrdiff_t offsetIdx = distance(smbTileList[i-1].startTime.begin(),
                                            find(smbTileList[i-1].startTime.begin(),
                                                 smbTileList[i-1].startTime.end(),
                                                 smbTileList[i].startTime.front()));

                smbTileList[i].begIdx=posIdx+1;
                smbTileList[i].endIdx=drt*3;
                smbTileList[i].len=smbTileList[i].endIdx - smbTileList[i].begIdx+1;
                //smbTileList[i].offset=smbTileList[i-1].countsIndex[offsetIdx];
                offsetTile += smbTileList[i-1].countsIndex[offsetIdx];
                smbTileList[i].offset=offsetTile;

                smbTileList[i].csBeg=smbTileList[i].countsIndex[smbTileList[i].begIdx];
                smbTileList[i].csEnd=smbTileList[i].countsIndex[smbTileList[i].endIdx+1];
                smbTileList[i].csLen=smbTileList[i].csEnd - smbTileList[i].csBeg;
                smbTileList[i].csIdxLen=smbTileList[i].len;
                if (centroidMode == true)
                {
                    smbTileList[i].mzBeg=smbTileList[i].countsIndex[smbTileList[i].begIdx];
                    smbTileList[i].mzEnd=smbTileList[i].countsIndex[smbTileList[i].endIdx+1];
                    smbTileList[i].mzLen=smbTileList[i].csEnd - smbTileList[i].csBeg;
                }
                else
                {
                    smbTileList[i].mzBeg=smbTileList[i].csBeg + smbTileList[i].begIdx;
                    smbTileList[i].mzEnd=smbTileList[i].csEnd + smbTileList[i].endIdx;
                    smbTileList[i].mzLen=smbTileList[i].mzEnd - smbTileList[i].mzBeg+1;
                }
            }

            smbTileList.back().begIdx=distance(smbTileList.back().startTime.begin(),
                                          find(smbTileList.back().startTime.begin(),
                                               smbTileList.back().startTime.end(),
                                               smbTileList[smbTileList.size()-2].startTime[smbTileList[smbTileList.size()-2].endIdx]))+1;

            smbTileList.back().endIdx=smbTileList.back().startTime.size()-1;
            smbTileList.back().len=smbTileList.back().endIdx - smbTileList.back().begIdx+1;

            smbTileList.back().csBeg=smbTileList.back().countsIndex[smbTileList.back().begIdx];
            smbTileList.back().csEnd=smbTileList.back().countsIndex[smbTileList.back().endIdx+1];
            smbTileList.back().csLen=smbTileList.back().csEnd - smbTileList.back().csBeg;
            smbTileList.back().csIdxLen=smbTileList.back().len+1;
            if (centroidMode == true)
            {
                smbTileList.back().mzBeg=smbTileList.back().countsIndex[smbTileList.back().begIdx];
                smbTileList.back().mzEnd=smbTileList.back().countsIndex[smbTileList.back().endIdx+1];
                smbTileList.back().mzLen=smbTileList.back().csEnd - smbTileList.back().csBeg;
            }
            else
            {
                smbTileList.back().mzBeg=smbTileList.back().csBeg + smbTileList.back().begIdx;
                smbTileList.back().mzEnd=smbTileList.back().csEnd + smbTileList.back().endIdx;
                smbTileList.back().mzLen=smbTileList.back().mzEnd - smbTileList.back().mzBeg+1;
            }

            //smbTileList.back().offset=smbTileList[smbTileList.size()-2].countsIndex[distance(smbTileList[smbTileList.size()-2].startTime.begin(),
            //                             find(smbTileList[smbTileList.size()-2].startTime.begin(),
            //                                 smbTileList[smbTileList.size()-2].startTime.end(),
            //                                 smbTileList.back().startTime.front()))];
            offsetTile += smbTileList[smbTileList.size()-2].countsIndex[distance(smbTileList[smbTileList.size()-2].startTime.begin(),
                                          find(smbTileList[smbTileList.size()-2].startTime.begin(),
                                               smbTileList[smbTileList.size()-2].startTime.end(),
                                               smbTileList.back().startTime.front()))];
            smbTileList.back().offset = offsetTile;

            cout<<"Boundaries for Subsection tiles:"<<endl;
            for (SmbTile i: smbTileList)
            {
                cout<<"BegIdx: "<<i.begIdx<<"    Val: "<<i.startTime[i.begIdx]<<"\tEndIdx: "<<i.endIdx<<"    Val: "<<i.startTime[i.endIdx];
                cout<<"    offset: "<<i.offset<<endl;
            }

            // Begin output of Tiles into one file.
            vector<double> binLocations;
            vector<fp> counts;
            vector<li> countsIndex;
            vector<fp>  exposures;
            vector<double> startTimes;
            vector<double> finishTimes;

            string fileNameOut = smbTileList.front().fileName.stem().string().substr(0,
                smbTileList.front().fileName.stem().string().find_last_of(".")) + ".smb";

            for (li i = 0; i < smbTileList.size(); ++i)
            {
                vector<double> _binLocations;
                vector<fp> _counts;
                vector<li> _countsIndex;
                vector<fp>  _exposures;
                vector<double> _startTimes;
                vector<double> _finishTimes;

                FileNetcdf smbTileIn(smbTileList[i].fileName.string());
                SmbTile *ptr = &smbTileList[i];

                //cout<<"Injecting File "<<i+1<<"/"<<smbTileList.size()<<": "<<smbTileList[i].fileName.string()<<"\r";
                cout<<"Injecting File "<<i+1<<"/"<<smbTileList.size()<<": "<<smbTileList[i].fileName.string()<<endl;

                size_t begptr = size_t(smbTileList[i].begIdx);
                size_t lenptr = size_t(smbTileList[i].len);
                smbTileIn.read_HypVecNC("startTimes", _startTimes, &begptr, &lenptr);
                smbTileIn.read_HypVecNC("finishTimes", _finishTimes, &begptr, &lenptr);
                smbTileIn.read_HypVecNC("exposures", _exposures, &begptr, &lenptr);
                begptr = size_t(smbTileList[i].begIdx);
                lenptr = size_t(smbTileList[i].csIdxLen);
                smbTileIn.read_HypVecNC("countsIndex", _countsIndex, &begptr, &lenptr);
                begptr = size_t(smbTileList[i].csBeg);
                lenptr = size_t(smbTileList[i].csLen);
                smbTileIn.read_HypVecNC("counts", _counts, &begptr, &lenptr);
                begptr = size_t(smbTileList[i].mzBeg);
                lenptr = size_t(smbTileList[i].mzLen);
                if (centroidMode == true)
                {
                    smbTileIn.read_HypVecNC("centroidLocations", _binLocations, &begptr,&lenptr);
                }
                else
                {
                    smbTileIn.read_HypVecNC("binLocations", _binLocations, &begptr,&lenptr);
                }

                smbTileIn.close();

                for (li x = 0; x < _countsIndex.size(); ++x)
                {
                    li offSet = smbTileList[i].offset;
                    _countsIndex[x] = _countsIndex[x]  + offSet;
                }

                /*
                li offSet=smbTileList[i].offset;
                transform(_countsIndex.begin(),_countsIndex.end(),_countsIndex.begin(),
                    [offSet](li x) -> li {return offSet + x; });
                */

                binLocations.insert(binLocations.end(),_binLocations.begin(),_binLocations.end());
                counts.insert(counts.end(),_counts.begin(),_counts.end());
                countsIndex.insert(countsIndex.end(),_countsIndex.begin(),_countsIndex.end());
                exposures.insert(exposures.end(),_exposures.begin(),_exposures.end());
                startTimes.insert(startTimes.end(),_startTimes.begin(),_startTimes.end());
                finishTimes.insert(finishTimes.end(),_finishTimes.begin(),_finishTimes.end());

            }


            cout<<"\nWirting out file: "<<fileNameOut<<endl;

            FileNetcdf fileSmbOut(fileNameOut,NC_NETCDF4);


            if (centroidMode == true)
            {
                fileSmbOut.write_VecNC("centroidLocations",binLocations,NC_DOUBLE);
            }
            else
            {
                fileSmbOut.write_VecNC("binLocations",binLocations,NC_DOUBLE);
            }
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

void createSmbTiles(int sections, const boost::filesystem::path &fileIn)
{
    create_directories(fileIn.stem());
    boost::filesystem::path pathOut;
    pathOut = fileIn.stem() / fileIn.stem();

    cout<<"Output file path: "<<pathOut.string()<<endl;

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
    auto dN = size_t(ceil(double(N)/double(sections)));
    //size_t dN = (N/tiles)+1;

    int tiles = 2*sections-1;

    vector<double> binLocations;
    vector<fp> counts;
    vector<li> countsIndex;
    vector<fp>  exposures;
    vector<double> startTimes;
    vector<double> finishTimes;

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
        startTimes.clear();
        finishTimes.clear();

        size_t beginIdx;
        size_t endIdx;
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

        li indexOffset=countsIndex.front();
        for (li i = 0; i < countsIndex.size(); ++i)
            countsIndex[i]=countsIndex[i]-indexOffset;

        fileSmbOut.write_VecNC("binLocations",binLocations,NC_DOUBLE);
        fileSmbOut.write_VecNC("counts",counts,NC_FLOAT);
        fileSmbOut.write_VecNC("countsIndex",countsIndex,NC_INT64);
        fileSmbOut.write_VecNC("exposures",exposures,NC_FLOAT);
        fileSmbOut.write_VecNC("finishTimes",finishTimes,NC_DOUBLE);
        fileSmbOut.write_VecNC("startTimes",startTimes,NC_DOUBLE);

        fileSmbOut.close();
    }
}
