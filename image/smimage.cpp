//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Bristol, UK
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

#include "imagecore.hpp"
#include <boost/program_options.hpp>

#include "../io/FileNetcdf.hpp"
#include "../io/VecMat.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv){

    string smbFileName, channel;
    //int mz_res=0;
    pair<double,double> mzInputRange;
    pair<double,double> rtInputRange;
    pair<ui,ui> xyview;
    bool rescaleExp = false;

    po::options_description usage("Usage: smimage [OPTION...] [SMB FILE]\n"
            "Example: image -m [200,300] -r [20,140] -o [3000,1000] --ms2 datafile.smb\n"
            "         image -o [,500] -m[,300] -r [12,] -f datafile.smj\n\n"
            "Options");

    usage.add_options()
        ("help,h", "Produce help message")
        ("file,f",
            po::value<string>(&smbFileName),
                "Filename - Load data from seaMass input file format (smj).")
        ("output-image,o",
            po::value<string>(),
                "Output range of generated SMG file (default values [1280,1024]). "
                "The argument takes the following format: [MZ Pixels,RT Pixels], "
                "omitting a value will result in the appropriate default values "
                ".i.e. [,Val] -> will give a range [1280,Val], "
                "[Val,] -> will give a range [Val,1024].")
        ("mz-range,m",
            po::value<string>(),
                "MZ range units (m/z Th) argument takes the following format: [Min Value,Max Value], "
                "omitting a value will result in the appropriate [Min,Max] value "
                ".i.e. [,Val] -> will give a range [Min,Val], "
                "[Val,] -> will give a range [Val,Max].")
        ("rt-range,r",
            po::value<string>(),
                "RT range units (s) argument takes the following format: [Min Value,Max Value], "
                "omitting a value will result in the appropriate [Min,Max] value "
                ".i.e. [,Val] -> will give a range [Min,Val], "
                "[Val,] -> will give a range [Val,Max].")
            ("rescale-exposure,e",
             "Rescale exposure gain control, default applyed, use this switch to rescale to origonal.");
    try
    {
        po::positional_options_description pod;
        pod.add("file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(usage).positional(pod).run(), vm);
        po::notify(vm);

        if(vm.count("help"))
        {
            cout<<usage<<endl;
            return 0;
        }
        if(vm.count("output-image"))
        {
            xyview=numStrBraket<ui,ui>(vm["output-image"].as<string>());
            if(xyview.first == -1) xyview.first = 1280;
            if(xyview.second == -1) xyview.second = 1024;
        }
        else
        {
            xyview = make_pair(1280,1024);
        }
        if(vm.count("mz-range"))
        {
            mzInputRange=numStrBraket<double,double>(vm["mz-range"].as<string>());
        }
        else
        {
            mzInputRange = make_pair(-1.0,-1.0);
        }
        if(vm.count("rt-range"))
        {
            rtInputRange=numStrBraket<double,double>(vm["rt-range"].as<string>());
            if(rtInputRange.first != -1) rtInputRange.first = rtInputRange.first;
            if(rtInputRange.second != -1) rtInputRange.second = rtInputRange.second;
        }
        else
        {
            rtInputRange = make_pair(-1.0,-1.0);
        }
        if(vm.count("rescale-exposure"))
        {
            rescaleExp=true;
        }
        if(vm.count("file"))
        {
            cout<<"Opening SMB file: "<<vm["file"].as<string>()<<endl;
        }
        else
        {
            cout<<usage<<endl;
            throw "Input SMB file was not give...";
        }
    }
    catch(exception& e)
    {
        cerr<<"error: " << e.what() <<endl;
        return 1;
    }
    catch(const char* msg)
    {
        cerr<<"error: "<<msg<<endl;
        return 1;
    }
    catch(...)
        {
        cerr<<"Exception of unknown type!\n";
    }

    string outFileName=smbFileName.substr(0,smbFileName.size()-4);
    FileNetcdf smbInput(smbFileName);

    cout<<"Load raw data from data files:"<<endl;
    MassSpecData raw;

    // Load Data from smb file
    cout<<"Load Start Time from SMB:"<<endl;
    int data2D = smbInput.searchGroup("startTimes");

    if(data2D == 0)
    {
        raw.N = 2;
        raw.rt.push_back(1);
        raw.rt.push_back(2);
        smbInput.readVector(raw.mz, "binLocations", 0);
        smbInput.readVector(raw.sc, "counts", 0);
        smbInput.readVector(raw.exp, "exposures", 0);
        raw.sci.push_back(0);
        raw.sci.push_back(lli(raw.mz.size())-1);
        xyview.second=1;
    }
    else
    {
        vector<double> finishTimes;
        vector<double> rtBuff;
        smbInput.readVector(raw.rt, "startTimes", 0);
        smbInput.readVector(finishTimes, "finishTimes", 0);
        for(size_t i =0; i < finishTimes.size(); ++i)
        {
            rtBuff.push_back(raw.rt[i]);
            rtBuff.push_back(finishTimes[i]);
        }
        raw.rt.clear();
        raw.rt=rtBuff;
        raw.N = raw.rt.size();
        smbInput.readVector(raw.mz, "binLocations", 0);
        smbInput.readVector(raw.sc, "counts", 0);
        smbInput.readVector(raw.sci, "countsIndex", 0);
        smbInput.readVector(raw.exp, "exposures", 0);
    }

    MassSpecImage imgBox(xyview);

    if(rescaleExp)
    {
        cout<<"Running Rescale!!"<<endl;
        for(size_t i=0; i < raw.exp.size(); ++i)
        {
            lli idx_beg=raw.sci[i];
            lli idx_end=raw.sci[i+1];

            raw.exp[i]=1.0/raw.exp[i];

            for(lli j=idx_beg; j <= idx_end; ++j)
            {
                raw.sc[j]=raw.sc[j]*raw.exp[i];
            }
        }
    }

    cout<<"Processing data to create new SMR raw image file..."<<endl;
    // Calculate MZ index
    raw.calMZi();
    // Find data limits in mass spec data.
    raw.calRange();

    raw.rti.resize(raw.N,-1);
    //for(size_t i=0; i < raw.N-1; ++i)
    //	raw.rti[i]=i;
    //for(size_t i=0; i < raw.N/2-1; ++i)
    for(size_t i=0; i < raw.N; i=i+2)
        raw.rti[i]=lli(i/2);

    // Clip Box if needed
    // Min MZ Value
    (mzInputRange.first < raw.mzRange.first || mzInputRange.first < 0 ) ?
            imgBox.mzRange.first=raw.mzRange.first:
            imgBox.mzRange.first=mzInputRange.first;
    // Max MZ Value
    (mzInputRange.second > raw.mzRange.second || mzInputRange.second < 0 ) ?
            imgBox.mzRange.second=raw.mzRange.second:
            imgBox.mzRange.second=mzInputRange.second;
    // Min RT Value
    (rtInputRange.first < raw.rtRange.first || rtInputRange.first < 0 ) ?
            imgBox.rtRange.first=raw.rtRange.first:
            imgBox.rtRange.first=rtInputRange.first;
    // Max RT Value
    (rtInputRange.second > raw.rtRange.second || rtInputRange.second < 0 ) ?
            imgBox.rtRange.second=raw.rtRange.second:
            imgBox.rtRange.second=rtInputRange.second;

    imgBox.dmz=genAxis(imgBox.mz, imgBox.mzRange.first, imgBox.mzRange.second);
    imgBox.drt=genAxis(imgBox.rt, imgBox.rtRange.first, imgBox.rtRange.second);

    // Find RT range index n raw data, as data has irregular DRT we have to scan manually.
    size_t rtidxBegin=0, rtidxEnd=0;
    for(size_t i = 0; i < raw.rt.size(); ++i)
    {
        if(imgBox.rt[0] < raw.rt[i])
        {
            rtidxBegin=i-1;
            break;
        }
    }
    for(size_t i = rtidxBegin; i < raw.rt.size(); ++i)
    {
        if(imgBox.rt[imgBox.xypxl.second-1] < raw.rt[i])
        {
            rtidxEnd=i;
            break;
        }
        else if(i == raw.rt.size()-1)
        {
            rtidxEnd=raw.rt.size()-1;
        }
    }

    if(data2D)
    {
        //for (size_t idx = rtidxBegin; idx <= rtidxEnd; ++idx) {
        for (size_t idx = rtidxBegin; idx <= rtidxEnd; idx=idx+2) {
            double scaleRT = 0.0;
            double drt = fabs(raw.rt[idx + 1] - raw.rt[idx]);
            double yn = (raw.rt[idx] - imgBox.rt[0]) / imgBox.drt;
            double yp = (raw.rt[idx + 1] - imgBox.rt[0]) / imgBox.drt;
            lli rowN = lli(floor(yn));
            lli rowP = lli(floor(yp));

            if (rowP > imgBox.xypxl.second - 2) rowP = imgBox.xypxl.second - 1;
            if (rowN >= imgBox.xypxl.second - 1) break;
            // Overlapping image 1st box on boundary
            if (rowN < 0 && rowP >= 0) {
                rowN = 0;
                for (lli row = rowN; row <= rowP; ++row) {
                    double boxBegin = 0.0;
                    double boxEnd = 0.0;
                    // Cut top overlap else its contained
                    (imgBox.rt[row] < raw.rt[idx]) ? boxBegin = raw.rt[idx]
                                                   : boxBegin = imgBox.rt[row];
                    // Cut bottom overlap else its contained
                    (imgBox.rt[row + 1] < raw.rt[idx + 1]) ? boxEnd = imgBox.rt[row + 1]
                                                           : boxEnd = raw.rt[idx + 1];
                    scaleRT = (boxEnd - boxBegin) / drt;
                    // Run MZ scan calculations
                    if (raw.rti[idx] >= 0)
                        scanMZ(raw, imgBox, raw.rti[idx], row, scaleRT);
                }
            }
                // Multiple smaller image boxes within single raw data scanline
            else if (rowN < rowP && rowN >= 0) {
                for (lli row = rowN; row <= rowP; ++row) {
                    double boxBegin = 0.0;
                    double boxEnd = 0.0;
                    // Cut top overlap else its contained
                    (imgBox.rt[row] < raw.rt[idx]) ? boxBegin = raw.rt[idx]
                                                   : boxBegin = imgBox.rt[row];
                    // Cut bottom overlap else its contained
                    (imgBox.rt[row + 1] < raw.rt[idx + 1]) ? boxEnd = imgBox.rt[row + 1]
                                                           : boxEnd = raw.rt[idx + 1];
                    scaleRT = (boxEnd - boxBegin) / drt;
                    // Run MZ scan calculations
                    if (raw.rti[idx] >= 0)
                        scanMZ(raw, imgBox, raw.rti[idx], row, scaleRT);
                }
            }
                // Raw data scanline contained within a Image box hence no scale needed.
            else if (rowN == rowP) {
                scaleRT = 1.0;
                // Run MZ scan calculations
                if (raw.rti[idx] >= 0) scanMZ(raw, imgBox, raw.rti[idx], rowN, scaleRT);
            }
        }
    }
    else
    {
        size_t idx = 0;
        double scaleRT = 1.0;
        //double drt = fabs(raw.rt[idx + 1] - raw.rt[idx]);
        double yn = (raw.rt[idx] - imgBox.rt[0]) / imgBox.drt;
        //double yp = (raw.rt[idx + 1] - imgBox.rt[0]) / imgBox.drt;
        lli rowN = lli(floor(yn));
        //lli rowP = lli(floor(yp));
        if (raw.rti[idx] >= 0) scanMZ(raw, imgBox, raw.rti[idx], rowN, scaleRT);
    }

    cout<<"Writing out data to seaMass image data file: "<<endl;
    cout<<"    "<<outFileName+".smr"<<endl;

    //Convert MZ and StartTimes to bin centres
    FileNetcdf smrFile(outFileName+".smr",NC_NETCDF4);

    if(data2D)
    {
        centreBin(imgBox.mz);
        centreBin(imgBox.rt);
        VecMat<float> scImage(imgBox.xypxl.second,imgBox.xypxl.first,imgBox.sc);
        smrFile.write_MatAxisNC("SpectrumCount",scImage,NC_FLOAT,
                                imgBox.mz,NC_DOUBLE,
                                imgBox.rt,NC_DOUBLE,
                                "SpectrumMZ","StartTime");
    }
    else
    {
        centreBin(imgBox.mz);
        smrFile.writeVector(imgBox.mz, "SpectrumMZ");
        smrFile.writeVector(imgBox.sc, "SpectrumCount");
    }


    return 0;
}
