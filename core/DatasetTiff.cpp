//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
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


#include "DatasetTiff.hpp"
#include <kernel.hpp>
#include <iomanip>
#include <boost/filesystem/convenience.hpp>

using namespace kernel;


DatasetTiff::DatasetTiff(const std::string& filePathIn, const std::string& filePathStemOut, Dataset::WriteType writeType) : fileIn_(0), fileOut_(0), finished_(false)
{
    TIFFSetWarningHandler(NULL);

    if (!filePathIn.empty())
    {
        cout << "Reading " << filePathIn << "..." << endl;
        fileIn_ = TIFFOpen(filePathIn.c_str(), "r");
        if (!fileIn_) throw runtime_error("");
    }

    if (filePathStemOut.empty())
    {
        string fileNameOut = boost::filesystem::path(filePathIn).stem().replace_extension("sml").string();
        fileOut_ = new FileNetcdf(filePathStemOut, NC_NETCDF4);
    }
    else
        fileOut_ = new FileNetcdf(filePathStemOut + (writeType == Dataset::WriteType::InputOutput ? ".sml" : ".smv"), NC_NETCDF4);
}


DatasetTiff::~DatasetTiff()
{
    if (fileIn_)
        TIFFClose(fileIn_);

    if (fileOut_)
        delete fileOut_;
}


bool DatasetTiff::read(Seamass::Input &input, std::string &id)
{
    input = Seamass::Input();

    if(finished_ == true)
        return false;

    uint32 width, height;
    int16 bps, spp;
    if (TIFFGetField(fileIn_, TIFFTAG_IMAGEWIDTH, &width) != 1) throw runtime_error("");
    if (TIFFGetField(fileIn_, TIFFTAG_IMAGELENGTH, &height) != 1) throw runtime_error("");
    if (TIFFGetField(fileIn_, TIFFTAG_BITSPERSAMPLE, &bps) != 1) throw runtime_error("");
    if (TIFFGetField(fileIn_, TIFFTAG_SAMPLESPERPIXEL, &spp) != 1) throw runtime_error("");
    if (bps != 32 || spp != 1) throw runtime_error("ERROR: Input is not a TIFF image in the expected format (single channel 32bit floating point).");

    ii m = height;
    //ii n = width;

    li scanPos = 0;
    input.countsIndex.push_back(scanPos);

    // Make space for image in memory
    float* scanSingle = (float*) _TIFFmalloc(TIFFScanlineSize(fileIn_));
    ii* pixelIdx = new ii[width];

    for (uint32 i = 0; i < m; ++i) {
        // Read image data allocating space for each line as we get it
        TIFFReadScanline(fileIn_,scanSingle,i);

        // make sparse
        ii scan_n = 0;
        ii jOut = 0;
        for(ii jIn = 0; jIn < width; jIn++)
        {
            if(scanSingle[jIn] > 0.0f)
            {
                scanSingle[jOut] = scanSingle[jIn];
                pixelIdx[jOut] = jIn;

                scan_n++;
                jOut++;
            }
        }

        input.countsIndex.push_back(scanPos);
        input.counts.insert(end(input.counts),scanSingle,scanSingle+scan_n);
        input.locations.insert(end(input.locations),pixelIdx,pixelIdx+scan_n);

        scanPos += scan_n;

        if ((i+1) % 100 == 0 || i+1 == m) cout << i+1 << "/" << m << " scanlines processed" << endl;
    }

    _TIFFfree(scanSingle);

    return finished_ = true;
}


void DatasetTiff::write(const Seamass::Input &input, const std::string &id)
{

    ii n = width;
    int matrixId = fileOut_->createGroup("xScale=0");
    fileOut_->writeAttribute(n, "n", "", matrixId);

    vector<long int> pixelIdx(input.locations.begin(),input.locations.end());

    int ijsId = fileOut_->write_VecNC("ijs",input.countsIndex, NC_LONG, matrixId, true);
    int jsId = fileOut_->write_VecNC("js", pixelIdx, NC_LONG, matrixId, true);
    int vsId = fileOut_->write_VecNC("vs", input.counts, NC_FLOAT, matrixId, true);

}


bool DatasetTiff::read(Seamass::Input &input, Seamass::Output &output, std::string &id)
{
    output = Seamass::Output();

    if (!read(input, id))
        return false;

    id = "";

    return finished_ = true;
}


void DatasetTiff::write(const Seamass::Input &input, const Seamass::Output &output, const std::string &id)
{
    write(input, id);

    // Write new Tiff data file after seamass processed...

}




