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
        fileIn_ = TIFFOpen(filePathIn.c_str(), "r");
        if (!fileIn_) throw runtime_error("");
    }

    if (!filePathStemOut.empty())
    {
        filePathSml_ = filePathStemOut + ".sml";
        string filePathOut = filePathIn + ".seamass.tiff";
        fileOut_ = TIFFOpen(filePathOut.c_str(), "w");
        if (!fileOut_) throw runtime_error("");
    }
}


DatasetTiff::~DatasetTiff()
{
    if (fileIn_)
        TIFFClose(fileIn_);

    if (fileOut_)
        TIFFClose(fileOut_);
}


bool DatasetTiff::read(std::string& filePathSml, std::string &id)
{
    if(finished_ == true)
        return false;
    
    filePathSml = filePathSml_;

    // Ranjeet todo: write into SML file rather than into Seamass::Input
    /*input = Seamass::Input();

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

    _TIFFfree(scanSingle);*/

    return finished_ = true;
}


void DatasetTiff::write(const std::string& filePathSml, const std::string &id)
{
    // Ranjeet todo: read from SML file rather than from Seamass::Input, writing to TIFF in fileOut_
    /*
    ii n = width;
    int matrixId = fileOut_->createGroup("xScale=0");
    fileOut_->writeAttribute(n, "n", "", matrixId);

    vector<long int> pixelIdx(input.locations.begin(),input.locations.end());

    int ijsId = fileOut_->write_VecNC("ijs",input.countsIndex, NC_LONG, matrixId, true);
    int jsId = fileOut_->write_VecNC("js", pixelIdx, NC_LONG, matrixId, true);
    int vsId = fileOut_->write_VecNC("vs", input.counts, NC_FLOAT, matrixId, true);*/
}


bool DatasetTiff::read(std::string& filePathSml, Seamass::Output &output, std::string &id)
{
    output = Seamass::Output();

    if (!read(filePathSml, id))
        return false;

    id = "";

    return finished_ = true;
}


void DatasetTiff::write(const std::string& filePathSml, const Seamass::Output &output, const std::string &id)
{
    // Ranjeet todo: write seamass output to TIFF in fileOut_
 }




