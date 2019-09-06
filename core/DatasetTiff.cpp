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
        string filePathOut = filePathStemOut + ".seamass.tiff";
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

    uint32 width, height;
    int16 bps, spp;
    if (TIFFGetField(fileIn_, TIFFTAG_IMAGEWIDTH, &width) != 1) throw runtime_error("");
    if (TIFFGetField(fileIn_, TIFFTAG_IMAGELENGTH, &height) != 1) throw runtime_error("");
    if (TIFFGetField(fileIn_, TIFFTAG_BITSPERSAMPLE, &bps) != 1) throw runtime_error("");
    if (TIFFGetField(fileIn_, TIFFTAG_SAMPLESPERPIXEL, &spp) != 1) throw runtime_error("");
    if (bps != 32 || spp != 1) throw runtime_error("ERROR: Input is not a TIFF image in the expected format (single channel 32bit floating point).");

    // output SML
    cout << "Writing " << filePathSml << "..." << endl;
    FileNetcdf sml(filePathSml, NC_NETCDF4);

    ii m = height;
    ii n = width;
    li zero;
    int matrixId = sml.createGroup("xScale=0");
    sml.writeAttribute(n, "n", "", matrixId);
    int ijsId = sml.write_VecNC("ijs", &zero, 1, NC_LONG, matrixId, true);
    int jsId = sml.write_VecNC("js", vector<long int>(), NC_LONG, matrixId, true);
    int vsId = sml.write_VecNC("vs", vector<float>(), NC_FLOAT, matrixId, true);
    li scanPos = 0;
    sml.update_VecNC(ijsId, 0, &scanPos, 1, matrixId);

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

        sml.update_VecNC(jsId, scanPos, pixelIdx, scan_n, matrixId);
        sml.update_VecNC(vsId, scanPos, scanSingle, scan_n, matrixId);

        scanPos += scan_n;

        sml.update_VecNC(ijsId, i+1, &scanPos, 1, matrixId);

        if ((i+1) % 100 == 0 || i+1 == m) cout << i+1 << "/" << m << " scanlines processed" << endl;
    }

    _TIFFfree(scanSingle);

    cout << "Completed" << endl;

    delete[] pixelIdx;
    cout << endl;

    return finished_ = true;
}


void DatasetTiff::write(const std::string& filePathSml, const std::string &id)
{
    FileNetcdf fileSml(filePathSml_,NC_NOWRITE);
    int samplePerPixel = 1;
    int bitsPerSample = 32;
    int bytesPerSample = 4;
    size_t level=1; // how deep to scan sml file for groups.

    MatrixSparse dataCsr;
    ii m;
    ii n;

    ii* ijs;
    ii* js;
    fp* vs;

    // Search for xScale data
    string grpId = fileSml.searchGroup<string>(level);

    fileSml.readMatrixSparseCsr(dataCsr, grpId);

    cout<<"What Scale: "<<grpId<<endl;

    m = dataCsr.m();
    n = dataCsr.n();
    ijs = dataCsr.ijs();
    js = dataCsr.js();
    vs = dataCsr.vs();

    vector<ii> idx;

    /*
    idx.resize(m);
    for (int i = 0; i < m; ++i)
        idx[i]=ijs[i];

    for (int i = 0; i < idx.size(); ++i)
        cout<<"idx["<<i<<"]: "<<idx[i]<<endl;
    */

    vector<float> image(n*m*samplePerPixel,0.0);

    char *metaData;
    if (TIFFGetField(fileIn_, TIFFTAG_IMAGEDESCRIPTION, &metaData) != 1) throw runtime_error("");
    //cout<<"Lets See if this works!!!"<<endl;
    //cout<<metaData<<endl;


    if (TIFFSetField(fileOut_, TIFFTAG_IMAGEWIDTH, n) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_IMAGELENGTH, m) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_XRESOLUTION, 1.0) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_YRESOLUTION, 1.0) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_BITSPERSAMPLE, bitsPerSample) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_SAMPLESPERPIXEL, samplePerPixel) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_ROWSPERSTRIP, 1) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_COMPRESSION, COMPRESSION_NONE) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_IMAGEDESCRIPTION, metaData) != 1) throw runtime_error("Error writing MetaData:");  // Put in original metaData.

    tsize_t lineBytes = samplePerPixel*n*bytesPerSample;     // length in memory of one row of pixel in the image.
    float *buf = NULL;        // buffer used to store the row of pixel information for writing to file

    // Allocating memory to store the pixels of current row
    if (TIFFScanlineSize(fileOut_) != lineBytes)
        buf = (float*)_TIFFmalloc(lineBytes);
    else
        buf = (float*)_TIFFmalloc(TIFFScanlineSize(fileOut_));

    // recreate image from CSR format...
    for (int i = 0; i < m-1; ++i)
    {
        ii idxBeg = ijs[i];
        ii idxEnd = ijs[i+1];
        for (ii j = idxBeg; j < idxEnd; ++j)
        {
            image[i*n + js[j]] = vs[j];
        }
    }

    // Now writing image to the file one strip at a time
    for (int i = 0; i < m; i++)
    {
        memcpy(buf, &image[i*n], n * sizeof(float));
        TIFFWriteScanline(fileOut_, buf, i, 0);
    }

    // Finally we destroy the buffer
    _TIFFfree(buf);
}


bool DatasetTiff::read(std::string& filePathSml, Seamass::Output &output, std::string &id)
{
    throw runtime_error("BUG: not yet implemented!");
}


void DatasetTiff::write(const std::string& filePathSml, const Seamass &seamass, const std::string &id)
{
    Seamass::ControlPoints controlPoints;
    seamass.getOutputControlPoints(controlPoints);
    
    cout << "going to write seaMass output with extent " << controlPoints.extent[0] << "x" << controlPoints.extent[1] << endl;

    ii m = controlPoints.extent[1];
    ii n = controlPoints.extent[0];

    int samplePerPixel = 1;
    int bitsPerSample = 32;
    int bytesPerSample = 4;

    // Read in metadata from the original TIFF file.
    char *metaData;
    if (TIFFGetField(fileIn_, TIFFTAG_IMAGEDESCRIPTION, &metaData) != 1) throw runtime_error("");
    //cout<<"Lets See if this works!!!"<<endl;
    //cout<<metaData<<endl;


    if (TIFFSetField(fileOut_, TIFFTAG_IMAGEWIDTH, n) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_IMAGELENGTH, m) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_XRESOLUTION, 1.0) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_YRESOLUTION, 1.0) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_BITSPERSAMPLE, bitsPerSample) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_SAMPLESPERPIXEL, samplePerPixel) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_ROWSPERSTRIP, 1) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_COMPRESSION, COMPRESSION_NONE) != 1) throw runtime_error("");
    if (TIFFSetField(fileOut_, TIFFTAG_IMAGEDESCRIPTION, metaData) != 1) throw runtime_error("Error writing MetaData:");  // Put in original metaData.

    tsize_t lineBytes = samplePerPixel*n*bytesPerSample;     // length in memory of one row of pixel in the image.
    float *buf = NULL;        // buffer used to store the row of pixel information for writing to file

    // Allocating memory to store the pixels of current row
    if (TIFFScanlineSize(fileOut_) != lineBytes)
        buf = (float*)_TIFFmalloc(lineBytes);
    else
        buf = (float*)_TIFFmalloc(TIFFScanlineSize(fileOut_));

    // Now writing image to the file one strip at a time
    for (int i = 0; i < m; i++)
    {
        memcpy(buf, &controlPoints.coeffs[i*n], n * sizeof(float));
        TIFFWriteScanline(fileOut_, buf, i, 0);
    }

    // Finally destroy the buffer
    _TIFFfree(buf);
}




