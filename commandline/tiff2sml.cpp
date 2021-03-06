//
// Original author: Ranjeet Bhamber <ranjeet.bhamber <a.t> bristol.ac.uk>
//
// Copyright (C) 2018  biospi Laboratory, University of Bristol, UK
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

#include <string>
#include <iostream>
#include <kernel.hpp>
#include "../io/FileNetcdf.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
#include <tiffio.h>

using namespace std;
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{

#ifdef NDEBUG
    try
#endif
    {
        string filePathIn;

        po::options_description general(
                "Usage\n"
                "-----\n"
                "tiff2sml [OPTIONS...] [tiff FILE]\n"
                "tiff2sml <file>"
        );

        general.add_options()
            ("help,h",
             "Produce this help message")
            ("file,f", po::value<string>(&filePathIn),
             "Input TIFF file.")
         ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "tiff2sml : Copyright (C) 2019 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;

        if(vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

 
        // Open Tiff file for reading and read tags
        TIFFSetWarningHandler(NULL);
        cout << "Reading " << filePathIn << "..." << endl;
        TIFF* imgTiff = TIFFOpen(filePathIn.c_str(), "r");
        if (!imgTiff) throw runtime_error("");
        
        uint32 width, height;
        int16 bps, spp;
        if (TIFFGetField(imgTiff, TIFFTAG_IMAGEWIDTH, &width) != 1) throw runtime_error("");
        if (TIFFGetField(imgTiff, TIFFTAG_IMAGELENGTH, &height) != 1) throw runtime_error("");
        if (TIFFGetField(imgTiff, TIFFTAG_BITSPERSAMPLE, &bps) != 1) throw runtime_error("");
        if (TIFFGetField(imgTiff, TIFFTAG_SAMPLESPERPIXEL, &spp) != 1) throw runtime_error("");
        if (bps != 32 || spp != 1) throw runtime_error("ERROR: Input is not a TIFF image in the expected format (single channel 32bit floating point).");
 
        // output SML
        string fileNameOut = boost::filesystem::path(filePathIn).stem().replace_extension("sml").string();
        cout << "Writing " << fileNameOut << "..." << endl;
        FileNetcdf sml(fileNameOut, NC_NETCDF4);
        
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
        float* scanSingle = (float*) _TIFFmalloc(TIFFScanlineSize(imgTiff));
        ii* pixelIdx = new ii[width];

        for (uint32 i = 0; i < m; ++i) {
            // Read image data allocating space for each line as we get it
            TIFFReadScanline(imgTiff,scanSingle,i);

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
        TIFFClose(imgTiff);
        
        cout << "Completed" << endl;

        delete[] pixelIdx;
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
