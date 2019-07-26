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
        string fileName;
        int reScale;

        po::options_description general(
                "Usage\n"
                "-----\n"
                "image2sml [OPTIONS...] [image FILE]\n"
                "image2sml <file>"
        );

        general.add_options()
                ("help,h",
                 "Produce this help message")
                ("file,f", po::value<string>(&fileName),
                 "Input image file in TIFF format.")
                ("bin_scale,b", po::value<int>(&reScale)->default_value(0),
                 "Output re-BIN resolution, currently = 0.")
                ;

        po::options_description desc;
        desc.add(general);

        po::positional_options_description pod;
        pod.add("file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(general).positional(pod).run(), vm);
        po::notify(vm);

        cout << endl;
        cout << "image2sml : Copyright (C) 2018 - biospi Laboratory, University of Bristol, UK" << endl;
        cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
        cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
        cout << endl;

        if(vm.count("help") || !vm.count("file"))
        {
            cout << desc << endl;
            return 0;
        }

        // libtiff variables
        uint32 width, height;
        short bitPerSample;
        tsize_t scanLength;

        li scanPos = 0;

        cout<<"Loading TIFF file: "<< fileName;

        // Open Tiff file for reading
        //TIFF*imgTiff = TIFFOpen(fileName.c_str(),"r");
        TIFF*imgTiff = TIFFOpen("test.tiff","r");
        if (!imgTiff) {
            cerr << "Failed to open image" << endl;
            exit(1);
        }

        // Read dimensions of image
        if (TIFFGetField(imgTiff,TIFFTAG_IMAGEWIDTH,&width) != 1) {
            cerr << "Failed to read width" << endl;
            exit(1);
        }
        if (TIFFGetField(imgTiff,TIFFTAG_IMAGELENGTH, &height) != 1) {
            cerr << "Failed to read height" << endl;
            exit(1);
        }
        if (TIFFGetField(imgTiff,TIFFTAG_BITSPERSAMPLE, &bitPerSample) != 1) {
            cerr << "Failed to read the number of bits per sample" << endl;
            exit(1);
        }
        scanLength = TIFFScanlineSize(imgTiff);

        cout << "Image dimensions (width x height): " << width << "x" << height << endl;
        cout << "Number of bits per sample: " << bitPerSample << endl;
        cout << "Line buffer length (bytes): " << scanLength << endl;

        // output SML
        string fileOut = boost::filesystem::path(fileName).stem().replace_extension("sml").string();
        FileNetcdf sml(fileOut, NC_NETCDF4);

        li zero;
        ostringstream oss; oss << "mzScale=" << reScale;
        int imageGroup_id = sml.createGroup(oss.str());
        int imageGroup_SpecIdx_id = sml.write_VecNC("spectrumIndex", &zero, 1, NC_LONG, imageGroup_id, true);
        int imageGroup_BinLoc_id = sml.write_VecNC("binLocations", vector<long int>(), NC_LONG, imageGroup_id, true);
        int imageGroup_BinInten_id = sml.write_VecNC("binIntensities", vector<float>(), NC_FLOAT, imageGroup_id, true);
        sml.update_VecNC(imageGroup_SpecIdx_id, 0, &scanPos, 1, imageGroup_id);

        // Make space for image in memory
        //float** singleScan = (float**)malloc(sizeof (float*)*height);
        float* scanSingle = new float[width];
        ii* pixelIdx = new ii[width];

        for (uint32 i = 0; i < height; ++i) {
            // Read image data allocating space for each line as we get it
            cout<<"Processing image scan: "<<height<<"/"<<i<<"\r";
            TIFFReadScanline(imgTiff,scanSingle,i);

            //cout<<"Scanline("<<i<<"): "<<scanSingle[0]<<", "<<scanSingle[1]<<", "<<scanSingle[2]
            //    <<", "<<scanSingle[3]<<", "<<scanSingle[4]<<", ... "<<scanSingle[width-1]<<endl;

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

            sml.update_VecNC(imageGroup_BinLoc_id, scanPos, pixelIdx, scan_n, imageGroup_id);
            sml.update_VecNC(imageGroup_BinInten_id, scanPos, scanSingle, scan_n, imageGroup_id);

            scanPos += scan_n;

            sml.update_VecNC(imageGroup_SpecIdx_id, i+1, &scanPos, 1, imageGroup_id);
        }

        TIFFClose(imgTiff);

        delete[] scanSingle;
        delete[] pixelIdx;
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