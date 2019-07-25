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

#include <vector>
#include <iostream>
#include <string>
#include <kernel.hpp>
#include "../io/FileNetcdf.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem/convenience.hpp>
//#include <boost/gil/gil_all.hpp>// This seems to be boost version dependent, v1.63
//#include <boost/gil.hpp>// This seems to be boost version dependent, v1.68
//#include <boost/gil/extension/io/jpeg.hpp>
//#include <boost/gil/extension/io/tiff.hpp>

#include <boost/mpl/vector.hpp>
//#include <boost/gil/extension/dynamic_image/any_image.hpp>
//#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <boost/gil/extension/dynamic_image/dynamic_image_all.hpp>

using namespace std;
namespace po = boost::program_options;
//using namespace boost::gil;


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


        // Load TIFF image using Boost GIL
        // any_image<> imgTiff;

        using namespace boost::gil;

        typedef boost::mpl::vector<gray8_image_t, rgb8_image_t, gray16_image_t, rgb16_image_t> my_images_t;

        any_image<my_images_t> dynamic_img;
        jpeg_read_image("test.jpg",dynamic_img);

        // Save the image upside down, preserving its native color space and channel depth
        jpeg_write_view("out-dynamic_image.jpg",flipped_up_down_view(const_view(dynamic_img)));




        // output SML
        cout<<"File name in: "<< fileName;
        string fileOut = boost::filesystem::path(fileName).stem().replace_extension("sml").string();
        FileNetcdf sml(fileOut, NC_NETCDF4);

        li zero;
        ostringstream oss; oss << "mzScale=" << reScale;
        int imageGroup_id = sml.createGroup(oss.str());
        int imageGroup_i_id = sml.write_VecNC("spectrumIndex", &zero, 1, NC_LONG, imageGroup_id, true);
        int imageGroup_j_id = sml.write_VecNC("binLocations", vector<long int>(), NC_LONG, imageGroup_id, true);
        int imageGroup_v_id = sml.write_VecNC("binIntensities", vector<float>(), NC_FLOAT, imageGroup_id, true);



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