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


#include "MatrixSparse.hpp"
#include <iomanip>
#include <sstream>
#if defined(_OPENMP)
  #include <omp.h>
#endif
#include <ippcore.h>
using namespace std;


namespace kernel {


static int debugLevel_ = 0;


void initKernel(int debugLevel)
{
    debugLevel_ = debugLevel;

    // Init IPP library
    ippInit();

    if (debugLevel % 10 >= 2)
    {
        // Get MKL library version info
        MKLVersion Version;
        mkl_get_version(&Version);
        cout << " Intel MKL " << Version.MajorVersion << "." << Version.MinorVersion << "." << Version.UpdateVersion;
        cout << " (" << Version.ProductStatus << " " << Version.Build << ")" << endl;
        cout << "  Platform               : " << Version.Platform << endl;
        cout << "  Processor optimization : " << Version.Processor << endl;

        // Get IPP library version info
        const IppLibraryVersion *lib;
        IppStatus status;
        Ipp64u mask, emask;
        lib = ippGetLibVersion();
        cout << " Intel " << lib->Name << " " << lib->Version << endl;

        // Get CPU features and features enabled with selected library level
        status = ippGetCpuFeatures(&mask, 0);
        if (ippStsNoErr == status)
        {
            emask = ippGetEnabledCpuFeatures();
            printf("  Features supported  by CPU by IPP\n");
            printf("   ippCPUID_MMX          : ");
            printf("%c\t%c\t", (mask & ippCPUID_MMX) ? 'Y' : 'N', (emask & ippCPUID_MMX) ? 'Y' : 'N');
            printf("Intel(R) Architecture MMX technology supported\n");
            printf("   ippCPUID_SSE          : ");
            printf("%c\t%c\t", (mask & ippCPUID_SSE) ? 'Y' : 'N', (emask & ippCPUID_SSE) ? 'Y' : 'N');
            printf("Intel(R) Streaming SIMD Extensions\n");
            printf("   ippCPUID_SSE2         : ");
            printf("%c\t%c\t", (mask & ippCPUID_SSE2) ? 'Y' : 'N', (emask & ippCPUID_SSE2) ? 'Y' : 'N');
            printf("Intel(R) Streaming SIMD Extensions 2\n");
            printf("   ippCPUID_SSE3         : ");
            printf("%c\t%c\t", (mask & ippCPUID_SSE3) ? 'Y' : 'N', (emask & ippCPUID_SSE3) ? 'Y' : 'N');
            printf("Intel(R) Streaming SIMD Extensions 3\n");
            printf("   ippCPUID_SSSE3        : ");
            printf("%c\t%c\t", (mask & ippCPUID_SSSE3) ? 'Y' : 'N', (emask & ippCPUID_SSSE3) ? 'Y' : 'N');
            printf("Intel(R) Supplemental Streaming SIMD Extensions 3\n");
            printf("   ippCPUID_MOVBE        : ");
            printf("%c\t%c\t", (mask & ippCPUID_MOVBE) ? 'Y' : 'N', (emask & ippCPUID_MOVBE) ? 'Y' : 'N');
            printf("The processor supports MOVBE instruction\n");
            printf("   ippCPUID_SSE41        : ");
            printf("%c\t%c\t", (mask & ippCPUID_SSE41) ? 'Y' : 'N', (emask & ippCPUID_SSE41) ? 'Y' : 'N');
            printf("Intel(R) Streaming SIMD Extensions 4.1\n");
            printf("   ippCPUID_SSE42        : ");
            printf("%c\t%c\t", (mask & ippCPUID_SSE42) ? 'Y' : 'N', (emask & ippCPUID_SSE42) ? 'Y' : 'N');
            printf("Intel(R) Streaming SIMD Extensions 4.2\n");
            printf("   ippCPUID_AVX          : ");
            printf("%c\t%c\t", (mask & ippCPUID_AVX) ? 'Y' : 'N', (emask & ippCPUID_AVX) ? 'Y' : 'N');
            printf("Intel(R) Advanced Vector Extensions instruction set\n");
            printf("   ippAVX_ENABLEDBYOS    : ");
            printf("%c\t%c\t", (mask & ippAVX_ENABLEDBYOS) ? 'Y' : 'N', (emask & ippAVX_ENABLEDBYOS) ? 'Y' : 'N');
            printf("The operating system supports Intel(R) AVX\n");
            printf("   ippCPUID_AES          : ");
            printf("%c\t%c\t", (mask & ippCPUID_AES) ? 'Y' : 'N', (emask & ippCPUID_AES) ? 'Y' : 'N');
            printf("Intel(R) AES instruction\n");
            printf("   ippCPUID_SHA          : ");
            printf("%c\t%c\t", (mask & ippCPUID_SHA) ? 'Y' : 'N', (emask & ippCPUID_SHA) ? 'Y' : 'N');
            printf("Intel(R) SHA new instructions\n");
            printf("   ippCPUID_CLMUL        : ");
            printf("%c\t%c\t", (mask & ippCPUID_CLMUL) ? 'Y' : 'N', (emask & ippCPUID_CLMUL) ? 'Y' : 'N');
            printf("PCLMULQDQ instruction\n");
            printf("   ippCPUID_RDRAND       : ");
            printf("%c\t%c\t", (mask & ippCPUID_RDRAND) ? 'Y' : 'N', (emask & ippCPUID_RDRAND) ? 'Y' : 'N');
            printf("Read Random Number instructions\n");
            printf("   ippCPUID_F16C         : ");
            printf("%c\t%c\t", (mask & ippCPUID_F16C) ? 'Y' : 'N', (emask & ippCPUID_F16C) ? 'Y' : 'N');
            printf("Float16 instructions\n");
            printf("   ippCPUID_AVX2         : ");
            printf("%c\t%c\t", (mask & ippCPUID_AVX2) ? 'Y' : 'N', (emask & ippCPUID_AVX2) ? 'Y' : 'N');
            printf("Intel(R) Advanced Vector Extensions 2 instruction set\n");
            printf("   ippCPUID_AVX512F      : ");
            printf("%c\t%c\t", (mask & ippCPUID_AVX512F) ? 'Y' : 'N', (emask & ippCPUID_AVX512F) ? 'Y' : 'N');
            printf("Intel(R) Advanced Vector Extensions 3.1 instruction set\n");
            printf("   ippCPUID_AVX512CD     : ");
            printf("%c\t%c\t", (mask & ippCPUID_AVX512CD) ? 'Y' : 'N', (emask & ippCPUID_AVX512CD) ? 'Y' : 'N');
            printf("Intel(R) Advanced Vector Extensions CD (Conflict Detection) instruction set\n");
            printf("   ippCPUID_AVX512ER     : ");
            printf("%c\t%c\t", (mask & ippCPUID_AVX512ER) ? 'Y' : 'N', (emask & ippCPUID_AVX512ER) ? 'Y' : 'N');
            printf("Intel(R) Advanced Vector Extensions ER instruction set\n");
            printf("   ippCPUID_ADCOX        : ");
            printf("%c\t%c\t", (mask & ippCPUID_ADCOX) ? 'Y' : 'N', (emask & ippCPUID_ADCOX) ? 'Y' : 'N');
            printf("ADCX and ADOX instructions\n");
            printf("   ippCPUID_RDSEED       : ");
            printf("%c\t%c\t", (mask & ippCPUID_RDSEED) ? 'Y' : 'N', (emask & ippCPUID_RDSEED) ? 'Y' : 'N');
            printf("The RDSEED instruction\n");
            printf("   ippCPUID_PREFETCHW    : ");
            printf("%c\t%c\t", (mask & ippCPUID_PREFETCHW) ? 'Y' : 'N', (emask & ippCPUID_PREFETCHW) ? 'Y' : 'N');
            printf("The PREFETCHW instruction\n");
            printf("   ippCPUID_KNC          : ");
            printf("%c\t%c\t", (mask & ippCPUID_KNC) ? 'Y' : 'N', (emask & ippCPUID_KNC) ? 'Y' : 'N');
            printf("Intel(R) Xeon Phi(TM) Coprocessor instruction set\n");
        }
    }

    // Thread Info
    if (debugLevel % 10 >= 1)
    {
        cout << " Config: " << 8 * sizeof(ii) << "bit MKL addressing, " << mkl_get_max_threads() << " MKL threads, ";
#if defined(_OPENMP)
        cout << omp_get_max_threads() << " OpenMP threads";
#else
        cout << "non-OpenMP build";
#endif
        cout << endl;
    }
}


int getDebugLevel()
{
    return debugLevel_;
}


static li id_ = -1;


li getId()
{
    return id_;
}


static double startTime_ = dsecnd();


void resetElapsedTime()
{
    startTime_ = dsecnd();
}


double getElapsedTime()
{
    return dsecnd() - startTime_;
}


li getUsedMemory()
{
    int allocatedBuffers;
    return mkl_mem_stat(&allocatedBuffers);
}


string getTimeStamp()
{
    ostringstream out;
    out << "[" << setw(9) << ++id_ << "," << fixed << internal << setw(9) << std::setprecision(3) << getElapsedTime() << "," << setw(9) << getUsedMemory()/1024.0/1024.0 << "] ";
    return out.str();
}


}