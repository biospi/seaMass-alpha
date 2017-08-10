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


#ifndef SEAMASS_CORE_SEAMASS_HPP
#define SEAMASS_CORE_SEAMASS_HPP


#include "BasisBspline.hpp"
#include "../asrl/OptimizerSrl.hpp"


/**
* SeamassCore fitting of a 1-d curve or multi-dimensional surface to the input spectr(um|a).
*/
class Seamass : public Subject
{
public:
    static void notice();

    struct Input {
        enum class Type { Binned, Sampled, Centroided } type;
        std::vector<fp> counts;
        std::vector<li> countsIndex;
        std::vector<double> locations;
        std::vector<double> startTimes;
        std::vector<double> finishTimes;
        std::vector<fp> exposures;
    };

    struct Output
    {
        std::vector<char> scale; // scale of finest basis functions, vector of size dimensions (i.e. 1 or 2)
        double shrinkage;        // shrinkage used
        double tolerance;        // tolerance used
        double peakFwhm;
        short chargeStates;

        std::vector<BasisBspline::GridInfo> gridInfos;
        std::vector<MatrixSparse> xs;
        std::vector<MatrixSparse> l2s;
        std::vector<MatrixSparse> l1l2s;

        /*std::vector<char> baselineScale; // scale of finest basis functions, vector of size dimensions (i.e. 1 or 2)
        std::vector<ii> baselineOffset;  // offset of finest basis functions
        std::vector<ii> baselineExtent;  // extent of finest basis functions
        double shrinkage;                // shrinkage used
        double tolerance;                // tolerance used
        std::vector< std::vector<char> > scales; // scales of each basis function for each dimension
        std::vector< std::vector<ii> > offsets; // offsets of each basis functions for each dimension
        std::vector<fp> weights;         // weight of each basis functions (i.e. xs)*/
    };

    struct ControlPoints {
        std::vector<fp> coeffs;
        std::vector<char> scale;
        std::vector<ii> offset;
        std::vector<ii> extent;
    };

    Seamass(Input& input, const std::vector<char>& scale, fp lambda, bool taperShrinkage, fp tolerance,
            double peakFwhm, short chargeStates);
    Seamass(Input& input, const Output& seed);
    virtual ~Seamass();

    bool step();
    ii getIteration() const;

    // get seaMass convolved input or restored output
    void getInput(Input &input, bool reconstruct = false) const;

    // get seaMass output (for smv file)
    void getOutput(Output& output, bool synthesize) const;

    // get restored 1D control points (i.e. per spectra) derived from seaMass output
    void getOutputControlPoints1d(ControlPoints& controlPoints, bool deconvolve, bool density) const;

    // get restored control points with dimension depending on input (i.e. 1D or 2D)
    void getOutputControlPoints(ControlPoints& controlPoints, bool deconvolve) const;

private:
    void init(Input& input, const std::vector<char>& scales, short chargeStates, bool seed);

    char dimensions_;
    std::vector<Basis*> bases_;
    std::vector<MatrixSparse> b_;

    Optimizer* innerOptimizer_;
    Optimizer* optimizer_;

    fp lambda_;
    fp lambdaStart_;
    bool taperShrinkage_;
    fp tolerance_;
    int iteration_;
    double peakFwhm_;
    short chargeStates_;
};


#endif

