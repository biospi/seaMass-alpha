module Seamass

using HDF5
using LibExpat

type MzmlSpectrum
  mzs
  intensities
end

function MzmlSpectrum(filename::AbstractString, spectrumID::Number)
    h5open(filename, "r") do file
        spectrumIndex = file["mzML_spectrumIndex"][spectrumID:spectrumID+1] + 1

        mzML = xp_parse(String(file["mzML"][spectrumIndex[1]:spectrumIndex[2]-1]))

        arrayLength = parse(Int, mzML["@defaultArrayLength"][1])
        mzsDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@externalDataset"][1]
        mzsOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@offset"][1]) + 1
        intensitiesDataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@externalDataset"][1]
        intensitiesOffset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@offset"][1]) + 1

        MzmlSpectrum(
            file[mzsDataset][mzsOffset:mzsOffset+arrayLength-1],
            file[intensitiesDataset][intensitiesOffset:intensitiesOffset+arrayLength-1]
        )
    end
end

end