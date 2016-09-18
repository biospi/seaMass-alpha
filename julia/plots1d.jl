using PyPlot
using HDF5
using LibExpat
using PyCall
@pyimport scipy.interpolate as interpolate


spectrumID = 1

filename = "/Volumes/C/s/swath/P02U_Swath_1_1D.mzMLb"
# load spectrum from mzMLb file
type mzMLb_spectrum
  mzs
  intensities
end
mzmlb_spectrum = h5open(filename, "r") do file
  spectrumIndex = file["mzML_spectrumIndex"][spectrumID:spectrumID+1] + 1

  mzML = xp_parse(utf8(file["mzML"][spectrumIndex[1]:spectrumIndex[2]-1]))

  arrayLength = parse(Int, mzML["@defaultArrayLength"][1])
  mzs_dataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@externalDataset"][1]
  mzs_offset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000514']/../binary"]["@offset"][1]) + 1
  intensities_dataset = mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@externalDataset"][1]
  intensities_offset = parse(Int, mzML["binaryDataArrayList/binaryDataArray/cvParam[@accession='MS:1000515']/../binary"]["@offset"][1]) + 1

  mzMLb_spectrum(
    file[mzs_dataset][mzs_offset:mzs_offset+arrayLength-1],
    file[intensities_dataset][intensities_offset:intensities_offset+arrayLength-1]
  )
end

figure()
axhline(color="grey")
plot(mzmlb_spectrum.mzs, mzmlb_spectrum.intensities)
xlabel("m/z (Th)")
ylabel("intensity")


# load spectrum from smi file
type SMI_Spectrum
  binEdges
  binCounts
  exposure
end
smi_spectrum = h5open("/Volumes/C/s/swath/P02U_Swath_1_1D.pos_617.55.smi", "r") do file
  spectrumIndex = [1 size(file["binCounts"])[1]+1]
  try
    spectrumIndex = file["spectrumIndex"][spectrumID:spectrumID+1] + 1
  end
  SMI_Spectrum(
    file["binEdges"][spectrumIndex[1]+spectrumID-1:spectrumIndex[2]+spectrumID-1],
    file["binCounts"][spectrumIndex[1]:spectrumIndex[2]-1],
    file["exposures"][spectrumID][1]
  )
end


figure()
plot(mzmlb_spectrum.mzs, mzmlb_spectrum.intensities)
step(smi_spectrum.binEdges, vcat(0.0, smi_spectrum.binCounts / smi_spectrum.exposure))
xlabel("m/z (Th)")
ylabel("intensity")


# load spectrum from smv file
type SMV_Spectrum
  binCounts
end
filename = "/Volumes/C/s/swath/P02U_Swath_1_1D.pos_617.55.smv"
smv_spectrum = h5open(filename, "r") do file
  spectrumIndex = [1 size(file["binCounts"])[1]+1]
  try
    spectrumIndex = file["spectrumIndex"][spectrumID:spectrumID+1] + 1
  end
  SMV_Spectrum(
    file["binCounts"][spectrumIndex[1]:spectrumIndex[2]-1],
  )
end


figure()
axhline(color="grey")
smi_spectrum_bin_widths = smi_spectrum.binEdges[2:length(smi_spectrum.binEdges)] - smi_spectrum.binEdges[1:length(smi_spectrum.binEdges)-1]
step(smi_spectrum.binEdges, vcat(0.0, smi_spectrum.binCounts ./ smi_spectrum_bin_widths ))
step(smi_spectrum.binEdges, vcat(0.0, (smi_spectrum.binCounts - smv_spectrum.binCounts) ./ smi_spectrum_bin_widths ))
xlabel("m/z (Th)")
ylabel("ion count per Th")


# load spectrum from smo file
type SMO_Spectrum
  controlPoints
  offset
  scale
end
filename = "/Volumes/C/s/swath/P02U_Swath_1_1D.pos_617.55.smo"
smo_spectrum = h5open(filename, "r") do file
  SMO_Spectrum(
    read(file, "controlPoints"),
    h5readattr(filename, "controlPoints")["offset"][1],
    h5readattr(filename, "controlPoints")["scale"][1]
  )
end


npoints = 10000
coefs = vcat(smo_spectrum.controlPoints, 0, 0, 0, 0) # no idea why I have to pad this
knots = -3:length(coefs)-4
tck = (knots, coefs, 3)
x = linspace(0, length(smo_spectrum.controlPoints)-3, npoints)
intensities_bspline = interpolate.splev(x, tck, ext=1)

mz_range = (smo_spectrum.offset + [0, length(smo_spectrum.controlPoints) - 3]) * 1.0033548378 / (60 * 2^smo_spectrum.scale)
mzs_bspline = linspace(mz_range[1], mz_range[2], npoints)

figure()
axhline(color="grey")
axvline(x=mz_range[1],color="grey")
axvline(x=mz_range[2],color="grey")
#step(smi_spectrum.binEdges, vcat(0.0, smi_spectrum.binCounts ./ smi_spectrum_bin_widths ))
#step(smi_spectrum.binEdges, vcat(0.0, smo_spectrum.outputBinCounts ./ smi_spectrum_bin_widths ))
plot(mzs_bspline, intensities_bspline ./ (mzs_bspline[2] - mzs_bspline[1]))




















figure()
scatter(-1:length(smo_spectrum.controlPoints)-2, smo_spectrum.controlPoints)



axhline(color="grey")
axvline(x=x[1],color="grey")
axvline(x=x[end],color="grey")
plot(x, outputBspline)
scatter(-1:length(smo_spectrum.controlPoints)-2, smo_spectrum.controlPoints)

figure()
axhline(color="grey")
axvline(x=x[1],color="grey")
axvline(x=x[end],color="grey")
plot(x, intensities_bspline)
scatter(-1:length(fcs)-2, fcs)















figure()

axhline(color="grey")
plot(input_spectrum[1], input_spectrum[2])


grr = join(h5read(mzMLb_filename, "mzML"))


#17_617.55
mzMLspecID, mzs, gs, fs, fcs, offset = h5open("P02U_Swath_1__682_690.smo", "r") do file
  mzMLspecID = file["51_0/2/SpectrumIndex"][specID+1][1]
  is = file["51_0/2/SpectrumCountIndex"][specID:specID+1]
  mzs = vec(file["51_0/2/SpectrumMZ"][is[1]+specID:is[2]+specID])
  gs = vcat(0, vec(file["51_0/2/SpectrumCount"][is[1]+1:is[2]]) ./ diff(mzs))
  fs = vcat(0, vec(file["51_0/2/5/0/-9/fs"][is[1]+1:is[2]]) ./ diff(mzs))
  fcs = vec(file["51_0/2/5/0/-9/fcs"][:, specID+1])
  offset = h5readattr("P02U_Swath_1__682_690.smo", "51_0/2/5/0/-9/fcs")["Offset"][1]
  mzMLspecID, mzs, gs, fs, fcs, offset
end




coefs = vcat(fcs, 0, 0, 0, 0) # no idea why I have to pad this
knots = -3:length(coefs)-4
tck = (knots, coefs, 3)
x = linspace(0, length(fcs)-3, npoints)
intensities_bspline = interpolate.splev(x, tck, ext=1)

figure()
axhline(color="grey")
axvline(x=x[1],color="grey")
axvline(x=x[end],color="grey")
plot(x, intensities_bspline)
scatter(-1:length(fcs)-2, fcs)


mz_range = (offset + [0, length(fcs) - 3]) * 1.0033548378 / (60 * 2^2)
mzs_bspline = linspace(mz_range[1], mz_range[2], npoints)

figure()
subplot(211)

axhline(color="grey")
axvline(x=mz_range[1],color="grey")
axvline(x=mz_range[2],color="grey")
plot(mzML_spectrum[1], mzML_spectrum[2])
plot(mzs_bspline, intensities_bspline * 30)

subplot(212)
step(mzs, gs)
step(mzs, fs)
#plot(mzs_bspline, intensities_bspline / (mzs_bspline[2] - mzs_bspline[1]))
plot(mzs_bspline, intensities_bspline * (length(fcs) - 3) / (mz_range[2]-mz_range[1]))
