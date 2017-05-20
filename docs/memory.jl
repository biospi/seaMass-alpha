using Plots
inspectdr()

f = open("../data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500/2.seamass/" *
         "P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.out");

ids = Int64[]
times = Float64[]
mems = Float64[]
for line in eachline(f)
    m = match(r"^\[ *(?<id>[0-9])+, *(?<time>[+-]?([0-9]*[.])?[0-9]+), *(?<mem>[+-]?([0-9]*[.])?[0-9]+)]", line)
    if m != nothing
        append!(ids, parse(Int64, m[:id]))
        append!(times, parse(Float64, m[:time]))
        append!(mems, parse(Float64, m[:mem]))
    end
end
close(f)

plt = plot(times, mems, label="Memory (Mb)")
gui(plt)
