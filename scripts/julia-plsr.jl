using NCDatasets
using Statistics

results_dir = "/home/ashiklom/projects/sbg-uncertainty/isofit/examples/py-hypertrace/output/ht-trait-analysis"

ncfiles = [f for f in readdir(results_dir, join=true) if occursin(r".nc$", f)]

ncf = ncfiles[1]

nc = NCDataset(ncf)

refl = nc["estimated_reflectance"][:,:,1:10,1]
minimum(refl[:,1,1])
maximum(refl[:,1,1])
quantile(refl[:,1,1], [0.025, 0.1, 0.5, 0.9, 0.975])
