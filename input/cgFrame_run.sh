export LD_LIBRARY_PATH=/work/lcvmm/crea/cgFrames_Final/armadillo-9.100.5:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/work/lcvmm/crea/cgFrames_Final/netcdf-cxx4-4.3.0/lib:$LD_LIBRARY_PATH

rm ${1}_cgf.fra ${1}_cgf.pfra
/work/lcvmm/crea/cgFrames_Final/bin/genFrames ${1}.pdb ""   ideal_bases_tsukuba.txt ${1}_cgf.fra ${1}_cgf.pfra

#rm ${1}_cgf_nc.fra ${1}_cgf_nc.pfra
#/work/lcvmm/crea/cgFrames_Final/bin/genFrames ${1}.nc ${1}.top ideal_bases_tsukuba.txt ${1}_cgf.fra ${1}_cgf.pfra
