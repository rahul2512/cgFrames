export LD_LIBRARY_PATH=/work/lcvmm/crea/cgFrames_Final/armadillo-9.100.5:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/work/lcvmm/crea/cgFrames_Final/netcdf-cxx4-4.3.0/lib:$LD_LIBRARY_PATH
echo "-----------------------------------------------------"
echo "huh" $1  $2 $3 $4
/work/lcvmm/crea/cgFrames_Final/bin/genFrames $1  $2 /work/lcvmm/crea/cgFrames_Final/input/ideal_bases_tsukuba.txt $3 $4
