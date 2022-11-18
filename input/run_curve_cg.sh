echo "set the number of bases in the curve"
echo "set the number of bases in the cgFrame"
#./curveplus_run.sh ${i}.nc ${i}.top ${i} ${i}

rm ${1}*fra ${1}_X.pdb  ${1}_B.pdb

./cgFrame_run.sh ${1}
./curveplus_run.sh ${1}.pdb ${1}.top ${1} ${1}

rm *lis *cda ${1}_X.pdb  ${1}_B.pdb
