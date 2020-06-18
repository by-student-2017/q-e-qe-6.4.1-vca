#!/bin/csh -f 

#set QEPATH = $HOME/q-e-qe-6.4.1/bin
set QEPATH = /usr/local/bin

setenv OMP_NUM_THREADS 1
set num_core = `grep 'core id' /proc/cpuinfo | sort -u | wc -l`

touch cell_param.txt etot.txt upf_list.txt
cd mix
set upf_list=`ls *.UPF`
cd ..
foreach upf_name ( ${upf_list} )
  cp ./mix/${upf_name} ./
  sed 's/mix.UPF/'${upf_name}'/g' tmp.vc-relax.in > run.vc-relax.in 
  mpirun -np ${num_core} ${QEPATH}/pw.x < run.vc-relax.in | tee run.vc-relax.out
  grep -A 3 "CELL_PARAMETERS" run.vc-relax.out | tail -4 >> cell_param.txt
  grep "!    total energy" run.vc-relax.out | tail -1 >> etot.txt
  echo ${upf_name} >> upf_list.txt
  rm ${upf_name}
end
cat cell_param.txt
cat etot.txt
cat upf_list.txt
