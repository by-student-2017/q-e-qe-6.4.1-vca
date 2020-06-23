#!/bin/csh -f 

set QEPATH = $HOME/q-e-qe-6.4.1

setenv OMP_NUM_THREADS 1
set num_core = `grep 'core id' /proc/cpuinfo | sort -u | wc -l`

mkdir mix
cd mix
cp ../1.UPF 1.UPF
cp ../2.UPF 2.UPF
cp ../x.dat x.dat_tmp
foreach x_value ( `cat x.dat_tmp` )
  echo ${x_value} > x.dat
  mpirun -np ${num_core} ${QEPATH}/upftools/virtual_v2_auto.x
  mv NewPseudo.UPF Xx"${x_value}".UPF
end
rm x.dat x.dat_tmp
mv 1.UPF Xx0.00000.UPF
mv 2.UPF Xx1.00000.UPF
set upf_list=`ls *.UPF`
cd ..

touch cell_param.txt etot.txt upf_list.txt
foreach upf_name ( ${upf_list} )
  cp ./mix/${upf_name} ./
  sed 's/mix.UPF/'${upf_name}'/g' tmp.vc-relax.in > run.vc-relax.in 
  sed 's/mix.UPF/'${upf_name}'/g' tmp.nscf.in > run.nscf.in
  mpirun -np ${num_core} ${QEPATH}/bin/pw.x < run.vc-relax.in | tee run.vc-relax.out
  mpirun -np ${num_core} ${QEPATH}/bin/pw.x < run.nscf.in > run.nscf.out
  mpirun -np ${num_core} ${QEPATH}/bin/projwfc.x < input.pr.in
  grep "Fermi" run.nscf.out | sed 's/the Fermi energy is//g'| sed 's/ev/0.0/g' > ef.txt
  gnuplot tdos.gp
  #mv tdos.png ${upf_name%UPF}_tdos.png
  set upf_1st_name = `echo ${upf_name} | sed -e 's/.UPF//'`
  mv tdos.png ${upf_1st_name}_tdos.png
  grep -A 3 "CELL_PARAMETERS" run.vc-relax.out | tail -4 >> cell_param.txt
  grep "!    total energy" run.vc-relax.out | tail -1 >> etot.txt
  echo ${upf_name} >> upf_list.txt
  rm ${upf_name}
end
cat cell_param.txt
cat etot.txt
cat upf_list.txt
