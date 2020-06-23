#!/bin/csh -f 

set QEPATH = $HOME/q-e-qe-6.4.1

setenv OMP_NUM_THREADS 1
set num_core = `grep 'core id' /proc/cpuinfo | sort -u | wc -l`

mkdir mix
mkdir plot
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
mv 1.UPF Xx1.00000.UPF
mv 2.UPF Xx0.00000.UPF
rm -f -r Xx1.00000.UPF
rm -f -r Xx0.00000.UPF
set upf_list=`ls *.UPF`
cd ..

set mass1UPF = `cat mass.dat | head -1`
set mass2UPF = `cat mass.dat | head -2 | tail -1`

touch cell_param.txt etot.txt upf_list.txt lat.dat lat_line.dat
foreach upf_name ( ${upf_list} )
  set upf_1st_name = `echo ${upf_name} | sed -e 's/.UPF//'`
  mkdir ${upf_1st_name}
  cd ${upf_1st_name}
  cp ../*.UPF ./
  cp ../mix/${upf_name} ./${upf_name}
  set x_value = `echo ${upf_name} | sed -e 's/Xx//' -e 's/.UPF//'`
  set mass = `echo "${x_value}*${mass1UPF}+(1.0-${x_value})*${mass2UPF}" | bc`
  sed 's/mix.UPF/'${upf_name}'/g' ../tmp.vc-relax.in > run.vc-relax.in 
  sed -i 's/mix_mass/'${mass}'/g' run.vc-relax.in
  sed 's/mix.UPF/'${upf_name}'/g' ../tmp.bands.in > run.bands.in
  sed -i 's/mix_mass/'${mass}'/g' run.bands.in
  mpirun -np ${num_core} ${QEPATH}/bin/pw.x < run.vc-relax.in | tee run.vc-relax.out
  mpirun -np ${num_core} ${QEPATH}/bin/pw.x < run.bands.in > run.bands.out
  mpirun -np ${num_core} ${QEPATH}/bin/bands.x < ../input.pp.in
  ../pwout2in.py run.vc-relax
  set n = `grep -n 'K_POINTS' run.vc-relax.new.in | sed 's/:.*//g'`
  awk -v line=${n} '{if(NR==line+1){print "  " $1*2 " " $2*2 " " $3*2 "  " $4 " " $5 " " $6}else{print $0}}' run.vc-relax.new.in > run.nscf.in
  sed -i 's/vc-relax/nscf/g' run.nscf.in
  sed -i 's/from_scratch/restart/g' run.nscf.in
  mpirun -np ${num_core} ${QEPATH}/bin/pw.x < run.nscf.in > run.nscf.out
  mpirun -np ${num_core} ${QEPATH}/bin/projwfc.x < ../input.pr.in
  grep "Fermi" run.nscf.out | sed 's/the Fermi energy is//g'| sed 's/ev/0.0/g' > ef.txt
  gnuplot ../band.gp
  gnuplot ../tdos.gp
  mv band.png ${upf_1st_name}_band.png
  cp ${upf_1st_name}_band.png ../plot/${upf_1st_name}_band.png
  mv tdos.png ${upf_1st_name}_tdos.png
  cp ${upf_1st_name}_tdos.png ../plot/${upf_1st_name}_tdos.png
  grep -A 3 "CELL_PARAMETERS" run.vc-relax.out | tail -4 >> ../cell_param.txt
  grep "!    total energy" run.vc-relax.out | tail -1 >> ../etot.txt
  #
  cp run.vc-relax.new.in run.scf.ph.in
  sed -i 's/vc-relax/scf/g' run.scf.ph.in
  mpirun -np ${num_core} ${QEPATH}/bin/pw.x < run.scf.ph.in
  mpirun -np ${num_core} ${QEPATH}/bin/ph.x < ../ph.in
  mpirun -np ${num_core} ${QEPATH}/bin/q2r.x < ../q2r.in
  mpirun -np ${num_core} ${QEPATH}/bin/matdyn.x < ../matdyn.in
  mpirun -np ${num_core} ${QEPATH}/bin/plotband.x < ../plotband.in > /dev/null
  gnuplot ../gnuplot.tmp
  mv "case.dispersions.png" ${upf_1st_name}.disp.png
  cp ${upf_1st_name}_disp.png ../plot/${upf_1st_name}_disp.png
  mpirun -np ${num_core} ${QEPATH}/bin/matdyn.x < ../phdos.in
  gnuplot ../gnuplot1.tmp
  mv "case.phdos.png" ${upf_1st_name}.phdos.png
  cp ${upf_1st_name}_phdos.png ../plot/${upf_1st_name}_phdos.png
  #
  echo ${upf_name} >> ../upf_list.txt
  set A = `grep "CELL_PARAMETERS" run.vc-relax.out | tail -1 | sed -e 's/.*(//' -e 's/alat=//' -e 's/).*//'`
  set a = `grep -A 1 "CELL_PARAMETERS" run.vc-relax.out | tail -1 | awk -v L="${A}" '{print L*($1^2+$2^2+$3^2)^0.5}'`
  echo "${x_value} ${a}" >> ../lat.dat
  cd ..
end
cat cell_param.txt
cat etot.txt
cat upf_list.txt

cat lat.dat | head -1 > lat_line.dat
cat lat.dat | tail -1 >> lat_line.dat
gnuplot lat.gp
cp lat.png ./plot/lat.png

