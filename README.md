# q-e-qe-6.4.1-vca


-----


settings


1. cd ~/q-e-qe-6.4.1/upftools


2. cp -r ~/q-e-qe-6.4.1-vca/upftools/* ./


3. make


-----


python setting


1. sudo apt install -y python


2. chmod +x pwout2in.py


-----


Run


1. cp Ga.pbe-dnl-rrkjus_psl.1.0.0.UPF 1.UPF


2. cp In.pbe-dn-rrkjus_psl.1.0.0.UPF 2.UPF


3. gedit x.dat


4. chmod +x run.csh


5. ./run.csh


-----


Note -1-


1.UPF = Xx1.00000.UPF, X=1.0


2.UPF = Xx0.00000.UPF, X=0.0


-----


Note -2-


x.dat


X.XXXXX


Y.YYYYY


Z.ZZZZZ


...


-----


Note -3-


mass.dat


First Line: mass of 1.UPF


Second Line: mass of 2.UPF


-----


# Note: virtual_v3.x and virtual_v3.f90 are development version. Please, you would develop them.


-----