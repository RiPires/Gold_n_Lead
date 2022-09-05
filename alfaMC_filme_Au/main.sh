#/bin/bash

gfortran  input_main.f main.f ulgeo.f ulsource.f ullib.f ulhistos.f AlfaMCLIB.f -o main.exe
chmod +x main.exe
time ./main.exe < material.in
rm -f main.exe
rm -f *.mod
gnuplot plot.gnu
