gfortran  input_main.f main.f ulgeo.f ulsource.f ullib.f ulhistos.f AlfaMCLIB.f -o main.exe
echo %time%
main.exe < material.in
echo %time%
del main.exe
del *.mod
