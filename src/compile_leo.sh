module purge
module load nvhpc/24.3
module load cuda/12.1
rm *.mod
rm mhit36
make clean
make
make

rm -rf output
mkdir output
rm *.o

#run the code
#./mhit36