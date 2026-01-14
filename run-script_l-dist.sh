for i in {1..99}
do
./a.out 10000 ../../LAMMPS_script/data/set1/Data_full_e0-10_S0-1_B-0.01_steps-100000.dump "$i"
done
