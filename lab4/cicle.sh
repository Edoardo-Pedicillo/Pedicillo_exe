./clean.sh
g++ -o MolDyn_NVE.exe MolDyn_NVE.cpp
for i in {1..5}
do
  echo "------------ SIMULATION RUN N ..... $i ------- "
   ./MolDyn_NVE.exe
done
