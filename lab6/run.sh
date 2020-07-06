./clean.sh

for ((k=50;k<201;k +=5))
do
    j=$(bc <<<"scale=2; $k / 100" )    
    sed -i '1 s/^.*$/'$j'/' input.dat
    ./Monte_Carlo_ISING_1D.exe

done
