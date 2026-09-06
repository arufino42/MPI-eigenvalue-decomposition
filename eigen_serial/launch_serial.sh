#! /bin/bash

for n in 144 216 288 576 648 1152
do
    echo "n=$n"
    sbatch -t 120 --export=n=$n serial.run
done