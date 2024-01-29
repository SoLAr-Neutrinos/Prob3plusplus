#! /bin/bash

# dm2_start=3.5e-5
# dm2_end=0.95e-4
# dm2_step=1e-5
dm2_start=7.4e-5
dm2_end=7.4e-5
dm2_step=1e-6

sin13_start=0.00
sin13_end=0.05
sin13_step=0.005
# sin13_start=0.021
# sin13_end=0.021
# sin13_step=0.01

sin12_start=0.265
sin12_end=0.355
sin12_step=0.005
# sin12_start=0.303
# sin12_end=0.303
# sin12_step=0.001

# Make a counter for the number of iterations
counter=1
for dm2 in $(seq $dm2_start $dm2_step $dm2_end)
do
    for sin12 in $(seq $sin12_start $sin12_step $sin12_end)
    do
        for sin13 in $(seq $sin13_start $sin13_step $sin13_end)
        do
            # Make a colorfull print to the terminal
            date
            echo -e "\e[0;36mIteration: $counter\e[0m"
            echo -e "\e[0;36mParameters: $sin12 $sin13 $dm2\e[0m"
            # Run the program with the parameters
            ./solar $sin12 $sin13 $dm2
            counter=$((counter+1))
        # echo $sin12
        done
    done
done

echo "Done!"