#!/bin/bash
source /users/students/r0754386/Documents/bmtk/bin/activate
module load mpi
cd /users/students/r0754386/Documents/bmtk/examples/v1
echo '_______________'
TOT=0
DONE=0
for config_file in $(find /users/students/r0754386/Documents/bmtk/examples/v1/exp2/config -type f -print)
do
    let TOT++
    output_path=${config_file/config/"output"}
    file_name="$(basename "$output_path")"
    file_name="${file_name:7:-5}"
    output_path="$(dirname "$output_path")"
    output_path="${output_path}/${file_name}/spikes.csv"
    if [[ -f $output_path ]]; then
        echo "$output_path already exists"
        let DONE++
    else 
        echo "$output_path does not exist"
    fi
done
echo "$DONE/$TOT"