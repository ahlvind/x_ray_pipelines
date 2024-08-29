#! /bin/bash
# This script creates a smaller image around the SN region that is used by wavdetect. The smaller image is used to minimise the runing time of wavdetect.

# Start the timer
start_time=$(date +%s)

input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_Galax.txt"

cd ../../Chandra/data
{
    read -r
    while read name obsid epoch Type dist expt efft dateSN dateXray ra dec nH
    do

        cd $obsid
        cd "repro"
        
        # Read the original line from the file
        original_line=$(cat "src_$name.reg")

        # Extract the last number and multiply it by x
        last_number=$(echo "$original_line" | awk -F, '{gsub(/[^0-9.]/, "", $3); print $3}')
        updated_last_number=$(awk "BEGIN {printf \"%.5f\", $last_number * 10}")

        # replace , with . for the decimal number
        updated_last_number2=$(echo "$updated_last_number" | tr ',' '.')
        
        # Update the original line with the new last number
        updated_line=$(echo "$original_line" | awk -v new_number="$updated_last_number2" 'BEGIN{OFS=FS=","}{$3=new_number; print}')
        
        # Extract the last two characters of the original line
        last_two_chars="${original_line: -2}"

        # Update the original line with the new last number and last two characters
        updated_line=$(echo "$original_line" | awk -v new_number="$updated_last_number2" -v last_chars="$last_two_chars" 'BEGIN{OFS=FS=","}{$3=new_number; print $0 last_chars}')

        echo "Old string: $original_line"
        echo "New string: $updated_line"

        dmcopy "image_0312_bin1.fits[sky="$updated_line"]" "image03_12_${name}.fits" clobber=yes
            
        
        cd ../..
        
    done
}< "$input"

# Stop the timer
end_time=$(date +%s)

# Calculate and display the elapsed time
elapsed_time=$((end_time - start_time))
echo "Script execution time: $elapsed_time seconds"

cd /Users/juliaahlvind/Documents/projekt_1/pipelines/chandra
