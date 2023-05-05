
input_folder=$1
sample=$2
output_folder=$3

echo Trimming sample: $sample
trim_galore $input_folder/$sample*.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $output_folder
echo Trimming $sample done.
