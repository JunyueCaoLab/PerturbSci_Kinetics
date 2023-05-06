input_folder=$1
sample=$2
output_folder=$3

split -l 50000000 $input_folder/${sample}.align --additional-suffix=.align $output_folder/${sample}
echo split alignment file done: $sample
