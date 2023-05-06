input_folder=$1
sample=$2
output_folder=$3
fa_file=$4

echo generate alignment file: $sample
java_zx="/zxu/tools/openjdk_v11/jdk-11.0.2/bin/java"
sam2tsv_zx="/zxu/tools/jvarkit/dist/sam2tsv.jar"
$java_zx -jar $sam2tsv_zx -R $fa_file $input_folder/$sample.sam \
> $output_folder/$sample.align
