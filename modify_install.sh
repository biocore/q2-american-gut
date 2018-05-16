# this script filters out unneeded packages in the installation file to cut down on environment size

declare -a arr=("q2-alignment" "q2-composition" "q2-cutadapt" "q2-dada2" "q2-deblur"
                "q2-demux" "q2-emperor" "q2-gneiss" "q2-longitudinal"
                "q2-vsearch" "q2-sample-classifier" "q2templates", "r-" "bioconductor")

for i in "${arr[@]}"
do 
    grep -v "^- $i" qiime2-latest-py35-linux-conda.yml > qiime2-install-mod.yml; 
    mv qiime2-install-mod.yml qiime2-latest-py35-linux-conda.yml
done
