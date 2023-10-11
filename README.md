# PhotobacteriumNanopore
# Scripts used for the genome assembly and analysis of Photobacterium spp. from ONT minION sequence data
# Nanopore genome assembly pipeline: Flye + Circlator + Medaka + Homopolish #
# you will need to setup and activate each conda environment before running each program


# change into dirctory with fastq. files
cd data/nanopore_photobact/fastqs


# first run filtlong on reads for QC #
# remove any reads shorter than 1 kbp and also exclude the worst 5% of reads #

# make new directory for filtered reads
mkdir filt

for i in *fastq.gz
do
filtlong --min_length 1000 --keep_percent 95 --split 1000 ${i} | gzip > filt1kb_95_split/${i%%_filt.fastq.gz}""
echo "${i}: Has been completed"
done





### run Flye for genome assmeblies ###

mkdir /data/nanopore_photobact/flye/

cd ~/data/nanopore_photobact/fastq/filt

for i in *fastq.gz
do
mkdir ~/data/nanopore_photobact/flye/${i%%.fastq.gz}""
flye --nano-raw ${i} --out-dir ~/data/SSI/flye/${i%%.fastq.gz}"" --threads 32
echo "${i}: Has been completed"
done





# run Circlator on assemblies #
cd data/nanopore_photobact/flye/filt

for i in *
do
circlator all --threads 24 ${i}/assembly.fasta ~/data/nanopore_photobact/fastq/filt/${i%%}".fastq.gz" ~/data/nanopore_photobact/flye/filt/circlator/${i}
echo "${i}: Has been completed"
done



# run medaka on circularized assemblies #

mkdir ~/data/nanopore_photobact/flye/filt/medaka

cd ~/data/nanopore_photobact/flye/filt/circlator

for i in barcode2*
do
medaka_consensus -i ~/data/nanopore_photobact/fastq/filt/${i%%}".fastq.gz" -d ${i}/assembly.fasta -o ~/data/nanopore_photobact/flye/filt/medaka/${i} -t 32
echo "${i}: Has been completed"
done


# run homopolisher on medaka consensus output #

mkdir ~/data/nanopore_photobact/flye/homopolish

cd ~/data/nanopore_photobact/flye/homopolish

#download bacteria DB# 
wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz
gunzip bacteria.msh.gz

cd ~/data/nanopore_photobact/flye/medaka/

#run for loop for homopolish - NOTE: this has a flag (-m) that indicates the flowcell chemistry -> need to lookup and use the R10 chemistry!!!
for i in *
do
python3 ~/homopolish/homopolish.py polish -a ${i}/consensus.fasta -g Photobacterium_leiognathi --download_contig_nums 20 -m R10.3.pkl -t 32 -o ~/data/nanopore_photobact/flye/filt/homopolish/${i}
echo "${i}: Has been completed"
done


# run ragtag on homopolished assemblies #

cd ~/data/nanopore_photobact/flye/filt/homopolish

# example for one sample:
ragtag_run.sh scaffold ~/data/nanopore_photobact/ref/Ik8.2.fa barcode08/*.fasta -o ~/data/SSI/nanopore_photobact/filt/ragtag/barcode08 -u -t 32




# run prokka on flye_homopolish90 draft assemblies #

cd ~/data/nanopore_photobact/flye/filt/final_asm
mkdir prokka


for i in *
do
prokka --kingdom Bacteria --outdir prokka/${i} --genus Photobacterium --prefix ${i} ${i} --cpus 32 --force
echo "${i}: Has been completed"
done




# run BUSCO on assemblies #

cd ~/data/nanopore_photobact/flye/filt/final_asm
mkdir busco

for i in *.fasta
do
busco -i ${i} -l vibrionales_odb10 -o busco/${i%%.fasta}"_busco" -m genome -c 32
echo "${i}: Has been completed"
done




# make softlinks for all of the .gffs to run in roary #

cd ~/data/nanopore_photobact/flye/filt/final_asm/prokka/

for i in *fasta
do
ln -s $i/*gff ~/data/nanopore_photobact/flye/filt/final_asm/prokka/gffs/${i%%.fasta.gff}".gff"
echo "${i}: Has been completed"
done




# run roary on annotated genomes #

mkdir ~/data/SSI/flye/filt/final_asm/roary

#first need to make soflinks to all of the .gff files in the prokka output directories - can do this in a for loop! (+ add some additional "outgroup" .gff files for tree
mkdir gffs
cd gffs

#example softlink: ln -s path_to_gff_file current_directory

#run w/out mafft and 95% genomes to be considered "core" gene 
roary -f out -cd 95 -e -p 32 -v gffs/*.gff

roary -f out2 -cd 95 -e -p 32 -v gffs/*.gff





# run IQTREE on full alignment file produced by roary #

cd ~/data/nanopore_photobact/roary/out

#use core_gene_alignment.aln output to build tree!!
#fist run iqtree with TESTONLY to determine the best model for the analysis 
iqtree -s core90.aln -m TESTONLY 
#best model predicted: GTR+F+I+G4 --> run with 1000 BS #

iqtree -s core90.aln -m GTR+F+I+G4 -B 1000 -alrt 1000 -T 32 -redo








# run fastANI with all strains #


cd ~/data/nanopore_photobact/flye/filt/final_asm/
mkdir ANI

# have to be in same directory as .fasta files (made softlinks to everything) #

for i in *fasta 
do 
ln -s ~/data/nanopore_photobact/flye/filt/final_asm/$i ~/data/nanopore_photobact/flye/filt/final_asm/ANI/$i
done

fastANI --ql strains.txt --rl refs.txt -o fastANI --matrix --threads 16

conda activate aniclustermap

ANIclustermap -i . -o cluster --cmap_colors white,yellow,red



