JOB_FOLDER="/mnt/lustre/home/epantaleo/altra/data/"
GenePredToSim=$JOB_FOLDER"/GenePredToSim"
lambdas=1.,.1:.1,1.
LOCUS=chr11:484400-497300
line=1,2
pK=0
nK=2
RL=46
OVERHANG=9
read_labels=ind1,ind2
NORM=1,1

sim_sam -r $read_labels -g $GenePredToSim -L $LOCUS -l $line -R $RL -M $OVERHANG -m $lambdas -o $JOB_FOLDER



GenePredRef=$JOB_FOLDER"ensGene_e62_hg19.GenePred.gz"
altra -L $LOCUS -o $JOB_FOLDER -r $JOB_FOLDER/ind1/,$JOB_FOLDER/ind2/ -R $read_labels -c $NORM -a $RL -O $OVERHANG -g $GenePredRef -d $pK -y $nK -p 2000 -n 1000 -q 1000 -e 0
