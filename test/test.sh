JOB_FOLDER="/mnt/lustre/home/epantaleo/altra/data/"
GenePredToSim=$JOB_FOLDER"/GenePredToSim"
lambdas=2.,.8:.8,2.
LOCUS=chr11:448267-491000
line=1,2
pK=2
nK=0
RL=46
OVERHANG=9
read_labels=ind1,ind2
NORM=1,1

sim_sam -r $read_labels -g $GenePredToSim -L $LOCUS -l $line -R $RL -M $OVERHANG -m $lambdas -o $JOB_FOLDER

GenePredRef=$JOB_FOLDER"ensGene_e62_hg19.GenePred.gz"

altra -L $LOCUS -o $JOB_FOLDER -r $JOB_FOLDER/ind1/,$JOB_FOLDER/ind2/ -R $read_labels -c $NORM -a $RL -O $OVERHANG -g $GenePredRef -d $pK -y $nK -p 4000 -n 1000 -q 1000 -e 0
