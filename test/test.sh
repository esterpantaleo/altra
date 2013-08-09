
# \file altra
#

# Copyright (C) 2012-2013 Ester Pantaleo


#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by


# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of


# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


#

function absdir() {
if [[ -n "$1" ]] ; then
    if [[ -d $1 ]] ; then
        pushd $1 2>&1 >/dev/null;
        echo $PWD;
        popd 2>&1 >/dev/null;
    else
        dname=`dirname $1`;
        pushd ${dname} 2>&1 >/dev/null;
        echo ${PWD};
        popd 2>&1 >/dev/null;
    fi
fi
}

BASEDIR=`dirname ${BASH_SOURCE-$0}`
BASEDIR=`absdir $BASEDIR`
JOB_FOLDER=${BASEDIR}"/../data/"
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

mkdir ${BASEDIR}"/../data/test"
${BASEDIR}"/../scripts/sim_sam" -r $read_labels -g $GenePredToSim -L $LOCUS -l $line -R $RL -M $OVERHANG -m $lambdas -o $JOB_FOLDER"test/"

GenePredRef=$JOB_FOLDER"ensGene_e62_hg19.GenePred.gz"

${BASEDIR}"/../scripts/altra" -L $LOCUS -o $JOB_FOLDER"test/" -r $JOB_FOLDER"test/ind1/",$JOB_FOLDER"test/ind2/" -R $read_labels -c $NORM -a $RL -O $OVERHANG -g $GenePredRef -d $pK -y $nK -p 4000 -n 1000 -q 1000 -e 0
