# \file utils.sh
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


#================== FUNCTION ============================================================== 
#         NAME: absdir                                                    
#  DESCRIPTION: This function returns the absolute path to the given path, i.e. the absolute 
#               path if the argument is a directory, or the absolute basename otherwise.   
#   ARGUMENT 1: path (a string)                          
#==========================================================================================
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


#-------------------------------------------------------------------------------------------
#   Get the relative path of this script, use that to load the                         
#   utilities, and finally use the absdir function from the utilities to                    
#   transform the relative path to an absolute path.                                  
#-------------------------------------------------------------------------------------------
BASEDIR=`dirname ${BASH_SOURCE-$0}`
BASEDIR=`absdir $BASEDIR`           
utils_awk=$BASEDIR"/utils.awk"


#============= FUNCTION =========================             
#         NAME:     try  
#================================================
function try() {
    "$@"
    status=$?
    if [ $status -ne 0 ]; then
        echo "error with $1"
    fi
    return $status
}

#============= FUNCTION =========================
#         NAME:     tic, tac, print_tic_tac
#================================================  
function tic() {
    _START_TIME=`date  +%s`
}

function tac() {
    _END_TIME=`date  +%s`
}

function print_tic_tac() {
    echo  "$((_END_TIME - _START_TIME)) secs."
}
 
#============= FUNCTION ==========================                                            
#         NAME:     regionToBed
#================================================= 
function regionToBed(){
    echo $1 | tr ":" "\t" | tr "-" "\t"
}

#============= FUNCTION ===========================
#         NAME:     my_genePredLineToBed
#  DESCRIPTION:     reads in a specific line of a  GenePred file and 
#                   converts it into Bed file
#   ARGUMENT 1:     GenePred line through stdin
#=================================================
function my_genePredLineToBed() {
    local _var
    read _var
    echo $_var |  awk '
     BEGIN {
       OFS="\t";
     }
     {
       name=$1
       chrom=$2
       strand=$3
       l = split($9, starts, ",");
       split($10, ends, ",");
       for (i=1; i<l; i++) 
           print chrom, starts[i], ends[i], name, 1, strand
     }'
}

#============= FUNCTION =======================================================================
#         NAME:     genePredView                                             
#        USAGE:     genePredView my_genePred_file +                               #print transcripts on positive strand
#                   genePredView my_genePred_file + which                         #print which lines in the file are on the positive strand
#                   genePredView my_genePred_file filter MIN_IN_LEN MIN_EX_LEN RL #modify the file if it has short introns, short exons
#                   genePredView my_genePred_file chr3:11133752-11135688 which    #print which lines in the file contain transcripts that overlap with the region chr3:11133752-11135688
#                   genePredView my_genePred_file chr3:11133752-11135688 compressed #print transcripts that overlap with the region chr3:11133752-11135688 given an gzipped genePred file
#                   genePredView my_genePred_file chr3:11133752-11135688 multiple_compressed #print transcripts that overlap with the region chr3:11133752-11135688 given a list of comma separated  gzipped genePred files
#                   genePredView my_genePred_file chr3:11133752-11135688          #print transcripts that overlap with the region chr3:11133752-11135688
#                   genePredView my_genePred_file chr3:11133752-11135688 cut
#==============================================================================================  
function genePredView(){
    if [[ $VERBOSE -eq 1 ]];then
        echo command line: "$FUNCNAME $@" >&2
        echo "Reading GenePredFile" >&2
    fi

    if [ -z $1 ] || [ -z $2 ] ;then
        echo "ERROR: incorrect number of arguments in function $FUNCNAME" >&2
        exit 1;
    fi
    local _GenePred=$1
    local _LOCUS=$2
    local _lines=$3

    case $_LOCUS in
	"+" | "-")
	    case "$_lines" in
		"which" )
		    cat $_GenePred | awk -v _sign="$_LOCUS" '                         
                       BEGIN{                                               
                          counter=0                       
                       }         
                       {                
                          counter++;   
                          if ($3 == _sign)                                
                             printf("%s,", counter)    
                       }'
		    ;;
		* ) 
		    cat $_GenePred | awk -v _sign="$_LOCUS" '
                       {     
                          if ($3 == _sign)                          
                             print     
                       }'
		    ;;
	    esac
	    ;;
	    
	"filter")
	    local _MIN_IN_LEN=$3
	    local _MIN_EX_LEN=$4
            local _RL=$5
	    cat ${_GenePred} | awk -v min_in_len="$_MIN_IN_LEN" -v read_length="$_RL" '{
               #if intron length <_MIN_IN_LEN merge exons     
               l = split($9, field9, ",");                                                       
               split($10, field10, ",");
               #remove transcripts that are shorter than the read length
               transcript_length=0 
               for (i=1; i < l; i++){
                  transcript_length += field10[i]
                  transcript_length -= field9[i]
               }
               if (transcript_length >= read_length){              
                  new_field9 = field9[1]",";                                 
                  new_field10 = "";                                             
                  counter = 0;                                                      
                  for (i = 1;i < l - 1; i ++){                                             
                     if (field9[i + 1] - field10[i] >= min_in_len){               
                        new_field9 = new_field9 field9[i + 1]",";          
                        new_field10 = new_field10 field10[i]",";              
                     } else                                                    
                     counter ++;                                              
                  }                                                          
                  new_field10 = new_field10 field10[l - 1]",";                       
                  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, $3, $4, $5, $6, $7, $8 - counter, new_field9, new_field10);       
                  for (i = 11; i < NF; i ++) 
                     printf("%s\t", $i); 
                  printf("%s\n", $NF);
               }                                                              
               }' | awk -v min_ex_len="${_MIN_EX_LEN}" '{ 
               #if exon length < _MIN_EX_LEN remove exon; if number of removed exons < number of exons then print else don t print
               l = split($9, field9, ",");                      
               split($10, field10, ",");                                     
               counter = 0;                 
               new_field9 = "";                                                 
               new_field10 = "";                                         
               for (i = 1; i < l; i ++){                                        
                  if (field10[i] - field9[i] >= min_ex_len){                       
                     new_field9 = new_field9 field9[i]",";                      
                     new_field10 = new_field10 field10[i]",";                 
                  } else                                                 
                  counter ++;                                           
               }
               if (counter < l - 1){                                                          
                  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, $3, $4, $5, $6, $7, $8 - counter, new_field9, new_field10);
                  for (i = 11; i < NF; i ++) 
                     printf("%s\t", $i); 
                  printf("%s\n", $NF);
               }                                                              
               }' | sort -k 9,9 -k 10,10 | awk '                   
               BEGIN{                                                   
               counter = 0    
               }                                                 
               {                                                      
               counter ++;                                                   
               if (counter > 1){                                           
                  if (!(col9 == $9 && col10 == $10))                           
                     print $0                                          
               } else                                                                            
                  print $0
               col9 = $9;                                                     
               col10 = $10;                                   
               }'
	    ;;
	* )
	    local _LOCUS_v=( `regionToBed $_LOCUS` )
            local _chr=${_LOCUS_v[0]}
	    local _locusStart=${_LOCUS_v[1]}
	    local _locusEnd=${_LOCUS_v[2]}
	    
	    case "$_lines" in
		"which" )
		    cat $_GenePred | awk \
                        -v _regStart="$_locusStart" -v _regEnd="$_locusEnd" -v _chrom="$_chr" ' 
                       BEGIN{counter = 0}                                   
                       {                                              
                          counter++;                                     
                          if ($2 == _chrom && (($4 >= _regStart && $4 < _regEnd)||($5 > _regStart && $5 <= _regEnd)||($4 <= _regStart  &&  $5 >= _regEnd)))     
                             printf("%s,",counter)               
                       }'
		    ;;
		"compressed" )
		    zcat $_GenePred | awk \
                            -v _regStart="$_locusStart" -v _regEnd="$_locusEnd" -v _chrom="$_chr" '                 
                        BEGIN{OFS = "\t"}
                        {                                             
                           if ($2 == _chrom && (($4 >= _regStart && $4 < _regEnd)||($5 > _regStart && $5 <= _regEnd)||($4 <= _regStart  &&  $5 >= _regEnd))) 
                              print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10              
                        }'
		    ;;
		"multiple_compressed" )
                    # extract from the genome-wide annotation files transcripts that       
                    # (partially) overlap with the locus (remove identical transcripts        
                    # in multiple annotations) 
                    local _i
		    
                    for _i in ${_GenePred//,/ }; do
                        zcat $_i | awk \
                            -v _regStart="$_locusStart" -v _regEnd="$_locusEnd" -v _chrom="$_chr" '
                            {                                                
                               if ($2 == _chrom && (($4 >= _regStart && $4 < _regEnd)||($5 > _regStart && $5 <= _regEnd)||($4 <= _regStart  &&  $5 >= _regEnd))) 
                                  print $0                                   
                            }'
                    done | sort -k 9,9 -k 10,10 | awk '                                          
                            BEGIN{OFS = "\t"; counter = 0}                        
                            {                                                  
                               counter ++;                                    
                               if (counter > 1){                                      
                                  if (!(col9 == $9 && col10 == $10))               
                                     print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10        
                               } else                                                     
                               print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10               
                               col9 = $9;                         
                               col10 = $10;         
                            }'
                    ;;
		"cut" )
		    local _BedOut=$_GenePred"_tmp"
                    regionToBed $_LOCUS > $_BedOut
		    
		    for (( _l=1; _l<=`cat $_GenePred | wc -l`; _l++ )); do
			cat $_GenePred | \
                        sed -n -e "$_l"p | \
                        my_genePredLineToBed | \
                        ${INTERSECT} -a ${_BedOut} -b stdin -wb | \
                        sort -k 2,2n -k 3,3n | \
                        awk -f $utils_awk --source '
                         BEGIN{OFS = "\t"}   
                         {start[NR] = $2; end[NR] = $3}
                         END{
                            print $7, $1, $9, start[1], end[NR], start[1], end[NR], NR, join(start, 1, NR, ","), join(end, 1, NR, ","), 0
                         }'
		    done
                    rm $_BedOut
		    ;;
		* )
		    cat $_GenePred | awk \
			-v _regStart="$_locusStart" -v _regEnd="$_locusEnd" -v _chrom="$_chr" '
                        {                                                           
                            if ($2 == _chrom && (($4 >= _regStart && $4 < _regEnd) || ($5 > _regStart && $5 <= _regEnd)))      
                               printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10)     
                        }'
		    ;;
	    esac
    esac
}

#============== FUNCTION ===========================       
#         NAME:      numCommonSS                               
#  DESCRIPTION:      given two comma separated lists (ending with a comma)
#                    returns how many elements the two lists share 
#   ARGUMENT 1:      list1                                        
#   ARGUMENT 2:      list2 
#===================================================     
function numCommonSS() {
    local list1=$1
    local list2=$2
    local l1
    local l2
  
    echo " " | awk -v l1="$list1" -v l2="$list2" '
        BEGIN{
        counter = 0;
        }
        {
        L1 = split(l1, list1, ",");
        L2 = split(l2, list2, ",");
        for (i = 1; i < L1; i ++) { 
            for (j = 1; j < L2; j ++) { 
                if (list1[i] == list2[j]) 
                    counter ++;
            }
        }
        }
        END{ 
        print counter
        }'
}

#============== FUNCTION ===========================                                  
#         NAME:      ListRemoveFirst                             
#  DESCRIPTION:      given a comma separated list
#                    returns a comma separated list without the first element           
#   ARGUMENT 1:      list1 (a comma separated list ending with comma)    
#      EXAMPLE:      ListRemoveFirst 1,2,3,4,5, > file.txt; cat file.txt 
#                    >2,3,4,5,                          
#=================================================== 
function ListRemoveFirst() {
    local list1=$1
    local l1
    echo " " | awk -v l1="$list1" '
        {
        L1 = split(l1, list1, ","); 
        for (i = 2; i < L1; i ++) 
            printf("%s, ", list1[i])
        }'
}

#============== FUNCTION ===========================                                 
#         NAME:      ListRemoveLast                                                       
#  DESCRIPTION:      given a comma separated list (ending with a comma) 
#                    returns a comma separated list without the last element       
#   ARGUMENT 1:      list1 (a comma separated list ending with comma)         
#      EXAMPLE:      ListRemoveLast 1,2,3,4,5, > file.txt; cat file.txt
#                    >1,2,3,4,
#===================================================   
function ListRemoveLast(){
    local list1=$1
    local l1
    echo " " | awk -v l1="$list1" '
        {
        L1 = split(l1, list1, ","); 
        for (i = 1; i < L1 - 1; i ++) 
            printf("%s, ", list1[i])
        }'
}

#============== FUNCTION ===========================                         
#         NAME:      select_lines_randomly                                     
#  DESCRIPTION:      select lines randomly from file $1 with repetition
#   ARGUMENT 1:      the input GenePred file      
#   ARGUMENT 2:      the log file
#         NOTE:      modifies nK and pK
#===================================================
function select_lines_randomly(){
    local _GenePredIn=$1
    local _pK=$2
    local _nK=$3

    local _l

    if [[ $VERBOSE -eq 1 ]]; then
        echo command line: "$FUNCNAME $@" >&2
	echo "Randomly selecting $_pK positive transcripts and $_nK negative transcripts from file $_GenePredIn" >&2
    fi

    if [ $_pK -ge 1 ]; then
    # sample pos lines randomly
	if [ -z `genePredView $_GenePredIn + which` ]; then
	    _pK=0
	else
	    _l=0
	    while (( $_l<$_pK )); do
		genePredView $_GenePredIn + | awk '       
               BEGIN{ 
                  srand() 
               }     
               { 
                  print rand() "\t" $0 
               }' | 
		sort -n |   # Sort numerically on first (random number) column          
		cut -f2- |  # Remove sorting column                      
		head -n $_pK >> $_GenePredIn"_pos_tmp"
		_l=`cat $_GenePredIn"_pos_tmp" | wc -l`
	    done 
	    cat $_GenePredIn"_pos_tmp" | head -n $_pK
	    rm $_GenePredIn"_pos_tmp"
	fi
    fi
    
    if [ $_nK -ge 1 ]; then
    # sample neg lines randomly                                          
	if [ -z `genePredView $_GenePredIn - which` ]; then
	    _nK=0
	else
	    _l=0
            while (( $_l<$_nK )); do
		genePredView $_GenePredIn - | awk '                      
               BEGIN{                                                         
                  srand()                                          
               }                                                         
               {                                      
                  print rand() "\t" $0                   
               }' | 
		sort -n |   # Sort numerically on first (random number) column    
		cut -f2- |  # Remove sorting column                      
		head -n $_nK >> $_GenePredIn"_neg_tmp"
		_l=`cat $_GenePredIn"_neg_tmp" | wc -l`
            done 
	    cat $_GenePredIn"_neg_tmp" | head -n $_nK
	    rm $_GenePredIn"_neg_tmp"
	fi
    fi
}
 
#============== FUNCTION ===========================                            
#         NAME:      write_SSfile  
#===================================================  
function write_SSfile(){
    if [[ $VERBOSE -eq 1 ]];then
        echo command line: "$FUNCNAME $@" >&2
        echo "Writing Splice Sites file" >&2
    fi

    if [ -z $1 ] || [ -z $2 ];then
	echo "ERROR: incorrect number of arguments for function $FUNCNAME" >&2;
    fi

    local _GenePredIn=$1
    local _JunctionsIn=$2

    #write SSfile
    #get junctions, start and end from genes annotated in $GenePredIn
    cat $_GenePredIn | awk '{
        l = split($9, field9, ",");
        split($10, field10, ",");
        for (ll = 1; ll < l; ll ++){
           printf("%s %s 3SS %s\n", $2, $3, field9[ll] + 1);
           printf("%s %s 5SS %s\n", $2, $3, field10[ll]);
        }
     }' | sort -k 2,2 -k 3,3 -k 4,4n | uniq > $_GenePredIn"_SSfile"
    cat $_JunctionsIn | awk '{printf("%s %s 5SS %s 1\n%s %s 3SS %s 1\n", $1, $4, $2, $1, $4, $3)}' |\
    sort -k 2,2 -k 3,3 -k 4,4n | uniq > $_JunctionsIn"_SSfile"
    
    _l=`cat $_JunctionsIn"_SSfile" | wc -l`
    if [ $_l -eq 0 ]; then
	cat $_GenePredIn"_SSfile" 
    else
	(for (( _i=1;_i<=$_l;_i++ ));do
	    _my_SS=`sed -n -e "$_i"p $_JunctionsIn"_SSfile" | awk '{for (i=1;i<=NF;i++) printf("%s,", $i)}'` 
            cat $_GenePredIn"_SSfile" | awk -v my_SS="$_my_SS" '{
               l = split(my_SS,my_SS_v,","); 
               if (!(my_SS_v[2] == $2 && 
                   my_SS_v[3] == $3 && 
                   my_SS_v[4] - $4 < 3 &&  
                   $4 - my_SS_v[4] < 3)) print $0
                }'  
        done
	cat $_JunctionsIn"_SSfile" )| sort -k 2,2 -k 3,3 -k 4,4n | uniq
    fi 
    rm $_GenePredIn"_SSfile"
    rm $_JunctionsIn"_SSfile"
}

#============== FUNCTION ===========================       
#         NAME:      write_lists                                                  
#  DESCRIPTION:      from the set of splice sites found by the mapping method, annotated or found 
#                    by FLLat (after filtering) write_lists defines a set of coordinates (or 
#                    equivalently intervals); each coordinate will correspond to a 
#                    positive splice site and/or a negative splice site. Note that splice sites 
#                    found by FLLat are saved as both positive and negative splice sites; 
#                    splice sites found by the mapping method have a directionality that can be
#                    inferred from the consensus splice sequence in the intron)
#       OUTPUT:      listCoordinate=locusStart,coor1,coor2,coor3,LocusEnd ordered and unique  
#===================================================    
function write_lists(){
    if [[ $VERBOSE -eq 1 ]];then
        echo command line: "$FUNCNAME $@" >&2
        echo "Writing list of Coordinates" >&2
    fi

    if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ];then
        echo "ERROR: incorrect number of arguments for function $FUNCNAME" >&2;
    fi

    local _LOCUS=$1
    local _locusStart=$2
    local _locusEnd=$3
    local _SSfile=$4

    listCoordinate=""
    local _alistpos5=`cat $_SSfile | grep -v " -" | grep "5SS" | awk '{printf("%s ",$4)}'`
    if [[ -z $_alistpos5 ]];then  
	_alistpos5="NA";
    else
	local _alistpos5_v=(`echo $_alistpos5" "$_locusEnd   | tr " " "\n" | sort -g -u`)
	listCoordinate=$listCoordinate" "`echo $_alistpos5" "$_locusEnd | awk '{for (i=1;i<=NF;i++) printf("%19.1f  ",$i+0.5)}'`
    fi
    
    local _alistpos3=`cat $_SSfile | grep -v " -" | grep "3SS" | awk '{printf("%s ",$4)}'`
    if [[ -z $_alistpos3 ]];then  
	_alistpos3="NA";
    else
	#using locusEnd as a 3' Splice Site
	local _alistpos3_v=(`echo $_alistpos3" "$_locusStart | tr " " "\n" | sort -u -g`)
	listCoordinate=$listCoordinate" "`echo $_alistpos3" "$_locusStart | awk '{for (i=1;i<=NF;i++) printf("%19.1f ",$i-0.5)}'`
    fi
    
    local _alistneg5=`cat $_SSfile | grep -v " +" | grep "5SS" | awk '{printf("%s ",$4)}'`
    if [[ -z $_alistneg5 ]];then  
	_alistneg5="NA";
    else
	local _alistneg5_v=(`echo $_alistneg5" "$_locusEnd | tr " " "\n" | sort -g -u`)
	listCoordinate=$listCoordinate" "`echo $_alistneg5" "$_locusEnd | awk '{for (i=1;i<=NF;i++) printf("%19.1f ",$i+0.5)}'`  
    fi
    
    local _alistneg3=`cat $_SSfile | grep -v " +" | grep "3SS" | awk '{printf("%s ",$4)}'`
    if [[ -z $_alistneg3 ]];then  
	_alistneg3="NA";
    else
	_alistneg3_v=(`echo $_alistneg3" "$_locusStart | tr " " "\n" | sort -g -u`)
	listCoordinate=$listCoordinate" "`echo $_alistneg3" "$_locusStart | awk '{for (i=1;i<=NF;i++) printf("%19.1f ",$i-0.5)}'`
    fi

    listCoordinate_v=(`echo $listCoordinate | tr " " "\n" | sort -g -u`)
    lC=${#listCoordinate_v[@]}
    listCoordinate=`echo ${listCoordinate_v[@]} | sed 's/ /,/g'`  

    if [[ "$_alistpos5" = "NA" ]];then
	listpos5="NA"
    else
	listpos5=""
	for _i in `echo ${_alistpos5_v[@]}`;do
	    _counter=0
	    for _j in `echo ${listCoordinate_v[@]}`;do
		_jj=`awk -v _jj="$_j" 'BEGIN{print _jj-0.5}'`
		if [[ $_i -eq $_jj ]];then
		    listpos5=$listpos5$_counter","
		    break;
		fi
		let _counter=$_counter+1
	    done
	done
    fi
    if [[ "$_alistpos3" = "NA" ]];then
        listpos3="NA"
    else
        listpos3=""
        for _i in `echo ${_alistpos3_v[@]}`;do
            _counter=0
            for _j in `echo ${listCoordinate_v[@]}`;do
                _jj=`awk -v _jj="$_j" 'BEGIN{print _jj+0.5}'`
                if [[ $_i -eq $_jj ]];then
                    listpos3=$listpos3$_counter","
		    break; 
		fi
                let _counter=$_counter+1
            done
        done
    fi
    if [[ "$_alistneg5" = "NA" ]];then
        listneg5="NA"
    else
        listneg5=""
        for _i in `echo ${_alistneg5_v[@]}`;do
            _counter=0
            for _j in `echo ${listCoordinate_v[@]}`;do
                _jj=`awk -v _jj="$_j" 'BEGIN{print _jj-0.5}'`
                if [[ $_i -eq $_jj ]];then
                    listneg5=$listneg5$_counter","
		    break;
                fi
                let _counter=$_counter+1 
            done
        done
    fi
    if [[ "$_alistneg3" = "NA" ]];then
        listneg3="NA"
    else
        listneg3=""
        for _i in `echo ${_alistneg3_v[@]}`;do
            _counter=0
            for _j in `echo ${listCoordinate_v[@]}`;do
                _jj=`awk -v _jj="$_j" 'BEGIN{print _jj+0.5}'`
                if [[ $_i -eq $_jj ]];then
                    listneg3=$listneg3$_counter","
		    break;
                fi
                let _counter=$_counter+1
            done
        done
    fi
}

#IMPORTANT: if the 3rd and 4th arguments are provided this function is strictly looking at reads that start after locus_start and end before locus_end             
function bam2JunctionBed6(){
    local _bamfile=$1
    local _chr=$2
    local _locus_start=$3
    local _locus_end=$4

    if [[ -z $_bamfile ]]; then
        echo "ERROR: wrong number of arguments for function $FUNCNAME" >&2
        exit 1
    fi

    if [[ -z $_locus_start ]] && [[ -z $_locus_end ]] && [[ -z $_chr ]]; then
        _locus=""
    elif [[ -z $_locus_start ]]; then
        echo "ERROR: wrong number of arguments for function $FUNCNAME" >&2
        exit 1
    elif [[ -z $_locus_end ]]; then
        echo "ERROR: wrong number of arguments for function $FUNCNAME" >&2
        exit 1
    elif [[ -z $_chr ]]; then
        echo "ERROR: wrong number of arguments for function $FUNCNAME" >&2
        exit 1
    else
        _locus=$_chr":"$_locus_start"-"$_locus_end
    fi
#flags 0x200 0x400 remove low quality and duplicates                           
    $SAMTOOLS view -F 0x200 -F 0x400 $_bamfile $_locus | \
        awk -v s="$_locus_start" -v n="$_locus_end" '            
        {if (!($6~/N/)) next; if (s!="") if (s>$4) next; if ($6~/[PSHX=*]/) next;   
             cigar_scalar=$6                                                  
             tophat_strand=$20            
             #get strand                                                  
             if (index(tophat_strand, "XS:A:")!=1){                       
                for (i=12;i<=NF;i++){
                   if (index($i, "XS:A:")==1){               
                      tophat_strand=$i                               
                      strand=substr(tophat_strand,length(tophat_strand),1)  
                      break                                                    
                   }                                     
                   strand="."                                                 
                }                                                             
             }else                                                          
             strand=substr(tophat_strand,length(tophat_strand),1)                                 
             gsub(/[MNID]/," & ",cigar_scalar)           
             c=split(cigar_scalar,cigar," ")              
             start=$4;                                
             for (i=2;i<=c;i=i+2){                                         
                if (cigar[i]=="M")              
                   start=start+cigar[i-1]-1                              
                else if (cigar[i]=="N"){                             
                   end=start+cigar[i-1]+1              
                   if (n=="")                               
                      print $3,start,end,strand     
                   else if (end<=n)                                        
                      print $3,start,end,strand                          
                   start=end                                 
                }else if (cigar[i]=="I")                                  
                   start=start-cigar[i-1]         
                else if (cigar[i]=="D")    
                   start=start+cigar[i-1]                            
             }                                              
        }' |  \
        sort -k1,1 -k2,2n -k3,3n -k4,4 | \
        uniq -c | \
        awk 'BEGIN{OFS="\t"}{print $2,$3,$4,".",$1,$5}'
}



#============== FUNCTION ============================================
#         NAME:      write_JunctionsIn 
#  DESCRIPTION:      extract junctions from bam files     
#====================================================================
function write_JunctionsIn(){
    if [[ $VERBOSE -eq 1 ]];then
        echo command line: "$FUNCNAME $@" >&2
        echo "Extracting splice sites from bam files" >&2
    fi

    if [ -z $1 ] || [ -z $2 ]; then
        echo "ERROR: incorrect number of arguments for function $FUNCNAME" >&2
        exit 1;
    fi
    local -a _readfiles_v=("${!1}")
    local _chr=$2
    local _locus_start=$3
    local _locus_end=$4

    local _j_tmp=`mktemp`
    local _j_2tmp=`mktemp`
    local _j_3tmp=`mktemp`
    local _j_4tmp=`mktemp`
    local _N=${#_readfiles_v[@]}

    for ((_i=0;_i<$_N;_i++)); do
        bam2JunctionBed6 ${_readfiles_v[$_i]} $_chr $_locus_start $_locus_end | \
            awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,$5}' > $_j_tmp
        if [ $_i -eq 0 ]; then
            mv $_j_tmp $_j_3tmp
        else
            for ((_ll=`cat $_j_3tmp | wc -l`;_ll>0;_ll--));do
                _firstline=`head -1 $_j_3tmp`
                _grepheader=`echo "$_firstline" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}'`
                _count=`grep "$_grepheader" $_j_tmp | awk '{print $5}'`
                if [ "$_count" = "" ]; then _count=0; fi
                echo $_firstline $_count >> $_j_2tmp
                grep -v "$_grepheader" $_j_tmp > $_j_4tmp; mv $_j_4tmp $_j_tmp
                tail -n +2 $_j_3tmp > $_j_4tmp; mv $_j_4tmp $_j_3tmp
            done
            if [ -s $_j_2tmp ]; then
                mv $_j_2tmp $_j_3tmp
            fi
            if [ -s $_j_tmp ];then
                cat $_j_tmp | awk -v I="$_i" '{          
                   printf("%s %s %s %s ", $1, $2, $3, $4);    
                   for (i = 1; i < I + 1; i ++)      
                      printf("0 ");   
                   printf("%s\n", $5); 
                }' >> $_j_3tmp
                rm $_j_tmp
            fi
        fi
    done
    echo chr 5ss 3ss strand ${_readfiles_v[@]}
    cat $_j_3tmp | awk -v n="$_N" '{        
       out = ""          
       for (i = NF - 4; i < n; i ++)          
          out = out" 0"                          
       printf("%s %s\n", $0, out)   
    }'
    if [ -e $_j_tmp ];then
        rm $_j_tmp
    fi
    rm $_j_3tmp
}

#============== FUNCTION =======================================================================
#         NAME:      filter_JunctionsIn                               
#  DESCRIPTION:      from the set of junctions in input file, 
#                    only take those junctions 
#                    whose relative expression is greater that _PERCENT_J (over all individuals) 
#                    and that are supported by at least _MINIMUM_J reads in at least one individual
#=============================================================================================== 
filter_JunctionsIn(){
    local _junctions=$1
    local _PERCENT_J=$2
    local _MINIMUM_J=$3

    if [[ $VERBOSE -eq 1 ]];then
        echo command line: "$FUNCNAME $@" >&2
        echo "Filtering junctions" >&2
    fi
    
    local _j_tmp=`mktemp`
    let wc_j=`cat $_junctions | wc -l`-1
    if [ $wc_j -ge 1 ];then
        tail -n $wc_j $_junctions | awk '{
           sum = 0; 
           for (i = 5; i <= NF; i ++) 
              sum = sum + $i; 
           print $0, sum
        }' > $_j_tmp
	sum_j=`cat $_j_tmp | awk 'BEGIN{sumj=0;}{sumj+=$NF}END{printf("%s",sumj);}'`
        cat $_j_tmp | awk -v sumj="$sum_j" -v wcj="$wc_j" -v pj="$_PERCENT_J" -v mj="$_MINIMUM_J" '{ 
           if ($NF > sumj * pj / (100 * wcj)){
              ind = 0;
              for (i = 5; i < NF; i ++){
                 if ($i >= mj) 
                    ind += 1;
              } 
              if (ind != 0) 
                 print $0
           }
           }' | sort -k 2,2n -k 3,3n 
    rm $_j_tmp
    fi
}

#============== FUNCTION ======================================================================
#         NAME:      bitset2GenePred                          
#  DESCRIPTION:      convert line _which_line of summary.txt file to GenePred file    
#   ARGUMENT 1:      file name                                                         
#   ARGUMENT 2:      line in file                                          
#   ARGUMENT 3:      chr                           
#   ARGUMENT 4:      listCoordinate                                   
#==============================================================================================
function bitset2GenePred() {
    _summary=$1
    _which_line=$2
    _chr=$3    
    _listCoordinate=$4
  
    if [[ ${VERBOSE} -eq 1 ]];then
        echo command line: "${FUNCNAME} $@" >&2
        echo "Converting line ${_which_line} in file ${_summary} from a bitset to a GenePred file" >&2
    fi

    sed -n -e "${_which_line}"p ${_summary} |  awk -v chrom="${_chr}" -v listCoordinate="${_listCoordinate}" -v ll="${_which_line}" -f ${utils_awk} --source '
    {
       #first field in summary is a set of transcripts in bitset format
       #convert first field into GenePred format
       bitset2GenePred($1, listCoordinate, chrom, ll)

       #check if also second field is a set of transcripts in bitset format
       split($2, transcript, ",")
       #if so convert second field into GenePred format
       if (transcript[1] == "-")  
          bitset2GenePred($2, listCoordinate, chrom, ll)
    }'
}

#============== FUNCTION ===========================                                           
#         NAME:     similarity_aux                                                                  
#  DESCRIPTION:     reads in a (sorted) Bed file and returns the total length                  
#                      of all the elements in the Bed file                                     
#===================================================                                           
similarity_aux(){
    local var
    while read var;do
          echo $var
    done | awk '                                                                               
           BEGIN{                                                                              
           i=0                                                                                 
           }                                                                                   
           {                                                                                   
           i++;                                                                                
           starts[i]=$2;                                                                       
           ends[i]=$3;                                                                         
           }                                                                                   
           END{                                                                                
           lI=0;                                                                               
           for (j=1; j<=length(starts); j++)                                                   
               lI+=ends[j]-starts[j];                                                          
           print lI;                                                                           
           }'
}


#============== FUNCTION ===========================                                        
#         NAME:      similarity                
#  DESCRIPTION:      If Argument 4 is provided: compute similarity between transcript at line line1 in      
#                    file GenePred1 and transcript at line line2 in GenePred2              
#                    If Argument 4 is not provided: compute similarity between transcript at line line1 in
#                    file GenePred1 and transcripts in GenePred2
#
#Similarity between annotated transcript $T$ and discovered transcript $\tilde{T}$:
#At the base level
#S_n=\frac{n_{T\bigcap \tilde{T}}}{n_T}
#S_p=\frac{TN}{TN-FP}
#where $TN$ is the number of bases in the region minus number of bases in $T$.
#\frac{Number of splice sites in $T$ and $\tilde{T}$}{Number of splice sites in $T$}.
#\frac{Number of splice sites common to $T$ and the internal part of $\tilde{T}$}{Number of splice sites in $T$-2}
#Sort 4 3 1 2
#   ARGUMENT 1:      GenPred1 (string)                                                          
#   ARGUMENT 2:      line1 (integer)                   
#   ARGUMENT 3:      GenePred2 (string)                  
#   ARGUMENT 4:      (optional) line2 (integer)                                           
#===================================================
function similarity(){
    local _GenePred1=$1
    local _line1=$2
    local _GenePred2=$3
    if [[ ! -z $4 ]]; then
	local _line2=$4
        #extract appropriate lines from _GenePred1 and _GenePred2
	sed -n -e "$_line1"p $_GenePred1 | my_genePredLineToBed > $_GenePred1"_bed1"
	sed -n -e "$_line2"p $_GenePred2 | my_genePredLineToBed > $_GenePred2"_bed2"

	#at the base level                                                                 
	#_LI is length of the intersection of the Bed files $GenePred1"_bed1" and $GenePred2"_bed2"        
	local _LI=`multiIntersectBed -i $_GenePred1"_bed1" $_GenePred2"_bed2" | grep "1,2" | similarity_aux`
	#_LE and _lE2 are the lengths (in bases) of the intervals in Bed files $_GenePred1"_bed1" and $_GenePred2"_bed2" respectively                                                                           
	local _LE=`cat $_GenePred1"_bed1" | similarity_aux`
	local _LE2=`cat $_GenePred2"_bed2" | similarity_aux`
	#_LEN is the length (in bases) of the union of $_GenePred1"_bed1" and $_GenePred2"_bed2"       
	#local _LEN=`cat $_GenePred1"_bed1" $_GenePred2"_bed2" | regionLengthBed`        
	
	#at the splice site level                                                  
	#_starts1 is the list of 3' splice sites including the start                     
	#_internalStarts2 is the list of 3' splice sites (excluding the start)            
	local _starts1=`sed -n -e "$_line1"p $_GenePred1 | awk '{print $10}'`
	local _starts2=`sed -n -e "$_line2"p $_GenePred2 | awk '{print $10}'`
	local _ends1=`sed -n -e "$_line1"p $_GenePred1 | awk '{print $11}'`
	local _ends2=`sed -n -e "$_line2"p $_GenePred2 | awk '{print $11}'`
	local _starts1_v=(`echo ${_starts1//,/ }`)
	local _starts2_v=(`echo ${_starts2//,/ }`)
	#_E2 is the number of exons in $GenePred2"_bed2"                                  
	local _E1=${#_starts1_v[@]}
	local _E2=${#_starts2_v[@]}
	#_TP is the number of common splice sites btw _GenePred1 and _GenePred2                              
	local _TP=`numCommonSS $_starts1 $_starts2`
	let _TP=$_TP+`numCommonSS $_ends1 $_ends2`
	#_internalTp is the number of common start sites 
	local _internalStarts1
	local _internalEnds1
	local _internalTP=0
	if [ $_E1 -ne 1 ];then
            _internalStarts1=`ListRemoveFirst $_starts1`
            _internalEnds1=`ListRemoveLast $_ends1`
            _internalTP=`numCommonSS $_starts2 $_internalStarts1`
            let _internalTP=$_internalTP+`numCommonSS $_ends2 $_internalEnds1`
	fi
	local _li
	local _le2
	local _le
	local _tp
	echo " " | awk -v _li="$_LI" -v _le="$_LE" -v _le2="$_LE2" -v _tp="$_TP" \
            -v _e2="$_E2" -v _internaltp="$_internalTP" '                                        
           {                                                              
           # TruePositive/(True Positive + FalseNegative) at the base level    
           BLSn=_li/_le2                                                                   
           # TrueNegative/(TruePositive + FalsePositive) or in the literature TruePositive/(TruePositive+FalsePositive)                                                                                        
           BLSp=_li/_le                                                                                 
           ELSn=_tp/(2*_e2)                                                                         
           elsn=_internaltp/(2*_e2-2);                                                               
           printf("%s %s %s %s",BLSn,BLSp,ELSn,elsn)                                
           }'
	rm $_GenePred1"_bed1"
	rm $_GenePred2"_bed2"
    else
	local _ll=`cat $_GenePred2 | wc -l `
	local _ff
	local _s
	for (( _ff=1; _ff<=$_ll; _ff++ ));do
	    _s=`similarity $_GenePred1 $_line1 $_GenePred2 $_ff`
	    echo $_ff $_s
	done | sort -k 4,4n -k 3,3n -k 1,1n -k 2,2n | tail -1
    fi
}

#============== FUNCTION ===========================                                   
#         NAME:      printSimilarity                                                   
#  DESCRIPTION:      print to sdout similarity between transcripts in GenePredOut              
#                    and transcripts in GenePredIn for each transcript in GenePredOut                       
#   ARGUMENT 1:      GenPredOut                                                            
#   ARGUMENT 2:      GenePredIn (the reference genePred file)                                  
#        USAGE:      printSimilarity GenepredOut GenePredIn > similarityFile
#===================================================                                           
function printSimilarity(){
    local _GenePredOut=$1
    local _GenePredIn=$2
    local _l=`cat $_GenePredOut | wc -l `
    local _ii
    for (( _ii=1; _ii<=$_l; _ii++ ));do
      similarity $_GenePredOut $_ii $_GenePredIn
    done
}

function version () {
    msg="$0 1.0\n"
    msg+="\n"
    msg+="Copyright (C) 2012 E. Pantaleo\n"
    msg+="License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
    msg+="This is free software; see the source for copying conditions.  There is NO\n"
    msg+="warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
    msg+="\n"
    msg+="Written by E. Pantaleo.\n"
    echo -e "$msg"
}

function help () {
    msg="altra: a Bayesian method for simultaneous transcript reconstruction and abundance estimation with RNA-Seq data in multiple samples.\n"
    msg+="\n"
    msg+="Usage: \n"
    msg+="`basename $0` -L <chr:locusStart-locusEnd> -r <bamfile1[,...,bamfileN]> -c <count1[,...,countN]> -a <int> -O <int>  -g <string> \n"
    msg+="Example: \n"
    msg+="`basename $0` -L chr13:19875805-19998012 -r ind1.bam,ind2.bam,ind3.bam,ind4.bam -c 7749699,11630601,10614370,35456718 -d 0 -y 5 -a 46 -O 9 -g annotations/knownGene_hg18.txt.gz,annotations/ensembl_4_30.gz \n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\t\tDisplay the help and exit.\n"
    msg+="  -V, --version\t\tOutput version information and exit.\n"
    msg+="  -v, --verbose\n"
    msg+="  -L, --locus\t\t<chr:locusStart-locusEnd>\tA string specifying the region: chromosome:locusStart-locusEnd.\n"
    msg+="  -o, --out\t\t<string>\t\t\tSet the path to the output folder. Default output folder is ./out_altra/\n"
    msg+="  -c, --totalnreads\t\tA comma separated string of integers specifying the total number of reads for sample.\n"
    msg+="  -r, --readfiles\t\tComma separated string of path to the bam files (in the future would like to replace this with a samplesheet with 3 columns specifying: sample ID (column 1), path to bam file containing reads for sample (column 2), total number of reads for sample (column 3); this file has N rows, one row for each sample. From the bam files altra will extract reads that start in the locus (more precisely that start between the start of the locus and the end of the locus minus the read length)).\n"
    msg+="  -a, --read_length\t<int>\t\t\t\tThe read length; all reads must have same length; altra will remove from the input GenePred file those transcripts that are shorter than the read length.\n"
    msg+="  -O, --overhang\t<int>\t\t\t\tThe overhang of the mapping method. Default value is 9.\n"
    
    msg+="\nOPTIONS ON THE TRANSCRIPT MODEL\n"
    msg+="By default altra simultaneously reconstructs the transcript model and estimates abundances (default option -G 1). In this case a (list) of (compressed) GenePred annotation files must be provided with option -g (see below for more details). Also in this case the number of transcripts in the transcript model must be provided using options -d and -y followed by the number of transcripts on the positive and negative strand, pK and nK, respectively: altra will (randomly) select pK and nK transcripts from the list of annotated transcripts in the locus.\nIf option -G 2 is set, altra will only estimate abundances for a fixed transcript model. In this case a GenePred file must be specified with option -f; option -F can be used to specify which lines in the provided GenePred file are to be used (if option -F is not provided, altra will use all the transcripts in the GenePred file).\n"
    msg+="  -G, --genepred_initialization\t<default=1/2>\t\tReconstruct the transcript model(1); assume the transcript model is known and only estimate abundances(2).\n"
    msg+="  -f, --genepred_file\t\t<string>\t\tThis option should only be specified when -G 2, i.e., when the gene model is known and is therefore not reconstructed by altra. GenePred file containing the input transcript model; note that altra will cut the annotated transcripts if they extend beyond the locus (specified by option -L).\n"
    msg+="  -g, --genepred_annotation\t<file1[...,fileG]\tThis option should only be specified when -G 1, i.e., when altra the gene model is unknown. A comma separated list of annotation files in (compressed gz) GenePred format; from these genome-wide annotation files altra will extract transcripts that (partially) overlap with the locus and will remove identical transcripts in multiple annotations. altra will truncate the annotated transcripts if they extend beyond the locus (the locus is specified by option -L).\n"
    msg+="  -s, --minimum_exon_length\t<int>\t\t\tMinimum exon length: for each transcript in the GenePred file, if the lenght of an exon is less than the minimum exon length, altra will remove that exon. Default value is 10. The minimum exon length cannot be smaller than the overhang provided with option -O. \n"
    msg+="  -S, --maximum_exon_length\t<int>\t\t\tMaximum exon length. Default value is 3000. If the annotated GenePred file contains an exon that is longer than the minimum exon length provided with option -s, then the value of this parameter will be progressively increased until it is larger than the annotated exon.\n"
    msg+="  -z, --minimum_intron_length\t<int>\t\t\tMinimum intron length: for each intron in the GenePred file, if the intron length is less than the minimum intron lenght (provided with option -z) altra will merge adjacent exons in the GenePred file. Default value is 50.\n"
    msg+="  -Z, --maximum_intron_length\t<int>\t\t\tMaximum intron length. Default value is 50000\n"
    msg+="  -d, --positiveK\t\t<int>\t\t\tNumber of positive transcripts pK; if no transcripts are present on the positive strand in the GenePred files specified with option -g altra will set the number of positive transcripts to 0. Default value is 3.\n"
    msg+="  -y, --negativeK\t\t<int>\t\t\tNumber of negative transcripts nK; if no transcripts are present on the negative strand in the GenePred files specified with option -g altra will set the number of negative transcripts to 0. Default value is 3.\n"
        
    msg+="\nOPTIONS ON SPLICE SITES\n"
    msg+="altra needs a list of positive splice sites to reconstruct transcripts on the positive strand and a list of negative splice sites to reconstruct transcripts on the negative strand (in altra a start and an end of a transcript are equivalent to a 3' and a 5' splice site, respectively).\n"
    msg+="altra will extract junctions from the bam files. Also, altra will extract annotated junctions from the annotations (GenePred files): if an annotated splice sites is less than 3 bases away from a splice site found in the data, altra will ignore it. Furthermore, altra will discover junctions using FLLat (option -e 0 will disable the use of FLLat). By using FLLat, altra will find potential 5' and 3' splice sites, and will save them both as positive and as negative 3' or 5' splice sits.\n"
    msg+="If pK or nK is 0 and after filtering altra finds junctions on the positive (negative) strand then altra will print a warning.\n"
    msg+="JUNCTIONS FILTERING\n"
    msg+="From the set of junctions found in spliced reads (in the bam files), altra will only take those junctions with relative expression greater than PERCENT_J (provided with option -J) (over all individuals) and that are supported by at least _MINIMUM_J (provided with option -M)reads in at least one individual.\n"
    msg+="From the set of  5' and 3' splice sites (on the positive strand and on the negative strand) extracted from the filtered GenePred file altra will remove splice sites that are less than 3 nucleotides away from a 3' (5') splice site (on the same strand) extracted from the mapped reads."
    msg+="The position of a potential splice site discovered by FLLat has a large error; for this reason, if a 3' (5') splice site discovered by FLLat is within less than DELTA_J (specified by option -D) nucleotides from a 3' (5') splice site extracted from the reads and/or annotated, altra will ignore it. Also if the distance between two adjacent 5' and 3' splice sites found by FLLat is less than the minimum exon length both splice sites are discarded.\n"
    msg+="\n"
    msg+="  -j, --scale_junction_count_by\t<int>\t\t\tMultiply the count of spliced reads by scale_junction_count_by; this factor is only used when reconstructing the transcript model and is set to 1 when abundances are inferred. Default value is 1.\n"
    msg+="  -J, --filter_junctions\t<double>\t\tFilter out junctions that have a frequency lower than PERCENT_J in all individuals (if the total number of junction reads -- summed over all N samples -- that map to a specific junction in a specific sample is less than a percentage PERCENTAGE_J of the total number of junction reads that map to all junctions in the specified locus then ignore that specific junction). Default value is 2.\n"
    msg+="  -M, --minimum_junctions\t<int>\t\t\tFilter out junctions that have less than MINIMUM_J counts in each individual for all individuals. Default value is 2.\n"
    msg+="  -e, --find_breakpoints\t<0/default=1>\t\t\tUse FLLat to find breakpoints (junctions) that are common to all samples (1); disable the use of FLLat- don't look for breakpoints (0).\n"
    msg+="  -D, --distance_of_redundant_breakpoint\t<int>\tIgnore a splice sites discovered by FLLat if it is less than DELTA_J nucleotides away from an annotated splice site. Default value is 12.\n"
        
    msg+="\nINTERNAL OPTIONS\n"
    msg+="  -p, --MC_STEPS\t\t<int>\t\t\tSet the number of MCMC steps.\n"
    msg+="  -n, --MC_BURNIN\t\t<int>\t\t\tSet the number of steps to be thrown as burn in.\n"
    msg+="  -t, --MC_THIN\t\t\t<int>\t\t\tSet the thinning of the chain.\n"
    msg+="  -q, --MC_EQ\t\t\t<int>\t\t\tSet the number of steps for final MCMC.\n"
    msg+="  -l, --list_loci_line\t\t<int>\t\n"
    msg+="  -b, --list_loci_file\t\t<string>\t\n"
    
    msg+="\n"
    msg+="\nINPUT AND OUTPUT FORMATS\n"
    msg+="GenePred format: see http://genome.ucsc.edu/FAQ/FAQformat#format9 for a description of the GenePred format\n"
    msg+="BED format: see http://genome.ucsc.edu/FAQ/FAQformat.html#format1\n"
    msg+="BAM format: see http://genome.ucsc.edu/FAQ/FAQformat.html#format5.1 Note: reads that map with junctions must have the XS attribute tag as in Tophat output.\n"
    msg+="\n"
    msg+="summary.txt This file contains the list of explored gene models. Specifically, each line refers to a gene model; the first field is the gene model in  a format internal to altra, the second field is the log-likelihood of the gene model, the third field is the posterior frequency of the gene model. Function bitset2GenePred in utils.sh can be used to convert a line of summary.txt into a GenePred file.\n"
    msg+="ExprOut This file contains a N by 2*(pK+nK+1) table with estimated posterior mean and standard deviation for lambda and epsilon. Even columns contain the posterior mean for each transcript, odd colums contain the posterior standard deviation for each transcipt, the first 2*(pK+nK) columns report lambda, the last two columns report epsilon. Different samples are reported on different lines.\n"    
    echo -e "$msg" 
}

function parseArgs () {
    TEMP=`getopt -o hVvo:r:c:a:O:f:F:g:s:S:z:Z:G:L:j:J:M:e:D:p:n:t:q:m:d:y:l:b: -l help,version,verbose: \
	-n "$0" -- "$@"`
    if [ $? != 0 ] ; then echo "ERROR: getopt failed" >&2 ; exit 1 ; fi
    eval set -- "$TEMP"
    while true; do
	case "$1" in
	    -h|--help) help; exit 0; shift;;
	    -V|--version) version; exit 0; shift;;
	    -v|--verbose) VERBOSE=1; shift;;
	    -o|--out) JOB_FOLDER=$2; shift 2;;
	    -r|--readfiles) readfiles=$2; shift 2;;
	    -c|--total_counts) C=$2; shift 2;;
            -a|--read_length) RL=$2; shift 2;;
	    -O|--overhang) OVERHANG=$2; shift 2;;
	    -f|--genepred_file) GenePredIn=$2; shift 2;;
	    -g|--genepred_annotation) GenePredReference=$2; shift 2;;
	    -G|--genepred_initialization) genePrediction=$2; shift 2;;
	    -s|--minimum_exon_length) MIN_EX_LENGTH=$2; shift 2;;
	    -S|--maximum_exon_length) MAX_EX_LENGTH=$2; shift 2;;
	    -z|--minimum_intron_length) MIN_IN_LENGTH=$2; shift 2;;
	    -Z|--maximum_intron_length) MAX_IN_LENGTH=$2; shift 2;;
	    -m|--both_strands) both_strands=$2; shift 2;;
            -d|--positiveK) pK=$2; shift 2;;
            -y|--negativeK) nK=$2; shift 2;;
	    -L|--region) LOCUS=$2; shift 2;;
	    -j|--scale_junction_count_by) scale_junction_count_by=$2; shift 2;;
	    -J|--filter_junctions) PERCENT_J=$2; shift 2;;
	    -M|--minimum_junctions) MINIMUM_J=$2; shift 2;;
	    -e|--find_breakpoints) breakpoints=$2; shift 2;;
	    -D|--distance_of_redundant_breakpoint) DELTA_J=$2; shift 2;;
	    -p|--MC_STEPS) MC_STEPS=$2; shift 2;;
	    -n|--MC_BURNIN) MC_BURNIN=$2; shift 2;;
	    -t|--MC_THIN) MC_THIN=$2; shift 2;;
	    -q|--MC_EQ) MC_EQ=$2; shift 2;;
	    -l|--list_loci_line) list_loci_line=$2; shift 2;;
	    -b|--list_loci_file) list_loci_file=$2; shift 2;;
	    --) shift; break;;
	    *) echo 1=$1 "ERROR: options parsing failed"; exit 1;;
	esac 
    done
    if [[ -z $LOCUS ]];then
	echo "ERROR: specify the LOCUS" >&2 
	help >&2 
	exit 1
    fi
    
    if [[ -z $RL ]];then
	echo "ERROR: specify the read length" >&2
	help >&2
	exit 1
    fi
    
    if [[ -z $DELTA_J ]]; then
	echo "ERROR: specify the distance DELTA_J (if a splice site discovered by FLLat is within less than DELTA_J nucleotides from an annotated splice site ignore it)" >&2
	help >&2
	exit 1
    fi
        
    if [[ ! -d $JOB_FOLDER ]]; then
	mkdir $JOB_FOLDER
    fi

    last=`echo $JOB_FOLDER | tail -c 2`
    if [[ "$last" != "/" ]]; then
	JOB_FOLDER=$JOB_FOLDER"/"
    fi
    echo "Writing output to $JOB_FOLDER" >&2
    
    if [[ -z $readfiles ]]; then
	echo "ERROR: specify path to BAM files using option -r" >&2
	help >&2
	exit 1
    fi

    if [[ -z $C ]]; then
	echo "ERROR: set value of C (sequencing depth) using option -c"
        help >&2
        exit 1
    fi

    if [ $genePrediction -eq 1 ]; then
	if [[ -z $GenePredReference ]]; then
            echo "ERROR: specify GenePred annotation files using option -g" >&2
	    exit 1
	fi
	if [[ ! -z $GenePredIn ]]; then
	    echo "ERROR: option -g is incompatible with (default) option -G 1" >&2
	    exit 1
	fi
	if [[ -z $pK || -z $nK ]]; then
	   echo "ERROR: specify the number of transcripts using options -d and -y" >&2
	   exit 1
	fi
    fi
 
    if [ $genePrediction -eq 2 ]; then
	if [[ -z $GenePredIn ]]; then
            echo "ERROR: option -G is set to 2; you must provide a GenePred file with option -f" >&2
            exit 1
	fi
    fi
}
