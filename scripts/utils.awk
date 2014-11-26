
# \file utils.awk
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

#####################################################################                              
#       definition of function join                                                
#       x[1]=1; x[2]=2; x[3]=3; join(x,1,3,",")                                               
#       returns 1,2,3,                                                               
#####################################################################                   
function join(array, start, end, sep){
    if (sep == "")
       sep = " "
    else if (sep == SUBSEP) # magic value                                                      
       sep = ""
    result = array[start] sep
    for (i = start + 1; i <= end; i++)
        result = result array[i] sep
    return result
}


#####################################################################                          
#       definition of function sam_start2vector_intervals                                      
#       example: start=$4 CIGAR=$6 in sam file                                                 
#                start=100 CIGAR=10M20N30M                                                     
#                sam_start2vector_intervals(start,CIGAR,read)                                  
#                read=(100,109,130,159)                                                        
#####################################################################                          
function sam_start_cigar2vector_intervals(CIGAR, start, read){                                   
        #split CIGAR                                                                           
        gsub("M", " ", CIGAR); gsub("N"," ",CIGAR);l=split(CIGAR,CIGAR_v," ");                    
                                                                                               
        delete read                                                                            
        read[1] = start                                                                        
        for (i = 1; i <= l; i ++){                                                             
	  if ((i % 2) == 0)                                                                  
	    const = -1;                                                                    
	  else                                                                               
	    const = 1;                                                                     
	  read[i + 1] = read[i] + CIGAR_v[i] - const;                                        
        }
}                                                                                              
#####################################################################   

#####################################################################                       
#       definition of function bitset2GenePred                                                
#####################################################################                          
function bitset2GenePred(bitsets, listCoordinate, chrom, ll){                              
  T = split(bitsets,transcript, ",") - 1;                                                    
  C = split(listCoordinate, l , ",");                                                          
  for (k = 2; k <= T; k ++){                                                                     
    counter = 0;                                                                                     
    for (c = 1; c < C; c ++){                                                                    
      if (substr(transcript[k], c, 1) == "0")                                        
	counter ++;                                                                     
    }                                                                                            
    if (counter == C - 1)                                                                             
      continue;                                                                             
    counter = 0;                                                           
    delete list;                                                                                   
    printf("altra.%s.%s.%s.%s\t%s\t%s\t", chrom, l[1], l[C], ll, chrom, transcript[1]);                                             
    if (substr(transcript[k], 1, 1) == "1"){                                              
      counter ++;                                                                          
      list[counter] = l[1] - 1 + 0.5;                                                        
    }                                                                                             
    for (c = 1; c < C - 1; c ++){                                                                 
      if (substr(transcript[k], c, 1) == "1" && substr(transcript[k], c + 1, 1) == "0"){ 
	counter ++;                                                                  
	list[counter] = l[c + 1] - 0.5;                                         
      }                                                                               
      else if(substr(transcript[k], c, 1) == "0" && substr(transcript[k], c + 1, 1) == "1" ){        
	counter ++;                                                                           
	list[counter] = l[c + 1] - 1 + 0.5;                                                   
      }                                                                                             
    }                                                                                    
    if(substr(transcript[k], C - 1, 1) == "1"){                                   
      counter ++;                                                                                    
      list[counter] = l[C] - 0.5;                                                             
    }                                                                             
    ll = length(list);                                                                              
    printf("%s\t%s\t%s\t%s\t%s\t", list[1], list[ll], list[1], list[ll], ll/2);             
    for (i = 1; i <= ll; i ++){                                                              
      printf("%s,",list[i]);                                                                   
      i ++;                                                                                   
    }                                                                                      
    printf("\t");                                                                 
    for (i = 1; i < ll; i ++){                                                              
      i ++;                                                                                       
      printf("%s,", list[i]);                                                                      
    }	   
    printf("\n");
  }                                                                                 
}                                                                                                
############################################################################ 

#####################################################################              
#       definition of function sam2bitset 
##################################################################### 
function sam2bitset(lC_v, C){
    sam_start_cigar2vector_intervals($6, $4, read)                                    
    l = length(read)                   

    #if read has no junctions 
    #if ($7=="="){split($1, field1, "["); split(field1[2], read_number, "]"); line=sprintf("%s",read_number[1]);for (i=2;i<7;i++) line=sprintf("%s\t%s",line, $i);line=sprintf("%s\t%s",line, "*");for (i=8;i<=NF;i++) line=sprintf("%s\t%s",line, $i); if ($4<$8) print line | "sort -k 1n,1n" > file1; else print line | "sort -k 1n,1n" > file2}}
    if (l == 2){                                                          
	for (i = 2; i <= C; i++){                                                                                                                 
	    if (read[1] < lC_v[i]){                                                 
		toprint = i - 2;                                                                                                                  
		for (j = i; j <= C; j ++){                                                                         
		    if (read[2] <= lC_v[j]){                                                                                    
			toprint = toprint","j - 1",";                                         
			break;                                                                                       
		    }                                                                                                     
		}                                                                                                                             
		break;                                                                                                                         
	    }                                                                                                                             
	}                                                    
        #if ($7 == "=") {print "="; print $1}                                                                      
	print "n";                                                                                                               
	print toprint;                                                                            
    }else{                                                                                                                                                
        #if read has junctions  
	for (i = 1; i <= C; i ++){                                                                                       
	    if (read[1] < lC_v[i]){                                                                              
                        toprint = i - 2","                                       
                        break;                       
	    }                                                                                         
	}                                                                                                               
	counter = 0;         
	for (ll = 2; ll < l - 1; ll += 2){                                                                                  
	    for (j = i; j <= C; j ++){                                                                                 
		if (read[ll] + 0.5 == lC_v[j]){                                                                 
		    counter ++;                                                                                                         
                    toprint = toprint j - 1","       
                    for (jj = j + 1; jj <= C; jj ++){                                                            
                        if (read[ll + 1] - 0.5 == lC_v[jj]){                                                 
                            toprint = toprint jj - 1","                                                            
                            counter ++;                                                                                      
                            break;                                                                          
                        }                                                                                                    
                    }                                                                                       
                    break;                                                                                                            
		}                                                                                                             
	    }                     
	}                                                                                                                                                                  
	if (counter != l - 2){                                                                      
            #if ($7 == "=") {print "="; print $1;}                                                       
	    print "z";                                                                                           
	}else{                                                             
	    for (jjj = jj; jjj <= C; jjj ++){                                                           
		if (read[l] < lC_v[jjj])                                                                                
		    break;                                                                       
	    }                                                                                       
	    if (jjj <= C)                                                                                                     
		toprint = toprint jjj - 1",";                                                                         
	    else                                                                                                                              
		toprint = toprint C - 1",";                                                                         
                      #if ($7 == "="){                                          
                      #   print "=";                                                                                            
                      #   print $1;                                                                               
                      #}                                                                                            
                      # if a field in the sam file contains the string $field=="XS:A:+" || $field=="XS:A:-"                             
                      # than the strand= substr($field, 6, 1)= "+" or "-" respectively                                                  
	    strand_is_present=0;                                                                                              
	    for (field=12; field<=20; field++)                                                                                 
		if (substr($field, 1, 4) == "XS:A"){                                                                         
		    print substr($field, 6, 1);                                                              
		    strand_is_present=1;
		    break;
		}                                              
	    if (strand_is_present==0){                                                    
		print "strand_is_not_available_in_sam_file";                                  
		next;                                                                                        
	    }                                                                                       
	    print toprint;
	}                                                                                                       
    }                                                                        
}
###############################################################################################
