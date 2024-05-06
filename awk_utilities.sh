#!/usr/bin/awk
# Universitat Potsdam
# Date: 2024-3-11
# for preparing the data for the visualization of the coverage or the length of the assembled unitigs from the pacbiohifi assembly. 
# test.cov is the coverage file coming from the Verkko assembly.
gem install youplot
# coverage 
for i in $(ls *.cov); \ 
           do cat $i | awk '{ print $2 }'; done
# length
for i in $(ls *.cov); \
             do cat $i | awk '{ print $2 }'; done
# filtering specific to the coverage values ( storing the coverage values in a hash value and then implementing the awk over the same)
coverage="value"
for i in $(for i in $(ls *.csv); \
      do cat $i | awk '{ print $2 }'; done); \ 
                do if [[ $i -ge "${coverage}" ]]; then echo $i; fi; done
length="value"
for i in $(for i in $(ls *.csv); \
      do cat $i | awk '{ print $3 }'; done); \ 
                do if [[ $i -ge "${length}" ]]; then echo $i; fi; done

# normalizing the coverage according to the length
coveragefile="name"
for i in $(cat "${coveragefile}" | awk '{ print $2 }'); \
                do for j in $(cat "${coveragefile}" | awk '{ print $3 }'); \
                                        do echo $i"\t"$j; done | awk '{ print $1*$2 }'
        
#plotting the length before filtering out the short unitigs
lengthselectionsort="variable"
for i in $(cat test.cov | awk '{ print $3 }'); \
                do if [[ $i -ge "${lengthselectionsort}" ]] then; \ 
                                        echo $i; fi; done | youplot barplot
# binning them according to the length filter and then making the sense of the assembled unitigs
 for i in $(cat test.cov | awk '{ print $3 }'); \
                do if [[ $i -ge "${lengthselectionsort}" ]] then; \ 
                                        echo $i; fi; done | youplot histogram  
                                        
# genome assembled following length filter and the filtered uitigs greater than 10000
cat test.cov | awk '$3 > 10000 { print $3 }' | gawk '{ sum += $1 }; \
                      END { print sum }' && cat test.cov | \
                                            awk '$3 > 10000 { print  $1"\t"$2"\t"$3 }'
