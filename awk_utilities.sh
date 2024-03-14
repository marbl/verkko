#!/usr/bin/awk
# Universitat Potsdam
# Date: 2024-3-11
# for preparing the data for the visualization of the coverage or the length of the assembled unitigs from the pacbiohifi assembly. 
gem install youplot
# coverage 
for i in $(ls -la *.csv); \ 
        do cat $i | awk '{ print $2 }'; done
# length
for i in $(ls -la *.csv); \
        do cat $i | awk '{ print $2 }'; done
# filtering specific to the coverage values ( storing the coverage values in a hash value and then implementing the awk over the same)
coverage="value"
for i in $(for i in $(ls -la *.csv); \
      do cat $i | awk '{ print $2 }'; done); \ 
                do if [[ $i -ge "${coverage}" ]]; then echo $i; fi; done
length="value"
for i in $(for i in $(ls -la *.csv); \
      do cat $i | awk '{ print $3 }'; done); \ 
                do if [[ $i -ge "${length}" ]]; then echo $i; fi; done

# normalizing the coverage according to the length
covergaefile="covergae"
for i in $(cat "${coverage}" | awk '{ print $2 }'); \ 
          do for j in $(cat "${coverage}" | awk ' { print $3 }'); \
                                        do expr ${i} * ${j}; done; done
#plotting the length right in the terminal after filtering out the short unitigs
# binning them according to the length filter and then making the sense of the assembled unitigs
lengthselectionsort="variable"
for i in $(cat test.cov | awk '{ print $3 }'); \
                do if [[ $i -ge "${lengthselectionsort}" ]] then; \ 
                                        echo $i; fi; done | youplot barplot

 for i in $(cat test.cov | awk '{ print $3 }'); \
                do if [[ $i -ge "${lengthselectionsort}" ]] then; \ 
                                        echo $i; fi; done | youplot histogram
