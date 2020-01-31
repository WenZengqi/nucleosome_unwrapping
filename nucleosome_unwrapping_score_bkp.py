#!/usr/bin/env python
#Usage: python nucleosome_unwrapping_score.py <fragment-file> <peak-file> <peak-size:m bp> <nus_output_file>
#   <fragment-file> = <chr,start,end, len>
#		Typically a file of paired illumina reads
#		strand is not used
#	<peak-file> = <chr, start, end, pkid,...>, sorted peak summit file, only the first 4 fileds were used, python script /data1/zqwen/script/print_bed_add_pkname.py can be used to add pkid_i to each peak of the sorted peak list. pk_id format is pkid_1.
#       
#	To calculate nuncleosome unwrapping score (nus) within peak-size region of each peak summit, the peak summits are firstly sloped by peak-size/2 on both direction. The 1 bp overlapped fragments with the slopped peaks are first found by intersectBed, but only the fragments within 0~200 bp and with middle points within the slopped peaks are counted during downstream analyses. Within each slopped peak, the read number of fragments within length 0 ~ x bp is counted as a, within x ~ 200 bp as counted as b. Then we calculate the nuncleosome unwrapping score of the slopped peak as nus = (a-b)/(a+b). Thus, -1 < nus < 1, the closer nus nears 1, the more the nucleosomes within the peak region unwrap; the closer nus nears -1, the less the nucleosomes within the peak region unwrap. nus may be affected by break point (x), but theoretically not by peak size (m) or sequence deapth. The nus of the same peak size (m) and break point (x) from different genome regions should be comparable. 
#  This python script is used to caculate break point (50 =< x <= 160) dependent nus at 1 bp resolution.

from __future__ import print_function
import re,sys,os,numpy

if len(sys.argv) != 5:
	print("python nucleosome_unwrapping_score.py <fragment-file> <peak-file> <peak-size:m bp> <nus_output_file>")
	sys.exit(1)

fm_num=int(os.popen('wc -l %s' % sys.argv[1]).read().split()[0])
pk_num=int(os.popen('wc -l %s' % sys.argv[2]).read().split()[0])
intFileCol1=len(os.popen('head -1 %s' % sys.argv[1]).read().split('\t')) #4
intFileCol2=len(os.popen('head -1 %s' % sys.argv[2]).read().split('\t')) 
pk_size = int(sys.argv[3])
pklist=os.popen("slopBed -i %s -g /data1/zqwen/genome/mm9.chrom.sizes -b %s" % (sys.argv[2], pk_size/2)).read().split('\n') ##extend the peak summit to pk_size bp
#bkp = int(sys.argv[4]) #break point
#fm_size = int(sys.argv[4]) #fragment length range, we use fm_size = 200 bp in this analysis.
nus_file = open(sys.argv[4],"wa")

uniqueIndex=str(os.getpid())
tmppeak='tmppeak'+uniqueIndex
fp=open(tmppeak,'w')
for peak in pklist:
	if peak!='':
		fp.writelines(peak+'\n')
	else:
		pass
fp.close()
interInfo=os.popen('intersectBed -a %s -b %s -wo' % (sys.argv[1],tmppeak)).read().split('\n')
os.system('rm -rf %s' % tmppeak)

#matrix to store all the 0-200bp DNA fragments within each peak (row)
count = numpy.zeros((pk_num, 200), dtype = float) 

for line in interInfo:
	if line!='':
		(fm_chr, fm_start, fm_end, fm_len, pk_chr, pk_start, pk_end, pk_id) = line.split('\t')[0:8]
		pk_mid = int(pk_start) + 1 + pk_size/2
		fm_mid = int(fm_start) + 1 + int(fm_len)/2
		dist = fm_mid - pk_mid 
		i = int(pk_id.split('_')[1]) - 1 #pk_id format: pkid_1
		j = int(fm_len) - 1
		if abs(dist) > pk_size/2:
			continue
		if j >= 0 and j <= 199: #we use 0~200 bp fragment length
			count[i,j] += 1
		else:
			continue

#nus matrix with 50 =< bkp <= 160, nus_bkp[i][j] is nus of peak_i with bkp = j+50
nus_bkp = numpy.zeros((pk_num, 112), dtype = float) 

for i in range(0, pk_num):
	nus_bkp[i,0] = sum(count[i,]) # store the total fragment for peak_i
	for k in range(1,112):
		bkp = k + 49
		a = sum(count[i,0:bkp])
		b = sum(count[i,bkp:200])
		if (a+b) == 0:
			nus_bkp[i,k] = 100 # we use 100 to indicate "NA"
		else:
			nus_bkp[i,k] = round((a-b)/(a+b), 8)
	
#print nus_file title

print('pkid\ttotal_fragment', end='\t', file=nus_file)
for k in range(1, 112):
	if k != 111:
		print("bkp_"+str(k+49), end='\t', file=nus_file)
	else:
		print("bkp_"+str(k+49), end='\r\n', file=nus_file)
		
#print nus_file body	
for i in range(0, pk_num):
	pkid = "pkid_" + str(i+1)
	print(pkid, end='\t', file=nus_file)
	for k in range(0, 112):
		if k != 111:
			print(nus_bkp[i,k], end='\t', file=nus_file)
		else:
			print(nus_bkp[i,k], end='\r\n', file=nus_file)
			
nus_file.close()
