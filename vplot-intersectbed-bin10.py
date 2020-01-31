#!/usr/bin/env python
#Usage: python vplot-intersectbed.py <fragment-file-name> <peak-file-name> <v-plot-file-name> <fraglength-file-name> <peak-size:m bp> <peak-bin: n bp> <fraglength-size: 0-j bp>
#   <fragment-file> = <chr,start,end, len>
#		Typically a file of paired illumina reads
#		strand is not used
#	<peak-file> = <chr,start,end,id,strand,...>, only the first 5 fileds were used
#		peak summit file
#       strand is used
#	To plot peak-size/2 upstream and downsteam of peak summits, the peak summits are firstly sloped by peak-size/2 on both direction, then find the 1 bp overlapped  fragments
#   
#	Computes distances from midpoints of data_file records to
#		point_file records for point_file records within
#		intervals around point_file midpoint +- pk_size

from __future__ import print_function
import re,sys,os,numpy

if len(sys.argv) != 8:
	print("python vplot-intersectbed.py <fragment-file-name> <peak-file-name> <v-plot-file-name> <fraglength-file-name> <peak-size:m bp> <peak-bin: n bp> <fraglength-size: 0-j bp>")
	sys.exit(1)

fm_num=int(os.popen('wc -l %s' % sys.argv[1]).read().split()[0])
pk_num=int(os.popen('wc -l %s' % sys.argv[2]).read().split()[0])
intFileCol1=len(os.popen('head -1 %s' % sys.argv[1]).read().split('\t')) #4
intFileCol2=len(os.popen('head -1 %s' % sys.argv[2]).read().split('\t')) 
pk_size = int(sys.argv[5])
sloplist=os.popen("slopBed -i %s -g /data1/zqwen/genome/mm9.chrom.sizes -b %s" % (sys.argv[2], pk_size/2)).read().split('\n') ##extend the peak summit to pk_size bp
vfile = open(sys.argv[3],"wa")
frac = open(sys.argv[4],"wa")
pk_bin = int(sys.argv[6])
fm_size = int(sys.argv[7])

uniqueIndex=str(os.getpid())
sloppeak='sloppeak'+uniqueIndex
fp=open(sloppeak,'w')
for peak in sloplist:
	if peak!='':
		fp.writelines(peak+'\n')
	else:
		pass
fp.close()
interInfo=os.popen('intersectBed -a %s -b %s -wo' % (sys.argv[1],sloppeak)).read().split('\n')
os.system('rm -rf %s' % sloppeak)

midlen='tmp'+uniqueIndex
tmp = open(midlen, "wa") #
for line in interInfo:
	if line!='':	
		(fm_chr, fm_start, fm_end, fm_len, pk_chr, pk_start, pk_end, pk_id, pk_strand) = line.split('\t')[0:9]
		pk_mid = int(pk_start) + 1 + pk_size/2
		fm_mid = int(fm_start) + 1 + int(fm_len)/2
		dist = fm_mid - pk_mid
		if abs(dist) <= pk_size/2:
			if pk_strand == '0' or pk_strand == '+': #'0' stands for '+' strand in homer
				print(pk_id, fm_len, dist, sep='\t', file=tmp) #
			elif pk_strand == '1' or pk_strand == '-':
				print(pk_id, fm_len, -dist, sep='\t', file=tmp) #
			else:
				print("Not correct strand information")
tmp.close()
tmp = open(midlen, "r")

#fm_size = 200
#pk_size = 1001,201 
count = numpy.zeros((fm_size, pk_size+1), dtype = 'int') #(fraglength, distance)
for line in tmp:
	#print(line, end="")
	lsplit = line.strip().split()
	i = int(lsplit[1]) - 1 #project fragmeng length to row of matrix
	j = int(lsplit[2]) + pk_size/2  #project distance to colum of matrix 
	if i < 0 or i > 199 :
		continue
	count[i][j] += 1

#print V-plot title
print('fraglength', end='\t', file=vfile)
for j in range(-pk_size/2, pk_size/2+1):
	if j != pk_size/2:
		print(j, end='\t', file=vfile)
	else:
		print(j, end='\r\n', file=vfile)
		
#print density file for V-plot	
i = fm_size
while i > 0 :
	print(i, end='\t', file=vfile)
	for j in range(pk_size+1):
		cpg = round(1000000000.0*count[i-1][j]/(fm_num*pk_num), 8) #cpm: centers per million fragments per peak
		if j != pk_size:
			print(100*cpg, end='\t', file=vfile)  #print element of a row separated by tab
		else:
			print(100*cpg, end='\r\n', file=vfile)
	i -= 1

xbins = pk_size/pk_bin #1000/10=100 #200/2

#print fractioned ratio file for meta-plot
count=numpy.transpose(count) #array(pk_size+1,fm_size)(distance, fraglength)
print('distance','80 bp','100 bp','120 bp','140 bp','168 bp',sep='\t', end='\r\n', file=frac)
for i in range(xbins-1):
	print(pk_bin*(i - xbins/2 + 1), end='\t', file=frac)
	len_sum = 0
	for k in range(pk_bin*i, pk_bin*(i + 2)): 
		for j in range(81):
			len_sum = len_sum + count[k][j]
	cpg = round(1000000000.0*len_sum/(2*fm_num*pk_num), 8) #via (i+2) and (2*fm_num*pk_num), cpg is smoothed by the next bin 
	print(cpg, end='\t', file=frac)  #print element of a row separated by tab
	len_sum = 0
	for k in range(pk_bin*i, pk_bin*(i + 2)): 
		for j in range(81,101):
			len_sum = len_sum + count[k][j]
	cpg = round(1000000000.0*len_sum/(2*fm_num*pk_num), 8)
	print(cpg, end='\t', file=frac)
	len_sum = 0
	for k in range(pk_bin*i, pk_bin*(i + 2)): 
		for j in range(101,121):
			len_sum = len_sum + count[k][j]
	cpg = round(1000000000.0*len_sum/(2*fm_num*pk_num), 8)
	print(cpg, end='\t', file=frac)
	len_sum = 0
	for k in range(pk_bin*i, pk_bin*(i + 2)):
		for j in range(121,141):
			len_sum = len_sum + count[k][j]
	cpg = round(1000000000.0*len_sum/(2*fm_num*pk_num), 8)
	print(cpg, end='\t', file=frac)
	len_sum = 0
	for k in range(pk_bin*i, pk_bin*(i + 2)):
		for j in range(141,169):
			len_sum = len_sum + count[k][j]
	cpg = round(1000000000.0*len_sum/(2*fm_num*pk_num), 8)
	print(cpg, end='\r\n', file=frac)

#print the last line of ratio file, as the previous lines are smoothed
i = xbins-1
print(pk_bin*(i - xbins/2 + 1), end='\t', file=frac)
len_sum = 0
for k in range(pk_bin*i, pk_bin*(i + 1)):
	for j in range(81):
		len_sum = len_sum + count[k][j]
cpg = round(1000000000.0*len_sum/(fm_num*pk_num), 8) 
print(cpg, end='\t', file=frac)  #print element of a row separated by tab
len_sum = 0
for k in range(pk_bin*i, pk_bin*(i + 1)):
	for j in range(81,101):
		len_sum = len_sum + count[k][j]
cpg = round(1000000000.0*len_sum/(fm_num*pk_num), 8)
print(cpg, end='\t', file=frac)
len_sum = 0
for k in range(pk_bin*i, pk_bin*(i + 1)):
	for j in range(101,121):
		len_sum = len_sum + count[k][j]
cpg = round(1000000000.0*len_sum/(fm_num*pk_num), 8)
print(cpg, end='\t', file=frac)
len_sum = 0
for k in range(pk_bin*i, pk_bin*(i + 1)):
	for j in range(121,141):
		len_sum = len_sum + count[k][j]
cpg = round(1000000000.0*len_sum/(fm_num*pk_num), 8)
print(cpg, end='\t', file=frac)
len_sum = 0
for k in range(pk_bin*i, pk_bin*(i + 1)):
	for j in range(141,169):
		len_sum = len_sum + count[k][j]
cpg = round(1000000000.0*len_sum/(fm_num*pk_num), 8)
print(cpg, end='\r\n', file=frac)	
	
tmp.close()
os.system('rm -rf %s' % midlen)
vfile.close()
frac.close()