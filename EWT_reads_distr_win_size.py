# -*- coding: utf-8 -*-
"""
Created on Tue Nov 08 16:38:08 2016

@author: lenovo
"""

#==============================================================================
# this is for getting reads distribution according to window size
#==============================================================================

from __future__ import division
#import numpy as np
#from scipy import stats
#import math
import os,glob,sys,re
from optparse import OptionParser
import subprocess as sp
from multiprocessing import Pool

parser=OptionParser()

parser.add_option(
	'-D',
	'--sorted-dir',
	dest='sortedDir',
	help='dir to sorted bam dir'
	)
	
parser.add_option(
	'-R',
	'--result-dir',
	dest='resDir',
	help='dir where intermediate files are stored'
	)

parser.add_option(
	'-P',
	'--cpu-nbr',
	dest='cpuNbr',
	help='cpu nbrs, default set to 1',
	default='1'
	)
	
(options,args)=parser.parse_args()

bedFile='/results/plugins/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed'

if not options.sortedDir or not options.resDir or not os.path.exists(bedFile):
	parser.print_help()
	sys.exit(1)

sortedDir=options.sortedDir
resDir=options.resDir
cpuNbr=int(options.cpuNbr)

def main():
	global sortedDir,resDir,cpuNbr
	winSizeLst=[10,20,30,40,50,60,70,80,90,100]
#	call mainExc
	
	headPat=re.compile(r'.*?/Home/(.*?)/plugin_out')
	fileHead=headPat.findall(sortedDir)[-1]
	
	bedDct=getBedInfo()
	os.chdir(sortedDir)
	sortedBams=glob.glob('*sorted.bam')
	
	pool=Pool(cpuNbr)
	for subBam in sortedBams:
		pool.apply_async(getSams,args=(sortedDir,subBam,bedDct,resDir))
	pool.close()
	pool.join()
	
	for winSize in winSizeLst:
		pool=Pool(cpuNbr)
		for sortedBam in sortedBams:
			pool.apply_async(calReads,args=(sortedDir,winSize,sortedBam,bedDct,resDir))
		pool.close()
		pool.join()
		
#		calReads(sortedDir,winSize,sortedBams[0],bedDct,resDir)
		
		os.chdir(resDir)
		pat=str(winSize)+'.nbr'
		winFiles=glob.glob('*'+pat)
		
		sp.call('touch %s' % (fileHead+'_'+pat+'.univ'),shell=True)
		for subfile in winFiles:
			catCmd='cat %s >> %s' % (subfile,fileHead+'_'+pat+'.univ')
			sp.call(catCmd,shell=True)
			os.remove(subfile)
		print(str(winSize)+' results merged...')
		
	os.chdir(resDir)
	samFiles=glob.glob('*.sam')
	for samfile in samFiles:
		os.remove(samfile)
	return

def getBedInfo():
	global bedFile
	bedDctGbi={}
	
	infh=open(bedFile,'r')
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
#		note that annotated linear[3] have overlaps, thus using linear[1] as key
		bedDctGbi[int(linear[1])]=[int(linear[2]),linear[3]]
	infh.close()

#	print bedDctGbi
	
	return bedDctGbi

def getSams(sortedDirGs,sortedBamGs,bedDctGs,resDirGs):
	os.chdir(sortedDirGs)
	
	print(sortedBamGs+' samtools starts...')
	
	tmpSam=resDirGs+'/'+sortedBamGs+'.sam'
	
	posLst=sorted(bedDctGs.keys())
	
	sp.call('touch %s' % (tmpSam),shell=True)
	for subStt in posLst:
#		use q20 for filtration of reads quality
		subEnd=bedDctGs[subStt][0]
		samCmd='samtools view %s -q 20 chrX:%s-%s >> %s' % (sortedBamGs,str(subStt),str(subEnd),tmpSam)
#		sp.call('samtools view %s -q 20 chrX:%s-%s' % (sortedBamCr,str(subStt),str(subEnd)),shell=True)
		sp.call(samCmd,shell=True)
	return

def calReads(sortedDirCr,winSizeCr,sortedBamCr,bedDctCr,resDirCr):
#	counting reads
#	put sam in resDir
#	one sorted bam is linked to a
	print(sortedBamCr+' winSize:'+str(winSizeCr)+' calReads starts...')
#==============================================================================
# 	os.chdir(sortedDirCr)
# 	posLst=sorted(bedDctCr.keys())
# #	print posLst
# #	posLst=sorted(posLst,key=lambda x:x[0])
# 	tmpSam=resDirCr+'/'+sortedBamCr+'.sam'
# 	
# 	sp.call('touch %s' % (tmpSam),shell=True)
# 	for subStt in posLst:
# #		use q20 for filtration of reads quality
# 		subEnd=bedDctCr[subStt][0]
# 		samCmd='samtools view %s -q 20 chrX:%s-%s >> %s' % (sortedBamCr,str(subStt),str(subEnd),tmpSam)
# #		sp.call('samtools view %s -q 20 chrX:%s-%s' % (sortedBamCr,str(subStt),str(subEnd)),shell=True)
# 		sp.call(samCmd,shell=True)
# #	passing in reads count in different gaps
#==============================================================================
	os.chdir(resDirCr)
	tmpSam=resDirCr+'/'+sortedBamCr+'.sam'
	readsNbrLst=bedReads(tmpSam,bedDctCr,winSizeCr)
	nbrInfoFile=sortedBamCr+'.'+str(winSizeCr)+'.nbr'
	outfh=open(nbrInfoFile,'w')
	for subNbr in readsNbrLst:
		outfh.write(str(subNbr)+'\n')
	outfh.close

	return
	
def bedReads(tgtSamBr,bedDctBr,winSizeBr):
	infh=open(tgtSamBr,'r')
	sttPosLst=[]
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
		sttPosLst.append(int(linear[3]))
	nbrLst=[]
	for keyP in bedDctBr.keys():
		endP=bedDctBr[keyP][0]+winSizeBr
		candStt=keyP
		candEnd=candStt+winSizeBr
		
		while candEnd <= endP:
			subLst=[i for i in sttPosLst if candStt <= i < candEnd]
			nbrLst.append(len(subLst))
			candStt=candEnd
			candEnd += winSizeBr
			
	return nbrLst

if __name__=='__main__':
	main()