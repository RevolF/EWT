# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 16:24:45 2016

@author: lenovo
"""

#==============================================================================
# use q30 as filtration threshold
# 2016.11.11:
#	we removed the event significance testing using one-side Z-method
#		this is considerred unnecessary if all event has a tail Pv larger than threshold
#		we shall see the results
#	shit like logic
#
# 2016.11.12:
#	a quick filtration for large deletion could be added, using percentage of reads change, and merge with gaoSize
#==============================================================================

from __future__ import division
import os,sys,glob,re
from optparse import OptionParser
import subprocess as sp
from multiprocessing import Pool
import numpy as np
from scipy import stats
import math

parser=OptionParser()

parser.add_option(
	'-D',
	'--sorted-dir',
	dest='sortedDir',
	help='dir to sorted bam, /dmd_project is automatically added'
	)
	
parser.add_option(
	'-R',
	'--result-dir',
	dest='resDir',
	help='dir where all results will be stored'
	)
	
parser.add_option(
	'-B',
	'--bed-file',
	dest='bedFile',
	help='bed file used, default set to: /results/plugins/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed',
	default='/results/plugins/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed'
	)
	
parser.add_option(
	'-P',
	'--cpu-nbr',
	dest='cpuNbr',
	help='cpu nbrs, default set to 1',
	default='1'
	)
	
(options,args)=parser.parse_args()

bedFile=options.bedFile
if not options.sortedDir or not options.resDir or not os.path.exists(options.bedFile):
	parser.print_help()
	sys.exit(1)
	
if not os.path.exists(options.resDir):
	print('file path not exists: %s'	%	(options.resDir))
	sys.exit(1)
	
bedFile=options.bedFile
sortedDir=options.sortedDir+'/dmd_project'
resDir=options.resDir
cpuNbr=int(options.cpuNbr)

################
winSize=50
fpr=0.05
# windows size is sensitive with fpr, change this after calculation with simulation
################

def main():
	global bedFile,sortedDir,resDir,cpuNbr,winSize,fpr
	
	headPat=re.compile(r'.*?/Home/(.*?)/plugin_out')
	fileHead=headPat.findall(sortedDir)[-1]
	
	if not os.path.exists(resDir+'/'+fileHead):
		os.mkdir(resDir+'/'+fileHead)
	
	bedDct=getBedInfo()
	
#	print(bedDct)
	
	os.chdir(sortedDir)
	sortedBams=glob.glob('*sorted.bam')
	
#==============================================================================
# 	pool=Pool(cpuNbr)
# 	for subsorted in sortedBams:
# #		getSams(sortedDirGs,sortedBamGs,bedDctGs,resDirGs,headGs)
# 		pool.apply_async(getSams,args=(sortedDir,subsorted,bedDct,resDir,fileHead))
# 	pool.close()
# 	pool.join()
#==============================================================================
#	AVOID multiprocessing in testing procedure
	
#	getSams(sortedDir,sortedBams[0],bedDct,resDir,fileHead)
	
	totalWinbr=getotalWinbr(bedDct,winSize)
	maxWinbr=getMaxwinbr(fpr,totalWinbr)
	
	fileHeadRef=resDir+'/'+fileHead+'/'
	os.chdir(resDir+'/'+fileHead+'/')
	
	samFiles=glob.glob('*.sam')
	
#==============================================================================
# 	pool=Pool(cpuNbr)
# 	for samFile in samFiles:
# #		note that directory has changed to 'resDir/fileHead/'
# 		pool.apply_async(mainExc,args=(bedDct,samFile,winSize,totalWinbr,maxWinbr,fpr,fileHeadRef,resDir))
# 	pool.close()
# 	pool.join()
#==============================================================================

#	avoid applying multiprocessing in the beginning
	print(samFiles[0]+' mainExc processing...')
	mainExc(bedDct,samFiles[0],winSize,totalWinbr,maxWinbr,fpr,fileHead,resDir)
	
	return
	
def mainExc(bedDctMe,samMe,winSizeMe,totalWinbrMe,maxWinbrMe,fprMe,fileHeadRefMe,resDirMe):
#	readsMat structure
#	stt end reads zScore lowerTail upperTail
#	note that mainExc is designed for parallel, taking each sam file as input
	print(resDirMe)
	print(resDirMe+'/'+fileHeadRefMe)
	os.chdir(resDirMe+'/'+fileHeadRefMe)
	
	readsCtLst=getReadsMat(bedDctMe,samMe,winSizeMe,resDirMe,fileHeadRefMe)
	readsMat=np.array(readsCtLst)
	
	meanDepth=np.mean(readsMat[:,2])
	stdDepth=np.std(readsMat[:,2])
	medianDepth=np.median(readsMat[:,2])
	
	readsMat[:,3]=np.round((readsMat[:,2]-meanDepth)/stdDepth,5)
	readsMat[:,4]=np.round(stats.norm.cdf(readsMat[:,3]),5)
	readsMat[:,5]=1-readsMat[:,4]
	
#	getEvtinfo(readsMatGe,totalWinbrGe,maxWinbrGe,winSizeGe,fprGe,meanDptGe,stdDptGe,mediDptGe)
#		return dupEvtMat,delEvtMat,dupFilMerged,delFilMerged
	
#	for info matrix:
#		stt end reads zScore lowerTail upperTail del/dupLength
#	for filMerged:
#		sttP endP sttIdx endIdx
	
	dupEvtMatMe,delEvtMatMe,dupFilMerged,delFilMerged = getEvtinfo(readsMat,totalWinbrMe,maxWinbrMe,winSizeMe,fprMe,meanDepth,stdDepth,medianDepth)
	
	dupEvtFile=fileHeadRefMe+samMe+'.dupEvtMat'
	delEvtFile=fileHeadRefMe+samMe+'.delEvtMat'
	dupMergeFile=fileHeadRefMe+samMe+'.dupMergeWin'
	delMergeFile=fileHeadRefMe+samMe+'.delMergeWin'
	
	dupEvtMatMe.savetxt(dupEvtFile,delimiter='\t')
	delEvtMatMe.savetxt(delEvtFile,delimiter='\t')
	dupFilMerged.savetxt(dupMergeFile,delimiter='\t')
	delFilMerged.savetxt(delMergeFile,delimiter='\t')
	
	return
	
def getBedInfo():
#	for main
	print('getBedInfo starts...')
	
	global bedFile
	bedDctGbi={}
	
	infh=open(bedFile,'r')
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
#		note that annotated linear[3] have overlaps, thus using linear[1] as key
		bedDctGbi[int(linear[1])]=[int(linear[2]),linear[3]]
	infh.close()
	return bedDctGbi
	
def getotalWinbr(bedDctGw,winSizeGw):
#	for main
	print('getotalWinbr starts...')
	winNbr=0
	for sttP in bedDctGw.keys():
		endP=bedDctGw[sttP][0]+winSizeGw
		
		candStt=sttP
		candEnd=sttP+winSizeGw
		
		while candEnd <= endP:
			candStt += winSizeGw
			candEnd += winSizeGw
			winNbr += 1
	return winNbr
	
def getMaxwinbr(fprGm,winNbrGm):
#	for main
	print('getMaxwinbr starts...')
	maxWinbr=0
	for i in range(1,winNbrGm):
		if math.pow((fprGm/(winNbrGm/i)),1/i)>0.5:
			maxWinbr=i
			break
	return maxWinbr
	
def getSams(sortedDirGs,sortedBamGs,bedDctGs,resDirGs,headGs):
#	for main

	print('%s getSams starts...'	%	(sortedBamGs))
	
	os.chdir(sortedDirGs)
#	print(sortedBamGs+' samtools starts...')
	
	tmpSam=resDirGs+'/'+headGs+'/'+sortedBamGs+'.sam'
	posLst=sorted(bedDctGs.keys())
	
#	print(posLst)
	
	sp.call('touch %s' % (tmpSam),shell=True)
	for subStt in posLst:
#		use q20 for filtration of reads quality
		subEnd=bedDctGs[subStt][0]
		samCmd='samtools view %s -q 30 chrX:%s-%s >> %s' % (sortedBamGs,str(subStt),str(subEnd),tmpSam)
		sp.call(samCmd,shell=True)
	return
	
def getReadsMat(bedDctGrm,samGrm,winSizeGrm,resDirGrm,fileHeadGrm):
	
	print('%s getReadsMat starts...'	%	(samGrm))
	os.chdir(resDirGrm+'/'+fileHeadGrm)
	
	infh=open(samGrm,'r')
	sttPosLst=[]
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
		sttPosLst.append(int(linear[3]))
	
	samReadsInfo=[]
	
	for sttP in sorted(bedDctGrm.keys()):
		endP=bedDctGrm[sttP][0]+winSizeGrm
		
		candStt=sttP
		candEnd=sttP+winSizeGrm
		
		while candEnd <= endP:
			readsCt=0
			for subPos in sttPosLst:
				if candStt <= subPos < candEnd:
					readsCt += 1
				elif subPos >=candEnd:
#					note that matrix size can not be modified in np, thus change shape in original construction
#					stt end reads zScore lowerTail upperTail
					samReadsInfo.append([candStt,candEnd,readsCt,0,0,0])
					break
			candStt=candEnd
			candEnd+=winSizeGrm
	
	return samReadsInfo
	
def getEvtinfo(readsMatGe,totalWinbrGe,maxWinbrGe,winSizeGe,fprGe,meanDptGe,stdDptGe,mediDptGe):
	
#	for each position in readsMatGe, get largest l for both dup and del conditions
#	thus returns a dup matrix and a del matrix
#	stt end reads zScore lowerTail upperTail
#	all returns with numpy object
	
#	dupEvtMat=dupEvt(readsMatGe,totalWinbrGe,fprGe,maxWinbrGe)
	
	delEvtMat=delEvt(readsMatGe,totalWinbrGe,fprGe,maxWinbrGe)
	
#	dupEvtMat=np.array(dupEvtMat)
	delEvtMat=np.array(delEvtMat)
	
#	dupMergedPos=mergeEvt(dupEvtMat)
	delMergedPos=mergeEvt(delEvtMat)
	
	sys.exit(1)
	
	dupMergedPos=np.array(dupMergedPos)
	delMergedPos=np.array(delMergedPos)
	
#	postFiltra(readsMatSz,oveMeanDptSz,oveStdDptSz,oveMediDptSz,mergedLstSz)
	dupFilMerged=postFiltra(dupEvtMat,meanDptGe,stdDptGe,mediDptGe,dupMergedPos)
	delFilMerged=postFiltra(delEvtMat,meanDptGe,stdDptGe,mediDptGe,delMergedPos)
	
	dupFilMerged=np.array(dupFilMerged)
	delFilMerged=np.array(delFilMerged)
	
	return dupEvtMat,delEvtMat,dupFilMerged,delFilMerged


def dupEvt(readsMatDe,totalWinbrDe,fprDe,maxWinbrDe):
#	dupEvt and delEvt returns a list in form of:
#	stt end reads zScore lowerTail upperTail del/dupLength
#	del and dup merge procedure would be carried out independently
	
	print('dupEvt starts...')

	dupEvtLst=[]
	for i in range(readsMatDe.shape[0]):
		curLen=0
		curThres=math.pow(fprDe/(totalWinbrDe/(curLen + 1)),1/(curLen + 1))
		
		print 'current Thres:' + str(curThres)
		
		curDupPv=[readsMatDe[i,5]]
		while np.max(curDupPv) < curThres and curLen < maxWinbrDe and (i+curLen)<readsMatDe.shape[0]:
			curLen += 1
			curThres=math.pow(fprDe/(totalWinbrDe/(curLen + 1)),1/(curLen + 1))

			print 'while curThres: ' + str(curThres)
			
			if (i+curLen)<readsMatDe.shape[0]:
				curDupPv.append(readsMatDe[i+curLen,5])
				print('sttdup: %s    curLen: %s'	%	(str(readsMatDe[i+curLen,0]),str(curLen)))

		print 'current rowNbr: ' + str(i)
		
		if curLen > 0:
			curLst=[]
			for subItem in readsMatDe[i,:]:
				curLst.append(subItem)
			curLst.append(curLen)
			dupEvtLst.append(curLst)
			
			print 'current Length: ' + str(curLen)
			
	return dupEvtLst
	
	
def delEvt(readsMatEl,totalWinbrEl,fprEl,maxWinbrEl):
	print('delEvt starts...')
	delEvtLst=[]
	for i in range(readsMatEl.shape[0]):
		curLen=0
		curThres=math.pow(fprEl/(totalWinbrEl/(curLen + 1)),1/(curLen + 1))
		
		print 'current Thres:' + str(curThres)
		
		curDelPv=[readsMatEl[i,4]]
		while np.max(curDelPv) < curThres and curLen < maxWinbrEl and (i+curLen)<readsMatEl.shape[0]:
			curLen += 1
			curThres=math.pow(fprEl/(totalWinbrEl/(curLen + 1)),1/(curLen + 1))
			
			print 'while curThres: ' + str(curThres)
			
			if (i+curLen)<readsMatEl.shape[0]:
				curDelPv.append(readsMatEl[i+curLen,4])
				print('sttdel: %s    curLen: %s'	%	(str(readsMatEl[i+curLen,0]),str(curLen)))

		print 'cur row nbr:' + str(i)
		
		if curLen>0:
			curLst=[]
			for subItem in readsMatEl[i,:]:
				curLst.append(subItem)
			curLst.append(curLen)
			delEvtLst.append(curLst)
			
			print 'curLen: ' + str(curLen)
			
	return delEvtLst


#	define a global gap size using for merge dup event
#	readsMat would be in form of:
#		stt end reads zScore lowerTail upperTail del/dupLength
def mergeEvt(evtMatMe):
#	this return a matrix containning sttP and endP after merging, alongside with readsMat for one-side Z-testing
	ind=0
#	define a gapThres here
	gapThres=500
	
	upper=evtMatMe.shape[0]
	print('upper: '+str(upper))
	
	if upper==0:
		return 0
	
	mergedLst=[]
	
	while ind < upper:
		curLen=1
		while (ind+curLen-1) < upper and (evtMatMe[ind+curLen,0]-evtMatMe[ind+curLen-1,1]) < gapThres:
			curLen += 1
#		sttP endP sttIdx endIdx
		mergedLst.append([evtMatMe[ind,0],evtMatMe[ind+curLen-1,1],ind,ind+curLen-1])
		ind += curLen
	return mergedLst
	

def postFiltra(readsMatSz,oveMeanDptSz,oveStdDptSz,oveMediDptSz,mergedLstSz):
#	note that oneSideZtest were carried out using its mean versus overall mean for 27K
#	stt end reads zScore lowerTail upperTail del/dupLength	
#	sttP endP sttIdx endIdx
	if mergedLstSz==0:
		return 0
	
	filtedMergeLst=[]
	
	#################################################
	#	define a threshold for one-side Z test
	#################################################
	zValThres=0.005
	
	lmediThres=0.75*oveMediDptSz
	umediThres=1.25*oveMediDptSz
	
	for ind in range(mergedLstSz.shape[0]):
		sttidx=mergedLstSz[ind,2]
		endidx=mergedLstSz[ind,3]
		
		gapMediDpt=np.median(readsMatSz[sttidx:(endidx+1),2])
		gapMeanDpt=np.mean(readsMatSz[sttidx:(endidx+1),2])
		
		zVal=np.round((gapMeanDpt-oveMeanDptSz)/oveStdDptSz,6)
		
		zlPvl=stats.norm.cdf(zVal)
		zuPvl=stats.norm.cdf(zVal)
		
		if zlPvl < zValThres or zuPvl < zValThres:
			if gapMediDpt > umediThres or gapMediDpt < lmediThres:
				filtedMergeLst.append(mergedLstSz[ind].tolist())
			
	return filtedMergeLst


if __name__=='__main__':
	main()
