# using pandas dataframe for rewriting EWT for boosting
# save GC window file for chrX in /home/ljzhang/data/hg19, saving resource for batch running
# using params for indicating regions of inner reference and targeted regions
# best practise with a bed file which would indicate targeted sequencing region for experiment protocals
# start from 
# current experiment design strategie is to targeted on DMD region and ref region
# while in original paper, L denotes the overall windows number, l is local window nbr in region A
# if not using FPR method, then duplication region would have a big biase
# try using a moderate method
# set L to be overall windows nbr, l to be local, no adjustment need to be made

from __future__ import division
from optparse import OptionParser
from scipy import stats
from multiprocessing import Pool
import math
import numpy as np
import pandas as pd
import os,re,glob,sys
import traceback
import subprocess as sp

# except Exception:
	# traceback.print_exc()
	# pass

parser=OptionParser()

# put multiprocessing to last
parser.add_option(
	'-D',
	'--raw-dir',
	dest='rawdir',
	help='dir containning sorted bams'
	)
	
parser.add_option(
	'-T',
	'--work-dir',
	dest='workdir',
	help='work dir for EWT method'
	)
	
parser.add_option(
	'-F',
	'--fasta-dir',
	dest='fadir',
	help='dir containing ucsc.hg19.chrX.fasta, default set to /home/ljzhang/data/hg19',
	default='/home/ljzhang/data/hg19'
	)
	
parser.add_option(
	'-R',
	'--reference-region',
	dest='refRegion',
	help='indicating targeted region for GC adjustment, format: chrX:31132344-33234673'
	)
	
parser.add_option(
	'-M',
	'--dmd-region',
	dest='dmdRegion',
	help='indicating whole DMD region for CNV detection, format: chrX:31132344-33234673'
	)
	
parser.add_option(
	'-P',
	'--fpr-value',
	dest='fpr',
	help='fpr value for overall false positive ratio, default set to 0.05',
	default='0.05'
	)
	
parser.add_option(
	'-W',
	'--windows-nbr',
	dest='windnbr',
	help='total windows nbr, typically chrX length / 100, default 155270560/100 = 1552705',
	default='1552705'
	)
	
(options,args)=parser.parse_args()

if not options.rawdir or not options.workdir or not options.refRegion or not options.dmdRegion:
	parser.print_help()
	sys.exit(1)
	
rawdir=options.rawdir
workdir=rawdir+'/'+options.workdir
refRegion=options.refRegion
dmdRegion=options.dmdRegion
fpr=float(options.fpr)
windnbr=int(options.windnbr)


def main():
	# 
	global rawdir,workdir,refRegion,dmdRegion,fadir
	getGCcontent()
	
def mainExc(filemE):
	global rawdir,workdir,refRegion,dmdRegion,fadir
	os.chdir(rawdir)
	
	dmdRegionspan=dmdRegion.lstrip('chrX:').split('-')
	refRegionspan=refRegion.lstrip('chrX:').split('-')
	
	samDmdCmd='samtools view %s %s -q 30 > %s'	%	(filemE,dmdRegion,workdir+'/'+filemE+'.dmdRegion.sam')
	try:
		sp.call(samDmdCmd,shell=True)
	except Exception:
		print traceback.print_exc()
		
	samRefCmd='samtools view %s %s -q 30 > %s'	%	(filemE,refRegion,workdir+'/'+filemE+'.refRegion.sam')
	try:
		sp.call(samRefCmd,shell=True)
	except Exception:
		print traceback.print_exc()
		
	# changed to working directory now
	os.chdir(workdir)
	
	# 'pos','gcNbr','depth'
	refDataFrame=mapReads(gcFile='ref_GC_content.dct',samFile=filemE+'.refRegion.sam')
	dmdDataFrame=mapReads(gcFile='dmd_GC_content.dct',samFile=filemE+'.dmdRegion.sam')
	
	# gcNbr => medianDepth
	refMedianDf=getMedianDf(dataframegM=refDataFrame)
	overallMedianDepth=np.median(refDataFrame.depth)
	overallMeanDepth=np.mean(refDataFrame.depth)
	
	# pos gcNbr depth medianDepth adjDepth zScore lowerCDF upperCDF
	# note that position nbr is sorted
	dmdDfreshaped=gcAdjust(mediandfGc=refMedianDf,dmdDfGc=dmdDataFrame,overallMedian=overallMedianDepth)
	
	dupLst=[]
	delLst=[]
	
	for rownbr in range(dmdDfreshaped.shape[0]):
		rownbr
		curLen=0
		thres=np.power(fpr/(windnbr/curLen),1/curLen)
		
		while dmdDfreshaped.ix[curLen-1+rownbr,'upperCDF']<thres:
			curLen+=1
			thres=np.power(fpr/(windnbr/curLen),1/curLen)
			
		# saving all possible results in one dataframe, prepare for future possible analysis
		if curLen>0:
			dupLst.append([rownbr,curLen])
	
		curDelLen=0
		thresDel=np.power(fpr/(windnbr/curDelLen),1/curDelLen)
		
		while dmdDfreshaped.ix[curDelLen-1+rownbr,'lowerCDF']<thresDel:
			curDelLen+=1
			thresDel=np.power(fpr/(windnbr/curDelLen),1/curDelLen)
			
		if curDelLen>0:
			delLst.append([rownbr,curDelLen])
			
	dupDataFrame=pd.DataFrame(dupLst,columns=['Pos','winSize'])
	delDataFrame=pd.DataFrame(delLst,columns=['Pos','winSize'])
	
	dupGapInfo=mergeGap(dupDataFrame)
	delGapInfo=mergeGap(delDataFrame)
	
	# write raw filter to txt
	dupGapInfo.to_csv('dupRawRegionInfo.txt',sep='\t',index=False)
	delGapInfo.to_csv('delRawRegionInfo.txt',sep='\t',index=False)
	
	dupFiltered=extraFileter(dfEf=dupGapInfo,overallMeamEf=overallMeanDepth,shapedRawInfoDf=dmdDfreshaped)
	delFiltered=extraFileter(dfEf=delGapInfo,overallMeamEf=overallMeanDepth,shapedRawInfoDf=dmdDfreshaped)
	
	dupFiltered.to_csv('dupFilteredInfo.txt',sep='\t',index=False)
	delFiltered.to_csv('delFilteredInfo.txt',sep='\t',index=False)


def getGCcontent():
	global rawdir,fadir,refRegion,dmdRegion
	
	os.chdir(fadir)
	referfa=fadir+'/ucsc.hg19.chrX.fasta'
	
	infh=open(referfa,'r')
	outfh=open('chrX_strip','w')
	# use 100bp window
	for line in infh.xreadlines():
		if line.startswith('>'):
			continue
		outfh.write(line.strip())
		
	infh.close()
	outfh.close()
	
	print('chrX_strip done!')
	
	dmdRegionspan=dmdRegion.lstrip('chrX:').split('-')
	refRegionspan=refRegion.lstrip('chrX:').split('-')
	
	infh=open('chrX_strip','r')
	# here for dmd and ref region with GC count
	dmdoutfh=open('dmd_GC_content.dct','w')
	refoutfh=open('ref_GC_content.dct','w')
	
	dmdStt=(int(dmdRegionspan[0])//100)*100+1
	dmdEnd=int(dmdRegionspan[1])
	
	# be careful for seek method
	infh.seek(dmdStt,0)
	while (dmdStt-100)<dmdEnd:
		# dmdoutfh.write(dmdStt+'\t')
		line=infh.read(100)
		line=line.upper()
		gCount=line.count('G')
		cCount=line.count('C')
		dmdoutfh.write(dmdStt+'\t'+str(gCount+cCount)+'\n')
		dmdStt+=100
	dmdoutfh.close()
	
	refStt=(int(refRegionspan[0])//100)*100+1
	refEnd=int(refRegionspan[1])
	
	infh.seek(refStt,0)
	while (refStt-100)<refEnd:
		# refoutfh.write(refStt+'\t')
		line=infh.read(100)
		line=line.upper()
		gCount=line.count('G')
		cCount=line.count('C')
		refoutfh.write(str(refStt)+'\t'+str(gCount+cCount)+'\n')
		refStt+=100
	refoutfh.close()
	
	print('GC for dmd and ref finished')
	
	os.remove('chrX_strip')
	os.chdir(rawdir)
	
# mapReads will be called twice, for dmd and ref
# for each time being called
# return file in workdir in format of DataFrame
def mapReads(gcFile,samFile):
	global fadir,rawdir,workdir
	
	gcDct=fadir+'/'+gcFile
	samInfo=rawdir+'/'+samFile
	
	gcDctinfh=open(gcDct,'r')
	samInfoinfh=open(samInfo,'r')
	
	gcDctlst=gcDctinfh.readlines()
	gcDctlst=[i.strip().split('\t') for i in gcDctlst]
	
	gcDataFrame=pd.DataFrame(gcDctlst,columns=['pos','gcNbr'])
	# gcDataFrame.ix[:,'readCount']=0
	
	gcDataFrame=gcDataFrame.astype(np.int64)
	
	saminfolst=[]
	# use pd.value_counts()
	for line in samInfoinfh.xreadlines():
		linear=line.strip().split('\t')
		poskey=(int(linear[3])//100)*100+1
		saminfolst.append(poskey)
		
	ct=pd.value_counts(saminfolst)
	
	saminfoDataFrame=pd.DataFrame({'pos':ct.index,'depth':ct.values},columns=['pos','depth'])
	
	dfDepthAdded=pd.merge(gcDataFrame,saminfoDataFrame,on='pos')
	
	return dfDepthAdded
		
# 'pos','gcNbr','depth'
def getMedianDf(dataframegM):
	uniqGcvalue=sorted(pd.unique(dataframegM.gcNbr))
	medianDct=[]
	for i in uniqGcvalue:
		medianDct.append([i,np.median(dataframegM[dataframegM.gcNbr==i].ix[:,'depth'])])
	
	medianDataFrame=pd.DataFrame(medianDct,columns=['gc','medianDepth'])
	return medianDataFrame
	
def gcAdjust(mediandfGc,dmdDfGc,overallMedian):
	# pos gcNbr depth medianDepth(ref)
	mergedDf=pd.merge(dmdDfGc,mediandfGc,on='gcNbr')
	
	mergedDf.ix[:,'adjDepth']=np.round(mergedDf.ix[:,'depth']*overallMedian/mergedDf.ix[:,'medianDepth'],3)
	adjustedMean=np.mean(mergedDf.adjDepth)
	adjustedStd=np.std(mergedDf.adjDepth)
	
	mergedDf.ix[:,'zScore']=np.round((mergedDf.ix[:,'adjDepth']-adjustedMean)/adjustedStd,4)
	
	mergedDf.ix[:,'lowerCDF']=np.round(stats.norm.cdf(mergedDf.ix[:,'zScore']),6)
	mergedDf.ix[:,'upperCDF']=np.round(1-mergedDf.ix[:,'lowerCDF'],6)
	
	return mergedDf
	
def mergeGap(canDf):
	posDiff=np.diff(canDf.ix[:,'Pos'])
	posGap=canDf[:-1][posDiff>500].ix[:,'Pos']
	posGap=posGap.values
	
	gapList=[]
	for i in range(posGap.size):
		if i==0:
			gapList.append([canDf.ix[0,'Pos'],posGap[0]])
			sttP=posGap[0]+100
		else:
			endP=posGap[i]
			gapList.append([sttP,endP])
			sttP=posGap[i]+100
	gapList.append([sttP,canDf[-1].ix[:,'Pos']])
	
	gapDf=pd.DataFrame(gapList,columns=['startP','endP'])
	
	gapDf.ix[:,'gapNbr']=(gapDf.ix[:,'startP']-gapDf.ix[:,'endP'])//100
	
	return gapDf

def extraFileter(dfEf,overallMeamEf,shapedRawInfoDf):
	zthres=stats.norm.ppf(math.pow(10,-6))
	zthresU=abs(zthres)
	# pos gcNbr depth medianDepth adjDepth zScore lowerCDF upperCDF
	# startP	endP	gapNbr
	filterLst=[]
	
	for i in range(dfEf.shape[0]):
		subDf=shapedRawInfoDf[(shapedRawInfoDf.pos>=dfEf.startP[i]) & (shapedRawInfoDf.pos<=dfEf.endP[i])]
		subDfMedian=np.median(subDf.depth)
		
		if subDfMedian>0.75*overallMeamEf or subDfMedian<1.25*overallMeamEf:
			subDfMean=np.mean(subDf.depth)
			subDfStd=np.std(subDf.depth)
			subDfzScore=(subDfMean-overallMeamEf)/(subDfStd/math.sqrt(dfEf.gapNbr[i]))
			
			subcdf=stats.norm.cdf(subDfzScore)
			if subcdf<zthres or subcdf>zthresU:
				filterLst.append(dfEf.ix[i].tolist())
			else:
				continue
		else:
			continue
			
	filterDf=pd.DataFrame(filterLst,columns=['startP','endP','gapNbr'])
	
	return filterDf
	
	
