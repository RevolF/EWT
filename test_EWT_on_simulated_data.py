# using subprocess.Popen and raw mapped and sorted bams as input
# for this version, bamdir should be directed to *.chrX.sorted.bam

from __future__ import division
import re,sys,os,glob
from optparse import OptionParser
from multiprocessing import Pool
import subprocess as sp
import numpy as np
import math
from scipy import stats

parser=OptionParser()

parser.add_option(
	'-D',
	'--sortedbam-dir',
	dest='bamdir',
	help='sorted bam dir'
	)
	
parser.add_option(
	'-M',
	'--processing-mode',
	dest='mode',
	help='processing mode, 1 for single 0 for filefolder batch, default set to 0',
	default='0'
	)
	
parser.add_option(
	'-T',
	'--work-dir',
	dest='workdir',
	help='work dir for GC results'
	)
	
parser.add_option(
	'-F',
	'--fasta-dir',
	dest='fadir',
	help='reference fasta folder,default set to /home/ljzhang/data/hg19/',
	default='/home/ljzhang/data/hg19/'
	)
	
parser.add_option(
	'-R',
	'--FPR-thres',
	dest='fpr',
	help='predefined false positive rate,default set to 0.05',
	default='0.05'
	)
	
parser.add_option(
	'-I',
	'--input-depthFile',
	dest='depthFile',
	help='simulated depthfile'
	)

(options,args)=parser.parse_args()

# if not options.bamdir or not options.workdir:
	## print(111)
	# parser.print_help()
	# sys.exit(1)
	
if not options.depthFile or not options.bamdir or not options.workdir:
	parser.print_help()
	sys.exit(1)
	
bamdir=options.bamdir
workdir=bamdir+'/'+options.workdir
referfa=options.fadir+'ucsc.hg19.chrX.fasta'
fpr=float(options.fpr)
depthFile=options.depthFile

if not os.path.exists(workdir):
	os.mkdir(workdir)

# a GC_content file was generated at: workdir+'chrX_GC_content'
	
def main():
	# print(bamdir)
	global bamdir,workdir,referfa,fpr
	
	os.chdir(bamdir)
	getGCcontent()
	# files=glob.glob('*.sorted.bam')
	# print(files)
	# pool=Pool(len(files))
	# pool.map(mainExc,files)
	# pool.close()
	# pool.join()
	
	mainExc(depthFile)
	
def mainExc(file):
	# samtoolpipecmd='samtools view %s chrX'	%	(file)
	# pipe=sp.Popen(samtoolpipecmd,stdout=sp.PIPE,shell=True)
	# for line in iter(pipe.stdout.readline,''):
	global bamdir,workdir,fpr
	
	os.chdir(bamdir)
	
	print('changed dir')
	
	infh=open(workdir+'/chrX_GC_content','r')
	gcDct={}
	
	# gcIndex=[]
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
		# pos => GC nbr per 100 bp, 
		gcDct[int(linear[0])]=[int(linear[1]),0]
		# gcIndex.append(int(linear[0]))
		# print(linear[0])
		
	infh.close()
	
	print('ori gcDct generated')
	
	## huge time consumming step if reading from pipe
	# samtoolsPipeCmd='samtools view %s chrX -q 30 -@ 8 > %s'	%	(file,file+'.chrX.Q30.sorted.sam')
	## pipe=sp.Popen(samtoolsPipeCmd,stdout=sp.PIPE,shell=True)
	# if not os.path.exists(file+'.chrX.Q30.sorted.sam'):
		# sp.call(samtoolsPipeCmd,shell=True)
	
	## for line in iter(pipe.stdout.readline,''):
	
	for line in open(file,'r').xreadlines():
		linear=line.strip().split('\t')
		
		arIndx=int(linear[1])//100
		arIndx=arIndx*100+1
		
		gcDct[arIndx][1]+=int(float(linear[2]))
		## print(gcDct[arIndx][1])
	
	# os.remove(file+'.chrX.Q30.sorted.sam')
	print('gcDct done')
	
	outfh=open(workdir+'/'+file+'.GC.windowed','w')
	
	keyar=gcDct.keys()
	keyar.sort()
	
	for key in keyar:
		outfh.write(str(key)+'\t'+'\t'.join([str(j) for j in gcDct[key]])+'\n')
		
	outfh.close()
	
	print('write out done 1')
	
	gcAdjustedDct=gcAdjust(dct=gcDct)
	del(gcDct)
	
	outfh=open(workdir+'/'+file+'.GC.windowed.adjusted_Info','w')
	keyar=gcAdjustedDct.keys()
	keyar.sort()
	
	# gcAdjustedDct:	pos => [gcNbr,depth,correctedDepth,zScore,lowerCdf,higherCdf]
	for key in keyar:
		outfh.write(str(key)+'\t'+'\t'.join([str(j) for j in gcAdjustedDct[key]])+'\n')
	
	outfh.close()
	
	print('write out done 2')
	
	allWinNbr=len(keyar)
	
	# maxWinsize=1
	# currentThres=round(math.pow((fpr/(allWinNbr/maxWinsize)),1/maxWinsize),4)
	# while currentThres<0.5:
		# maxWinsize+=1
		# currentThres=round(math.pow((fpr/(allWinNbr/maxWinsize)),1/maxWinsize),4)
	
	# pos => windowSize
	overallDup={}
	overallDel={}
	
	for rownbr in range(allWinNbr):
		currentPos=keyar[rownbr]
		currentNbr=rownbr
		
		currentDupLst=[]
		currentDupUpper=gcAdjustedDct[currentPos][5]
		
		locLen=1
		currentThres=round(math.pow((fpr/(allWinNbr/locLen)),1/locLen),4)
		
		while currentDupUpper<currentThres:
			locLen+=1
			currentDupLst.append(keyar[currentNbr])
			currentThres=round(math.pow((fpr/(allWinNbr/locLen)),1/locLen),4)
			currentNbr+=1
			currentDupUpper=gcAdjustedDct[keyar[currentNbr]][5]
		
		if currentDupLst:
			# pos => [positions after pos]
			overallDup[currentPos]=currentDupLst
		
		
		currentNbr=rownbr
		currentDelLst=[]
		currentDelLower=gcAdjustedDct[currentPos][4]
		
		locLen=1
		currentThres=round(math.pow((fpr/(allWinNbr/locLen)),1/locLen),4)
		
		while currentDelLower<currentThres:
			locLen+=1
			currentDelLst.append(keyar[currentNbr])
			currentThres=round(math.pow((fpr/(allWinNbr/locLen)),1/locLen),4)
			currentNbr+=1
			currentDelLower=gcAdjustedDct[keyar[currentNbr]][4]
			
		if currentDelLst:
			# pos => [positions after pos]
			overallDel[currentPos]=currentDelLst
		
	print('overalls done')
		
	# gcAdjustedDct:	pos => [gcNbr,depth,correctedDepth,zScore,lowerCdf,higherCdf]
	overallAdjDepth=[ar[2] for ar in gcAdjustedDct.itervalues()]
	overallAdjDepthMean=np.array(overallAdjDepth).mean()
	del(overallAdjDepth)
	
	for key in overallDup.iterkeys():
		tmpLst=[gcAdjustedDct[i][2] for i in overallDup[key]]
		tmpMedian=getMedian(tmpLst)
		
		if 0.75*overallAdjDepthMean<tmpMedian<1.25*overallAdjDepthMean:
			del overallDup[key]
	
	for key in overallDel.iterkeys():
		tmpLst=[gcAdjustedDct[i][2] for i in overallDel[key]]
		tmpMedian=getMedian(tmpLst)
		
		if 0.75*overallAdjDepthMean<tmpMedian<1.25*overallAdjDepthMean:
			del overallDel[key]
		
	print('two fors done')
	
	overallDup=MergeGaps(dct=overallDup)
	overallDel=MergeGaps(dct=overallDel)
	
	dupKeyar=overallDup.keys()
	dupKeyar.sort()
	# gcAdjustedDct:	pos => [gcNbr,depth,correctedDepth,zScore,lowerCdf,higherCdf]
	for dupKey in dupKeyar:
		winNbr=len(overallDup[dupKey])
		tmpDpthLst=[gcAdjustedDct[i][2] for i in overallDup[dupKey]]
		tmpMean=np.array(tmpDpthLst).mean()
		tmpZsc=abs(tmpMean-overallAdjDepthMean)/math.sqrt(winNbr)
		sigLvl=stats.norm.cdf(tmpZsc)
		sigLvl=min([sigLvl,1-sigLvl])
		
		if sigLvl>0.000001:
			del overallDup[dupKey]
			
	delKeyar=overallDel.keys()
	delKeyar.sort()
	# gcAdjustedDct:	pos => [gcNbr,depth,correctedDepth,zScore,lowerCdf,higherCdf]
	for delKey in delKeyar:
		winNbr=len(overallDel[delKey])
		tmpDpthLst=[gcAdjustedDct[i][2] for i in overallDel[delKey]]
		tmpMean=np.array(tmpDpthLst).mean()
		tmpZsc=abs(tmpMean-overallAdjDepthMean)/math.sqrt(winNbr)
		sigLvl=stats.norm.cdf(tmpZsc)
		sigLvl=min([sigLvl,1-sigLvl])
		
		if sigLvl>0.000001:
			del overallDel[delKey]
	
	print('two key processing done')
	
	os.chdir(workdir)
	
	overallDelSortedKeyar=sorted(overallDel,key=lambda x:overallDel[x],reverse=True)
	delWinOutfh=open(file+'.delWin.sum','w')
	delDetailOutfh=open(file+'.delDetail.info','w')
	for delSortedKey in overallDelSortedKeyar:
		delWinOutfh.write(str(delSortedKey)+'\t'+str(len(overallDel[delSortedKey]))+'\t')
		delWinOutfh.write('\t'.join([str(j) for j in overallDel[delSortedKey]]))
		delWinOutfh.write('\n')
		
		for subPos in overallDel[delSortedKey]:
			delDetailOutfh.write(str(delSortedKey)+'\t'+str(len(overallDel[delSortedKey]))+'\t')
			delDetailOutfh.write(str(subPos)+'\t'+'\t'.join([str(j) for j in gcAdjustedDct[subPos]])+'\n')
			
	overallDupSortedKeyar=sorted(overallDup,key=lambda x:overallDup[x],reverse=True)
	dupWinOutfh=open(file+'.dupWin.sum','w')
	dupDetailOutfh=open(file+'.dupDetail.info','w')
	for dupSortedKey in overallDupSortedKeyar:
		dupWinOutfh.write(str(dupSortedKey)+'\t'+str(len(overallDup[dupSortedKey]))+'\t')
		dupWinOutfh.write('\t'.join([str(j) for j in overallDup[dupSortedKey]]))
		dupWinOutfh.write('\n')
		
		for subPos in overallDup[dupSortedKey]:
			dupDetailOutfh.write(str(dupSortedKey)+'\t'+str(len(overallDup[dupSortedKey]))+'\t')
			dupDetailOutfh.write(str(subPos)+'\t'+'\t'.join([str(j) for j in gcAdjustedDct[subPos]])+'\n')
			
	return


# a GC_content file was generated at: workdir+'chrX_GC_content'
def getGCcontent():
	global bamdir,workdir,referfa,fpr
	
	infh=open(referfa,'r')
	outfh=open(workdir+'/chrX_strip','w')
	# use 100bp window
	for line in infh.xreadlines():
		if line.startswith('>'):
			continue
		outfh.write(line.strip())
		
	infh.close()
	outfh.close()
	
	print('chrX_strip done!')
	
	infh=open(workdir+'/chrX_strip','r')
	outfh=open(workdir+'/chrX_GC_content','w')
	ind='a'
	pos='1'
	while ind.isalpha():
		ind=infh.read(100)
		outfh.write(pos+'\t')
		ind=ind.upper()
		gCount=ind.count('G')
		cCount=ind.count('C')
		outfh.write(str(gCount+cCount)+'\n')
		pos=str(int(pos)+100)
	infh.close()
	outfh.close()
	
	print('GC content done')
	
	os.remove(workdir+'/chrX_strip')

def gcAdjust(dct):
	# subGCdct: GC nbr per 100bp => depth array
	subGCdct={}
	univDepth=[]
	
	for value in dct.itervalues():
		if value[1]<1:
			continue
		if subGCdct.has_key(value[0]):
			subGCdct[value[0]].append(value[1])
		else:
			subGCdct[value[0]]=[value[1]]
		univDepth.append(value[1])
	# gcMedianDct: GC nbr per 100bp => depth median
	gcMedianDct={}
	
	# too many 0 values which wouold make median value 0
	# for key,value in subGCdct.iteritems():
		# print(str(key)+'\t=>\t')
		# print(value)
	
	# gcMedianDct gcNbr => medianDepth
	for key in subGCdct.iterkeys():
		gcMedianDct[key]=getMedian(subGCdct[key])
		print(str(key)+'\t=>\t'+str(gcMedianDct[key]))
		
	del subGCdct
	
	univMedian=getMedian(univDepth)
	univDepthList=[]
	for key,value in dct.iteritems():
		# if gcMedianDct[value[0]]==0:
			# adjustedDepth=0
		adjustedDepth=0
		if gcMedianDct.has_key(value[0]):
			adjustedDepth=round(value[1]*univMedian/gcMedianDct[value[0]])
			# print adjustedDepth
		univDepthList.append(adjustedDepth)
		dct[key].append(adjustedDepth)
	
	univDepthList=np.array(univDepthList)
	
	univDepthMean=univDepthList.mean()
	univDepthStd=univDepthList.std()
	
	nt=map(lambda x:dct[x].append(round((dct[x][2]-univDepthMean)/univDepthStd,3)),dct.keys())
	# get cdf for different z values
	# append lower cdf
	nt=map(lambda x:dct[x].append(round(stats.norm.cdf(dct[x][3]),7)),dct.keys())
	# append higher cdf
	nt=map(lambda x:dct[x].append(round(1-stats.norm.cdf(dct[x][3]),7)),dct.keys())

	return dct
	
def getMedian(lst):
	lst.sort()
	size=len(lst)
	
	if size%2 == 0:
		median=(lst[size//2]+lst[size//2-1])/2
	else:
		median=lst[(size-1)//2]
	
	return median
	
def MergeGaps(dct):
	keyar=np.array(dct.keys())
	diffar=np.diff(keyar)
	oddvalue=keyar[:-1][diffar>100]
	oddindex=[keyar.index(i) for i in oddvalue]
	oddindex.insert(0,0)
	oddindex.insert(len(oddindex),len(keyar)-1)
	
	newDct={}
	for i in range(len(oddindex)):
		sttIndx=oddindex[i]
		endIndx=oddindex[i+1]
		
		intv=[keyar[i] for i in range(sttIndx,endIndx+1)]
		tmpPosLst=[]
		nt=[tmpPosLst.extend(dct[i]) for i in intv]
		
		newDct[keyar[sttIndx]]=tmpPosLst
		
	return newDct
	
	
if __name__ == '__main__':
	main()


