# we generate depth file from norm(150,100)
# duplication from norm(250,100)
# deletion from norm(20,10)
# depth were rounded and <0 were adjusted to 0


from __future__ import division
import numpy as np
from scipy import stats
import pandas as pd
from optparse import OptionParser
import os,sys

parser=OptionParser()

parser.add_option(
	'-I',
	'--input-template',
	dest='posTpl',
	help='input position template'
	)
	
parser.add_option(
	'-G',
	'--gap-intv',
	dest='gap',
	help='gap intv on chrX for duplication or deletion, eg: "31132344-33234673"'
	)
	
parser.add_option(
	'-T',
	'--variation-type',
	dest='type',
	help='variation type, dup for duplication or del for deletion'
	)
	
(options,args)=parser.parse_args()

if not options.posTpl or not options.gap or not options.type:
	parser.print_help()
	sys.exit(1)
	
posTplfile=options.posTpl
gap=options.gap
type=options.type

def main():
	global posTplfile,gap,type
	
	# print(gap)
	# print(type)
	# print(posTplfile)
	
	infh=open(posTplfile,'r')
	posLst=[int(i.strip()) for i in infh.readlines()]
	
	posLen=len(posLst)
	
	depthDf=pd.DataFrame({'Depth':np.round(stats.norm.rvs(150,100,posLen)),'Pos':posLst,'chr':'chrX'},columns=['chr','Pos','Depth'])
	
	intvLst=gap.split('-')
	rowNbr=depthDf.loc[(depthDf['Pos']>int(intvLst[0])) & (depthDf['Pos']<int(intvLst[1]))].shape[0]

	if type=='dup':
		depthDf.loc[(depthDf['Pos']>int(intvLst[0])) & (depthDf['Pos']<int(intvLst[1])),('Depth')]=np.round(stats.norm.rvs(250,100,rowNbr))
		# depthDf.loc[(depthDf['Pos']>int(intvLst[0])) & (depthDf['Pos']<int(intvLst[1])),('Depth')]=500
		
	elif type=='del':
		depthDf.loc[(depthDf['Pos']>int(intvLst[0])) & (depthDf['Pos']<int(intvLst[1])),('Depth')]=np.round(stats.norm.rvs(20,10,rowNbr))
		# depthDf.loc[(depthDf['Pos']>int(intvLst[0])) & (depthDf['Pos']<int(intvLst[1])),('Depth')]=0
		
	depthDf.loc[depthDf['Depth']<0,'Depth']=0
	
	depthDf.to_csv(posTplfile+'.simulatedDepthFile.'+type,header=False,index=False,sep='\t')
	
	
if __name__=='__main__':
	main()

