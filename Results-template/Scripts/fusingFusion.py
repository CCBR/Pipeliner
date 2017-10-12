#Last updated 10/12/17
import numpy as np
import pandas as pd 
#import seaborn as sns
#import matplotlib.pyplot as plt
#from scipy import stats 
import sys, os
#from matplotlib.backends.backend_pdf import PdfPages

fusionCatcher=pd.read_csv(sys.argv[1],sep='\t')
starFusion=pd.read_csv(sys.argv[2],sep='\t')

if os.path.isfile(sys.argv[3]):
	oncoFuseSF=pd.read_csv(sys.argv[3],sep='\t')
else:	
	oncoFuseSF=pd.DataFrame(columns=['5_FPG_GENE_NAME','3_FPG_GENE_NAME','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','5_DOMAINS_RETAINED','3_DOMAINS_RETAINED','5_DOMAINS_BROKEN','3_DOMAINS_BROKEN','5_PII_RETAINED','3_PII_RETAINED'])

if os.path.isfile(sys.argv[4]):
        oncoFuseFC=pd.read_csv(sys.argv[4],sep='\t')
else: 
	oncoFuseFC=pd.DataFrame(columns=['5_FPG_GENE_NAME','3_FPG_GENE_NAME','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','5_DOMAINS_RETAINED','3_DOMAINS_RETAINED','5_DOMAINS_BROKEN','3_DOMAINS_BROKEN','5_PII_RETAINED','3_PII_RETAINED'])

if os.path.isfile(sys.argv[5]):
        detectorSF=pd.read_csv(sys.argv[5],sep='\t')
else: 
        detectorFC=pd.DataFrame(columns=['detG1','detG2'])

if os.path.isfile(sys.argv[6]):
        detectorFC=pd.read_csv(sys.argv[6],sep='\t')
else: 
        detectorFC=pd.DataFrame(columns=['detG1','detG1'])
print(sys.argv[1].split('/')[1])


if fusionCatcher.shape[0]>0:
	fc=fusionCatcher[['Gene_1_symbol(5end_fusion_partner)','Gene_2_symbol(3end_fusion_partner)', 'Gene_1_id(5end_fusion_partner)','Gene_2_id(3end_fusion_partner)', 'Fusion_description', 'Counts_of_common_mapping_reads','Spanning_pairs','Spanning_unique_reads','Fusion_sequence','Predicted_effect','Predicted_fused_proteins']].drop_duplicates()  
	fc.columns = ['Gene1','Gene2','Gene1 ID', 'Gene2 ID', 'Fusion_description', 'Counts_of_common_mapping_reads','Spanning_pairs','Spanning_unique_reads','Fusion_sequence','Predicted_effect','Predicted_fused_proteins']
	fc = fc.astype(str)
	fc=fc.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID']).agg('|'.join)
	fc.reset_index(inplace=True)  
else:
	fc=fusionCatcher[['Gene_1_symbol(5end_fusion_partner)','Gene_2_symbol(3end_fusion_partner)', 'Gene_1_id(5end_fusion_partner)','Gene_2_id(3end_fusion_partner)', 'Fusion_description', 'Counts_of_common_mapping_reads','Spanning_pairs','Spanning_unique_reads','Fusion_sequence']].drop_duplicates()
        fc.columns = ['Gene1','Gene2','Gene1 ID', 'Gene2 ID', 'Fusion_description', 'Counts_of_common_mapping_reads','Spanning_pairs','Spanning_unique_reads','Fusion_sequence']
        fc = fc.astype(str)
        fc=fc.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID']).agg('|'.join)
        fc.reset_index(inplace=True)

if detectorFC.shape[0] <1:
    fc['Fusion Inspector']="False"
else:
    detectorFC['detG1'], detectorFC['detG2'] = detectorFC['#FusionName'].str.split('--', 1).str
    detFC=detectorFC[['detG1','detG2']].drop_duplicates()
    fc=pd.merge(fc, detFC, how='left', left_on=['Gene1','Gene2'], right_on = ['detG1','detG2'])    
    fc['Fusion Inspector']=(fc['Gene1']==fc['detG1']) & (fc['Gene2']==fc['detG2'])
    fc.drop(['detG1', 'detG2'], axis=1, inplace=True)

offc=oncoFuseFC[['5_FPG_GENE_NAME','3_FPG_GENE_NAME','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','5_DOMAINS_RETAINED','3_DOMAINS_RETAINED','5_DOMAINS_BROKEN','3_DOMAINS_BROKEN','5_PII_RETAINED','3_PII_RETAINED']].drop_duplicates()
offc.columns=['Gene1','Gene2','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','5_DOMAINS_RETAINED','3_DOMAINS_RETAINED','5_DOMAINS_BROKEN','3_DOMAINS_BROKEN','5_PII_RETAINED','3_PII_RETAINED']
offc = offc.astype(str)
offc=offc.groupby(['Gene1','Gene2']).agg('|'.join)
offc.reset_index(inplace=True)  

if oncoFuseFC.shape[0] >0:
    combFC=pd.merge(fc, offc, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
    combFC['Pred']= 'Catch'
else:
    combFC=fc
    combFC['Pred']= 'Catch'

print(combFC.head())

starFusion['Gene1'], starFusion['Gene2'] = starFusion['#FusionName'].str.split('--', 1).str
starFusion['Gene1 ID']=starFusion.LeftGene.apply(lambda x: x.split('^')[1]).apply(lambda x: x.split('.')[0])
starFusion['Gene2 ID']=starFusion.RightGene.apply(lambda x: x.split('^')[1]).apply(lambda x: x.split('.')[0])

sf=starFusion[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','SpliceType','SpanningFragCount','LeftBreakpoint','RightBreakpoint']].drop_duplicates()  
sf = sf.astype(str)
sf=sf.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID']).agg('|'.join)
sf.reset_index(inplace=True)  

if detectorSF.shape[0] <1:
    sf['Fusion Inspector']='False'
else:
    detectorSF['detG1'], detectorSF['detG2'] = detectorSF['#FusionName'].str.split('--', 1).str
    detSF=detectorSF[['detG1','detG2']].drop_duplicates()
    sf=pd.merge(sf, detSF, how='left', left_on=['Gene1','Gene2'], right_on = ['detG1','detG2'])    
    sf['Fusion Inspector']=(sf['Gene1']==sf['detG1']) & (sf['Gene2']==sf['detG2'])
    sf.drop(['detG1', 'detG2'], axis=1, inplace=True)


ofsf=oncoFuseSF[['5_FPG_GENE_NAME','3_FPG_GENE_NAME','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','5_DOMAINS_RETAINED','3_DOMAINS_RETAINED','5_DOMAINS_BROKEN','3_DOMAINS_BROKEN','5_PII_RETAINED','3_PII_RETAINED']].drop_duplicates()
ofsf.columns=['Gene1','Gene2','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','5_DOMAINS_RETAINED','3_DOMAINS_RETAINED','5_DOMAINS_BROKEN','3_DOMAINS_BROKEN','5_PII_RETAINED','3_PII_RETAINED']
ofsf = ofsf.astype(str)
ofsf=ofsf.groupby(['Gene1','Gene2']).agg('|'.join)
ofsf.reset_index(inplace=True)  

if oncoFuseSF.shape[0] >0:
    combSF=pd.merge(sf, ofsf, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
    combSF['Pred']= 'Star'
else:
    combSF=sf
    combSF['Pred']= 'Star'

if combFC.shape[0]<1:
	
	combFC = pd.DataFrame(columns=['Gene1','Gene2','Gene1 ID', 'Gene2 ID', 'Fusion_description', 'Counts_of_common_mapping_reads','Spanning_pairs','Spanning_unique_reads','Fusion_sequence','Pred','Fusion Inspector'])

#	combFC=newcombFC
 
if detectorFC.shape[0] < 1 and oncoFuseFC.shape[0] <1: 
    selectFC=combFC[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Pred','Fusion Inspector']]
elif detectorFC.shape[0] < 1: 
    selectFC=combFC[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','Pred','Fusion Inspector']]
elif oncoFuseFC.shape[0] <1:
    selectFC=combFC[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Fusion Inspector','Pred']]
else: 
    selectFC=combFC[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Fusion Inspector','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','Pred']]

if detectorSF.shape[0] < 1 and oncoFuseSF.shape[0] <1: 
    selectSF=combSF[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Pred','Fusion Inspector']]
elif detectorSF.shape[0] < 1: 
    selectSF=combSF[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','Pred','Fusion Inspector']]
elif oncoFuseSF.shape[0] <1:
    selectSF=combSF[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Fusion Inspector','Pred']]
else:
    selectSF=combSF[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Fusion Inspector','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN','Pred']]

x=pd.concat([selectFC,selectSF], axis=0)
x = x.astype(str)
x=x.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID']).agg('|'.join)
x.reset_index(inplace=True)
conditions = [    
    (x['Pred'] == 'Catch|Star'),
    (x['Pred'] == 'Star'),
    (x['Pred'] == 'Catch')]
choices = ['Both', 'Star', 'Catcher']
x['Predicted By'] = np.select(conditions, choices)
#x.drop(['Pred'], axis=1, inplace=True)
x['Sample'] = sys.argv[1].split('/')[1]
x.to_csv(path_or_buf=sys.argv[7],sep='\t',index=False)

combFC.drop(['Pred'], axis=1, inplace=True)
combSF.drop(['Pred'], axis=1, inplace=True)
combFC['Sample'] = sys.argv[1].split('/')[1]
combSF['Sample'] = sys.argv[1].split('/')[1]

combFC.to_csv(path_or_buf=sys.argv[8],sep='\t')
combSF.to_csv(path_or_buf=sys.argv[9],sep='\t')


ofsf=oncoFuseSF[['5_FPG_GENE_NAME','3_FPG_GENE_NAME','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN']].drop_duplicates()
ofsf.columns=['Gene1','Gene2','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN']
offc=oncoFuseFC[['5_FPG_GENE_NAME','3_FPG_GENE_NAME','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN']].drop_duplicates()
offc.columns=['Gene1','Gene2','P_VAL_CORR','DRIVER_PROB','EXPRESSION_GAIN']

if oncoFuseFC.shape[0] >0 or oncoFuseSF.shape[0] >0: 
	of=pd.concat([ofsf,offc],axis=0)
	of.groupby(['Gene1','Gene2'])
	of.drop_duplicates()#, sort=False)['DRIVER_PROB'].max()#.agg('|'.join)
	#of.reset_index(inplace=True)
	x=x[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Fusion Inspector','Predicted By','Sample','Pred']].drop_duplicates()	
	y=pd.merge(x, of, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
	y.fillna('notFound', inplace=True)
	y=y.drop_duplicates()
else:
	of=pd.concat([ofsf,offc],axis=0)
        of.groupby(['Gene1','Gene2'])
        of.drop_duplicates()#, sort=False)['DRIVER_PROB'].max()#.agg('|'.join)
        #of.reset_index(inplace=True)
        x=x[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Fusion Inspector','Predicted By','Sample','Pred']].drop_duplicates()
        y=pd.merge(x, of, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
        y.fillna('notFound', inplace=True)
        y=y.drop_duplicates()

	
y.to_csv(path_or_buf=sys.argv[10],sep='\t')
