import numpy as np
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1],sep='\t')

count=df.groupby(['Gene1','Gene2']).size().reset_index(name='Count')
catch = df.groupby(['Gene1', 'Gene2'])['Pred'].apply(lambda x: x[x.str.contains('Catch')].count())
star= df.groupby(['Gene1', 'Gene2'])['Pred'].apply(lambda x: x[x.str.contains('Star')].count())
det= df.groupby(['Gene1', 'Gene2'])['Fusion Inspector'].apply(lambda x: x[x.str.contains('True')].count())

driver =df.groupby(['Gene1','Gene2'], sort=False)['DRIVER_PROB'].max()
corr=df.groupby(['Gene1','Gene2'], sort=False)['P_VAL_CORR'].min()
exp=df.groupby(['Gene1','Gene2'], sort=False)['EXPRESSION_GAIN'].max()


count['Star Fusion'] = star.reset_index(inplace=False).iloc[:,2]
count['Fusion Catcher'] = catch.reset_index(inplace=False).iloc[:,2]
count['Fusion Inspector'] = det.reset_index(inplace=False).iloc[:,2]
count['DRIVER_PROB'] = driver.reset_index(inplace=False).iloc[:,2]
count['P_VAL_CORR'] = corr.reset_index(inplace=False).iloc[:,2]
count['EXPRESSION_GAIN'] = exp.reset_index(inplace=False).iloc[:,2]


wanted= df[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Sample']]

final=pd.merge(count, wanted, how='inner', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])

#final=pd.concat(['count','wanted'], axis=0)
final=final.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID', 'Count','Star Fusion','Fusion Catcher','Fusion Inspector','DRIVER_PROB','P_VAL_CORR','EXPRESSION_GAIN']).agg(','.join)
final.reset_index(inplace=True)

final.to_csv(path_or_buf='../summary/summaryCounts.tsv',sep='\t',index=False)



#fusionNames= pd.merge(fusions, fusions2, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
#fusionNames=fusionNames.groupby(['Gene1','Gene2','Count','Fusion Inspector','Predicted By']).agg(','.join)
#fusionNames.reset_index(inplace=True)
#fusionNames= fusionNames.sort_values('Count',ascending=False)



