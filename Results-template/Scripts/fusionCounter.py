import numpy as np
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1],sep='\t')

df1= df[['Gene1','Gene2','Sample','Pred','Fusion Inspector','Gene1 ID','Gene2 ID']].drop_duplicates()

count=df1.groupby(['Gene1','Gene2']).size().reset_index(name='Count')
catch = df1.groupby(['Gene1', 'Gene2'])['Pred'].apply(lambda x: x[x.str.contains('Catch')].count()).reset_index(name='Fusion Catcher')
star= df1.groupby(['Gene1', 'Gene2'])['Pred'].apply(lambda x: x[x.str.contains('Star')].count()).reset_index(name='Star Fusion')
det= df1.groupby(['Gene1', 'Gene2'])['Fusion Inspector'].apply(lambda x: x[x.str.contains('True')].count()).reset_index(name='Fusion Inspector')

driver =df.groupby(['Gene1','Gene2'], sort=False)['DRIVER_PROB'].max().reset_index(name='DRIVER_PROB')
corr=df.groupby(['Gene1','Gene2'], sort=False)['P_VAL_CORR'].min().reset_index(name='P_VAL_CORR')
exp=df.groupby(['Gene1','Gene2'], sort=False)['EXPRESSION_GAIN'].max().reset_index(name='EXPRESSION_GAIN')


count=pd.merge(count, catch, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
count=pd.merge(count, star, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
count=pd.merge(count, det, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
count=pd.merge(count, driver, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
count=pd.merge(count, corr, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])
count=pd.merge(count, exp, how='left', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])


wanted= df1[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Sample']]

final=pd.merge(count, wanted, how='inner', left_on=['Gene1','Gene2'], right_on = ['Gene1','Gene2'])

#final=pd.concat(['count','wanted'], axis=0)
final=final.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID', 'Count','Star Fusion','Fusion Catcher','Fusion Inspector','DRIVER_PROB','P_VAL_CORR','EXPRESSION_GAIN']).agg(','.join)
final.reset_index(inplace=True)

final.to_csv(path_or_buf='../summary/summaryCounts.tsv',sep='\t',index=False)






