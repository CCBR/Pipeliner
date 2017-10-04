#version1
import numpy as np
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1],sep='\t')

df['Sample']=df['Sample'].apply(str)
df['Gene1']=df['Gene1'].str.upper()
df['Gene2']=df['Gene2'].str.upper()

df1= df[['Gene1','Gene2','Sample','Pred','Fusion Inspector','Gene1 ID','Gene2 ID']].drop_duplicates()

df2=df[['Gene1 ID','Gene2 ID','Sample']].drop_duplicates()

count=df2.groupby(['Gene1 ID','Gene2 ID']).size().reset_index(name='Count')
catch = df1.groupby(['Gene1 ID', 'Gene2 ID'])['Pred'].apply(lambda x: x[x.str.contains('Catch')].count()).reset_index(name='Fusion Catcher')
star= df1.groupby(['Gene1 ID', 'Gene2 ID'])['Pred'].apply(lambda x: x[x.str.contains('Star')].count()).reset_index(name='Star Fusion')
det= df1.groupby(['Gene1 ID', 'Gene2 ID'])['Fusion Inspector'].apply(lambda x: x[x.str.contains('True')].count()).reset_index(name='Fusion Inspector')

driver =df.groupby(['Gene1 ID','Gene2 ID'], sort=False)['DRIVER_PROB'].max().reset_index(name='DRIVER_PROB')
corr=df.groupby(['Gene1 ID','Gene2 ID'], sort=False)['P_VAL_CORR'].min().reset_index(name='P_VAL_CORR')
exp=df.groupby(['Gene1 ID','Gene2 ID'], sort=False)['EXPRESSION_GAIN'].max().reset_index(name='EXPRESSION_GAIN')


count=pd.merge(count, catch, how='left', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])
count=pd.merge(count, star, how='left', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])
count=pd.merge(count, det, how='left', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])
count=pd.merge(count, driver, how='left', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])
count=pd.merge(count, corr, how='left', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])
count=pd.merge(count, exp, how='left', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])


wanted= df1[['Gene1','Gene2','Gene1 ID', 'Gene2 ID','Sample']].drop_duplicates()

final=pd.merge(count, wanted, how='inner', left_on=['Gene1 ID','Gene2 ID'], right_on = ['Gene1 ID','Gene2 ID'])

final=final.groupby(['Gene1','Gene2','Gene1 ID', 'Gene2 ID', 'Count','Star Fusion','Fusion Catcher','Fusion Inspector','DRIVER_PROB','P_VAL_CORR','EXPRESSION_GAIN']).agg(';'.join)

final.reset_index(inplace=True)


final.to_csv(path_or_buf='../summary/summaryCounts.tsv',sep='\t',index=False)




