import glob
import pandas as pd
import sys

#path =r'/home/abdelmaksoudaa/fusion_test/' # use your path
path = sys.argv[1]
allFiles = glob.glob(path + "*")
frame = pd.DataFrame()
list_ = []
for file_ in allFiles:
    df = pd.read_csv(file_,index_col=None, header=0,error_bad_lines=False,sep='\t')
    list_.append(df)
frame = pd.concat(list_)

frame.to_csv(path_or_buf='../all/easySum.txt',sep='\t')
