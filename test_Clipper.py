import numpy as np
import pandas as pd
import sys
base_dir = "dir/to/Clipper"
#base_dir = "/Users/gexinzhou/Dropbox/clipper/Github/Clipper_python/"
sys.path.append(base_dir)
import Clipper


exp_d = pd.read_csv('./data/exp_d.csv',sep='\t',index_col=0).values.tolist()
back_d = pd.read_csv('./data/back_d.csv', sep='\t',index_col=0).values.tolist()
exp_e = pd.read_csv('./data/exp_e.csv', sep='\t',index_col=0).values.tolist()
back_e = pd.read_csv('./data/back_e.csv', sep='\t',index_col=0).values.tolist()

print(f"...starting the first test on sample differential analysis")
re1 = Clipper.clipper(score_exp=exp_d, score_back=back_d, analysis="differential", FDR=[0.01, 0.05, 0.1])
trueid = np.arange(2000)
discoveries = pd.DataFrame(re1["discoveries"][0])
print(f"...the proportion of left-out discoveries: {np.sum(~discoveries.isin(trueid).values)/len(trueid)}")
print()

print(f"...starting the second test on sample enrichment analysis")
re2 = Clipper.clipper(score_exp=exp_e, score_back=back_e, analysis="enrichment", FDR=[0.01, 0.05, 0.1])
trueid = np.arange(1000)
discoveries = pd.DataFrame(re2["discoveries"][0])
print(f"...the accuracy of discoveries: {np.sum(discoveries.isin(trueid).values)/len(trueid)}")
print()