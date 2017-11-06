import pandas as pd
import sys

f = sys.argv[1]
out = sys.argv[2]
df = pd.read_csv(f, delimiter='\t', index_col=False, header=None)
df=df.reindex(columns=[1,0])
df.to_csv(out, sep='\t', index=False, header=False)

