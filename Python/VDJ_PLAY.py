# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:59:06 2020

@author: yjk12
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
filtered = pd.read_csv('BM2_filtered_contig_annotations.csv')
filtered=filtered[['barcode','umis']]
filtered=filtered.groupby('barcode').sum()

#filtered=filtered.set_index('barcode')
filtered=filtered.sort_values('umis')
#filtered.style.background_gradient(cmap='Blues')
plt.pcolor(filtered)
plt.yticks(np.arange(0.5, len(filtered.index), 1), filtered.index)
plt.xticks(np.arange(0.5, len(filtered.columns), 1), filtered.columns)
plt.show()

filtered.to_csv('sort_vdj.csv')