##path1 is the gem file file.paths, pasth2 is the output file path
###gem file should be decompress,please run "gzid -d *gem.gz",if file is txt file,could direcy running
####Author Jinpei Lin
####Times 2024.07

path1 = "/00gem/"
path2 = "/01RNA_image/"

###import env
import os
import glob
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import pickle
import loompy as lp
import os
import glob
import scipy
import scipy.sparse as sp
import cv2
from scipy.sparse import csr_matrix

dtype = {
    "geneID": "category",  # geneID
    "x": np.uint32,  # x
    "y": np.uint32,  # y
    # Multiple different names for total counts
    "MIDCounts": np.uint16,
    "MIDCount": np.uint16,
    "UMICount": np.uint16,
    "UMICounts": np.uint16,
}

rename = {
    "MIDCounts": "total",
    "MIDCount": "total",
    "UMICount": "total",
    "UMICounts": "total",
}


####extract filename
for key in [i [:-4] for i in os.listdir(path1)]:
    mtx_path=path1+key+".gem"
    data = pd.read_csv(
        mtx_path,
        sep="\t",
        dtype=dtype,
        comment="#",
    ).rename(columns=rename)
    y, x = data["x"].values, data["y"].values
    x_max, y_max = x.max(), y.max()
    shape = (x_max + 1, y_max + 1)
    x_min, y_min = data["x"].min(), data["y"].min()

    X = csr_matrix((data["total"].values, (x, y)), shape=shape, dtype=np.uint16)[y_min:, x_min:]

    fil = path2+key+"_RNA.tif"

    cv2.imwrite(fil,X.A.astype(np.uint8))