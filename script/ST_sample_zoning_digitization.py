import os
import matplotlib.pyplot as plt
import spateo as st
import dynamo as dyn
import scanpy as sc
import numpy as np
import plotly.express as px
import pandas as pd
import cv2


#####sample spilt
path1 = "./03sample_spilt/01h5ad/"
img="./03sample_spilt/00mask/"
outdir = "./03sample_spilt/02meta/"
for chip in [i [:-5] for i in os.listdir(path1)]:
    adata = sc.read_h5ad(path1 +chip+".h5ad")
    print(chip)
    adata.obsm["spatial"]=adata.obs[["x","y"]].to_numpy()
    adata.uns["__type"] = "UMI"
    img1 = img+chip+"/"
    adata.obs['sample'] = "filter"
    for key in [i [:-4] for i in os.listdir(img1)]:
        smaple = cv2.imread(img1+key+".tif",cv2.IMREAD_GRAYSCALE)
        for i in range(len(adata)):
            try:
                if smaple[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
                    adata.obs['sample'][i] = key
            except:
                pass
    st.pl.space(adata, color=['sample'], pointsize=0.2)
    adata.obs.to_csv(outdir + chip + "_smaplespilt_meta.csv")

#####sample spilt also use for eyeball area regions segmentation
####division spilt
%matplotlib inline
path1 = "./05divi_area/05sp_h5ad/"
tls = "./05divi_area/03mask/"
outdir = "./05divi_area/06divi_meta/"
for key in [i [:-5] for i in os.listdir(path1)]:
    adata = sc.read_h5ad(path1 +key+".h5ad")
    adata.obsm["spatial"]=adata.obs[["x","y"]].to_numpy()
    adata.uns["__type"] = "UMI"
    ###AD
    adata.obs['area'] = "filter"
    for s in [i [:-4] for i in sorted(os.listdir(tls+key+"/"))]:
        mask = cv2.imread(tls+key+"/"+s+".tif",cv2.IMREAD_GRAYSCALE)
        for i in range(len(adata)):
            try:
                if mask[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
                    adata.obs['area'][i] = s
            except:
                pass
    st.pl.space(adata, color=['area'], pointsize=0.2)
    adata.obs.to_csv(outdir + key + "_divi_meta.csv")
    
    
####layer
####digitization

import os
import matplotlib.pyplot as plt
import spateo as st
#import dynamo as dyn
import scanpy as sc
import numpy as np
import plotly.express as px
import pandas as pd
import cv2
#st.configuration.set_pub_style_mpltex()

%matplotlib inline
path1 = "./08digitization/01h5ad/"
outdir = "./08digitization/02meta/"
key="RIR1_1"#####use anterior or posterior
adata = sc.read_h5ad(path1 +key+".h5ad")
adata.obsm["spatial"]=adata.obs[["x","y"]].to_numpy()
adata.uns["__type"] = "UMI"
# Extract an area of interest
adata.obsm['spatial_bin20'] = adata.obsm['spatial']//25
cluster_label_image_lowres = st.dd.gen_cluster_image(adata, bin_size=1, spatial_key="spatial_bin20", cluster_key='dgdivi', show=False)
cluster_label_list = np.unique(adata[adata.obs['dgdivi'].isin(["1","2"]), :].obs["cluster_img_label"])
contours, cluster_image_close, cluster_image_contour = st.dd.extract_cluster_contours(cluster_label_image_lowres, 
                                                                                      cluster_label_list, bin_size=1, k_size=15, show=False)

px.imshow(cluster_image_contour, width=500, height=500)
import matplotlib.pyplot as plt
viridis_copy = plt.cm.viridis
plt.cm.register_cmap(name='my_viridis', cmap=viridis_copy)
# User input to specify a gridding direction
pnt_XY = (381,243)###maxcolum maxlayer
pnt_Xy = (369,239)###maxcolum minlayer 
pnt_xY = (146,189)##mincolum maxlayer
pnt_xy = (164,192)##mincolum minlayer 

# Digitize the area of interest
st.dd.digitize(
    adata=adata,
    ctrs=contours,
    ctr_idx=0,
    pnt_xy=pnt_xy,
    pnt_xY=pnt_xY,
    pnt_Xy=pnt_Xy,
    pnt_XY=pnt_XY,
    spatial_key="spatial_bin20"
)# Visualize digitized layers and columns
# Visualize digitized layers and columns
st.pl.space(adata, color=['digital_layer'], pointsize=0.2)
plt.savefig(outdir+key+"_digitization.pdf")
adata.obs.to_csv(outdir + key + "_digitization.csv")