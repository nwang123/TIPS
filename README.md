TIPS
===
a novel pathway-guided transcriptome-wide association studies method to reveal biological processes underlying complex polygenic traits

Installation
===
To install the `TIPS` package, you will first need to install `devtools` package and then execute the following code: 
```
#install.packages("devtools")
library(devtools)
install_github("nwang123/TIPS")
```
Usage
===========
The following help page will also provide quick references for TIPS package and the example command lines:
```
library(TIPS)
```

Data
===========
There are two data: heart data and brain data. Since the size of the GWAS dataset is too large, we split it to chunks in order to save computation time. There are limitation of files size in Github, so only part of the data is uploaded. If you are interested in the complete data file, please refer to <wangneng7877@gmail.com>.

```
library(TIPS)
data(y_gene_heart_chu1)
data(w1_heart_chu1)
data(z_heart)
```
Since there are .zip file for w2_heart_chu1, the best way to import it is to download the .zip file and unzip it. Then you can import the data in R.
