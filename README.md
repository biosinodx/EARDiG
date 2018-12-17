# EARDiG
Enrichment in ARD-associated GWAS (EARDiG)
Version 1.0.0
Updated date: 2018.12.17
## Authors and License
Author: Xiao Dong
Email: biosinodx@gmail.com, xiao.dong@einstein.yu.edu
Licensed under the GNU Affero General Public License version 3 or later
## Usage
Rscript EARDiG.R workingdir projectname inputgenelist.txt backgroundgenelist.txt repeattime diseaseclass
workingdir: e.g. ./
Projectname: anyname, e.g. test
Inputgenelist.txt: file name and its root of an input gene list, should be tab limited file with header, the first colunmn is the gene names of interest.
Backgroundgenelist.txt: file name and its root of an input background gene list, should be tab limited file with header, the first colunmn is the gene names of interest, 2rd colunmn gene start, 3rd colunmn gene end; (e.g. By Xiao Dong, for protein coding genes, as in the depository "./resource/background_genelist.txt")
Repeattime: recommend 2000
Diseaseclass: file name and its root of a predefined disease - gene classes; (e.g. by Simon C Johnson as in the depository "./resource/diseasecat_simon_agingcell_2015.RData")

