download all files 
unzip -x alloutput.zip   (directory output: all clustering outputs ) 
unzip -x allpathway.zip  (directory pathway:  pathway predictions ) 
unzip -x allmolfunc.zip  (directory molfunction: molecular function predictions ) 
unzip -x allgenelst.zip  (directories: genelst  genelstwithinbase ) 
unzip -x allrandtestseed.zip (directory: randseed  ) 
unzip -x allrandtestshuffle.zip (directory: randshuffle ) 

For pathway computation, see pathway/readme.txt
For GO function computation, see molfunction/readme.txt
For random number seed change test see randseed
For random shuffle test see randshuffle

For disease list see "all.lst"
For 2000 gene base vector see "ALL_enrichmoa_2000.lst"
For disease-gene list as predicted by LeMeDISCO, see genelst/*.txt 
disease genes within the 2000 base vetor "ALL_enrichmoa_2000.lst" are in genelstwithinbase
For Sankey input data:               see data/sankey_test.txt
For all disease cluster assignments" see data/all_dis_kmean_sortname_allclass.tsv

Steps to reproduce all 3 levels of  clustering: 
Running:
OS requirement:
need R package installed, fortran compile, python, tcsh shell.  

1
./scan_kmeanfea.job  generate features 
-> all_dis_kmean_fea.txt

before doing k-mean clustering & calculating Silhouette score (SC), compile the Fortran code using f77 or f95 by typing
f77 -o  silscore_avedis silscore_avedis.f. If you do not have f77 or f95 on your machine, you need to install it.

2
then doing clusering by  (default clustering number was set from 2 to 50)
./scan_loopk_alldis.job 
->dis_k.dat
1st column is cluster #
2nd column is between_SS / total_SS  (%)
3rd is summed within cluster sum of squares
last column  is SC score.

the clustering outputs are in  output/pred_dis_${i}.out  where i is the number of clusters

if in case 2 to 50 do not give you a local  peak of Silhouette score, modify the script by set "npp" to a larger value

from dis_k.dat, we see k=2 is the first peak. 

3. 
./scan_getcls.job  output/pred_dis_2.out   all_dis_kmean_fea.txt  >  all_dis_kmean_cls_2.txt 
to get the label for each disease ( represented by DOIDs ) at LEVEL 1
At this level, cluster  # 1 has only 416 and majority are infectious diseases  
whereas # 2 has 3110 members (~90% of total) and has majority of all other clinical classes, we need to further cluster 
#2 to examine subclasses.  

4.
get features  for C2 (label "2"):
./scan_getclsfea.job  output/pred_dis_2.out   all_dis_kmean_fea.txt  2  > all_dis_c2_kmean_fea.txt 
do clustering at LEVEL 2:
./scan_loopk_dis_c2.job
->dis_c2_k.dat      optimal # of clusters is k=9
and output/pred_dis_c2_${k}.out     k= number of clusters

we find that k=9 has the peak SC score. 
5. 
./scan_getcls.job  output/pred_dis_c2_9.out   all_dis_c2_kmean_fea.txt  >  all_dis_c2_kmean_cls_9.txt
to get the cluster members at LEVEL 2 
We find that after this step, cluster #9  still has  large number of memebers (1405, ~1/2 of total) 
larger than those of the rest    

6: doing LEVEL 3 clustering for LEVEL 2 cluster # 9
get features:
./scan_getclsfea.job  output/pred_dis_c2_9.out   all_dis_c2_kmean_fea.txt  9  > all_dis_c2-9_kmean_fea.txt
doing clustering and getting SC scores:
./scan_loopk_dis_c2-9.job 
-> dis_c2-9_k.dat  -> k=17 has peak SC. 
-> outputs at: output/output/pred_dis_c2-9_${k}.out  where k=cluster number
get cluster members:
./scan_getcls.job  output/pred_dis_c2-9_17.out   all_dis_c2-9_kmean_fea.txt  >  all_dis_c2-9_kmean_cls_17.txt 
->  all_dis_c2-9_kmean_cls_17.txt

