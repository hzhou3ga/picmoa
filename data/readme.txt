all disease's classification: all_dis_kmean_sortname_allclass.tsv
underline data for Sankey diagram : sankey_test.txt
top 500 non-redundant disease pairs within each cluster ranked by Euclidean distance (the shorter the better):  all_top500pair.txt
randomly paired disease 2 from other clusters for each of the disease 1 in all_top500pair.txt :  all_top500pair_rand.txt

To verify that all_top500pair_rand.txt is randomly generated, run
./scan_pairwisedis_red_top500_rand_kn.job > tmp.txt
tmp.txt should be the same as all_top500pair_rand.txt 

While we cannot provide the ChatGPT script that has private key, we can provide here the ChatGPT output here:
For pairs within clusters:  allrelateanswer.zip
For random pairs : allrandanswer.zip
