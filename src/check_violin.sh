time paste *.fisher_p.txt.upperlim.txt > all.signal.txt

ls *.fisher_p.txt.upperlim.txt | awk -F '.' '{print $1"."$2}' > signal_list.txt

time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.txt signal_list.txt all.signal.od.violin.png


time paste *.pknorm.txt > all.signal.pkn.txt

time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.pkn.txt signal_list.txt all.signal.pkn.violin.png



time paste *.pkn.txt.z1.txt > all.signal.pkn_z1_all.txt

time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.pkn_z1_all.txt signal_list.txt all.signal.pkn_z1_all.violin.png


time paste *.pkz.txt > all.signal.pkz_violin.txt

time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.pkz_violin.txt signal_list.txt all.signal.pkz_violin.violin.png





time paste *.atacrep.fisher_p.txt.signorm.txt.marknorm.txt > all.signal.txt

ls *.atacrep.fisher_p.txt.signorm.txt.marknorm.txt | awk -F '.' '{print $1"."$2}' > signal_list.txt

time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.txt signal_list.txt all.signal.violin.png

time paste *.atacrep.fisher_p.txt.signorm.txt.marknorm.txt > all.signal.txt

ls *.atacrep.fisher_p.txt.signorm.txt.marknorm.txt | awk -F '.' '{print $1"."$2}' > signal_list.txt

time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.txt signal_list.txt all.signal.violin.png


/storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim/get_pkn.sh
