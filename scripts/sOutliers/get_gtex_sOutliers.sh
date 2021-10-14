
## Download Outlier callset from v8 of GTEx 
wget https://storage.googleapis.com/gtex_analysis_v8/outlier_calls/GTEx_v8_outlier_calls.zip 

## unzip outlier zip file 
unzip GTEx_v8_outlier_calls.zip

## remove all other files but the sOutlier file
rm gtexV8.aseOutlier.stats.globalOutliers.removed.txt.gz
rm gtexV8.eOutlier.stats.globalOutliers.removed.txt.gz
rm GTEx_v8_outlier_calls.zip
rm README.md


## Move sOutlier file to the data dir
mv gtexV8.sOutlier.stats.globalOutliers.removed.txt.gz ../../data

