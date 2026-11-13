# q2-classo vegetation classification on corrected top-15 Atacama ASVs
qiime classo classify \
  --i-features atacama-top15-train-53-clr.qza \
  --m-y-file vegetation-train-53-metadata-q2.tsv \
  --m-y-column vegetation \
  --p-no-path \
  --p-cv \
  --p-cv-one-se \
  --p-cv-seed 1 \
  --p-cv-subsets 5 \
  --p-cv--nlam 100 \
  --p-cv-lamin 0.001 \
  --p-no-stabsel \
  --p-no-lamfixed \
  --p-no-huber \
  --p-intercept \
  --o-result atacama-top15-vegetation-classification-cv.qza

qiime classo predict \
  --i-features atacama-top15-test-53-clr.qza \
  --i-problem atacama-top15-vegetation-classification-cv.qza \
  --o-predictions atacama-top15-vegetation-classification-cv-test-predictions.qza
