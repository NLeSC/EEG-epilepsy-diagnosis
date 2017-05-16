# script to check whether standard deviations are different between countries
# statistically speaking this was a poor effort, because we are applying a lot of t-tests.
# however, this provided a quick impression (no signifiant differences between controls and epilepsy)
rm(list=ls())
graphics.off()
for (country in c("ni","gb")) { # loop through countries
  # load features per country to get labels (LAB object)
  if (country == "ni") {
    path = "/media/windows-share/EEG/EEGs_Nigeria_cleaned"
    shareddrive = "/media/windows-share/EEG/"
    load(paste0(shareddrive,"features_and_bestmodels/features/features_nigeria_4.RData"))
  } else {
    path = "/media/windows-share/EEG/EEGs_Guinea-Bissau_cleaned"
    shareddrive = "/media/windows-share/EEG/"
    metadatafile = paste0(shareddrive,"/merged-meta-data_nigeria.csv")
    load(paste0(shareddrive,"features_and_bestmodels/features/features_ginneabissau_4.RData"))
  }
  print(country)
  fn_full = dir(path,full.names = TRUE)
  fn = dir(path,full.names = FALSE)
  collect_mean = collect_sd = matrix(NA,length(fn_full),14)
  diagnosis = rep(NA,length(fn_full))
  for (i in 1:length(fn_full)) { # looped through cleaned EE data
    if (i/10 == round(i/10)) cat(paste0(i," "))
    x = read.csv(fn_full[i])
    collect_mean[i,] = apply(x[,2:15],2,mean) # <= maybe need to add weighting by dur variable here
    collect_sd[i,] = apply(x[,2:15],2,mean)
    diagnosis[i] = LAB$diagnosis[which(rownames(LAB) == fn[i])]
  }
  # print(paste0(mean(collect_sd[which(diagnosis == 1),]),mean(collect_sd[which(diagnosis == 2),])))
  # Is sd significantly different between diagnosis 1 and 2?
  cnt = 0
  for (k in 1:14) { # loop through 14 channels to compare different diagnosis
   tt = t.test(collect_sd[which(diagnosis == 1),k],collect_sd[which(diagnosis == 2),k])
   if (tt$p.value < 0.05) cnt = cnt + 1
  }
  print(cnt)
}