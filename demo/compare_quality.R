rm(list=ls())
# compare quality of data between two countries
# I am using the counts of quality scores generated from the preprocessing stage
G = read.csv("/media/windows-share/EEG/quality_overview_EEGs_gb_cleaned.csv")
N = read.csv("/media/windows-share/EEG/quality_overview_EEGs_ni_cleaned.csv")
cut = which(is.na(N$id) == TRUE)
if (length(cut) > 0) N = N[-cut,]
cut = which(is.na(G$id) == TRUE)
if (length(cut) > 0) G = G[-cut,]
# normalize scores

G$totalqc = G$closed_qc1 + G$closed_qc2 + G$closed_qc3 + G$closed_qc4
G$closed_qc1 = G$closed_qc1 / G$totalqc
G$closed_qc2 = G$closed_qc2 / G$totalqc
G$closed_qc3 = G$closed_qc3 / G$totalqc
G$closed_qc4 = G$closed_qc4 / G$totalqc

N$totalqc = N$closed_qc1 + N$closed_qc2 + N$closed_qc3 + N$closed_qc4
N$closed_qc1 = N$closed_qc1 / N$totalqc
N$closed_qc2 = N$closed_qc2 / N$totalqc
N$closed_qc3 = N$closed_qc3 / N$totalqc
N$closed_qc4 = N$closed_qc4 / N$totalqc

N$totalqu = N$closed_quality0 + N$closed_quality1 + N$closed_quality2
N$closed_quality0 = N$closed_quality0 / N$totalqu
N$closed_quality1 = N$closed_quality1 / N$totalqu
N$closed_quality2 = N$closed_quality2 / N$totalqu

G$totalqu = G$closed_quality0 + G$closed_quality1 + G$closed_quality2
G$closed_quality0 = G$closed_quality0 / G$totalqu
G$closed_quality1 = G$closed_quality1 / G$totalqu
G$closed_quality2 = G$closed_quality2 / G$totalqu

# tt_qc1 = t.test(N$closed_qc1,G$closed_qc1,paired = FALSE,conf.level = 0.05)
# tt_qc2 = t.test(N$closed_qc2,G$closed_qc2,paired = FALSE,conf.level = 0.05)
# tt_qc3 = t.test(N$closed_qc3,G$closed_qc3,paired = FALSE,conf.level = 0.05)
# tt_qc4 = t.test(N$closed_qc4,G$closed_qc4,paired = FALSE,conf.level = 0.05)
print("quality defined as qc score")
NN = c(mean(N$closed_qc1),mean(N$closed_qc2),mean(N$closed_qc3),mean(N$closed_qc4))
NNsd = c(sd(N$closed_qc1),sd(N$closed_qc2),sd(N$closed_qc3),sd(N$closed_qc4))
NN = round(NN,digits=2)
NNsd = round(NNsd,digits=2)
GG = c(mean(G$closed_qc1),mean(G$closed_qc2),mean(G$closed_qc3),mean(G$closed_qc4))
GGsd = c(sd(G$closed_qc1),sd(G$closed_qc2),sd(G$closed_qc3),sd(G$closed_qc4))
GG = round(GG,digits=2)
GGsd = round(GGsd,digits=2)

print(paste0(GG,collapse = " "))
print(paste0(GGsd,collapse = " "))
print(paste0(NN,collapse = " "))
print(paste0(NNsd,collapse = " "))

# quality
print("quality defined as artifacts detected")
NNu = c(mean(G$closed_quality0),mean(G$closed_quality1),mean(G$closed_quality2))
NNu = round(NNu,digits=2)
GGu = c(mean(N$closed_quality0),mean(N$closed_quality1),mean(N$closed_quality2))
GGu = round(GGu,digits=2)

print(paste0(GGu,collapse = " "))
print(paste0(NNu,collapse = " "))


# correlations with age?
NIMETA = read.csv("/media/windows-share/EEG/input/merged-meta-data_nigeria.csv")
GBMETA = read.csv("/media/windows-share/EEG/input/subject.id_with_meta-info__anonymized.csv")

G2 = merge(G,GBMETA,by.x="id",by.y="subject.id")
N2 = merge(N,NIMETA,by.x="id",by.y="subject.id")
