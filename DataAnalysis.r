setwd(file.path("Desktop","ParallelApplication","HW1"))
library(stringi,stringr)

mmul1 <- read.table("blockingnaive.txt", header = FALSE, sep = "", dec = ".")
Mulnaive <- data.frame(c(mmul1[,c(2,4,6)]))
Mulnaive2 <- cbind("Naive Code",Mulnaive)
colnames(Mulnaive2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")


setwd(file.path("Desktop","ParallelApplication","HW1"))

library(stringi,stringr)

mmul1 <- read.table("blockingnaive.txt", header = FALSE, sep = "", dec = ".")
Mulnaive <- data.frame(c(mmul1[,c(2,4,6)]))
Mulnaive2 <- cbind("Naive Code",Mulnaive)
colnames(Mulnaive2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked1 <- read.table("blocked2_32_256.txt", header = FALSE, sep = "", dec = ".")
blocked2_32_256 <- data.frame(c(blocked1[,c(2,4,6)]))
blocked2_32_256_2 <- cbind("Blocked 2_32_256", blocked2_32_256)
colnames(blocked2_32_256_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked2 <- read.table("blocked2_64_256.txt", header = FALSE, sep = "", dec = ".")
blocked2_64_256 <- data.frame(c(blocked2[,c(2,4,6)]))
blocked2_64_256_2 <- cbind("Blocked 2_64_256", blocked2_64_256)
colnames(blocked2_64_256_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked3 <- read.table("blocked4_32_256.txt", header = FALSE, sep = "", dec = ".")
blocked4_32_256 <- data.frame(c(blocked3[,c(2,4,6)]))
blocked4_32_256_2 <- cbind("Blocked 4_32_256", blocked4_32_256)
colnames(blocked4_32_256_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked4 <- read.table("blocked4_32_512.txt", header = FALSE, sep = "", dec = ".")
blocked4_32_512 <- data.frame(c(blocked4[,c(2,4,6)]))
blocked4_32_512_2 <- cbind("Blocked 4_32_512", blocked4_32_512)
colnames(blocked4_32_512_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked5 <- read.table("blocked4_64_512.txt", header = FALSE, sep = "", dec = ".")
blocked4_64_256 <- data.frame(c(blocked5[,c(2,4,6)]))
blocked4_64_256_2 <- cbind("Blocked 4_64_256", blocked4_64_256)
colnames(blocked4_64_256_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked6 <- read.table("blocked4_64_512.txt", header = FALSE, sep = "", dec = ".")
blocked4_64_512 <- data.frame(c(blocked6[,c(2,4,6)]))
blocked4_64_512_2 <- cbind("Blocked 4_64_512", blocked2_64_256)
colnames(blocked4_64_512_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked7 <- read.table("blocked4_80_256.txt", header = FALSE, sep = "", dec = ".")
blocked4_80_256 <- data.frame(c(blocked6[,c(2,4,6)]))
blocked4_80_256_2 <- cbind("Blocked 4_80_256", blocked4_80_256)
colnames(blocked4_80_256_2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")

blocked8 <- read.table("blocked4_80_512.txt", header = FALSE, sep = "", dec = ".")
Mulblas<- data.frame(c(mmul3[,c(2,4,6)]))
Mulblas2 <- cbind("BLAS", Mulblas)
colnames(Mulblas2) <- c("Method","Matrix_Dimension","MFlops_s","Percentage")


library(tidyr)
blocks <- rbind(blocked2_32_256_2,blocked2_64_256_2,blocked4_32_256_2,blocked4_32_512_2,
                blocked4_64_256_2,blocked4_64_512_2,blocked4_80_256_2,blocked4_80_512_2)
multotal <- rbind(Mulnaive2,blocked4_80_256_2,Mulblas2)

library(ggplot2);library(grid);library(gridExtra);library(plyr);library(dplyr);library(scales);
TTT1 <- ggplot(blocks, aes(x = Matrix_Dimension, y = Percentage, colour=Method, group=Method))
TTT2 <- TTT1 + geom_point() + geom_line();
TTT2
dev.copy2pdf(out.type="pdf", file = "BlockTuning.pdf" ) ## Copy my plot to a pdf file
dev.off() ## closing the pdf device!.
TTT3 <- ggplot(multotal, aes(x = Matrix_Dimension, y = Percentage, colour=Method, group=Method))
TTT4 <- TTT3 + geom_point() + geom_line();
TTT4
dev.copy2pdf(out.type="pdf", file = "BlockedNaiveBLAS.pdf" ) ## Copy my plot to a pdf file
dev.off() ## closing the pdf device!.


grid.arrange(TTTTTT2,TTTTTT4,TTTTTT6,TTTTTT8, ncol = 2)

xval=rbind(4,17,4,4,8,8,13,17,8,13);
yval=rbind(5,5,11,16,16,11,11,11,5,5);
