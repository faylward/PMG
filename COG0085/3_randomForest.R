#usr/bin/R

############################################
######### Random Forest Method #############
############################################

library(randomForest)
#x <- read.table(file="COG0085_data_table.txt", row.names=1, sep="\t", header=T)
x <- read.table(file="protein_annotations.txt", row.names=1, sep="\t", header=T)
test_data <- x[,3:47]

rf <- randomForest(category~., data=test_data, ntree=100, proximity=TRUE)
predicted <- as.character(rf$predicted)
final <- cbind(x, predicted)

write.table(final, file="rf_predictions.txt", sep="\t", quote=F)


fp <- dim(subset(final, final$predicted=="BH" & final$category == "NH"))[1]
fn <- dim(subset(final, final$predicted=="NH" & final$category == "BH"))[1]
tp <- dim(subset(final, final$predicted=="BH" & final$category == "BH"))[1]
tn <- dim(subset(final, final$predicted=="NH" & final$category == "NH"))[1]

sensitivity <- 100*(tp / (tp + fn))
specificity <- 100*(tn / (tn + fp))

specificity
sensitivity

## plot
library(ggplot2)
set <- final
bh <- set[which(set$predicted=="BH" & set$category == "BH"),]
fn <- set[which(set$predicted=="NH" & set$category == "BH"),]
fp <- set[which(set$predicted=="BH" & set$category == "NH"),]
nh <- set[which(set$predicted=="NH" & set$category == "NH"),]

a <- geom_point(data=bh, aes(x=COG0085_score, y=qlen), fill="blue", colour="blue", alpha=0.5)
b <- geom_point(data=fp, aes(x=COG0085_score, y=qlen), fill="red", colour="red", alpha=0.5)
c <- geom_point(data=fn, aes(x=COG0085_score, y=qlen), fill="green", colour="green4", alpha=0.5)
d <- geom_point(data=nh, aes(x=COG0085_score, y=qlen), fill="grey", colour="grey", alpha=0.5)

e <- geom_vline(xintercept=540)
f <- theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), axis.text.y = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), panel.grid.major = element_line(colour="grey80"), panel.grid.major.x = element_line(colour = "grey80"), panel.grid.major.y = element_line(colour="grey80"))

jpeg(file="COG0085_rf.jpg", height=6, width=6, units="in", quality=100, res=400)
ggplot() + d + a + b + c + e + f
dev.off()



############################################
######### Traditional Approach #############
############################################
library(ggplot2)
#set <- read.table(file="merged_COG0085_raw.txt", row.names=1, sep="\t", header=T, quote="")
set <- read.table(file="merged_COG0085_raw2.txt", row.names=1, sep="\t", header=T, quote="")

passed <- subset(set, set$score >= 540)
removed <- subset(set, set$score < 540)

tp <- sum(passed$category == "PRIMARY")
fp <- sum(passed$category != "PRIMARY")

tn <- sum(removed$category != "PRIMARY")
fn <- sum(removed$category == "PRIMARY" | removed$category == "SECO" | removed$category == "MAIN")

sensitivity <- 100*(tp / (tp + fn))
specificity <- 100*(tn / (tn + fp))

specificity
sensitivity

tp <- subset(set, set$category == "PRIMARY" & set$score >= 540)
fp <- subset(set, set$category != "PRIMARY" & set$score >= 540)
tn <- subset(set, set$category != "PRIMARY" & set$score < 540)
#fn <- subset(set, (set$category == "PRIMARY" | set$category == "SECO") & set$score < 540)
fn <- subset(set, (set$category == "PRIMARY" & set$score < 540 & set$length > 800) | (set$score < 540 & (set$category == "MAIN" | set$category == "SECO")))

a <- geom_point(data=tp, aes(x=score, y=length), fill="blue", colour="blue", alpha=0.5)
b <- geom_point(data=fp, aes(x=score, y=length), fill="red", colour="red", alpha=0.5)
c <- geom_point(data=fn, aes(x=score, y=length), fill="green", colour="green4", alpha=0.5)
d <- geom_point(data=tn, aes(x=score, y=length), fill="grey", colour="grey", alpha=0.5)

e <- geom_vline(xintercept=540)
f <- theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), axis.text.y = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), panel.grid.major = element_line(colour="grey80"), panel.grid.major.x = element_line(colour = "grey80"), panel.grid.major.y = element_line(colour="grey80"))

jpeg(file="COG0085_traditional.jpg", height=6, width=6, units="in", quality=100, res=400)
ggplot() + d + a + b + c + e + f
dev.off()













test_data <- x[,4:13]
rf <- randomForest(class~., data=test_data, ntree=100, proximity=TRUE)

table(predict(rf), test_data$class)

