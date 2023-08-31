## This code creates plots to look at different thresholds for 
## assigning guides to cells
tyler_filter <- function(x) { mean(x) + 6*sd(x) }
abovethresh <- function(x) { length(which(x>9))}
mytable2<-metadata2[,16:739]
row_ratio1 <- apply(mytable2, 1, tyler_filter)
tyler_filter2 <- function(x) { length(which(x>(mean(x) + 6*sd(x))))>0 & length(which(x>(mean(x) + 6*sd(x))))<6 & max(x)>=10 }
# Tyler's filter
row_ratio2 <- apply(mytable2, 1, tyler_filter2)
row_ratio3<-apply(mytable2,1,abovethresh)
abovethresh2 <- function(x) {length(which(x>9))>0 & length(which(x>9))<6}
# Noam's filter
row_ratio4<-apply(mytable2,1,abovethresh2)
# inds is where Noam and Tyler disagree; 33,022
inds<-which(row_ratio2!=row_ratio4)
# i is where Tyler doesn't call it but Noam does; 3,179
i<-which(row_ratio2[inds]==FALSE)
# j is where Noam doesn't call it but Tyler does; 29,843
j<-which(row_ratio4[inds]==FALSE)

myind<-inds[j[15]]; lll<-mytable2[myind,]; plot(as.numeric(lll), ylab = "nUMI"); abline(h=row_ratio1[myind], col="red"); abline(h=10, col="blue");


myind<-inds[i[1]]; lll<-mytable2[myind,];
p1<-ggplot(data.frame(x = 1:length(lll), y = as.numeric(lll)), aes(x = x, y = y)) +
geom_point() +
geom_hline(yintercept = row_ratio1[myind], color = "red") + geom_hline(yintercept = 10, color = "blue") + theme(axis.title.x = element_blank(),
axis.title.y = element_blank())

myind<-inds[i[2]]; lll<-mytable2[myind,];
p2<-ggplot(data.frame(x = 1:length(lll), y = as.numeric(lll)), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = row_ratio1[myind], color = "red") + geom_hline(yintercept = 10, color = "blue") + theme(axis.title.x = element_blank(),
                                                                                                                  axis.title.y = element_blank())
myind<-inds[i[3]]; lll<-mytable2[myind,];
p3<-ggplot(data.frame(x = 1:length(lll), y = as.numeric(lll)), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = row_ratio1[myind], color = "red") + geom_hline(yintercept = 10, color = "blue") + theme(axis.title.x = element_blank(),
                                                                                                                  axis.title.y = element_blank())
myind<-inds[i[4]]; lll<-mytable2[myind,];
p4<-ggplot(data.frame(x = 1:length(lll), y = as.numeric(lll)), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = row_ratio1[myind], color = "red") + geom_hline(yintercept = 10, color = "blue") + theme(axis.title.x = element_blank(),
                                                                                                                  axis.title.y = element_blank())
myind<-inds[i[5]]; lll<-mytable2[myind,];
p5<-ggplot(data.frame(x = 1:length(lll), y = as.numeric(lll)), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = row_ratio1[myind], color = "red") + geom_hline(yintercept = 10, color = "blue") + theme(axis.title.x = element_blank(),
                                                                                                                  axis.title.y = element_blank())
myind<-inds[i[6]]; lll<-mytable2[myind,];
p6<-ggplot(data.frame(x = 1:length(lll), y = as.numeric(lll)), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = row_ratio1[myind], color = "red") + geom_hline(yintercept = 10, color = "blue") + theme(axis.title.x = element_blank(),
                                                                                                                  axis.title.y = element_blank())
smooth_values <- function(data, window_size) {
  n <- nrow(data)
  first_pos <- data$position[1] - ( window_size / 2 )
  last_pos <- data$position[n] + ( window_size / 2 )
  smoothed_positions <- seq(first_pos, last_pos, by = 1)
  smoothed_values <- numeric(length(smoothed_positions))
  for (i in 1:length(smoothed_positions)) {
    
    start_pos <- smoothed_positions[i] - (window_size / 2)
    end_pos <- smoothed_positions[i] + (window_size / 2)
    # Filter the data within the window
    window_data <- data[data$position >= start_pos & data$position <= end_pos,]
    # Calculate the average per guide within the window
    smoothed_values[i] <- sum(window_data$value);#/length(which(data$position >= start_pos & data$position <= end_pos))
  }
  return(data.frame(positions = smoothed_positions, values = smoothed_values))
}
