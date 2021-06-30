data<-read.table("output.homoplasies.txt")

min_big_cluster_size<-5
min_other_occurences<-2
min_other_samples<-10

homoplastic_test<-tapply(data$V3, data$V1,function(x){
  if (
    (max(x) >= min_big_cluster_size) &&
    (length(x)-1 >= min_other_occurences) &&
    (sum(x)-max(x) >= min_other_samples)
    ){
      TRUE
    } else {
      FALSE
    }
})

homoplasic_pos<-as.numeric(names(homoplastic_test[which(homoplastic_test)]))

write.table(homoplasic_pos,"selected_homoplastic_positions.txt",row.names = F,col.names = F,quote = F)
