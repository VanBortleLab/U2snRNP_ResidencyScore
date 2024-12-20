#Loading needed libraries
library(dplyr)
library(stringr)

#########read in read counts (3'of introns) extracted from ENCODE eCLIP processed data
df <- read.delim("sample.bed")
df[,c(7:ncol(df))] <- df[,c(7:ncol(df))]+1
df <- df %>% mutate(element = str_extract(Index, "(?<=_)[^_/]+(?=/)"))
df.1 <- df %>%
  group_by(element) %>%
  mutate(group_id = element) %>%
  ungroup() %>%
  select(-group_id) 

df.1 <- df.1 %>%
  group_by(Index, Strand) %>%
  mutate(bin = if_else(Strand == "+", row_number(), 101 - row_number())) %>%
  ungroup()

#########Extract maximum count from 50bp downstream of the 3' intron-exon junction
max_values <- df.1 %>%
  filter(bin >= 45 & bin <= 50) %>% #5 bins down 
  group_by(Index) %>%
  summarize(across(matches("^(GSM|SF3)"), max, na.rm = TRUE))

colnames(max_values) <- c("Index", paste(colnames(max_values)[2:18], "_max_p"))



#########Calculate local lambda#########
#########Extract and round the mean of counts 450bp downstream of the 3' intron-exon junction
lambda_dos <- df.1 %>%
  filter(bin >= 1 & bin < 45) %>%
  group_by(Index) %>%
  summarize(across(matches("^(GSM|SF3)"), mean, na.rm = TRUE))
lambda_dos.2 <- lambda_dos %>%
  mutate(across(matches("^(GSM|SF3)"), ~ round(.x, 0)))
colnames(lambda_dos.2) <- c("Index", paste(colnames(lambda_dos.2)[2:ncol(lambda_dos.2)], "_lambda_dos"))

#########Extract and round the mean of counts 500bp upstream of the 3' intron-exon junction
lambda_ups <- df.1 %>%
  filter(bin > 50 & bin <= 100) %>%
  group_by(Index) %>%
  summarize(across(matches("^(GSM|SF3)"), mean, na.rm = TRUE))
lambda_ups.2 <- lambda_ups %>%
  mutate(across(matches("^(GSM|SF3)"), ~ round(.x, 0)))
colnames(lambda_ups.2) <- c("Index", paste(colnames(lambda_ups.2)[2:ncol(lambda_ups.2)], "_lambda_ups"))

#########Extract and round the mean of counts upstream and downstream of the 3' intron-exon junction
lambda_whole <- df.1 %>%
  filter((bin >= 1 & bin < 45) | (bin > 50 & bin <= 100)) %>%
  summarize(across(matches("^(GSM|SF3)"), ~ round(mean(.x, na.rm = TRUE), 0)))
#lambda_whole <- lambda_whole[,c(1:17)]
df.3 <- df.1 %>%
  left_join(max_values, by = "Index") %>%
  left_join(lambda_dos.2, by = "Index") %>%
  left_join(lambda_ups.2, by = "Index")
lambda_whole_expanded <- lambda_whole[rep(1, nrow(df.3)), ]
colnames(lambda_whole_expanded) <- c( paste(colnames(lambda_whole_expanded), "_lambda_whole"))

df.3 <- cbind(df.3,lambda_whole_expanded)
df.4 <- df.3 %>%
  filter(bin == 1) %>%
  group_by(Index) %>%
  ungroup()


#######Calculating the U2 occupancy at 3'(p.value)########
p.pois <- matrix(nrow = nrow(df.4), ncol = ncol(df)-5 )
p.pois[,ncol(p.pois)-1] <- df.4$Index
p.pois[,ncol(p.pois)] <- df.4$element
colnames(p.pois)<- c(colnames(df.4)[7:(ncol(df)-1)],"Index","element")
for (i in (ncol(df)+2):(ncol(df)+10)){
  print(i)
  for (k in 1:nrow(df.4)){
    print(k)
    max.p = as.numeric(df.4[k,i])
    lambda.d <- as.numeric(df.4[k, i+(ncol(df)-7)])
    lambda.u <- as.numeric(df.4[k, i+(ncol(df)-7)*2])
    lambda.w <- as.numeric(df.4[k, i+(ncol(df)-7)*3])
    max.lambda = pmax(lambda.d, lambda.u,lambda.w)
    p.value = ppois(max.p, max.lambda, lower.tail = F)
    p.pois[k,i-17] <- p.value
  }
}
p.pois <- as.data.frame(p.pois)
p.pois <- p.pois %>%
  mutate_at(vars(1:(ncol(p.pois)-2)), as.numeric)
