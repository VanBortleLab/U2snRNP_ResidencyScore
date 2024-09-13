#Loading needed libraries
library(dplyr)


#########read in read counts (3'/5' of introns) extracted from ENCODE eCLIP, process data
df <- read.delim("intron_3p_read.txt") #read in intron_5p_read.txt for calculating score at 5'
df[,c(7:18)] <- df[,c(7:18)]+1
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

#########Extract maximum count from a ± 70 bp range flanking 3'/5' of the intron
max_values <- df.1 %>%
  filter(bin >= 43 & bin <= 57) %>% #7 bins down and 7 bins up
  group_by(Index) %>%
  summarize(across(matches("^(GSM|SF3)"), max, na.rm = TRUE))

colnames(max_values) <- c("Index", paste(colnames(max_values)[2:18], "_max_p"))
df.2 <- df.1 %>%
  left_join(max_values, by = "Index")



#########Extract and round the mean of counts 430bp downstream of the intron 3'/5' region
lambda_dos <- df.1 %>%
  filter(bin >= 1 & bin < 43) %>%
  group_by(Index) %>%
  summarize(across(matches("^(GSM|SF3)"), mean, na.rm = TRUE))
lambda_dos.2 <- lambda_dos %>%
  mutate(across(matches("^(GSM|SF3)"), ~ round(.x, 0)))
colnames(lambda_dos.2) <- c("Index", paste(colnames(lambda_dos.2)[2:18], "_lambda_dos"))

#########Extract and round the mean of counts 430bp upstream of the intron 3'/5' region
lambda_ups <- df.1 %>%
  filter(bin > 57 & bin <= 100) %>%
  group_by(Index) %>%
  summarize(across(matches("^(GSM|SF3)"), mean, na.rm = TRUE))
lambda_ups.2 <- lambda_ups %>%
  mutate(across(matches("^(GSM|SF3)"), ~ round(.x, 0)))
colnames(lambda_ups.2) <- c("Index", paste(colnames(lambda_ups.2)[2:18], "_lambda_ups"))

#########Extract and round the mean of counts 430bp upstream and downstream of the intron 3'/5' region
lambda_whole <- df.2 %>%
  filter((bin >= 1 & bin < 43) | (bin > 57 & bin <= 100)) %>%
  summarize(across(matches("^(GSM|SF3)"), ~ round(mean(.x, na.rm = TRUE), 0)))
lambda_whole <- lambda_whole[,c(1:17)]
df.3 <- df.2 %>%
  left_join(lambda_dos.2, by = "Index") %>%
  left_join(lambda_ups.2, by = "Index")
lambda_whole_expanded <- lambda_whole[rep(1, nrow(df.3)), ]
colnames(lambda_whole_expanded) <- c( paste(colnames(lambda_whole_expanded)[1:17], "_lambda_whole"))
df.3 <- df.2 %>%
  left_join(lambda_dos.2, by = "Index") %>%
  left_join(lambda_ups.2, by = "Index")
df.3 <- cbind(df.3,lambda_whole_expanded)
df.4 <- df.3 %>%
  filter(bin == 1) %>%
  group_by(Index) %>%
  ungroup()


#######Calculating the U2 occupancy at 3'/5' (p.value)
p.pois <- matrix(nrow = nrow(df.4), ncol = 20)
p.pois[,18] <- df.4$Index
p.pois[,19] <- df.4$element
p.pois[,20] <- df.4$response
colnames(p.pois)<- c(colnames(df.4)[7:23],"Index","element","response")
for (i in 27:43){
  print(i)
  for (k in 1:nrow(df.4)){
    print(k)
    max.p = as.numeric(df.4[k,i])
    lambda.d <- as.numeric(df.4[k, i+17])
    lambda.u <- as.numeric(df.4[k, i+34])
    lambda.w <- as.numeric(df.4[k, i+51])
    max.lambda = pmax(lambda.d, lambda.u,lambda.w)
    p.value = ppois(max.p, max.lambda, lower.tail = F)
    p.pois[k,i-26] <- p.value
  }
}
p.pois <- as.data.frame(p.pois)
p.pois <- p.pois %>%
  mutate_at(vars(1:17), as.numeric)
