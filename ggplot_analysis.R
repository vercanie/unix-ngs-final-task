# Read data
setwd("~/projects/final/")
final_data <- read_tsv('data/02-data-with-DP-is-INDEL.tsv', col_names=T) 
colnames(final_data)

#Change first colname from #CHROM to CHROM
colnames(final_data)[1] <- "CHROM"
head(final_data)

# Task 1 - over whole genome

  ggplot(final_data, aes(QUAL)) +
  geom_histogram(binwidth=1) +
  # facet_wrap(~CHROM, ncol = 1) +
  ggtitle("Overall PHRED quality") +
  ylab("Count of variants") +
  xlab("PHRED quality")

# With filtered outliers
  final_data %>% 
    filter(QUAL < 999) %>% 
    ggplot(aes(QUAL)) +
    geom_histogram(binwidth=1) +
    # facet_wrap(~CHROM, ncol = 1) +
    ggtitle("Overall PHRED quality without outliers") +
    ylab("Count of variants") +
    xlab("PHRED quality")

# Task 1 - separated by chromosome

  final_data %>% ggplot(aes(factor(CHROM), QUAL)) + geom_boxplot()

# That looks a bit weird...
# There are two problems: 
# 1) outliers at value 999
# 2) chromosomes with "_random" in name
# I'll do the following:
# Ad 1) scale_y_log10() - as suggested
# Ad 2) keep normal chromosomes and store other chromosomes in "Other"

levels(factor(final_data$CHROM))

final_data %>% 
  mutate(CHROM_recode = case_when(grepl("_random",final_data$CHROM)==T ~ "Random",
                                  TRUE ~ CHROM)) -> final_data_2

final_data_2 %>% 
  group_by(factor(CHROM_recode)) %>% 
  tally()

final_data_2 %>% 
  ggplot(aes(factor(CHROM_recode), QUAL)) + 
  geom_boxplot() + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("PHRED quality for each chromosome") +
  ylab("PHRED quality") +
  xlab("Chromosome")

# With outliers filtered out
final_data_2 %>% 
  filter(QUAL < 999) %>% 
  ggplot(aes(factor(CHROM_recode), QUAL)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("PHRED quality for each chromosome without 999 outliers") +
  ylab("PHRED quality") +
  xlab("Chromosome")

# Task 2 - distribution of read depth whole genome
ggplot(final_data, aes(DP)) +
  geom_histogram(binwidth=1) +
  # facet_wrap(~CHROM, ncol = 1) +
  ggtitle("Overall read depth") +
  ylab("Count of variants") +
  xlab("Depth")

# Task 2 - distribution of read depth per chromosome
final_data_2 %>% 
  filter(!is.na(DP)) %>% 
  ggplot(aes(factor(CHROM_recode), DP)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Read depth for each chromosome") +
  ylab("Read depth") +
  xlab("Chromosome")

ggplot(final_data_2 %>% filter(!is.na(DP)), aes(DP)) +
  geom_histogram(binwidth=1) +
  facet_wrap(~CHROM_recode, ncol = 3) +
  ggtitle("Read depth") +
  ylab("Count of variants") +
  xlab("")

# Task 3. Distribution of PHRED qualities INDELS vs. SNPs
final_data_2 %>% 
  filter(!is.na(QUAL)) %>% 
  ggplot(aes(factor(ISINDEL), QUAL)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("PHRED quality for indels and SNPs") +
  ylab("PHRED quality") +
  xlab("Variant type")

# Task 3. Distribution of PHRED qualities INDELS vs. SNPs - each chromosome
final_data_2 %>% 
  filter(!is.na(QUAL)) %>% 
  ggplot(aes(factor(ISINDEL), QUAL)) + 
  facet_wrap(~CHROM_recode, ncol = 8) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("PHRED quality for indels and SNPs - each chromosome") +
  ylab("PHRED quality") +
  xlab("Variant type")

final_data_2 %>% 
  filter(!is.na(QUAL) & QUAL < 999) %>% 
  ggplot(aes(factor(ISINDEL), QUAL)) + 
  facet_wrap(~CHROM_recode, ncol = 8) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("PHRED quality for indels and SNPs - each chromosome, 999 outliers filtered out") +
  ylab("PHRED quality") +
  xlab("Variant type")

# Task 4. Distribution of read depth (DP) qualities INDELS vs. SNPs
final_data_2 %>% 
  filter(!is.na(DP)) %>% 
  ggplot(aes(factor(ISINDEL), DP)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Read depth for indels and SNPs") +
  ylab("Read depth") +
  xlab("Variant type")

# Task 4. Distribution of read depth (DP) qualities INDELS vs. SNPs - each chromosome
final_data_2 %>% 
  filter(!is.na(DP)) %>% 
  ggplot(aes(factor(ISINDEL), DP)) + 
  facet_wrap(~CHROM_recode, ncol = 8) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Read depth for indels and SNPs - each chromosome") +
  ylab("Read depth") +
  xlab("Variant type")

# Task 9. Correlation between PHRED and DP
library(stats)
# Correlation using non-parametric Spearman test
cor.test(final_data_2$QUAL,final_data_2$DP, method = "spearman")

final_data_2 %>% 
  filter(!is.na(DP)) %>% 
  ggplot(aes(QUAL,DP)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Read depth vs PHRED quality") +
  ylab("PHRED quality") +
  xlab("Read depth")



# Task 10. Correlation between PHRED and allele frequency
cor.test(final_data_2$QUAL,final_data_2$AF, method = "spearman")

final_data_2 %>% 
  filter(!is.na(QUAL)) %>% 
  ggplot(aes(QUAL,AF)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Allele frequency vs PHRED quality") +
  ylab("Allele frequency") +
  xlab("PHRED quality")

# Task 11. Correlation between read depth (DP) and allele frequency
cor.test(final_data_2$DP,final_data_2$AF, method = "spearman")

final_data_2 %>% 
  filter(!is.na(QUAL)) %>% 
  ggplot(aes(QUAL,AF)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Allele frequency vs PHRED quality") +
  ylab("Allele frequency") +
  xlab("PHRED quality")
