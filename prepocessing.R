#檔案大用fread() 開啟
#已完成 genotype calling 的 Axiom SNP microarray 資料

install.packages("data.table")
file_path <- "/Users/sheena/Desktop/summerintern/10_genotyping_data.txt"


#680865*77
#680,865 個 SNP／probeset * 10 位受試者的genotype欄位，加上67個SNP註解與QC欄位
#每一列是一個 SNP，每一欄是樣本 genotype 或 SNP 資訊

genotype_data <- fread(file_path,
  sep = "\t",
  header = TRUE,
  skip = "probeset_id",
  quote = "",
  na.strings = c("NoCall", "---", "", "NA"),
  showProgress = TRUE)
dim(genotype_data)
head(genotype_data)

str(
  genotype_data[
    ,
    .(
      FLD,
      HomFLD,
      HetSO,
      HomRO,
      n_AA,
      n_AB,
      n_BB,
      n_NC,
      hemizygous,
      specialSNP_chr,
      gender_metrics,
      ConversionType
    )
  ]
)
unique(genotype_data$ConversionType)



#10位受試者genotype
sample_columns <- grep("\\.CEL_call_code$",names(genotype_data),value = TRUE)
sample_columns
#10個人的snp id
genotype_data[1:10,c("probeset_id", sample_columns),with = FALSE]

#Step1:整理genotype與annotation----
#10位受試者genotype:sample_columns
head(genotype_data)

#step2:檢查 genotype 與基本 QC----
#2-1.檢查10 位受試者是否只有 AA／AB／BB／NoCall =>"AA" "AB" "BB"----
genotype_value<-sort(unique(unlist(genotype_data[,..sample_columns],use.names = FALSE)))
#genotype 數量
genotype_count<-lapply(sample_columns,function(x){table(genotype_data[[x]],useNA = 'ifany')})

#call rate:(AA+AB+BB)/total
sample_names <- sub("_\\(Axiom_TPM\\)_.*\\.CEL_call_code$","",sample_columns)
sample_names
#建立total snp
sample_qc <- data.table(sample = sample_names,total_snps = nrow(genotype_data))
head(sample_qc)

#先找出nocall
sample_qc[, no_call_count := sapply(sample_columns,function(x) {sum(is.na(genotype_data[[x]]) |genotype_data[[x]] == "NoCall")})]
#計算called genotype
sample_qc[, called_count := total_snps - no_call_count]
#計算call rate
sample_qc[, call_rate := called_count / total_snps]
##計算missing rate
sample_qc[, missing_rate := no_call_count / total_snps]

sample_qc[, call_rate_percent := round(call_rate * 100,2)]
sample_qc[, missing_rate_percent := round(missing_rate * 100,2)]
sample_qc

#plot
ggplot(sample_qc,
  aes(x = reorder(sample, call_rate_percent),y = call_rate_percent)) +geom_col() +coord_flip() +
  labs(x = "sample",y = "sall rate(%)",title = "sample-level genotype call rate") +theme_minimal()


#2-2每個snp缺失比例----
#SNP missing rate =NoCall 人數/10
snp_missing_count <- rowSums(is.na(genotype_data[, ..sample_columns]) |genotype_data[, ..sample_columns] == "NoCall")
snp_qc_10samples <- data.table(probeset_id = genotype_data$probeset_id,
                               missing_count_10 = snp_missing_count,
                               missing_rate_10 = snp_missing_count / length(sample_columns))
head(snp_qc_10samples)
#缺失比例分佈
table(snp_qc_10samples$missing_count_10)

round(prop.table(table(snp_qc_10samples$missing_count_10)) * 100,2)

ggplot(snp_qc_10samples,aes(x = missing_rate_10)) +geom_histogram(binwidth = 0.1,boundary = 0) +
  labs(x = "SNP missing rate among 10 samples",y = "number of SNP",title = "SNP-level missing") +theme_minimal()

#發現13,404個SNP剛好缺失6人 進一步探討----
snps_missing_6 <- which(snp_qc_10samples$missing_count_10 == 6)
#看13,404個SNP在10位樣本中的缺失比例
missing_pattern_6 <- sapply(sample_columns,function(x) {mean(is.na(genotype_data[[x]][snps_missing_6]))})
names(missing_pattern_6) <- sample_names
missing_pattern_6



#2-3檢查chromosome,position----
chromosome_position_qc <- data.table(variable = c("hg19_chromosome","hg19_position","Chr_id","Start"),
  missing_count = c(sum(is.na(genotype_data$hg19_chromosome) |genotype_data$hg19_chromosome == ""),
    sum(is.na(genotype_data$hg19_position)),
    sum(is.na(genotype_data$Chr_id) |genotype_data$Chr_id == ""),
    sum(is.na(genotype_data$Start))))

table(genotype_data$hg19_chromosome,useNA = "ifany")
#總和=total snp
sum(table(genotype_data$hg19_chromosome, useNA = "ifany"))
#NA:85
missing_chr_snps <- genotype_data[is.na(hg19_chromosome) | hg19_chromosome == ""]
missing_chr_snps[,.(probeset_id,dbSNP_RS_ID,Chr_id,Start,hg19_chromosome,hg19_position)]
table(is.na(missing_chr_snps$Chr_id) |missing_chr_snps$Chr_id == "",useNA = "ifany")
#補值
#hg19_chromosome:NA->利用Chr_id
genotype_data[,final_chromosome := fifelse(!is.na(hg19_chromosome) & hg19_chromosome != "",hg19_chromosome,Chr_id)]
#hg19_position:NA->利用Start
genotype_data[,final_position := fifelse(!is.na(hg19_position),hg19_position,Start)]
sum(is.na(genotype_data$final_chromosome) |genotype_data$final_chromosome == "")
#新增一個欄位position_source 如果原本hg19是NA就會用Chr_id,start去補值 而position_source就會顯示Chr_id_Start
genotype_data[,position_source := fifelse(!is.na(hg19_chromosome) & !is.na(hg19_position),"hg19","Chr_id_Start")]
table(genotype_data$position_source)


#2-4檢查重複的snp----
#重複 probeset_id ->0
sum(duplicated(genotype_data$probeset_id))#0


#重複 rsID ->0
valid_rsid <- !is.na(genotype_data$dbSNP_RS_ID) &
  genotype_data$dbSNP_RS_ID != "" &
  genotype_data$dbSNP_RS_ID != "---"

sum(duplicated(genotype_data$dbSNP_RS_ID[valid_rsid]))#0


#重複 chromosome + position->4143
position_key <- paste(genotype_data$hg19_chromosome,genotype_data$hg19_position,sep = ":")
sum(duplicated(position_key))

duplicated_position <- genotype_data[duplicated(position_key) |duplicated(position_key, fromLast = TRUE),.(
    probeset_id,
    dbSNP_RS_ID,
    hg19_chromosome,
    hg19_position,
    Allele_A,
    Allele_B,
    BestProbeset,
    BestandRecommended)]
head(duplicated_position)

#2-5檢查原始240人batch-levelQC分布----
batch_qc_columns <- c(
  "CR",
  "FLD",
  "HomFLD",
  "HetSO",
  "HomRO",
  "n_AA",
  "n_AB",
  "n_BB",
  "n_NC",
  "MinorAlleleFrequency",
  "H.W.p-Value"
)

batch_qc_columns <- batch_qc_columns[
  batch_qc_columns %in% names(genotype_data)
]

summary(genotype_data[, ..batch_qc_columns])

#檢查n_AA+n_AB+n_BB+n_NC 是否等於240
str(genotype_data[, .(n_AA, n_AB, n_BB, n_NC)])
#檢查資料型態
class(genotype_data$n_AA)
class(genotype_data$n_AB)
class(genotype_data$n_BB)
class(genotype_data$n_NC)
#建立數值版本（不覆蓋）
n_AA_num <- as.numeric(genotype_data$n_AA)
n_AB_num <- as.numeric(genotype_data$n_AB)
n_BB_num <- as.numeric(genotype_data$n_BB)
n_NC_num <- as.numeric(genotype_data$n_NC)
#n_bb n_cc皆有非純文字的註解
bad_n_BB <- genotype_data[is.na(suppressWarnings(as.numeric(n_BB))) & !is.na(n_BB),.N,by = n_BB]
bad_n_BB
bad_n_NC <- genotype_data[is.na(suppressWarnings(as.numeric(n_NC))) & !is.na(n_NC),.N,by = n_NC]
bad_n_NC


which(names(genotype_data) %in%
        c(
          "n_AA",
          "n_AB",
          "n_BB",
          "n_NC",
          "hemizygous",
          "specialSNP_chr",
          "gender_metrics"
        ))
names(genotype_data)[45:65]

#2-6檢查 CR 分布----
summary(genotype_data$CR)
ggplot(
  genotype_data,
  aes(x = CR)
) +
  geom_histogram(bins = 50) +
  labs(
    x = "Batch-level SNP call rate",
    y = "Number of SNPs",
    title = "Distribution of SNP call rate in original batch"
  ) +
  theme_minimal()

#2-7檢查 ConversionType----
table(genotype_data$ConversionType,useNA = "ifany")
round(prop.table(table(
      genotype_data$ConversionType,
      useNA = "ifany")) * 100,2)


genotype_data[
  ConversionType == "0",
  .(
    probeset_id,
    n_AA,
    n_AB,
    n_BB,
    n_NC,
    hemizygous,
    specialSNP_chr,
    gender_metrics,
    ConversionType,
    BestProbeset,
    BestandRecommended,
    HomHet
  )
][1:10]


genotype_data[
  ConversionType == "PolyHighResolution",
  .(
    probeset_id,
    n_AA,
    n_AB,
    n_BB,
    n_NC,
    hemizygous,
    specialSNP_chr,
    gender_metrics,
    ConversionType,
    BestProbeset,
    BestandRecommended,
    HomHet
  )
][1:10]





#2-8檢查 BestProbeset 與 BestandRecommended----
table(
  genotype_data$BestProbeset,
  useNA = "ifany"
)

table(
  genotype_data$BestandRecommended,
  useNA = "ifany"
)

#2-9檢查 FLD、HomFLD、HetSO、HomRO----
summary(
  genotype_data[
    ,
    .(
      FLD,
      HomFLD,
      HetSO,
      HomRO
    )
  ]
)
ggplot(genotype_data, aes(x = FLD)) +
  geom_histogram(bins = 50) +
  theme_minimal()

ggplot(genotype_data, aes(x = HomFLD)) +
  geom_histogram(bins = 50) +
  theme_minimal()


ggplot(genotype_data, aes(x = HetSO)) +
  geom_histogram(bins = 50) +
  theme_minimal()

ggplot(genotype_data, aes(x = HomRO)) +
  geom_histogram(bins = 50) +
  theme_minimal()