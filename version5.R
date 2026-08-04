#version5------------------------------
#分兩層：
#第一層nucleotide+hilbert curve排列 用彩色上色
#第二層:ld chrmosome

#載入套件與Version 5需要的資料
library(data.table)
snp_v2_snv<-readRDS("/Users/sheena/Desktop/summerintern/snp_v2_snv.rds")
snp_v4_ld10<-readRDS("/Users/sheena/Desktop/summerintern/snp_v4_ld10.rds")

#確認兩個物件的列數與欄位
dim(snp_v2_snv) #631696 保留snv mt
dim(snp_v4_ld10) #617808  刪除mt ,y chrmosome

#snp_v5:保留chr1-22, mt 排除y----
#N=618563
snp_v5<-copy(snp_v2_snv[Chr_id!="Y"])
#先照染色體排序
chr_levels_v5<-c(as.character(1:22),"X","MT")
#再依據位置排序
snp_v5[,chr_order_v5:=match(Chr_id,chr_levels_v5)]
setorder(snp_v5,chr_order_v5,Start,probeset_id)

#QC----
#number of chr
snp_v5[,.N,by=Chr_id][order(match(Chr_id,chr_levels_v5))]
#already exclude chr y
sum(snp_v5$Chr_id=="Y")
#確認10個nucleotide genotype欄位仍存在
sample_nucleotide_cols_v5<-grep("_nucleotide$",names(snp_v5),value=TRUE)
length(sample_nucleotide_cols_v5)
sample_nucleotide_cols_v5

#把snp_v4_ld10的ld加回 ld只會在常染色體中（chr1-22）
#建立檢索
snp_v5[,row_id_v5:=.I]
#擷取PLINK明確認定屬於LD block的SNP
ld_member_map_v5<-unique(snp_v4_ld10[!is.na(LD_block_id),.(probeset_id,LD_block_id,LD_block_start,LD_block_end,LD_block_kb)],by="probeset_id")
#依probeset_id將LD block member資訊加入Version 5主表
snp_v5[ld_member_map_v5,on="probeset_id",`:=`(LD_block_id=i.LD_block_id,LD_block_start=i.LD_block_start,LD_block_end=i.LD_block_end,LD_block_kb=i.LD_block_kb)]


#將核染色體SNP表示為Start到Start的單點區間，MT不參與LD區間判定
snp_points_v5<-snp_v5[Chr_id!="MT",.(row_id_v5,Chr_id,point_start=as.integer(Start),point_end=as.integer(Start))]

#整理632個LD block的染色體與物理範圍
ld_ranges_v5<-ld_interval_map10[,.(Chr_id,LD_interval_id=LD_block_id,LD_interval_start=as.integer(BP1),LD_interval_end=as.integer(BP2))]

#設定區間鍵值後，以foverlaps找出落在每個LD block範圍內的SNP
setkey(snp_points_v5,Chr_id,point_start,point_end)
setkey(ld_ranges_v5,Chr_id,LD_interval_start,LD_interval_end)
ld_interval_hits_v5<-foverlaps(snp_points_v5,ld_ranges_v5,by.x=c("Chr_id","point_start","point_end"),by.y=c("Chr_id","LD_interval_start","LD_interval_end"),type="within",nomatch=0L)[,.(row_id_v5,LD_interval_id,LD_interval_start,LD_interval_end)]

#依row_id將LD區間資訊接回Version 5主表
snp_v5[ld_interval_hits_v5,on="row_id_v5",`:=`(LD_interval_id=i.LD_interval_id,LD_interval_start=i.LD_interval_start,LD_interval_end=i.LD_interval_end)]

#建立LD member、LD interval與MT三個0/1結構mask
snp_v5[,`:=`(LD_member_mask=as.integer(!is.na(LD_block_id)),LD_interval_mask=as.integer(!is.na(LD_interval_id)),MT_mask=as.integer(Chr_id=="MT"))]
#查看Version 5整體結構標記數量
snp_v5[,.(total_SNP=.N,LD_member_N=sum(LD_member_mask),LD_interval_N=sum(LD_interval_mask),MT_N=sum(MT_mask))]

#確認MT全部標記為MT，且不會被標記為LD
snp_v5[Chr_id=="MT",.(SNP_N=.N,LD_member_N=sum(LD_member_mask),LD_interval_N=sum(LD_interval_mask),MT_N=sum(MT_mask))]

#確認每個SNP最多只落入一個LD區間
anyDuplicated(ld_interval_hits_v5$row_id_v5)

#hilbert curve----
#step1: order of chr
chr_levels_v5<-c(as.character(1:22),"X","MT")
#step2:計算每條染色體的snp數量
head(snp_v5)
hilbert_snp<-snp_v5[,.(SNP_N=.N),by=Chr_id]

#Hilbert curve容量為4的k次方，影像邊長為2的k次方
hilbert_snp[,hilbert_order:=ceiling(log(SNP_N)/log(4))]
hilbert_snp[,tile_side:=2^hilbert_order]
hilbert_snp[,tile_capacity:=tile_side^2]


#計算每條染色體需要補多少padding及實際使用率
hilbert_snp[,padding_N:=tile_capacity-SNP_N]
hilbert_snp[,usage_rate:=round(SNP_N/tile_capacity*100,2)]

#依1–22、X、MT排序
hilbert_snp[,chr_order_v5:=match(Chr_id,chr_levels_v5)]
setorder(hilbert_snp,chr_order_v5)
hilbert_snp

#每條染色體的Hilbert容量都必須大於或等於SNP數
all(hilbert_snp$tile_capacity>=hilbert_snp$SNP_N)

#不應出現無法辨認的染色體
sum(is.na(hilbert_snp$chr_order_v5))

#確認共有1–22、X與MT，共24個區域
nrow(hilbert_snp)
#儲存已加入LD與MT mask的Version 5主表
saveRDS(snp_v5,"/Users/sheena/Desktop/snp_v5_hilbert_base.rds")

#儲存每條染色體的Hilbert tile設定
saveRDS(hilbert_snp,"/Users/sheena/Desktop/hilbert_snp.rds")

#會根據 hilbert_order，替每條染色體的每個 SNP 
#建立:1.hilbert_index_chr染色體內的一為順序
#2.hilbert_row_chr：Hilbert方形中的列座標
#3.hilbert_col_chr：Hilbert方形中的欄座標
#將Hilbert curve的一維index轉成二維x、y座標
hilbert_d2xy<-function(n,d){
  #n必須是2的次方，例如32、128或256
  if(length(n)!=1L||n<1L||bitwAnd(as.integer(n),as.integer(n-1L))!=0L){
    stop("n必須是2的次方")
  }
  
  d<-as.integer(d)
  
  #index必須介於0與n平方減1之間
  if(any(d<0L|d>=n*n)){
    stop("Hilbert index超出tile容量")
  }
  
  x<-integer(length(d))
  y<-integer(length(d))
  t<-d
  s<-1L
  
  while(s<n){
    rx<-bitwAnd(bitwShiftR(t,1L),1L)
    ry<-bitwXor(bitwAnd(t,1L),rx)
    
    rotate_index<-which(ry==0L)
    
    if(length(rotate_index)>0L){
      flip_index<-rotate_index[rx[rotate_index]==1L]
      
      if(length(flip_index)>0L){
        x[flip_index]<-s-1L-x[flip_index]
        y[flip_index]<-s-1L-y[flip_index]
      }
      
      temporary_x<-x[rotate_index]
      x[rotate_index]<-y[rotate_index]
      y[rotate_index]<-temporary_x
    }
    
    x<-x+s*rx
    y<-y+s*ry
    t<-t%/%4L
    s<-s*2L
  }
  
  data.table(x=x,y=y)
}


# step2排序：Start最小的SNP → Hilbert index 0
#下一個SNP → Hilbert index 1
#下一個SNP → Hilbert index 2
#若chr_order_v5已存在，此行會重新確認其內容
snp_v5[,chr_order_v5:=match(Chr_id,chr_levels_v5)]

#同一染色體內依Start排序；同位置時再依probeset_id排序
setorder(snp_v5,chr_order_v5,Start,probeset_id)


#step3:把每條染色體的 tile 設定接回 SNP 主表
names(hilbert_snp)
hilbert_snp
#加回snp_v5
#擷取每條染色體的Hilbert tile設定
hilbert_tile_map_v5<-hilbert_snp[,.(Chr_id,hilbert_order,tile_side,tile_capacity)]

#將Hilbert order、邊長及容量接回每個SNP
snp_v5[hilbert_tile_map_v5,on="Chr_id",`:=`(
  hilbert_order=i.hilbert_order,
  tile_side=i.tile_side,
  tile_capacity=i.tile_capacity
)]

#step4:建立每條染色體內的 Hilbert index
#每條染色體重新從0開始編號,Chr2 也會重新從 0 開始，而不會接在 Chr1 後面
snp_v5[,hilbert_index_chr:=seq_len(.N)-1L,by=Chr_id]
snp_v5[,.(Chr_id,Start,probeset_id,hilbert_index_chr)][
  ,head(.SD,3),by=Chr_id
]

#Step 5：計算每個 SNP 的 Hilbert row 和 column
#建立Hilbert座標欄位
snp_v5[,`:=`(
  hilbert_row_chr=NA_integer_,
  hilbert_col_chr=NA_integer_
)]
#逐條染色體把Hilbert index轉成row與column
for(chr_now in chr_levels_v5){
  
  #找到目前染色體在snp_v5中的列位置
  row_index_now<-which(snp_v5$Chr_id==chr_now)
  
  #取得目前染色體的tile邊長
  tile_side_now<-unique(snp_v5$tile_side[row_index_now])
  
  #確認每條染色體只有一種tile大小
  stopifnot(length(tile_side_now)==1L)
  
  #將Hilbert index轉成0-based的x、y座標
  coordinate_now<-hilbert_d2xy(
    n=tile_side_now,
    d=snp_v5$hilbert_index_chr[row_index_now]
  )
  
  #轉成R使用的1-based row與column
  set(snp_v5,i=row_index_now,j="hilbert_row_chr",value=coordinate_now$y+1L)
  set(snp_v5,i=row_index_now,j="hilbert_col_chr",value=coordinate_now$x+1L)
}
snp_v5[,.(Chr_id,Start,probeset_id,
          hilbert_index_chr,
          hilbert_row_chr,
          hilbert_col_chr,
          tile_side)][1:20]

#step6:確認連續 SNP 沿著 Hilbert curve 相鄰:Hilbert curve 最重要的特性是，相鄰 index 的格子應該共享一條邊。
#檢查連續兩個Hilbert index的曼哈頓距離是否皆為1
hilbert_adjacency_qc_v5<-snp_v5[
  order(chr_order_v5,hilbert_index_chr),
  .(
    consecutive_pair_N=.N-1L,
    all_consecutive_adjacent=all(
      abs(diff(hilbert_row_chr))+
        abs(diff(hilbert_col_chr))==1L
    )
  ),
  by=Chr_id
]

hilbert_adjacency_qc_v5[
  order(match(Chr_id,chr_levels_v5))
]
#儲存已建立Hilbert座標的Version 5資料
saveRDS(snp_v5,"/Users/sheena/Desktop/summerintern/snp_v5_hilbert_coordinates.rds")

#儲存Hilbert coordinate QC結果
saveRDS(hilbert_coordinate_qc_v5,"/Users/sheena/Desktop/summerintern/hilbert_coordinate_qc_v5.rds")

#test-chr6
#確認Hilbert座標欄位已存在
stopifnot(all(c("hilbert_row_chr","hilbert_col_chr","tile_side")%in%names(snp_v5)))

#選擇Chr6作為測試染色體
chr_test_v5<-"6"
snp_chr_test_v5<-copy(snp_v5[Chr_id==chr_test_v5])

#抓取10個nucleotide genotype欄位
sample_nucleotide_cols_v5<-grep("_nucleotide$",names(snp_v5),value=TRUE)

#先使用第一位樣本
sample_col_test_v5<-sample_nucleotide_cols_v5[1]
sample_col_test_v5

#確認Chr6只有一種tile大小->256*256
tile_side_test_v5<-unique(snp_chr_test_v5$tile_side)
tile_side_test_v5

#將A/G與G/A等價基因型統一為相同格式
normalize_nucleotide_v5<-function(x){
  x<-as.character(x)
  x[is.na(x)]<-"NoCall"
  genotype_row<-grepl("^[ACGT]/[ACGT]$",x)
  
  x[genotype_row]<-vapply(
    strsplit(x[genotype_row],"/",fixed=TRUE),
    function(z)paste(sort(z),collapse="/"),
    character(1)
  )
  
  x
}

#整理目前樣本在Chr6的基因型
snp_chr_test_v5[,genotype_normalized:=normalize_nucleotide_v5(get(sample_col_test_v5))]
sort(unique(snp_chr_test_v5$genotype_normalized))
#建立nucleotide genotype數值編碼
#0:padding
genotype_code_map_v5<-c(
  "A/A"=1L, 
  "A/C"=2L,
  "A/G"=3L,
  "A/T"=4L,
  "C/C"=5L,
  "C/G"=6L,
  "C/T"=7L,
  "G/G"=8L,
  "G/T"=9L,
  "T/T"=10L,
  "NoCall"=11L
)

#將文字基因型轉成數值
snp_chr_test_v5[,genotype_code:=unname(genotype_code_map_v5[genotype_normalized])]
#建立空白矩陣；0代表padding
#genotype
genotype_matrix_test_v5<-matrix(
  0L,
  nrow=tile_side_test_v5,
  ncol=tile_side_test_v5
)
#valid_snp
valid_snp_matrix_test_v5<-matrix(
  0L,
  nrow=tile_side_test_v5,
  ncol=tile_side_test_v5
)
#ld
ld_interval_matrix_test_v5<-matrix(
  0L,
  nrow=tile_side_test_v5,
  ncol=tile_side_test_v5
)
#mt
mt_matrix_test_v5<-matrix(
  0L,
  nrow=tile_side_test_v5,
  ncol=tile_side_test_v5
)

#建立每個SNP對應的row與column
coordinate_test_v5<-cbind(
  snp_chr_test_v5$hilbert_row_chr,
  snp_chr_test_v5$hilbert_col_chr
)

#將genotype與三種structure mask放入相同Hilbert座標
genotype_matrix_test_v5[coordinate_test_v5]<-
  snp_chr_test_v5$genotype_code

valid_snp_matrix_test_v5[coordinate_test_v5]<-1L#建立真實 SNP 的位置標記 1 = 這個位置有真實SNP/0 = padding，沒有SNP

#每個snp LD_interval_mask = 1 代表該 SNP 的物理位置落在某個 LD block 的起點與終點之間/= 0 不在任何 LD block interval 中
ld_interval_matrix_test_v5[coordinate_test_v5]<-
  snp_chr_test_v5$LD_interval_mask


#MT_mask ：1 = 這個SNP屬於粒線體染色體MT/ 0 = 不是MT SNP
mt_matrix_test_v5[coordinate_test_v5]<-snp_chr_test_v5$MT_mask



#灰階圖
#0：padding
#1：真實SNP但不在LD interval
#2：位於LD interval
#建立人類視覺化使用的structure矩陣
structure_display_test_v5<-
  valid_snp_matrix_test_v5+
  ld_interval_matrix_test_v5

table(as.vector(structure_display_test_v5))


#設定genotype顏色；第1個顏色對應code 0的padding
genotype_palette_v5<-c(
  "black",       #0 Padding
  "red",         #1 A/A
  "purple",      #2 A/C
  "yellow",      #3 A/G
  "salmon",      #4 A/T
  "blue",        #5 C/C
  "cyan",        #6 C/G
  "lightblue",   #7 C/T
  "green",       #8 G/G
  "yellowgreen", #9 G/T
  "orange",      #10 T/T
  "grey"         #11 NoCall
)

#structure灰階：padding、一般SNP、LD interval
structure_palette_v5<-c("black","grey40","white")

#將矩陣轉成影像時，讓row 1顯示在最上方
plot_hilbert_matrix_v5<-function(mat,palette,main_text){
  image(
    t(mat[nrow(mat):1,,drop=FALSE]),
    col=palette,
    breaks=seq(-0.5,length(palette)-0.5,by=1),
    axes=FALSE,
    xlab="",
    ylab="",
    asp=1,
    main=main_text
  )
}

#並排查看同一Hilbert座標下的兩種資訊
par(mfrow=c(1,2),mar=c(1,1,3,1))

plot_hilbert_matrix_v5(
  genotype_matrix_test_v5,
  genotype_palette_v5,
  paste0(chr_test_v5," — Genotype\n",sample_col_test_v5)
)

plot_hilbert_matrix_v5(
  structure_display_test_v5,
  structure_palette_v5,
  paste0(chr_test_v5," — LD structure")
)

par(mfrow=c(1,1))


#建立所有染色體
#Step 1 設定染色體順序
chr_levels_v5<-c(as.character(1:22),"X","MT")
#10 個人的nucleotide 欄位
sample_nucleotide_cols_v5<-grep("_nucleotide$",names(snp_v5),value = TRUE)

#step2確認genotype編碼
#統一A/G與G/A等價基因型
normalize_nucleotide_v5<-function(x){
  x<-as.character(x)
  x[is.na(x)]<-"NoCall"
  genotype_row<-grepl("^[ACGT]/[ACGT]$",x)
  x[genotype_row]<-vapply(strsplit(x[genotype_row],"/",fixed=TRUE),
                          function(z)paste(sort(z),collapse="/"),
                          character(1))
  x
}

#0保留給padding；真實基因型使用1–11
genotype_code_map_v5<-c(
  "A/A"=1L,"A/C"=2L,"A/G"=3L,"A/T"=4L,
  "C/C"=5L,"C/G"=6L,"C/T"=7L,
  "G/G"=8L,"G/T"=9L,"T/T"=10L,
  "NoCall"=11L
)

#step3建立單一個染色體矩陣
#指定一條染色體和一位樣本後，把該染色體的 genotype、LD、MT 與有效 SNP 位置，建立成彼此座標完全一致的矩陣。
#輸入染色體與樣本欄位，輸出對齊的5種矩陣
build_chr_matrices_v5<-function(chr_now,sample_col_now){
  dt<-copy(snp_v5[Chr_id==chr_now]) #擷取指定染色體
  if(nrow(dt)==0L)stop(paste0("找不到Chr",chr_now))
  
  tile_side_now<-unique(dt$tile_side) #前面依 SNP 數量算出的 Hilbert square 邊長
  if(length(tile_side_now)!=1L)stop(paste0("Chr",chr_now,"有多種tile_side"))
  
  #姐取指定染色體後新增 genotype_normalized /genotype_code
  dt[,genotype_normalized:=normalize_nucleotide_v5(get(sample_col_now))] #normalize_nucleotide_v5()這個函數會把 genotype 統一:A/G=G/A
  dt[,genotype_code:=unname(genotype_code_map_v5[genotype_normalized])]#把 genotype 文字轉成數字
  
  #檢查是否有未定義基因型
  if(anyNA(dt$genotype_code)){
    bad_type<-unique(dt[is.na(genotype_code),genotype_normalized])
    stop(paste0("Chr",chr_now,"有未定義基因型：",paste(bad_type,collapse=", ")))
  }
  
  #建立Hilbert座標
  coordinate_now<-cbind(dt$hilbert_row_chr,dt$hilbert_col_chr)
  
  #建立空白矩陣；0代表padding
  genotype_matrix<-matrix(0L,tile_side_now,tile_side_now)
  valid_matrix<-matrix(0L,tile_side_now,tile_side_now)
  ld_interval_matrix<-matrix(0L,tile_side_now,tile_side_now)
  ld_member_matrix<-matrix(0L,tile_side_now,tile_side_now)
  mt_matrix<-matrix(0L,tile_side_now,tile_side_now)
  
  #將同一個SNP的資訊放到相同Hilbert座標
  genotype_matrix[coordinate_now]<-dt$genotype_code
  valid_matrix[coordinate_now]<-1L
  ld_interval_matrix[coordinate_now]<-dt$LD_interval_mask
  ld_member_matrix[coordinate_now]<-dt$LD_member_mask
  mt_matrix[coordinate_now]<-dt$MT_mask
  
  list(
    genotype=genotype_matrix,
    valid=valid_matrix,
    ld_interval=ld_interval_matrix,
    ld_member=ld_member_matrix,
    mt=mt_matrix
  )
}

#step4 批次建立所有染色體矩陣
#建立儲存不同染色體矩陣的list
genotype_matrix_list_v5<-list()
valid_snp_matrix_list_v5<-list()
ld_interval_matrix_list_v5<-list()
ld_member_matrix_list_v5<-list()
mt_matrix_list_v5<-list()

#逐條染色體建立矩陣
for(chr_now in chr_levels_v5){
  result_now<-build_chr_matrices_v5(chr_now,sample_col_v5)
  
  genotype_matrix_list_v5[[chr_now]]<-result_now$genotype
  valid_snp_matrix_list_v5[[chr_now]]<-result_now$valid
  ld_interval_matrix_list_v5[[chr_now]]<-result_now$ld_interval
  ld_member_matrix_list_v5[[chr_now]]<-result_now$ld_member
  mt_matrix_list_v5[[chr_now]]<-result_now$mt
}

#每個list都有24個元素
length(genotype_matrix_list_v5)
names(genotype_matrix_list_v5)

#step5#比較矩陣像素數與原始SNP資料
chr_matrix_qc_v5<-rbindlist(lapply(chr_levels_v5,function(chr_now){
  dt<-snp_v5[Chr_id==chr_now]
  
  data.table(
    Chr_id=chr_now,
    SNP_N=nrow(dt),
    genotype_pixel_N=sum(genotype_matrix_list_v5[[chr_now]]>0),
    valid_pixel_N=sum(valid_snp_matrix_list_v5[[chr_now]]),
    LD_interval_pixel_N=sum(ld_interval_matrix_list_v5[[chr_now]]),
    LD_interval_expected=sum(dt$LD_interval_mask),
    LD_member_pixel_N=sum(ld_member_matrix_list_v5[[chr_now]]),
    LD_member_expected=sum(dt$LD_member_mask),
    MT_pixel_N=sum(mt_matrix_list_v5[[chr_now]]),
    MT_expected=sum(dt$MT_mask)
  )
}))

#確認每一項都與原始資料一致
chr_matrix_qc_v5[,all_match:=
                   SNP_N==genotype_pixel_N&
                   SNP_N==valid_pixel_N&
                   LD_interval_pixel_N==LD_interval_expected&
                   LD_member_pixel_N==LD_member_expected&
                   MT_pixel_N==MT_expected
]

chr_matrix_qc_v5

#保存所有染色體的共同結構矩陣
structure_matrices_v5<-list(
  valid=valid_snp_matrix_list_v5,
  ld_interval=ld_interval_matrix_list_v5,
  ld_member=ld_member_matrix_list_v5,
  mt=mt_matrix_list_v5
)

saveRDS(structure_matrices_v5,"/Users/sheena/Desktop/summerintern/structure_matrices_v5.rds")
#保存第一位樣本的24條染色體genotype矩陣
saveRDS(genotype_matrix_list_v5,paste0("/Users/sheena/Desktop/summerintern/",sample_name_v5,"_genotype_matrices_v5.rds"))

all(chr_matrix_qc_v5$all_match)

#image
#建立Version 5 genotype matrix輸出資料夾
output_dir_v5<-"/Users/sheena/Desktop/summerintern/genotype_matrices_v5"
dir.create(output_dir_v5,showWarnings=FALSE,recursive=TRUE)

#建立指定樣本、指定染色體的genotype matrix
#只會建立genotype matrix 其他的都相同
build_chr_genotype_v5<-function(chr_now,sample_col_now){
  dt<-copy(snp_v5[Chr_id==chr_now])
  if(nrow(dt)==0L)stop(paste0("找不到Chr",chr_now))
  
  tile_side_now<-unique(dt$tile_side)
  if(length(tile_side_now)!=1L)stop(paste0("Chr",chr_now,"的tile_side不唯一"))
  
  genotype_normalized<-normalize_nucleotide_v5(dt[[sample_col_now]])
  genotype_code<-unname(genotype_code_map_v5[genotype_normalized])
  
  if(anyNA(genotype_code)){
    bad_type<-unique(genotype_normalized[is.na(genotype_code)])
    stop(paste0(sample_col_now,"在Chr",chr_now,"有未定義基因型：",paste(bad_type,collapse=", ")))
  }
  
  coordinate_now<-cbind(dt$hilbert_row_chr,dt$hilbert_col_chr)
  genotype_matrix<-matrix(0L,nrow=tile_side_now,ncol=tile_side_now)
  genotype_matrix[coordinate_now]<-genotype_code
  
  genotype_matrix
}

#建立儲存QC結果的list
genotype_qc_list_v5<-list()

#逐位樣本建立24條染色體的genotype matrices
for(sample_col_now in sample_nucleotide_cols_v5){
  sample_name_now<-sub("_\\(Axiom_TPM\\).*","",sample_col_now)
  genotype_list_now<-setNames(vector("list",length(chr_levels_v5)),chr_levels_v5)
  
  for(chr_now in chr_levels_v5){
    genotype_list_now[[chr_now]]<-build_chr_genotype_v5(chr_now,sample_col_now)
    
    expected_snp_n<-snp_v5[Chr_id==chr_now,.N]
    observed_pixel_n<-sum(genotype_list_now[[chr_now]]>0)
    
    genotype_qc_list_v5[[length(genotype_qc_list_v5)+1L]]<-data.table(
      sample_id=sample_name_now,
      Chr_id=chr_now,
      expected_SNP_N=expected_snp_n,
      genotype_pixel_N=observed_pixel_n,
      pixel_match=expected_snp_n==observed_pixel_n
    )
  }
  
  #每位樣本儲存成獨立RDS檔
  saveRDS(
    genotype_list_now,
    file.path(output_dir_v5,paste0(sample_name_now,"_genotype_matrices_v5.rds"))
  )
  
  message("完成：",sample_name_now)
}

#step4整理
#合併所有樣本與染色體QC結果
genotype_qc_v5<-rbindlist(genotype_qc_list_v5)

#確認10人乘以24個染色體區域，共240列
dim(genotype_qc_v5)

#確認每張矩陣的真實genotype pixel數都等於原始SNP數
all(genotype_qc_v5$pixel_match)

#查看不一致的結果，正常應為0列
genotype_qc_v5[pixel_match==FALSE]

#列出已完成的genotype matrix檔案
genotype_files_v5<-list.files(
  output_dir_v5,
  pattern="_genotype_matrices_v5\\.rds$",
  full.names=TRUE
)

length(genotype_files_v5)
basename(genotype_files_v5)

saveRDS(genotype_qc_v5,"/Users/sheena/Desktop/summerintern/genotype_matrix_qc_v5.rds")
fwrite(genotype_qc_v5,"/Users/sheena/Desktop/summerintern/genotype_matrix_qc_v5.csv")

#24 個染色體 Hilbert tiles 組成一張固定位置的 whole-genome chromosome atlas
#step1
#固定24個染色體的排列順序
chr_levels_v5<-c(as.character(1:22),"X","MT")

#每列放6條染色體，共4列
atlas_ncol_v5<-6L
atlas_nrow_v5<-4L

#每個格子使用最大Hilbert tile邊長，染色體間保留8 pixels
cell_side_v5<-max(snp_v5$tile_side)
gap_v5<-8L

#建立每條染色體在atlas中的列與欄
atlas_layout_v5<-data.table(
  Chr_id=chr_levels_v5,
  atlas_index=seq_along(chr_levels_v5)
)

atlas_layout_v5[,grid_row:=(atlas_index-1L)%/%atlas_ncol_v5+1L]
atlas_layout_v5[,grid_col:=(atlas_index-1L)%%atlas_ncol_v5+1L]

#取得各染色體原本的Hilbert tile大小
tile_size_map_v5<-unique(snp_v5[,.(Chr_id,tile_side)])
atlas_layout_v5[tile_size_map_v5,on="Chr_id",tile_side:=i.tile_side]

#step2 計算每條染色體位置
#較小的tile放在256×256格子的中央，不進行放大
atlas_layout_v5[,row_start:=
                  (grid_row-1L)*(cell_side_v5+gap_v5)+
                  1L+(cell_side_v5-tile_side)%/%2L]

atlas_layout_v5[,col_start:=
                  (grid_col-1L)*(cell_side_v5+gap_v5)+
                  1L+(cell_side_v5-tile_side)%/%2L]

atlas_layout_v5[,row_end:=row_start+tile_side-1L]
atlas_layout_v5[,col_end:=col_start+tile_side-1L]

#計算整張atlas的高度與寬度
atlas_height_v5<-atlas_nrow_v5*cell_side_v5+
  (atlas_nrow_v5-1L)*gap_v5

atlas_width_v5<-atlas_ncol_v5*cell_side_v5+
  (atlas_ncol_v5-1L)*gap_v5

atlas_layout_v5
c(atlas_height=atlas_height_v5,atlas_width=atlas_width_v5)

#step3
#把24個染色體矩陣放入固定atlas位置
build_genome_atlas_v5<-function(tile_list,fill_value=0L){
  atlas_matrix<-matrix(
    fill_value,
    nrow=atlas_height_v5,
    ncol=atlas_width_v5
  )
  
  for(chr_now in chr_levels_v5){
    location_now<-atlas_layout_v5[Chr_id==chr_now]
    tile_now<-tile_list[[chr_now]]
    
    #確認矩陣大小符合配置表
    if(is.null(tile_now))stop(paste0("缺少Chr",chr_now,"的矩陣"))
    if(nrow(tile_now)!=location_now$tile_side||
       ncol(tile_now)!=location_now$tile_side){
      stop(paste0("Chr",chr_now,"矩陣大小不正確"))
    }
    
    #將該染色體tile放入atlas
    atlas_matrix[
      location_now$row_start:location_now$row_end,
      location_now$col_start:location_now$col_end
    ]<-tile_now
  }
  
  atlas_matrix
}

#step4
#讀取DM-002的24條染色體genotype矩陣
DM002_genotype_v5<-readRDS(
  "/Users/sheena/Desktop/summerintern/genotype_matrices_v5/DM-002_genotype_matrices_v5.rds"
)

#讀取所有樣本共同使用的結構矩陣
structure_matrices_v5<-readRDS(
  "/Users/sheena/Desktop/summerintern/structure_matrices_v5.rds"
)

#step5 建立whole genome altas
#DM-002個人genotype atlas
genotype_atlas_test_v5<-
  build_genome_atlas_v5(DM002_genotype_v5)

#所有人共同的結構atlas
valid_atlas_v5<-
  build_genome_atlas_v5(structure_matrices_v5$valid)

ld_interval_atlas_v5<-
  build_genome_atlas_v5(structure_matrices_v5$ld_interval)

ld_member_atlas_v5<-
  build_genome_atlas_v5(structure_matrices_v5$ld_member)

mt_atlas_v5<-
  build_genome_atlas_v5(structure_matrices_v5$mt)

dim(genotype_atlas_test_v5)
dim(valid_atlas_v5)
dim(ld_interval_atlas_v5)
dim(ld_member_atlas_v5)
dim(mt_atlas_v5)


#Genotype非padding像素數應等於全部SNP數
sum(genotype_atlas_test_v5>0)==nrow(snp_v5)

#Valid SNP像素數應等於全部SNP數
sum(valid_atlas_v5)==nrow(snp_v5)

#LD interval像素數應與原始資料一致
sum(ld_interval_atlas_v5)==sum(snp_v5$LD_interval_mask)

#LD member像素數應與原始資料一致
sum(ld_member_atlas_v5)==sum(snp_v5$LD_member_mask)

#MT像素數應與原始資料一致
sum(mt_atlas_v5)==sum(snp_v5$MT_mask)

#step7視覺化
#0=padding、1=一般SNP、2=LD interval、3=MT
structure_display_atlas_v5<-
  valid_atlas_v5+
  ld_interval_atlas_v5+
  2L*mt_atlas_v5

table(as.vector(structure_display_atlas_v5))

par(mfrow=c(1,2),mar=c(1,1,3,1))

plot_hilbert_matrix_v5(
  genotype_atlas_test_v5,
  genotype_palette_v5,
  "DM-002 — Whole-genome genotype atlas"
)

plot_hilbert_matrix_v5(
  structure_display_atlas_v5,
  structure_palette_v5,
  "Whole-genome LD structure atlas"
)

par(mfrow=c(1,1))

#查看hilbert 的規則----
#依Hilbert index取得Chr6前256個SNP
chr6_path_test<-snp_v5[
  Chr_id=="6"
][order(hilbert_index_chr)][1:256]

plot(
  chr6_path_test$hilbert_col_chr,
  -chr6_path_test$hilbert_row_chr,
  type="l",
  asp=1,
  axes=FALSE,
  xlab="Hilbert column",
  ylab="Hilbert row",
  main="Chr6 — First 256 SNPs on Hilbert curve"
)

points(
  chr6_path_test$hilbert_col_chr,
  -chr6_path_test$hilbert_row_chr,
  pch=16,
  cex=0.3
)
text(
  chr6_path_test$hilbert_col_chr[1:20],
  -chr6_path_test$hilbert_row_chr[1:20],
  labels=chr6_path_test$hilbert_index_chr[1:20],
  pos=3,
  cex=0.6
)
#每條染色體的設定
snp_v5[
  ,
  .(
    SNP_N=.N,
    hilbert_order=unique(hilbert_order),
    tile_side=unique(tile_side),
    tile_capacity=unique(tile_capacity)
  ),
  by=Chr_id
][order(match(Chr_id,c(as.character(1:22),"X","MT")))]

#建立4×4 Hilbert curve的16個座標
hilbert_path_4x4<-hilbert_d2xy(n=4,d=0:15)

#加入Hilbert index及轉成R使用的1-based座標
hilbert_path_4x4[,`:=`(
  hilbert_index=0:15,
  row=y+1L,
  col=x+1L
)]

hilbert_path_4x4[,.(hilbert_index,row,col)]

#畫出4×4 Hilbert curve完整路徑
plot(
  hilbert_path_4x4$col,
  -hilbert_path_4x4$row,
  type="o",
  asp=1,
  axes=FALSE,
  xlab="Hilbert column",
  ylab="Hilbert row",
  main="4 × 4 Hilbert curve"
)

#在每個位置標上Hilbert index
text(
  hilbert_path_4x4$col,
  -hilbert_path_4x4$row,
  labels=hilbert_path_4x4$hilbert_index,
  pos=3,
  cex=0.9
)

#建立8×8 Hilbert curve，共64個位置
hilbert_path_8x8<-hilbert_d2xy(n=8,d=0:63)

hilbert_path_8x8[,`:=`(
  hilbert_index=0:63,
  row=y+1L,
  col=x+1L
)]

plot(
  hilbert_path_8x8$col,
  -hilbert_path_8x8$row,
  type="l",
  asp=1,
  axes=FALSE,
  xlab="Hilbert column",
  ylab="Hilbert row",
  main="8 × 8 Hilbert curve"
)

points(
  hilbert_path_8x8$col,
  -hilbert_path_8x8$row,
  pch=16,
  cex=0.4
)

#只標示前16個index，避免文字太擁擠
text(
  hilbert_path_8x8$col[1:16],
  -hilbert_path_8x8$row[1:16],
  labels=hilbert_path_8x8$hilbert_index[1:16],
  pos=3,
  cex=0.7
)
#取得Chr6並依Hilbert index排序
chr6_path_v5<-snp_v5[
  Chr_id=="6"
][order(hilbert_index_chr)]

#先顯示前256個SNP
chr6_path_first256<-chr6_path_v5[1:256]

plot(
  chr6_path_first256$hilbert_col_chr,
  -chr6_path_first256$hilbert_row_chr,
  type="l",
  asp=1,
  axes=FALSE,
  xlab="Hilbert column",
  ylab="Hilbert row",
  main="Chr6 — First 256 SNPs along Hilbert curve"
)

points(
  chr6_path_first256$hilbert_col_chr,
  -chr6_path_first256$hilbert_row_chr,
  pch=16,
  cex=0.25
)

#標記前20個Hilbert index
text(
  chr6_path_first256$hilbert_col_chr[1:20],
  -chr6_path_first256$hilbert_row_chr[1:20],
  labels=chr6_path_first256$hilbert_index_chr[1:20],
  pos=3,
  cex=0.6
)

#可重複使用的染色體路徑函數----
#視覺化指定染色體前n個SNP的Hilbert路徑
plot_chr_hilbert_path_v5<-function(chr_now,n_show=256L,label_n=20L){
  dt<-snp_v5[Chr_id==chr_now][order(hilbert_index_chr)]
  
  if(nrow(dt)==0L){
    stop(paste0("找不到Chr",chr_now))
  }
  
  n_show<-min(as.integer(n_show),nrow(dt))
  label_n<-min(as.integer(label_n),n_show)
  dt_show<-dt[1:n_show]
  
  plot(
    dt_show$hilbert_col_chr,
    -dt_show$hilbert_row_chr,
    type="l",
    asp=1,
    axes=FALSE,
    xlab="Hilbert column",
    ylab="Hilbert row",
    main=paste0("Chr",chr_now,
                " — First ",n_show,
                " SNPs on Hilbert curve")
  )
  
  points(
    dt_show$hilbert_col_chr,
    -dt_show$hilbert_row_chr,
    pch=16,
    cex=0.25
  )
  
  if(label_n>0L){
    text(
      dt_show$hilbert_col_chr[1:label_n],
      -dt_show$hilbert_row_chr[1:label_n],
      labels=dt_show$hilbert_index_chr[1:label_n],
      pos=3,
      cex=0.6
    )
  }
  
  #標示開始與目前顯示範圍的終點
  points(
    dt_show$hilbert_col_chr[1],
    -dt_show$hilbert_row_chr[1],
    pch=16,
    cex=1
  )
  
  text(
    dt_show$hilbert_col_chr[1],
    -dt_show$hilbert_row_chr[1],
    labels="Start",
    pos=2
  )
  
  points(
    dt_show$hilbert_col_chr[n_show],
    -dt_show$hilbert_row_chr[n_show],
    pch=17,
    cex=1
  )
  
  text(
    dt_show$hilbert_col_chr[n_show],
    -dt_show$hilbert_row_chr[n_show],
    labels="End",
    pos=4
  )
}

#Chr1前256個SNP
plot_chr_hilbert_path_v5("1",256,20)

#Chr6前1024個SNP
plot_chr_hilbert_path_v5("6",1024,20)

#MT前256個SNP
plot_chr_hilbert_path_v5("MT",256,20)

#顯示每條染色體在whole-genome atlas中的前256個Hilbert位置
plot_atlas_hilbert_paths_v5<-function(n_show_each=256L){
  plot(
    NA,
    xlim=c(1,atlas_width_v5),
    ylim=c(atlas_height_v5,1),
    asp=1,
    axes=FALSE,
    xlab="Whole-genome atlas column",
    ylab="Whole-genome atlas row",
    main=paste0(
      "Chromosome-specific Hilbert paths\n",
      "First ",n_show_each," SNPs per chromosome"
    )
  )
  
  for(chr_now in chr_levels_v5){
    dt<-snp_v5[
      Chr_id==chr_now
    ][order(hilbert_index_chr)]
    
    n_now<-min(as.integer(n_show_each),nrow(dt))
    dt_show<-dt[1:n_now]
    
    #畫該染色體Hilbert路徑
    lines(
      dt_show$atlas_col,
      dt_show$atlas_row
    )
    
    #標示該染色體名稱
    location_now<-atlas_layout_v5[Chr_id==chr_now]
    
    text(
      x=(location_now$col_start+location_now$col_end)/2,
      y=location_now$row_start-3,
      labels=paste0("Chr",chr_now),
      cex=0.7
    )
    
    #標示每條染色體的起點
    points(
      dt_show$atlas_col[1],
      dt_show$atlas_row[1],
      pch=16,
      cex=0.5
    )
  }
}
plot_atlas_hilbert_paths_v5(n_show_each=256)