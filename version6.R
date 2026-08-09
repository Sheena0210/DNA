#version6

#snp_v2_snv:631696
library(data.table)
snp_v2_snv<-readRDS("/Users/sheenazhengmei/Desktop/summerintern/snp_v2_snv.rds")

#snp_v6_del_y先排除chromosome y(13133):618563
snp_v6_del_y<-copy(snp_v2_snv[Chr_id!="Y"])


#Version 6A----
#step1:依 Start 建立 1 Mb windows----
#設定每個window大小為1,000,000 bp
window_size_bp_v6A<-1000000L

snp_v6_del_y[,Start_num:=as.numeric(Start)]
snp_v6_del_y[,chr_order_v6A:=match(Chr_id,c(as.character(1:22),"X","MT"))]

#計算每個snp所屬的window
snp_v6_del_y[,window_id_v6A:=
               as.integer(floor((Start_num-1)/window_size_bp_v6A)+1L)]

#起點終點
snp_v6_del_y[,window_bp_start_v6A:=(window_id_v6A-1L)*window_size_bp_v6A+1L]
snp_v6_del_y[,window_bp_end_v6A:=window_id_v6A*window_size_bp_v6A]

#window name
snp_v6_del_y[,window_key_v6A:=paste0("Chr",Chr_id,"_W",window_id_v6A)]
#排序
setorder(snp_v6_del_y,chr_order_v6A,Start_num,probeset_id)


#step2:統計每個window的snp----
window_snp_count_v6A<-snp_v6_del_y[,.(SNP_n=.N,observed_start_bp=min(Start_num),
                                      observed_end_bp=max(Start_num),
                                      window_bp_start_v6A=first(window_bp_start_v6A),
                                      window_bp_end_v6A=first(window_bp_end_v6A)),
                                   by=.(chr_order_v6A,Chr_id,window_id_v6A,window_key_v6A)]
#依染色體與window順序排列
setorder(window_snp_count_v6A,chr_order_v6A,window_id_v6A)

#查看前20個window
window_snp_count_v6A[1:20]

#全基因組中單一1 Mb window的最大SNP數:3551->hilbert curve:64*64
max_snp_per_window_v6A<-max(window_snp_count_v6A$SNP_n)
#列出所有並列最大SNP數的window
max_snp_window_v6A<-window_snp_count_v6A[SNP_n==max_snp_per_window_v6A]

max_snp_window_v6A

#每條染色體內SNP數最多的1 Mb window
max_window_each_chr_v6A<-window_snp_count_v6A[
  SNP_n==max(SNP_n),
  by=.(chr_order_v6A,Chr_id)
]

setorder(max_window_each_chr_v6A,chr_order_v6A)
max_window_each_chr_v6A


sum(window_snp_count_v6A$SNP_n)
nrow(snp_v6_del_y)

sum(window_snp_count_v6A$SNP_n)==nrow(snp_v6_del_y)

summary(window_snp_count_v6A$SNP_n)

max_snp_per_window_v6A


#step3 64 × 64 Hilbert template----

#設定Hilbert方格邊長與容量
hilbert_side_v6A<-64L
hilbert_capacity_v6A<-hilbert_side_v6A^2

#Version 6A Step 3：從頭建立Hilbert curve函數----------------------------

#旋轉或翻轉Hilbert子區塊
hilbert_rotate_v6A<-function(side,x,y,rx,ry){
  if(ry==0L){
    if(rx==1L){
      x<-side-1L-x
      y<-side-1L-y
    }
    temp<-x
    x<-y
    y<-temp
  }
  list(x=x,y=y)
}

#將一維Hilbert index轉成二維x、y座標
hilbert_d2xy_v6A<-function(side,index){
  
  #檢查side必須是2的次方，例如8、16、32、64
  if(side<=0L||bitwAnd(side,side-1L)!=0L){
    stop("side必須是2的次方")
  }
  
  #檢查index是否位於方格容量範圍內
  if(any(index<0L)||any(index>=side^2)){
    stop("Hilbert index超出方格容量")
  }
  
  index<-as.integer(index)
  x_result<-integer(length(index))
  y_result<-integer(length(index))
  
  for(i in seq_along(index)){
    x<-0L
    y<-0L
    temp_index<-index[i]
    current_side<-1L
    
    while(current_side<side){
      
      #決定目前Hilbert子區塊的方向
      rx<-bitwAnd(bitwShiftR(temp_index,1L),1L)
      ry<-bitwAnd(bitwXor(temp_index,rx),1L)
      
      #依方向旋轉或翻轉座標
      rotated<-hilbert_rotate_v6A(
        side=current_side,
        x=x,
        y=y,
        rx=rx,
        ry=ry
      )
      
      x<-rotated$x+current_side*rx
      y<-rotated$y+current_side*ry
      
      #每次處理Hilbert index的兩個bit
      temp_index<-bitwShiftR(temp_index,2L)
      current_side<-current_side*2L
    }
    
    x_result[i]<-x
    y_result[i]<-y
  }
  
  data.table(
    hilbert_index_v6A=index,
    x=x_result,
    y=y_result
  )
}

#設定所有1 Mb windows共用的Hilbert方格
hilbert_side_v6A<-64L
hilbert_capacity_v6A<-hilbert_side_v6A^2

#確認容量足以容納最大window的3551個SNP
hilbert_capacity_v6A
hilbert_capacity_v6A>=max_snp_per_window_v6A

#建立index 0–4095的共同Hilbert座標
hilbert_template_v6A<-hilbert_d2xy_v6A(
  side=hilbert_side_v6A,
  index=0:(hilbert_capacity_v6A-1L)
)

#將0-based座標轉為R matrix使用的1-based row與column
hilbert_template_v6A[,`:=`(
  patch_row_v6A=y+1L,
  patch_col_v6A=x+1L
)]

#只保留後續需要的欄位
hilbert_template_v6A<-hilbert_template_v6A[
  ,
  .(
    hilbert_index_v6A,
    patch_row_v6A,
    patch_col_v6A
  )
]

#step4每個window內排序並配置Hilbert座標----
#再次固定排序：染色體→ window→ Start→ probeset_id
setorder(snp_v6_del_y,chr_order_v6A,window_id_v6A,
  Start_num,probeset_id)
#每個window內，依Start由小到大建立SNP順序
snp_v6_del_y[,snp_order_within_window_v6A:=seq_len(.N),
             by=.(Chr_id,window_id_v6A)]

#Hilbert index從0開始
snp_v6_del_y[,hilbert_index_v6A:=
               snp_order_within_window_v6A-1L]
#依共同Hilbert template取得64×64 patch內的row與column
snp_v6_del_y[,patch_row_v6A:=
               hilbert_template_v6A$patch_row_v6A[
                 match(hilbert_index_v6A,
                       hilbert_template_v6A$hilbert_index_v6A)
               ]]

snp_v6_del_y[,patch_col_v6A:=
               hilbert_template_v6A$patch_col_v6A[
                 match(hilbert_index_v6A,
                       hilbert_template_v6A$hilbert_index_v6A)
               ]]



#start5建立每格snp的chromosome----

#資料中最大的window編號，決定整張atlas有幾個window欄位
max_window_id_v6A<-max(snp_v6_del_y$window_id_v6A,na.rm=TRUE)

#確認共有幾條染色體：1–22、X、MT，預期為24
n_chr_v6A<-uniqueN(snp_v6_del_y$chr_order_v6A)

#計算整張atlas的長與寬
atlas_height_v6A<-n_chr_v6A*hilbert_side_v6A
atlas_width_v6A<-max_window_id_v6A*hilbert_side_v6A

max_window_id_v6A
n_chr_v6A
atlas_height_v6A
atlas_width_v6A

#染色體決定垂直位置，patch_row決定染色體區塊內的位置
snp_v6_del_y[,atlas_row_v6A:=
               (chr_order_v6A-1L)*hilbert_side_v6A+patch_row_v6A]

#window決定水平位置，patch_col決定window區塊內的位置
snp_v6_del_y[,atlas_col_v6A:=
               (window_id_v6A-1L)*hilbert_side_v6A+patch_col_v6A]

#記錄每條染色體在整張atlas中的row範圍
chr_atlas_map_v6A<-unique(
  snp_v6_del_y[,.(chr_order_v6A,Chr_id)]
)

chr_atlas_map_v6A[,`:=`(
  atlas_row_start_v6A=(chr_order_v6A-1L)*hilbert_side_v6A+1L,
  atlas_row_end_v6A=chr_order_v6A*hilbert_side_v6A
)]

setorder(chr_atlas_map_v6A,chr_order_v6A)
chr_atlas_map_v6A

snp_v6_del_y[
  Chr_id=="1"&window_id_v6A<=2,
  .(
    probeset_id,
    Start_num,
    window_id_v6A,
    hilbert_index_v6A,
    patch_row_v6A,
    patch_col_v6A,
    atlas_row_v6A,
    atlas_col_v6A
  )
][1:30]


#step6 nucleotide genotype 合併成七類----
#0 = padding（下一步建立）
#1 = AA/TT
#2 = CC/GG
#3 = AG/CT
#4 = AC/GT
#5 = AT
#6 = CG
#7 = NoCall


#抓取10位受試者的nucleotide genotype欄位
sample_nucleotide_cols_v6A<-grep("_nucleotide$",names(snp_v6_del_y),value=TRUE)

#整理受試者名稱
sample_ids_v6A<-sub("_nucleotide$","",sample_nucleotide_cols_v6A)
sample_ids_v6A<-sub("_\\(Axiom_TPM\\).*","",sample_ids_v6A)

#確認欄位數與受試者名稱
length(sample_nucleotide_cols_v6A)
data.table(
  sample_id=sample_ids_v6A,
  nucleotide_col=sample_nucleotide_cols_v6A
)

#將diploid nucleotide genotype合併成七類
collapse_nucleotide_v6A<-function(x){
  x<-toupper(gsub("\\s+","",as.character(x)))
  x<-gsub("\\|","/",x)
  
  fcase(
    x%chin%c("A/A","AA","T/T","TT"),"AA/TT",
    x%chin%c("C/C","CC","G/G","GG"),"CC/GG",
    x%chin%c("A/G","G/A","AG","GA","C/T","T/C","CT","TC"),"AG/CT",
    x%chin%c("A/C","C/A","AC","CA","G/T","T/G","GT","TG"),"AC/GT",
    x%chin%c("A/T","T/A","AT","TA"),"AT",
    x%chin%c("C/G","G/C","CG","GC"),"CG",
    default="NoCall"
  )
}

#設定七類genotype的數值編碼；0保留給padding
genotype_code_map_v6A<-c(
  "AA/TT"=1L,
  "CC/GG"=2L,
  "AG/CT"=3L,
  "AC/GT"=4L,
  "AT"=5L,
  "CG"=6L,
  "NoCall"=7L
)

#逐位受試者新增七類文字欄位與數值欄位
for(i in seq_along(sample_nucleotide_cols_v6A)){
  genotype7_now<-collapse_nucleotide_v6A(
    snp_v6_del_y[[sample_nucleotide_cols_v6A[i]]]
  )
  
  text_col_now<-paste0(sample_ids_v6A[i],"_genotype7_v6A")
  code_col_now<-paste0(sample_ids_v6A[i],"_genotype_code_v6A")
  
  set(snp_v6_del_y,j=text_col_now,value=genotype7_now)
  set(snp_v6_del_y,j=code_col_now,
      value=unname(genotype_code_map_v6A[genotype7_now]))
}

#記錄後續會使用的欄位名稱
genotype7_cols_v6A<-paste0(sample_ids_v6A,"_genotype7_v6A")
genotype_code_cols_v6A<-paste0(sample_ids_v6A,"_genotype_code_v6A")


genotype7_qc_v6A<-rbindlist(lapply(seq_along(sample_ids_v6A),function(i){
  snp_v6_del_y[,.(SNP_n=.N),
               by=.(genotype7=get(genotype7_cols_v6A[i]))][
                 ,sample_id:=sample_ids_v6A[i]
               ][,.(sample_id,genotype7,SNP_n)]
}))

setorder(genotype7_qc_v6A,sample_id,genotype7)
genotype7_qc_v6A

#step7:建立matrix----


#選擇第一位受試者
sample_index_v6A<-1L
sample_id_now_v6A<-sample_ids_v6A[sample_index_v6A]
code_col_now_v6A<-genotype_code_cols_v6A[sample_index_v6A]

#建立全為0的atlas；0代表padding
atlas_matrix_test_v6A<-matrix(as.raw(0L),
  nrow=atlas_height_v6A,
  ncol=atlas_width_v6A
)

#將每個真實SNP的genotype code填入共同座標
atlas_matrix_test_v6A[
  cbind(
    snp_v6_del_y$atlas_row_v6A,
    snp_v6_del_y$atlas_col_v6A
  )
]<-as.raw(snp_v6_del_y[[code_col_now_v6A]])

dim(atlas_matrix_test_v6A)
atlas_height_v6A
atlas_width_v6A

#非padding pixel數應等於SNP總數
real_snp_pixel_n_v6A<-sum(
  atlas_matrix_test_v6A!=as.raw(0L)
)

real_snp_pixel_n_v6A
nrow(snp_v6_del_y)

real_snp_pixel_n_v6A==nrow(snp_v6_del_y)

#padding 數量
#整張atlas總pixel數
total_pixel_n_v6A<-atlas_height_v6A*atlas_width_v6A

#沒有真實SNP的位置皆為padding
padding_pixel_n_v6A<-total_pixel_n_v6A-nrow(snp_v6_del_y)

data.table(
  total_pixel_n=total_pixel_n_v6A,
  real_SNP_pixel_n=nrow(snp_v6_del_y),
  padding_pixel_n=padding_pixel_n_v6A,
  padding_rate=padding_pixel_n_v6A/total_pixel_n_v6A
)

#依SNP座標從atlas中取回genotype code
code_readback_v6A<-as.integer(
  atlas_matrix_test_v6A[
    cbind(
      snp_v6_del_y$atlas_row_v6A,
      snp_v6_del_y$atlas_col_v6A
    )
  ]
)

#和原始genotype code比較
all(
  code_readback_v6A==
    snp_v6_del_y[[code_col_now_v6A]]
)

#原始SNP資料中的七類數量
table(snp_v6_del_y[[code_col_now_v6A]])

#從atlas真實SNP座標讀回後的七類數量
table(code_readback_v6A)

#Chr6的atlas row範圍
chr6_row_start_v6A<-(6L-1L)*hilbert_side_v6A+1L

chr6_row_end_v6A<-6L*hilbert_side_v6A

#Window 32的atlas column範圍
window32_col_start_v6A<-(32L-1L)*hilbert_side_v6A+1L

window32_col_end_v6A<-32L*hilbert_side_v6A

#擷取Chr6_W32的64×64 Hilbert patch
chr6_w32_patch_test_v6A<-atlas_matrix_test_v6A[
  chr6_row_start_v6A:chr6_row_end_v6A,
  window32_col_start_v6A:window32_col_end_v6A
]

dim(chr6_w32_patch_test_v6A) #64*64=4096(Chr6_W32 snp:3551 padding:545)
table(as.integer(chr6_w32_patch_test_v6A))


#step8----
#Chr6_W32


#指定要查看的染色體與window
view_chr_v6A<-"6"
view_window_v6A<-32L

#找出該染色體在atlas中的順序
view_chr_order_v6A<-unique(snp_v6_del_y[Chr_id==view_chr_v6A,chr_order_v6A])

#計算該patch在完整atlas中的row與column範圍
view_row_start_v6A<-(view_chr_order_v6A-1L)*hilbert_side_v6A+1L
view_row_end_v6A<-view_chr_order_v6A*hilbert_side_v6A
view_col_start_v6A<-(view_window_v6A-1L)*hilbert_side_v6A+1L
view_col_end_v6A<-view_window_v6A*hilbert_side_v6A

#擷取64×64 patch
patch_test_v6A<-atlas_matrix_test_v6A[
  view_row_start_v6A:view_row_end_v6A,
  view_col_start_v6A:view_col_end_v6A
]

dim(patch_test_v6A)

#查看0–7各類pixel數量
patch_code_count_v6A<-table(factor(as.integer(patch_test_v6A),levels=0:7))

patch_code_count_v6A

#真實SNP與padding數量
sum(patch_test_v6A!=as.raw(0L))
sum(patch_test_v6A==as.raw(0L))

#0為padding；1–7為七類genotype
genotype_color_v6A<-c(
  "white",     #0 Padding
  "#E41A1C",   #1 AA/TT
  "#377EB8",   #2 CC/GG
  "#4DAF4A",   #3 AG/CT
  "#984EA3",   #4 AC/GT
  "#FF7F00",   #5 AT
  "#A65628",   #6 CG
  "grey40"     #7 NoCall
)

genotype_label_v6A<-c(
  "Padding",
  "AA/TT",
  "CC/GG",
  "AG/CT",
  "AC/GT",
  "AT",
  "CG",
  "NoCall"
)

#轉為一般integer matrix
patch_plot_v6A<-matrix(
  as.integer(patch_test_v6A),
  nrow=hilbert_side_v6A,
  ncol=hilbert_side_v6A
)

#轉置並上下翻轉，使matrix第1列顯示在圖片上方
patch_plot_display_v6A<-t(
  patch_plot_v6A[nrow(patch_plot_v6A):1,]
)

#畫出64×64 Hilbert genotype patch
par(mar=c(2,2,4,8),xpd=TRUE)
image(
  x=1:hilbert_side_v6A,
  y=1:hilbert_side_v6A,
  z=patch_plot_display_v6A,
  col=genotype_color_v6A,
  breaks=seq(-0.5,7.5,by=1),
  axes=FALSE,
  asp=1,
  useRaster=TRUE,
  main=paste0(sample_id_now_v6A,"－Chr6 Window 32")
)
box()
legend(
  "right",
  inset=c(-0.38,0),
  legend=genotype_label_v6A,
  fill=genotype_color_v6A,
  border=NA,
  bty="n",
  cex=0.85
)

snp_v6_del_y[
  Chr_id=="6"&window_id_v6A==32,
  .(
    probeset_id,
    Start_num,
    snp_order_within_window_v6A,
    hilbert_index_v6A,
    patch_row_v6A,
    patch_col_v6A,
    genotype=get(genotype7_cols_v6A[1]),
    genotype_code=get(genotype_code_cols_v6A[1])
  )
][1:20]

chr6_w32_order_check_v6A<-snp_v6_del_y[
  Chr_id=="6"&window_id_v6A==32
]

all(diff(chr6_w32_order_check_v6A$Start_num)>=0)
all(
  chr6_w32_order_check_v6A$hilbert_index_v6A==
    0:(nrow(chr6_w32_order_check_v6A)-1L)
)

#step9 10人的matrix----

#所有SNP共同使用的atlas座標
atlas_pixel_index_v6A<-cbind(
  snp_v6_del_y$atlas_row_v6A,
  snp_v6_del_y$atlas_col_v6A
)

#建立存放10張atlas的list
atlas_matrix_all_v6A<-setNames(vector("list",length(sample_ids_v6A)),sample_ids_v6A)


#逐位受試者建立raw matrix
for(i in seq_along(sample_ids_v6A)){
  atlas_now<-matrix(
    as.raw(0L),
    nrow=atlas_height_v6A,
    ncol=atlas_width_v6A
  )
  
  atlas_now[atlas_pixel_index_v6A]<-as.raw(snp_v6_del_y[[genotype_code_cols_v6A[i]]])
  
  atlas_matrix_all_v6A[[i]]<-atlas_now
  cat(i,"/",length(sample_ids_v6A),sample_ids_v6A[i],"完成\n")
}

#每張atlas 7種genotype數量
atlas_code_count_v6A<-rbindlist(lapply(seq_along(sample_ids_v6A),function(i){
  code_count<-table(
    factor(
      as.integer(atlas_matrix_all_v6A[[i]]),
      levels=0:7
    )
  )
  
  data.table(
    sample_id=sample_ids_v6A[i],
    genotype_code=0:7,
    pixel_n=as.integer(code_count)
  )
}))

atlas_code_count_v6A
#add genotype_class name
atlas_code_count_v6A[,genotype_class:=genotype_label_v6A[genotype_code+1L]]
setcolorder(atlas_code_count_v6A,c("sample_id","genotype_code","genotype_class","pixel_n"))

saveRDS(atlas_matrix_all_v6A,file="atlas_matrix_all_v6A.rds",compress="xz")

#step10----
#Version 6A Step 10：繪製單一受試者完整atlas（最終修正版）----------------

#選擇第一位受試者
view_sample_id_v6A<-sample_ids_v6A[1]

#取出該受試者的raw atlas matrix
atlas_view_raw_v6A<-atlas_matrix_all_v6A[[view_sample_id_v6A]]

#將raw matrix轉成0–7的integer matrix
atlas_view_plot_v6A<-matrix(
  as.integer(atlas_view_raw_v6A),
  nrow=atlas_height_v6A,
  ncol=atlas_width_v6A
)

#設定padding與七類genotype顏色
genotype_color_v6A<-c(
  "white",     #0 Padding
  "#E41A1C",   #1 AA/TT
  "#377EB8",   #2 CC/GG
  "#4DAF4A",   #3 AG/CT
  "#984EA3",   #4 AC/GT
  "#FF7F00",   #5 AT
  "#A65628",   #6 CG
  "grey40"     #7 NoCall
)

#設定圖例名稱
genotype_label_v6A<-c(
  "Padding",
  "AA/TT",
  "CC/GG",
  "AG/CT",
  "AC/GT",
  "AT",
  "CG",
  "NoCall"
)

#建立染色體順序與標籤位置
chr_label_v6A<-unique(
  snp_v6_del_y[,.(chr_order_v6A,Chr_id)]
)

setorder(chr_label_v6A,chr_order_v6A)

#因為Chr1放在圖片最上方，所以需要反向計算Y軸位置
chr_label_v6A[,display_y_v6A:=
                atlas_height_v6A-
                (chr_order_v6A-0.5)*hilbert_side_v6A
]

#建立X軸window刻度：第1個window後，每20 Mb標示一次
window_tick_id_v6A<-sort(unique(c(
  1L,
  seq(20L,max_window_id_v6A,by=20L)
)))

#若最大window不是20的倍數，將最後一個window也加入
window_tick_id_v6A<-sort(unique(c(
  window_tick_id_v6A,
  max_window_id_v6A
)))

#將window編號換算成atlas中的pixel位置
window_tick_x_v6A<-
  (window_tick_id_v6A-0.5)*hilbert_side_v6A

#輸出檔案位置
atlas_png_file_v6A<-file.path(
  getwd(),
  paste0(
    view_sample_id_v6A,
    "_Version6A_full_atlas_v3.png"
  )
)

#設定圖片解析度
atlas_png_res_v6A<-120

#設定四周留白，單位為inch
plot_margin_bottom_v6A<-0.55
plot_margin_left_v6A<-0.80
plot_margin_top_v6A<-1.10
plot_margin_right_v6A<-0.15

#依atlas實際pixel數計算輸出圖片大小
atlas_png_width_v6A<-as.integer(
  atlas_width_v6A+
    (plot_margin_left_v6A+plot_margin_right_v6A)*atlas_png_res_v6A
)

atlas_png_height_v6A<-as.integer(
  atlas_height_v6A+
    (plot_margin_bottom_v6A+plot_margin_top_v6A)*atlas_png_res_v6A
)

#開啟PNG繪圖裝置
png(
  filename=atlas_png_file_v6A,
  width=atlas_png_width_v6A,
  height=atlas_png_height_v6A,
  units="px",
  res=atlas_png_res_v6A
)

#固定四周留白，xpd=NA允許圖例畫在黑框上方
par(
  mai=c(
    plot_margin_bottom_v6A,
    plot_margin_left_v6A,
    plot_margin_top_v6A,
    plot_margin_right_v6A
  ),
  xpd=NA
)

#將matrix上下反轉並轉置
#使Chr1顯示於最上方，ChrMT顯示於最下方
atlas_display_v6A<-t(
  atlas_view_plot_v6A[
    atlas_height_v6A:1L,
    ,
    drop=FALSE
  ]
)

#繪製完整atlas
image(
  x=seq(0.5,atlas_width_v6A-0.5,by=1),
  y=seq(0.5,atlas_height_v6A-0.5,by=1),
  z=atlas_display_v6A,
  col=genotype_color_v6A,
  breaks=seq(-0.5,7.5,by=1),
  xlim=c(0,atlas_width_v6A),
  ylim=c(0,atlas_height_v6A),
  axes=FALSE,
  xlab="",
  ylab="",
  xaxs="i",
  yaxs="i",
  asp=1,
  useRaster=TRUE
)

#畫每條染色體之間的水平分隔線
abline(
  h=seq(
    hilbert_side_v6A,
    atlas_height_v6A-hilbert_side_v6A,
    by=hilbert_side_v6A
  ),
  col="grey65",
  lwd=0.4
)

#畫每個1 Mb window之間的垂直分隔線
abline(
  v=seq(
    hilbert_side_v6A,
    atlas_width_v6A-hilbert_side_v6A,
    by=hilbert_side_v6A
  ),
  col="grey88",
  lwd=0.2
)

#畫X軸，只顯示window編號
axis(
  side=1,
  at=window_tick_x_v6A,
  labels=window_tick_id_v6A,
  cex.axis=0.65,
  las=1,
  lwd=0.7,
  lwd.ticks=0.7
)

#畫Y軸，只顯示24條染色體名稱
axis(
  side=2,
  at=chr_label_v6A$display_y_v6A,
  labels=paste0("Chr",chr_label_v6A$Chr_id),
  las=1,
  cex.axis=0.65,
  lwd=0.7,
  lwd.ticks=0.7
)

#加入X軸與Y軸標題
mtext(
  "1-Mb genomic window",
  side=1,
  line=2.3,
  cex=0.8
)

mtext(
  "Chromosome",
  side=2,
  line=3.5,
  cex=0.8
)

#畫atlas黑色外框
rect(
  xleft=0,
  ybottom=0,
  xright=atlas_width_v6A,
  ytop=atlas_height_v6A,
  border="black",
  lwd=1
)

#橫向圖例，放在黑框上方
legend(
  x=atlas_width_v6A/2,
  y=atlas_height_v6A+hilbert_side_v6A*0.40,
  legend=genotype_label_v6A,
  fill=genotype_color_v6A,
  border=NA,
  bty="n",
  horiz=TRUE,
  xjust=0.5,
  yjust=0.5,
  cex=0.75,
  xpd=NA
)

#加入主標題
mtext(
  paste0(
    view_sample_id_v6A,
    " - Version 6A chromosome × window Hilbert atlas"
  ),
  side=3,
  line=4.2,
  cex=1,
  font=2
)

#完成並關閉圖片
dev.off()

#確認圖片成功輸出
atlas_png_file_v6A
file.exists(atlas_png_file_v6A)
file.info(atlas_png_file_v6A)$size/1024^2

#Mac直接開啟圖片
system(paste("open",shQuote(atlas_png_file_v6A)))

#刪除僅供繪圖使用的大型matrix
rm(
  atlas_view_plot_v6A,
  atlas_display_v6A
)

gc()

#step11. `10個人----


#建立輸出資料夾
atlas_output_dir_v6A<-"Version6A_full_atlas_v3"
dir.create(atlas_output_dir_v6A,showWarnings=FALSE,recursive=TRUE)

#建立染色體標籤位置
chr_label_v6A<-unique(
  snp_v6_del_y[,.(chr_order_v6A,Chr_id)]
)
setorder(chr_label_v6A,chr_order_v6A)

#Chr1顯示在上方，ChrMT顯示在下方
chr_label_v6A[,display_y_v6A:=
                atlas_height_v6A-
                (chr_order_v6A-0.5)*hilbert_side_v6A
]

#建立X軸刻度：Window 1、每20個window，以及最後一個window
window_tick_id_v6A<-sort(unique(c(
  1L,
  seq(20L,max_window_id_v6A,by=20L),
  max_window_id_v6A
)))

#將window編號轉成atlas pixel位置
window_tick_x_v6A<-
  (window_tick_id_v6A-0.5)*hilbert_side_v6A

#設定輸出解析度與四周留白
atlas_png_res_v6A<-120
plot_margin_bottom_v6A<-0.55
plot_margin_left_v6A<-0.80
plot_margin_top_v6A<-1.10
plot_margin_right_v6A<-0.15

#依atlas實際尺寸決定PNG大小
atlas_png_width_v6A<-as.integer(
  atlas_width_v6A+
    (plot_margin_left_v6A+plot_margin_right_v6A)*atlas_png_res_v6A
)

atlas_png_height_v6A<-as.integer(
  atlas_height_v6A+
    (plot_margin_bottom_v6A+plot_margin_top_v6A)*atlas_png_res_v6A
)
#建立單一受試者完整atlas繪圖函數
plot_full_atlas_v6A<-function(sample_id,atlas_raw){
  
  #將raw matrix轉成0–7 integer matrix
  atlas_plot<-matrix(
    as.integer(atlas_raw),
    nrow=atlas_height_v6A,
    ncol=atlas_width_v6A
  )
  
  #上下反轉並轉置，使Chr1在最上方、ChrMT在最下方
  atlas_display<-t(
    atlas_plot[
      atlas_height_v6A:1L,
      ,
      drop=FALSE
    ]
  )
  
  #設定輸出檔名
  output_file<-file.path(
    atlas_output_dir_v6A,
    paste0(sample_id,"_Version6A_full_atlas_v3.png")
  )
  
  #開啟PNG繪圖裝置
  png(
    filename=output_file,
    width=atlas_png_width_v6A,
    height=atlas_png_height_v6A,
    units="px",
    res=atlas_png_res_v6A
  )
  
  #設定四周留白，允許圖例畫在黑框外
  par(
    mai=c(
      plot_margin_bottom_v6A,
      plot_margin_left_v6A,
      plot_margin_top_v6A,
      plot_margin_right_v6A
    ),
    xpd=NA
  )
  
  #繪製完整atlas
  image(
    x=seq(0.5,atlas_width_v6A-0.5,by=1),
    y=seq(0.5,atlas_height_v6A-0.5,by=1),
    z=atlas_display,
    col=genotype_color_v6A,
    breaks=seq(-0.5,7.5,by=1),
    xlim=c(0,atlas_width_v6A),
    ylim=c(0,atlas_height_v6A),
    axes=FALSE,
    xlab="",
    ylab="",
    xaxs="i",
    yaxs="i",
    asp=1,
    useRaster=TRUE
  )
  
  #染色體之間的水平分隔線
  abline(
    h=seq(
      hilbert_side_v6A,
      atlas_height_v6A-hilbert_side_v6A,
      by=hilbert_side_v6A
    ),
    col="grey65",
    lwd=0.4
  )
  
  #1 Mb window之間的垂直分隔線
  abline(
    v=seq(
      hilbert_side_v6A,
      atlas_width_v6A-hilbert_side_v6A,
      by=hilbert_side_v6A
    ),
    col="grey88",
    lwd=0.2
  )
  
  #X軸：1 Mb window編號
  axis(
    side=1,
    at=window_tick_x_v6A,
    labels=window_tick_id_v6A,
    cex.axis=0.65,
    las=1,
    lwd=0.7,
    lwd.ticks=0.7
  )
  
  #Y軸：染色體名稱
  axis(
    side=2,
    at=chr_label_v6A$display_y_v6A,
    labels=paste0("Chr",chr_label_v6A$Chr_id),
    las=1,
    cex.axis=0.65,
    lwd=0.7,
    lwd.ticks=0.7
  )
  
  #座標軸標題
  mtext(
    "1-Mb genomic window",
    side=1,
    line=2.3,
    cex=0.8
  )
  
  mtext(
    "Chromosome",
    side=2,
    line=3.5,
    cex=0.8
  )
  
  #完整atlas黑色外框
  rect(
    xleft=0,
    ybottom=0,
    xright=atlas_width_v6A,
    ytop=atlas_height_v6A,
    border="black",
    lwd=1
  )
  
  #橫向圖例放在黑框上方
  legend(
    x=atlas_width_v6A/2,
    y=atlas_height_v6A+hilbert_side_v6A*0.40,
    legend=genotype_label_v6A,
    fill=genotype_color_v6A,
    border=NA,
    bty="n",
    horiz=TRUE,
    xjust=0.5,
    yjust=0.5,
    cex=0.75,
    xpd=NA
  )
  
  #主標題
  mtext(
    paste0(
      sample_id,
      " - Version 6A chromosome × window Hilbert atlas"
    ),
    side=3,
    line=4.2,
    cex=1,
    font=2
  )
  
  #關閉PNG裝置
  dev.off()
  
  #刪除本次繪圖暫存matrix
  rm(atlas_plot,atlas_display)
  invisible(gc())
  
  return(output_file)
}
#建立儲存輸出檔案路徑的向量
atlas_png_files_v6A<-setNames(
  character(length(sample_ids_v6A)),
  sample_ids_v6A
)

#逐位受試者輸出完整atlas
for(i in seq_along(sample_ids_v6A)){
  
  sample_id_now<-sample_ids_v6A[i]
  
  atlas_png_files_v6A[i]<-plot_full_atlas_v6A(
    sample_id=sample_id_now,
    atlas_raw=atlas_matrix_all_v6A[[sample_id_now]]
  )
  
  cat(
    i,"/",length(sample_ids_v6A),
    sample_id_now,
    "完成\n"
  )
}
#整理10張圖片的輸出資訊
atlas_png_qc_v6A<-data.table(
  sample_id=sample_ids_v6A,
  file_path=unname(atlas_png_files_v6A),
  file_exists=file.exists(atlas_png_files_v6A),
  file_size_MB=round(
    file.info(atlas_png_files_v6A)$size/1024^2,
    3
  )
)

atlas_png_qc_v6A

system(paste("open",shQuote(atlas_output_dir_v6A)))



#step12 SNP-to-pixel 對照表----

#建立真實SNP的pixel對照表；padding不屬於任何SNP，因此不會出現在這張表
pixel_snp_map_v6A<-snp_v6_del_y[,.(
  
  #原始SNP資訊
  probeset_id,
  dbSNP_RS_ID,
  Chr_id,
  Start=Start_num,
  Allele_A,
  Allele_B,
  
  #1 Mb window資訊
  chr_order_v6A,
  window_id_v6A,
  window_key_v6A,
  window_bp_start_v6A,
  window_bp_end_v6A,
  
  #window內排列資訊
  snp_order_within_window_v6A,
  hilbert_index_v6A,
  patch_row_v6A,
  patch_col_v6A,
  
  #完整atlas座標
  atlas_row_v6A,
  atlas_col_v6A
)]

#固定對照表順序
setorder(
  pixel_snp_map_v6A,
  chr_order_v6A,
  window_id_v6A,
  snp_order_within_window_v6A
)

#查看前20筆
pixel_snp_map_v6A[1:20]
saveRDS(
  pixel_snp_map_v6A,
  file="pixel_snp_map_v6A.rds",
  compress="xz"
)




#Version 6A：視覺化Hilbert curve路徑順序----------------


#複製Hilbert template；轉成繪圖座標
hilbert_path_plot_v6A<-copy(hilbert_template_v6A)
hilbert_path_plot_v6A[,plot_x:=patch_col_v6A]
hilbert_path_plot_v6A[,plot_y:=hilbert_side_v6A-patch_row_v6A+1L]

#前64步
n_show_v6A<-64L
first_path_v6A<-hilbert_path_plot_v6A[1:n_show_v6A]

#起點終點座標
start_x<-hilbert_path_plot_v6A$plot_x[1]
start_y<-hilbert_path_plot_v6A$plot_y[1]
end_x<-hilbert_path_plot_v6A$plot_x[nrow(hilbert_path_plot_v6A)]
end_y<-hilbert_path_plot_v6A$plot_y[nrow(hilbert_path_plot_v6A)]

png("Hilbert_path_v6A_fix.png",width=1800,height=900,res=150)
par(mfrow=c(1,2),mar=c(4,4,4,2))

#左圖：完整Hilbert path
plot(hilbert_path_plot_v6A$plot_x,hilbert_path_plot_v6A$plot_y,
     type="n",
     xlim=c(1,hilbert_side_v6A+4),
     ylim=c(1,hilbert_side_v6A+2),
     xlab="Column",
     ylab="Row",
     main=paste0("Version 6A Hilbert path (",hilbert_side_v6A,"×",hilbert_side_v6A,")"),
     asp=1,
     xaxt="n",
     yaxt="n")

abline(v=seq(0.5,hilbert_side_v6A+0.5,by=1),col="grey92",lwd=0.3)
abline(h=seq(0.5,hilbert_side_v6A+0.5,by=1),col="grey92",lwd=0.3)
axis(1,at=seq(1,hilbert_side_v6A,by=8))
axis(2,at=seq(1,hilbert_side_v6A,by=8),las=1)

lines(hilbert_path_plot_v6A$plot_x,hilbert_path_plot_v6A$plot_y,col="grey35",lwd=1)
points(hilbert_path_plot_v6A$plot_x,hilbert_path_plot_v6A$plot_y,pch=16,cex=0.12,col="grey35")

#起點終點
points(start_x,start_y,pch=16,cex=1.6,col="red")
points(end_x,end_y,pch=16,cex=1.6,col="blue")

text(start_x+1.5,start_y+0.5,labels="Start",col="red",cex=0.9,xpd=NA,adj=c(0,0.5))
text(end_x-1.5,end_y+0.5,labels="End",col="blue",cex=0.9,xpd=NA,adj=c(1,0.5))
box()

#右圖：前64步放大
plot(first_path_v6A$plot_x,first_path_v6A$plot_y,
     type="n",
     xlim=c(min(first_path_v6A$plot_x)-1,max(first_path_v6A$plot_x)+2),
     ylim=c(min(first_path_v6A$plot_y)-1,max(first_path_v6A$plot_y)+1.5),
     xlab="Column",
     ylab="Row",
     main=paste0("First ",n_show_v6A," Hilbert steps"),
     asp=1)

abline(v=seq(floor(min(first_path_v6A$plot_x))-1,ceiling(max(first_path_v6A$plot_x))+1,by=1),col="grey92",lwd=0.3)
abline(h=seq(floor(min(first_path_v6A$plot_y))-1,ceiling(max(first_path_v6A$plot_y))+1,by=1),col="grey92",lwd=0.3)

lines(first_path_v6A$plot_x,first_path_v6A$plot_y,col="grey35",lwd=1.2)
points(first_path_v6A$plot_x,first_path_v6A$plot_y,pch=16,cex=0.8,col="grey35")

#每一步的index
text(first_path_v6A$plot_x,first_path_v6A$plot_y,
     labels=first_path_v6A$hilbert_index_v6A,
     pos=3,
     cex=0.5)

#起點終點
points(first_path_v6A$plot_x[1],first_path_v6A$plot_y[1],pch=16,cex=1.6,col="red")
points(first_path_v6A$plot_x[nrow(first_path_v6A)],first_path_v6A$plot_y[nrow(first_path_v6A)],pch=16,cex=1.6,col="blue")

text(first_path_v6A$plot_x[1]-0.4,first_path_v6A$plot_y[1]-0.2,labels="Start",col="red",cex=0.9,xpd=NA,adj=c(0,1))
text(first_path_v6A$plot_x[nrow(first_path_v6A)]-0.4,first_path_v6A$plot_y[nrow(first_path_v6A)]-0.2,
     labels="End",col="blue",cex=0.9,xpd=NA,adj=c(0,1))
box()
dev.off()

system("open Hilbert_path_v6A_fix.png")

#Version 6A Step 13：視覺化實際window在Hilbert curve中的填入方式----------------

#選擇SNP最多的window：Chr6_W32
view_window_key_v6A<-"Chr6_W32"

#抓出該window的SNP資料
view_window_snp_v6A<-snp_v6_del_y[
  window_key_v6A==view_window_key_v6A
]

#確認SNP數
nrow(view_window_snp_v6A)

#複製完整64×64 Hilbert template
hilbert_window_view_v6A<-copy(hilbert_template_v6A)

#0-based Hilbert座標轉成畫圖座標
hilbert_window_view_v6A[,plot_x:=patch_col_v6A]
hilbert_window_view_v6A[,plot_y:=hilbert_side_v6A-patch_row_v6A+1L]

#標示每個Hilbert位置是真實SNP還是padding
hilbert_window_view_v6A[,pixel_type:=
                          ifelse(
                            hilbert_index_v6A<nrow(view_window_snp_v6A),
                            "SNP",
                            "Padding"
                          )
]

table(hilbert_window_view_v6A$pixel_type)


#輸出圖片
png(
  "Chr6_W32_Hilbert_fill_v6A.png",
  width=1000,
  height=1000,
  res=150
)

par(mar=c(4,4,4,2))

#建立空白64×64畫布
plot(
  hilbert_window_view_v6A$plot_x,
  hilbert_window_view_v6A$plot_y,
  type="n",
  xlim=c(0.5,64.5),
  ylim=c(0.5,64.5),
  xlab="Column",
  ylab="Row",
  main="Chr6_W32: SNP filling along the 64×64 Hilbert curve",
  asp=1,
  xaxt="n",
  yaxt="n"
)

#格線
abline(
  v=seq(0.5,64.5,by=1),
  h=seq(0.5,64.5,by=1),
  col="grey93",
  lwd=0.3
)

axis(1,at=seq(1,64,by=8))
axis(2,at=seq(1,64,by=8),las=1)

#完整Hilbert路徑
lines(
  hilbert_window_view_v6A$plot_x,
  hilbert_window_view_v6A$plot_y,
  col="grey75",
  lwd=0.8
)

#真實SNP位置
points(
  hilbert_window_view_v6A[pixel_type=="SNP",plot_x],
  hilbert_window_view_v6A[pixel_type=="SNP",plot_y],
  pch=16,
  cex=0.25,
  col="black"
)

#padding位置
points(
  hilbert_window_view_v6A[pixel_type=="Padding",plot_x],
  hilbert_window_view_v6A[pixel_type=="Padding",plot_y],
  pch=16,
  cex=0.25,
  col="grey80"
)

#起點
points(
  hilbert_window_view_v6A$plot_x[1],
  hilbert_window_view_v6A$plot_y[1],
  pch=16,
  cex=1.3,
  col="red"
)

text(
  hilbert_window_view_v6A$plot_x[1]+1,
  hilbert_window_view_v6A$plot_y[1],
  "Start",
  col="red",
  pos=4
)

#最後一個真實SNP
last_snp_index_v6A<-nrow(view_window_snp_v6A)

points(
  hilbert_window_view_v6A$plot_x[last_snp_index_v6A],
  hilbert_window_view_v6A$plot_y[last_snp_index_v6A],
  pch=16,
  cex=1.3,
  col="blue"
)

text(
  hilbert_window_view_v6A$plot_x[last_snp_index_v6A]-1,
  hilbert_window_view_v6A$plot_y[last_snp_index_v6A],
  "Last SNP",
  col="blue",
  pos=2
)

legend(
  "bottomright",
  legend=c("SNP","Padding","Start","Last SNP"),
  pch=16,
  col=c("black","grey80","red","blue"),
  bty="n"
)

box()
dev.off()

#Mac直接開啟
system("open Chr6_W32_Hilbert_fill_v6A.png")



















