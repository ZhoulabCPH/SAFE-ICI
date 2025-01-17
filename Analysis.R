drug_reac = read.csv(paste(path, '/Original_data/drug_reac.csv', sep=''))

# length(unique(drug_reac[drug_reac$drugname_standard %in% c('PD-1/PD-L1', 'Combined', 'CTLA-4'), ]$primaryid))

ICI_clinical = drug_reac[drug_reac$drugname_standard %in% c('PD-1/PD-L1', 'Combined'), ]
ICI_clinical[ICI_clinical$cancer_type == '', 'cancer_type'] = 'Other'
ICI_clinical$event_dt_year = as.integer(ICI_clinical$event_dt_new/10000)

ICI_clinical[ICI_clinical$event_dt_flag=='' & !is.na(ICI_clinical$start_dt), 'duration'] = 
  as.numeric(difftime(as.Date(as.character(ICI_clinical[ICI_clinical$event_dt_flag=='' & !is.na(ICI_clinical$start_dt),]$event_dt_new), format='%Y%m%d'), 
                      as.Date(as.character(ICI_clinical[ICI_clinical$event_dt_flag=='' & !is.na(ICI_clinical$start_dt),]$start_dt), format='%Y%m%d'), units = 'days'))

ICI_clinical[ICI_clinical$duration <= 0 | is.na(ICI_clinical$duration), 'duration'] = NA

ICI_clinical$event_dt_year = factor(ICI_clinical$event_dt_year, levels = seq(2014,2022,1))
ICI_clinical$cancer_type = factor(ICI_clinical$cancer_type, levels = c('NSCLC', 'MEL', 'RCC', 'Other'))
ICI_clinical$outc_cod = factor(ICI_clinical$outc_cod, levels = c('DE'))
ICI_clinical$occr_country = factor(ICI_clinical$occr_country, levels = c('US', 'non_US'))
ICI_clinical$drugname_standard = factor(ICI_clinical$drugname_standard, levels = c('PD-1/PD-L1', 'Combined'))

ICI_clinical[ICI_clinical$outc_cod %in% c('DE'), 'outc_cod_num'] = 1
ICI_clinical[is.na(ICI_clinical$outc_cod_num), 'outc_cod_num'] = 0
ICI_clinical$Drug = paste(ICI_clinical$drugname_standard, ICI_clinical$cancer_type, sep='_')
ICI_clinical$Reac = ICI_clinical$PT_standard
ICI_clinical_Backup = ICI_clinical

Fatal_rate = ICI_clinical %>% group_by(Drug,Reac) %>% summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))

ICI_clinical = ICI_clinical[!duplicated(ICI_clinical$primaryid),]
ICI_clinical$age_num = round(ICI_clinical$age_num)
ICI_clinical_Backup$age_num = round(ICI_clinical_Backup$age_num)
#########################################################################################################
# Table 1: clinical characteristic

label(ICI_clinical$sex) <- "Sex"
label(ICI_clinical$age_num) <- "Age (years)"
label(ICI_clinical$age) <- "Age (class)"
label(ICI_clinical$cancer_type) <- "Cancer types"
label(ICI_clinical$duration) = 'Time-to-Onset (days)'
label(ICI_clinical$outc_cod) <- 'Outcome-Death'
label(ICI_clinical$occr_country) = 'Region reporting'
label(ICI_clinical$event_dt_year) = 'Reporting year'

print(ks.test(ICI_clinical$age_num, rnorm(100, mean(ICI_clinical$age_num), sd(ICI_clinical$age_num))))
duration = ICI_clinical[!is.na(ICI_clinical$duration), ]$duration
print(ks.test(duration, rnorm(100, mean(duration), sd(duration))))

# Paper version
strata <- c(list(`Total` = ICI_clinical[ICI_clinical$Drug %in% ICIs_item[1:8],]),
            split(ICI_clinical[ICI_clinical$Drug %in% ICIs_item[1:8],],
                  ICI_clinical[ICI_clinical$Drug %in% ICIs_item[1:8],]$drugname_standard),
            split(ICI_clinical[ICI_clinical$Drug %in% ICIs_item[1:8] & ICI_clinical$outc_cod_num == 1,],
                  ICI_clinical[ICI_clinical$Drug %in% ICIs_item[1:8] & ICI_clinical$outc_cod_num == 1,]$drugname_standard)
)

labels <- list(variables = list(age_num = render.varlabel(ICI_clinical$age_num),
                                age = render.varlabel(ICI_clinical$age),
                                sex = render.varlabel(ICI_clinical$sex), 
                                cancer_type = render.varlabel(ICI_clinical$cancer_type),
                                duration = render.varlabel(ICI_clinical$duration),
                                outc_cod = render.varlabel(ICI_clinical$outc_cod),
                                occr_country = render.varlabel(ICI_clinical$occr_country),
                                event_dt_year = render.varlabel(ICI_clinical$event_dt_year)
                                )
               )

caption = "Clinical characteristic of xxxxxxxxxx"  
footnote = "*All irAEs were classifed according to the United States Health and Human Services Common Terminology Criteria for Adverse Events (CTCAE) v.5.0 with grade ≥ 3 considered to be severe (High grade)"
table1(strata, 
       labels, 
       render.continuous = c(. = "Mean (SD)",
                             . = "Median [Q1, Q3]"),
       caption = caption, 
       footnote = footnote)

# Chinese version
strata <- c(list(`Total` = ICI_clinical),
            split(ICI_clinical, ICI_clinical$drugname_standard))

labels <- list(variables = list(age_num = render.varlabel(ICI_clinical$age_num),
                                age = render.varlabel(ICI_clinical$age),
                                sex = render.varlabel(ICI_clinical$sex), 
                                cancer_type = render.varlabel(ICI_clinical$cancer_type),
                                duration = render.varlabel(ICI_clinical$duration),
                                outc_cod = render.varlabel(ICI_clinical$outc_cod),
                                occr_country = render.varlabel(ICI_clinical$occr_country),
                                event_dt_year = render.varlabel(ICI_clinical$event_dt_year)
)
)

caption = "Clinical characteristic of xxxxxxxxxx"  
footnote = "*All irAEs were classifed according to the United States Health and Human Services Common Terminology Criteria for Adverse Events (CTCAE) v.5.0 with grade ≥ 3 considered to be severe (High grade)"
table1(strata, 
       labels, 
       render.continuous = c(. = "Mean (SD)",
                             . = "Median [Q1, Q3]"),
       caption = caption, 
       footnote = footnote)

print(wilcox.test(ICI_clinical[!is.na(ICI_clinical$duration) & ICI_clinical$drugname_standard == 'Combined', ]$duration,
                  ICI_clinical[!is.na(ICI_clinical$duration) & ICI_clinical$drugname_standard == 'PD-1/PD-L1', ]$duration))

print(chisq.test(matrix(c(unname(table(ICI_clinical[ICI_clinical$drugname_standard == 'PD-1/PD-L1', ]$outc_cod_num)),
                          unname(table(ICI_clinical[ICI_clinical$drugname_standard == 'Combined', ]$outc_cod_num))), ncol = 2))$p.value)


#########################################################################################################

data = read.csv(paste(path, '/Data/ror_PT_all_cancer/split_month/2022.csv', sep=''), sep=',')
data = data[data$Drug %in% ICIs_item[1:8] & data$ROR > 3 & data$EBGM > 3, ]
print(cor(data[,c('ROR', 'EBGM')], method = 'pearson'))

data_MedDRA = merge(data, ATC_MedDRA_embed, by = 'Reac')
k = data_MedDRA %>% group_by(short_SOC) %>% summarise(count = n())
data_MedDRA$short_SOC = factor(data_MedDRA$short_SOC, levels = k[order(-k$count), ]$short_SOC)
data_MedDRA$Drug = factor(data_MedDRA$Drug, levels = rev(ICIs_item))
data_MedDRA = merge(data_MedDRA, Fatal_rate, by = c('Drug', 'Reac'))

data_MedDRA$ICI = str_split_fixed(data_MedDRA$Drug, '_', 2)[,1]
data_MedDRA$Cancer = str_split_fixed(data_MedDRA$Drug, '_', 2)[,2]


SOC_count = data_MedDRA %>% group_by(Drug, short_SOC) %>% summarise(count = n())
SOC_count = rbind(SOC_count,
                      expand.grid(Drug = unique(SOC_count$Drug), short_SOC = unique(SOC_count$short_SOC), count = 0))
SOC_count = SOC_count[!duplicated(SOC_count[, c('Drug', 'short_SOC')]),]

count_plot = function(drug, color){
  ggplot(SOC_count[SOC_count$Drug == drug,], aes(x = short_SOC, y = count, label = count))+
    geom_bar(stat = 'identity', width = 0.8, position = position_dodge(), fill = color, alpha = 0.9)+
    geom_text(color = 'black', angle=0, vjust=1.5, size = 3)+
    annotate("text", x = 15, y = 18, label = paste('Total of', dim(data[data$Drug == drug,])[1], 'alert signals'))+
    geom_hline(yintercept = 10, linetype = 'dashed', color = 'darkgrey')+
    scale_y_continuous(limits = c(0,20), expand = c(0,0))+
    theme_bw()+
    theme_classic()+
    theme(
      axis.ticks.length.x.bottom = unit(0.1,'cm'),
      axis.ticks.length.y.left = unit(0.1,'cm'),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(color='black')
    )+
    labs(x='',y= 'No. of alert signals')
}

p1 = count_plot(drug = ICIs_item[1], color = '#d96098') 
p2 = count_plot(drug = ICIs_item[2], color = '#325288') 
p3 = count_plot(drug = ICIs_item[3], color = '#24a19c') 
p4 = count_plot(drug = ICIs_item[4], color = '#d3d1d1') 

p5 = count_plot(drug = ICIs_item[5], color = '#d96098') 
p6 = count_plot(drug = ICIs_item[6], color = '#325288') 
p7 = count_plot(drug = ICIs_item[7], color = '#24a19c') 
p8 = count_plot(drug = ICIs_item[8], color = '#d3d1d1') 

pdf('毕业_SOC_count.pdf', height = 4, width = 8)
print(p1+p5+p2+p6)
print(p3+p7+p4+p8)
dev.off()

#########################################################################################################
# table 1

Reac_term_out = data.frame()
Reac_term = ATC_MedDRA_embed[ATC_MedDRA_embed$Reac %in% data$Reac, ]
for (i in unique(Reac_term$short_SOC)){
  temp = Reac_term[Reac_term$short_SOC == i,]
  Reac_term_out = rbind(Reac_term_out, data.frame(short_SOC = i, SOC = temp$SOC[1], PT = paste(str_to_sentence(tolower(temp$Reac)), collapse = ', ')))
}
SourceData = createWorkbook()
SourceData_info(fig_data = Reac_term_out,
                colnames = colnames(Reac_term_out),
                cols = 1:ncol(Reac_term_out),
                sheet = 'sheet',
                widths = c(10, 20, 40))

saveWorkbook(SourceData, '毕业_Table_soc_pt_term.xlsx', overwrite = TRUE)

# table 2
table_data = data
table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]

table_data = merge(table_data, Fatal_rate, by = c('Drug', 'Reac'))
table_data$Reac = str_to_title(tolower(table_data$Reac))
table_data$ROR = round(table_data$ROR,1)
table_data$EBGM = round(table_data$EBGM,1)

out_1 = reshape2::dcast(table_data[(table_data$ROR+table_data$EBGM)/2 >=50, ], 
                      Reac+Cancer~ICI, value.var = c('Report_count'))
out_2 = reshape2::dcast(table_data[(table_data$ROR+table_data$EBGM)/2 >=50, ], 
                        Reac+Cancer~ICI, value.var = c('EBGM'))
out_3 = reshape2::dcast(table_data[(table_data$ROR+table_data$EBGM)/2 >=50, ], 
                        Reac+Cancer~ICI, value.var = c('ROR'))
out = left_join(left_join(out_1, out_2, by = c('Reac', 'Cancer')), out_3, by = c('Reac', 'Cancer'))

out[is.na(out)] = '--'

SourceData = createWorkbook()
SourceData_info(fig_data = out,
                colnames = colnames(out),
                cols = 1:ncol(out),
                sheet = 'Table_3_2',
                widths = c(20, rep(10, ncol(out)-1)))

saveWorkbook(SourceData, '毕业_Table_high_ror_ebgm.xlsx', overwrite = TRUE)


out = table_data[,c(12,13,2, 4:6, 16)]
out = out[order(out$ICI, out$Cancer, -out$N),]

SourceData = createWorkbook()
SourceData_info(fig_data = out,
                colnames = colnames(out),
                cols = 1:ncol(out),
                sheet = 'Table_3_2',
                widths = rep(10, ncol(out)))

saveWorkbook(SourceData, '毕业_Table_ror_ebgm_all.xlsx', overwrite = TRUE)




# table 3
Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration<=90  & !is.na(ICI_clinical_Backup$duration), ] %>% 
  group_by(Drug,Reac) %>% 
  summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))

table_data = data
table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]

table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
table_data$Reac = str_to_title(tolower(table_data$Reac))
table_data$ROR = round(table_data$ROR,1)
table_data$EBGM = round(table_data$EBGM,1)

table_data$fatal = paste(table_data$Fatal_count, " (", round(table_data$Fatal_rate, 1), '%)',sep='')
table_data =  table_data[table_data$Fatal_rate >=20 , ]

w = table_data %>% group_by(Reac)%>% summarise(n=n())
table_data = table_data[table_data$Reac %in% w[w$n >= 4, ]$Reac, ]

out_1 = reshape2::dcast(table_data, 
                        Reac+Cancer~ICI, value.var = c('fatal'))
out_2 = reshape2::dcast(table_data, 
                        Reac+Cancer~ICI, value.var = c('EBGM'))
out_3 = reshape2::dcast(table_data, 
                        Reac+Cancer~ICI, value.var = c('ROR'))
out = left_join(left_join(out_1, out_2, by = c('Reac', 'Cancer')), out_3, by = c('Reac', 'Cancer'))

out[is.na(out)] = '--'

out$Cancer = factor(out$Cancer, levels = c('NSCLC', 'MEL', 'RCC', 'Other'))
out = out[order(out$Cancer), ]

out$Reac = factor(out$Reac, levels = w[order(-w$n), ]$Reac)
out = out[order(out$Reac), ]

SourceData = createWorkbook()
SourceData_info(fig_data = out,
                colnames = colnames(out),
                cols = 1:ncol(out),
                sheet = 'Table_3_3',
                widths = c(20, rep(10, ncol(out)-1)))

saveWorkbook(SourceData, '毕业_Table_high_fatality.xlsx', overwrite = TRUE)

SOC_ROR_EBGM_fatal = function(soc){
  Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration<=90, ] %>% 
    group_by(Drug,Reac) %>% 
    summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))
  
  table_data = data
  table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
  table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]
  
  table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
  table_data = merge(table_data, ATC_MedDRA_embed, by = 'Reac')
  
  table_data$short_SOC = factor(table_data$short_SOC, levels = levels(data_MedDRA$short_SOC))
  table_data$Drug = factor(table_data$Drug, levels = rev(ICIs_item))

  temp = table_data[table_data$short_SOC == soc, ]
  temp[temp$ROR >80, 'ROR'] = 80
  temp[temp$EBGM >80, 'EBGM'] = 80
  temp[temp$Fatal_rate >80, 'Fatal_rate'] = 80
  
  temp$Reac = str_to_title(tolower(temp$Reac))
  
  p = ggplot()+
    geom_point(data = temp[temp$ROR + temp$EBGM  <= 100,], aes(x = Reac, y = Drug, size = (ROR+EBGM)/2, fill = Fatal_rate-20), color='black', shape = 21)+
    geom_point(data = temp[temp$ROR + temp$EBGM > 100,], aes(x = Reac, y = Drug, size = (ROR+EBGM)/2, fill = Fatal_rate-20), color='black', shape = 22)+
    scale_fill_gradient2(low = 'blue', mid = '#f2e6e6', high = 'red', limits = c(0, 80)-20, breaks = seq(0, 80, 20)-20, labels = str(seq(0, 80, 20)))+
    scale_size(limits = c(0, 80), breaks = c(10, 20, 30, 50, 80))+
    ggtitle(soc)+
    theme_bw()+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=0, hjust = 1)
    )+
    labs(x = '', y = '')
  p
}

pdf('毕业_SOC_ROR_EBGM_fatal_1.pdf', height = 2.5, width = 8)
for(i in levels(data_MedDRA$short_SOC))
  print(SOC_ROR_EBGM_fatal(soc = i))
dev.off()

SOC_ROR_EBGM_fatal = function(soc){
  Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration<=90, ] %>% 
    group_by(Drug,Reac) %>% 
    summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))
  
  table_data = data
  table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
  table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]
  
  table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
  table_data = merge(table_data, ATC_MedDRA_embed, by = 'Reac')
  
  table_data$short_SOC = factor(table_data$short_SOC, levels = levels(data_MedDRA$short_SOC))
  table_data$Drug = factor(table_data$Drug, levels = rev(ICIs_item))
  
  temp = table_data[table_data$short_SOC == soc, ]
  temp[temp$ROR >80, 'ROR'] = 80
  temp[temp$EBGM >80, 'EBGM'] = 80
  temp[temp$Fatal_rate >80, 'Fatal_rate'] = 80
  
  temp$Reac = str_to_title(tolower(temp$Reac))
  
  p = ggplot()+
    geom_point(data = temp[temp$ROR + temp$EBGM  <= 100,], aes(x = Reac, y = Drug, size = (ROR+EBGM)/2, fill = Fatal_rate-20), color='black', shape = 21)+
    geom_point(data = temp[temp$ROR + temp$EBGM > 100,], aes(x = Reac, y = Drug, size = (ROR+EBGM)/2, fill = Fatal_rate-20), color='black', shape = 22)+
    scale_fill_gradient2(low = 'blue', mid = '#f2e6e6', high = 'red', limits = c(0, 80)-20, breaks = seq(0, 80, 20)-20, labels = str(seq(0, 80, 20)))+
    scale_size(limits = c(0, 80), breaks = c(10, 20, 30, 50, 80))+
    ggtitle(soc)+
    theme_bw()+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=0, hjust = 1),
          axis.text.y = element_blank()
    )+
    labs(x = '', y = '')
  p
}

p1 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[1])
p2 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[2])
p3 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[3])
p4 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[4])
p5 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[5])
p6 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[6])
p7 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[7])
p8 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[8])
p9 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[9])
p10 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[10])
p11 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[11])
p12 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[12])
p13 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[13])
p14 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[14])
p15 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[15])
p16 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[16])
p17 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[17])
p18 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[18])


pdf('毕业_SOC_ROR_EBGM_fatal.pdf', height = 3.5, width = 14)

print(p1+p2+p16+p18 + plot_layout(widths = c(3.2,2.1, 0.8,0.3),guides='collect') & theme(legend.position='bottom'))
print(p3+p4+p7+p15 +plot_layout(widths = c(2,2,2,1.4),guides='collect') & theme(legend.position='bottom'))
print(p5+p6+p8 +p13+p15+plot_layout(widths = c(3.2,2.7,1.9, 1.2,1.2),guides='collect') & theme(legend.position='bottom'))
print(p9+p10+p11+p12+p14+p17 + plot_layout(widths = c(1.8,2,1.2,1.5, 1.1,0.9),guides='collect') & theme(legend.position='bottom'))

dev.off()

SOC_ROR_EBGM_fatal = function(soc){
  Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration<=90, ] %>% 
    group_by(Drug,Reac) %>% 
    summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))
  
  table_data = data
  table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
  table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]
  
  table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
  table_data = merge(table_data, ATC_MedDRA_embed, by = 'Reac')
  
  table_data$short_SOC = factor(table_data$short_SOC, levels = levels(data_MedDRA$short_SOC))
  table_data$Drug = factor(table_data$Drug, levels = rev(ICIs_item))
  
  temp = table_data[table_data$short_SOC == soc, ]
  temp[temp$ROR >80, 'ROR'] = 80
  temp[temp$EBGM >80, 'EBGM'] = 80
  temp[temp$Fatal_rate >80, 'Fatal_rate'] = 80
  
  temp$Reac = str_to_title(tolower(temp$Reac))
  
  p = ggplot()+
    geom_point(data = temp[temp$ROR + temp$EBGM  <= 100,], aes(x = Reac, y = Drug, size = (ROR+EBGM)/2, fill = Fatal_rate-20), color='black', shape = 21)+
    geom_point(data = temp[temp$ROR + temp$EBGM > 100,], aes(x = Reac, y = Drug, size = (ROR+EBGM)/2, fill = Fatal_rate-20), color='black', shape = 22)+
    scale_fill_gradient(low = '#f2e6e6', high = 'red', limits = c(0, 80)-20, breaks = seq(0, 80, 20))+
    scale_size(limits = c(0, 80), breaks = c(10, 20, 30, 50, 80))+
    ggtitle(soc)+
    theme_bw()+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=90, hjust = 1),
          axis.text.y = element_blank()
    )+
    labs(x = '', y = '')
  p
}

p1 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[1])
p2 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[2])
p3 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[3])
p4 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[4])
p5 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[5])
p6 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[6])
p7 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[7])
p8 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[8])
p9 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[9])
p10 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[10])
p11 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[11])
p12 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[12])
p13 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[13])
p14 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[14])
p15 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[15])
p16 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[16])
p17 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[17])
p18 = SOC_ROR_EBGM_fatal(soc = levels(data_MedDRA$short_SOC)[18])


pdf('毕业_SOC_ROR_EBGM_fatal_2.pdf', height = 3.5, width = 14)

print(p1+p2+p16+p18 + plot_layout(widths = c(3.2,2.1, 0.8,0.3),guides='collect') & theme(legend.position='bottom'))
print(p3+p4+p7+p15 +plot_layout(widths = c(2,2,2,1.4),guides='collect') & theme(legend.position='bottom'))
print(p5+p6+p8 +p13+p15+plot_layout(widths = c(3.2,2.7,1.9, 1.2,1.2),guides='collect') & theme(legend.position='bottom'))
print(p9+p10+p11+p12+p14+p17 + plot_layout(widths = c(1.8,2,1.2,1.5, 1.1,0.9),guides='collect') & theme(legend.position='bottom'))

dev.off()

Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration<=90, ] %>% 
  group_by(Drug,Reac) %>% 
  summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))

table_data = data
table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]

table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
table_data = merge(table_data, ATC_MedDRA_embed, by = 'Reac')

temp = table_data

temp[temp$ROR >80, 'ROR'] = 80
temp[temp$EBGM >80, 'EBGM'] = 80
temp[temp$Fatal_rate >80, 'Fatal_rate'] = 80

temp$Reac = str_to_title(tolower(temp$Reac))

w = temp[temp$Fatal_rate >= 20,] %>% group_by(short_SOC, Reac) %>% summarise(n = n())

###################################################################################################################################

# check = ICI_clinical[!is.na(ICI_clinical$duration), ] %>% group_by(Drug, duration) %>% summarise(n = n(), fatal_n = sum(outc_cod_num))
# check = check[order(check$Drug, check$duration),]
# 
# fatal_out = data.frame()
# for (drug in ICIs_item[1:8]){
#   ready = check[check$Drug == drug, ]
#   fatal_out = rbind(fatal_out, data.frame(Drug = drug, duration = ready$duration,Fatal_rate = cumsum(ready$fatal_n)/cumsum(ready$n)))
# }
# 
# ggplot()+
#   geom_line(data = fatal_out, aes(x = duration, y = Fatal_rate*100, group = Drug))+
#   scale_x_continuous(limits = c(0, 120/7), breaks = seq(0, 120, 14)/7, labels = seq(0, 120, 14)/7)

input = data.frame()
Fatal_rate_2 = data.frame()  
for (soc in levels(data_MedDRA$short_SOC)[1:18]){
  time_data_ori = merge(ICI_clinical_Backup, data[,c('Drug', 'Reac')], by = c('Drug', 'Reac'))
  temp = time_data_ori[time_data_ori$Reac %in% ATC_MedDRA_embed[ATC_MedDRA_embed$short_SOC == soc, ]$Reac, ]
  temp = ICI_clinical[ICI_clinical$primaryid %in% temp$primaryid, ]
  Fatal_rate_2 = rbind(Fatal_rate_2 , cbind(short_SOC = soc, temp %>% group_by(Drug) %>% summarise(Fatal_rate = sum(outc_cod_num)/n() * 100)))
  for (drug in unique(temp$Drug)){
    line = density(temp[!is.na(temp$duration) & temp$Drug == drug,]$duration/7, kernel="gaussian", n=1000)
    input = rbind(input,data.frame(Drug = drug, short_SOC = soc, x = line$x, y = line$y, max_y = max(line$y)))
  }
}
input = merge(input, Fatal_rate_2, by = c('Drug', 'short_SOC'))

######
Fatal_rate_2$short_SOC = factor(Fatal_rate_2$short_SOC, levels = rev(levels(data_MedDRA$short_SOC)))
Fatal_rate_2$Drug = factor(Fatal_rate_2$Drug, levels = ICIs_item)

pdf('毕业_SOC_Time_to_Onset_fatal_value.pdf', height = 6, width = 5)

p = ggplot()+
  geom_point(data = Fatal_rate_2, 
             aes(x = Drug, y = short_SOC, size = 0, fill = Fatal_rate-30), 
             color='black', shape = 21)+
  geom_point(data = Fatal_rate_2[Fatal_rate_2$Fatal_rate < 20,], 
             aes(x = Drug, y = short_SOC, size = 0.5, fill = Fatal_rate-30), 
             color='black', shape = 21)+
  geom_point(data = Fatal_rate_2[Fatal_rate_2$Fatal_rate < 40 & Fatal_rate_2$Fatal_rate >= 20,], 
             aes(x = Drug, y = short_SOC, size = 1, fill = Fatal_rate-30), 
             color='black', shape = 22)+
  geom_point(data = Fatal_rate_2[Fatal_rate_2$Fatal_rate >= 40,], 
             aes(x = Drug, y = short_SOC, size = 3, fill = Fatal_rate-30), 
             color='black', shape = 22)+
  scale_fill_gradient2(low = '#527853', mid = '#f9e8d9', high = '#ee7214', 
                       limits = c(0, 62)-30, breaks = seq(0, 62, 20) - 30, labels = seq(0, 62, 20))+
  theme_bw()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle=90, hjust = 1)
  )+
  labs(x = '', y = '')
print(p)

dev.off()

######

Fatal_duration = function(drug, max){
  ggplot()+
    geom_vline(xintercept = 90/7, linetype ='dotdash', color = 'darkgrey')+
    geom_line(data = input[input$Drug == drug, ], aes(x = x, y = y*Fatal_rate, group = short_SOC, color = short_SOC))+
    scale_x_continuous(limits = c(0, 120/7), breaks = seq(0, 120, 14)/7, labels = seq(0, 120, 14)/7)+
    scale_y_continuous(limits = c(0, max*7), expand = c(0,0))+
    
    scale_color_manual(values = c("#ED0000", "#FFCD00", "#f79486", "#f9a828", "#FFDC91",'#fbfad3', "#a3de83","#d5def5", "#ADE2D0", "#00A087",
                                  "#0387B1", "#925E9F", '#8594e4', '#6643b5', "#3C5488", "#6F99AD", "#A9A9A9",'#494953'),
                       breaks = levels(data_MedDRA$short_SOC)[1:18])+

    ggtitle(drug)+
    theme_bw()+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),
          # axis.text.x = element_text(angle=0, hjust = 1),
          panel.grid.minor = element_blank()
    )+
    labs(x = 'Time-to-Onset (weeks)', y = 'Density*(Outcome to Death (%))')
}

pdf('毕业_SOC_Time_to_Onset_fatal.pdf', height = 8, width = 12)
p1 = Fatal_duration(ICIs_item[1], max = 1.1)
p2 = Fatal_duration(ICIs_item[2], max = 1.1)
p3 = Fatal_duration(ICIs_item[3], max = 1.1)
p4 = Fatal_duration(ICIs_item[4], max = 1.1)

p5 = Fatal_duration(ICIs_item[5], max = 1.1)
p6 = Fatal_duration(ICIs_item[6], max = 1.1)
p7 = Fatal_duration(ICIs_item[7], max = 1.1)
p8 = Fatal_duration(ICIs_item[8], max = 1.1)

print((p1|p2)/(p3|p4))
print((p5|p6)/(p7|p8))

p5 = Fatal_duration(ICIs_item[5], max = 2.6)
p6 = Fatal_duration(ICIs_item[6], max = 2.6)
p7 = Fatal_duration(ICIs_item[7], max = 2.6)
p8 = Fatal_duration(ICIs_item[8], max = 2.6)
print((p5|p6)/(p7|p8))

dev.off()


w=input[input$Drug == ICIs_item[1], ]
w$x[which.max(w$y*w$Fatal_rate)]*7


#######################################################################################################
time_data_ori = merge(ICI_clinical_Backup, data[,c('Drug', 'Reac')], by = c('Drug', 'Reac'))

label_data = data.frame()
pdf('毕业_SOC_Time_to_Onset.pdf', height = 5, width = 12)
for (soc in levels(data_MedDRA$short_SOC)[1:17]){
  temp = time_data_ori[time_data_ori$Reac %in% ATC_MedDRA_embed[ATC_MedDRA_embed$short_SOC == soc, ]$Reac, ]
  temp = ICI_clinical[ICI_clinical$primaryid %in% temp$primaryid, ]
  temp$status = 1
  temp$Drug = factor(temp$Drug, levels = ICIs_item)
  
  fit <- survfit(Surv(duration, status) ~ Drug, data = temp[temp$Drug %in% ICIs_item[1:4], ])
  label_data_1 = temp[temp$Drug %in% ICIs_item[1:4] & !is.na(temp$duration), ] %>% 
    group_by(Drug) %>% 
    summarise(short_SOC = soc,
              log_rank_p = survdiff(Surv(duration, status) ~ Drug, data = temp[temp$Drug %in% ICIs_item[1:4], ])$pvalue, 
              median = median(duration), y = 0.1, median_str = paste('D-', median(duration), sep=''))
  
  p = ggsurvplot(fit,
                 pval = TRUE, conf.int = TRUE,conf.int.alpha = 0.1,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 linetype = "solid",
                 surv.median.line = "hv",
                 ggtheme = theme_bw()+
                   theme_classic2()+
                   theme(axis.text = element_text(color = 'black'),
                         panel.grid = element_blank(),
                         legend.position = c(0.1, 0.8),
                         legend.justification = c(0, 1),
                         legend.box.just = "left"
                   ),
                 # palette = c('#ff4273', '#3e92a3', '#ff9d76',  '#d3d1d1'),
                 palette = c('#d96098', '#325288', '#24a19c',  '#d3d1d1'),
                 xlim = c(0, 90),
                 break.time.by = 14,
                 pval.coord = c(0, 0.2),
                 pval.prefix = "Log-rank test"
  )
  
  a = p$plot + 
    scale_x_continuous(breaks = seq(14, 90, 14), labels = seq(14, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    geom_text(data = label_data_1, aes(x = median, y = y, label = median_str, color = paste('Drug=', Drug, sep='')), size = 3, alpha=1)+
    ggtitle(paste("PD-1/PD-L1 for", soc))+
    theme(axis.text = element_text(color = 'black'),
          panel.grid = element_blank(),
          legend.position = c(0.8, 1),
          legend.justification = c(0, 1),
          legend.box.just = "left"
    )+
    labs(y = '1 - Onset probability')
  
  b = p$table + 
    scale_x_continuous(limits = c(0, 90), breaks = seq(14, 90, 14), labels = seq(14, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    theme(axis.text = element_text(color = 'black'),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none'
    )
  
  fit <- survfit(Surv(duration, status) ~ Drug, data = temp[temp$Drug %in% ICIs_item[5:8], ])
  label_data_2 = temp[temp$Drug %in% ICIs_item[5:8] & !is.na(temp$duration), ] %>% 
    group_by(Drug) %>% 
    summarise(short_SOC = soc,
              log_rank_p = survdiff(Surv(duration, status) ~ Drug, data = temp[temp$Drug %in% ICIs_item[5:8], ])$pvalue, 
              median = median(duration), y = 0.1, median_str = paste('D-', median(duration), sep=''))
  
  p = ggsurvplot(fit,
                 pval = TRUE, conf.int = TRUE,conf.int.alpha = 0.1,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 linetype = "solid",
                 surv.median.line = "hv",
                 ggtheme = theme_bw()+
                   theme_classic2()+
                   theme(axis.text = element_text(color = 'black'),
                         panel.grid = element_blank(),
                         legend.position = c(0.1, 0.8),
                         legend.justification = c(0, 1),
                         legend.box.just = "left"
                   ),
                 # palette = c('#ff4273', '#3e92a3', '#ff9d76',  '#d3d1d1'),
                 palette = c('#d96098', '#325288', '#24a19c',  '#d3d1d1'),
                 xlim = c(0, 90),
                 break.time.by = 14,
                 pval.coord = c(0, 0.2),
                 pval.prefix = "Log-rank test"
  )
  
  c = p$plot + 
    scale_x_continuous(limits = c(0, 90), breaks = seq(14, 90, 14), labels = seq(14, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    geom_text(data = label_data_2, aes(x = median, y = y, label = median_str, color = paste('Drug=', Drug, sep='')), size = 3, alpha=1)+
    ggtitle(paste("Combined for", soc))+
    theme(axis.text = element_text(color = 'black'),
          panel.grid = element_blank(),
          legend.position = c(0.8, 1),
          legend.justification = c(0, 1),
          legend.box.just = "left"
    )+
    labs(y = '1 - Onset probability')
  
  d = p$table + 
    scale_x_continuous(breaks = seq(14, 90, 14), labels = seq(14, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    theme(axis.text = element_text(color = 'black'),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none'
    )
  print((a /b + plot_layout(heights = c(5,2),guides='collect') & theme(legend.position='right'))|
    (c /d + plot_layout(heights = c(5,2),guides='collect') & theme(legend.position='right')))
  
  label_data = rbind(label_data, rbind(label_data_1, label_data_2))

}
dev.off()


label_data$Drug = factor(label_data$Drug, levels = ICIs_item[c(1,5,2,6,3,7,4,8)])
for (i in 1:4){
  label_data[label_data$Drug == ICIs_item[i], 'Drug_num'] = (i-1)*2 +1
}
for (i in 5:8){
  label_data[label_data$Drug == ICIs_item[i], 'Drug_num'] = (i-4)*2
}
  
pdf('毕业_SOC_Time_to_Onset_collected.pdf', height = 12, width = 12)
p = ggplot()+
  # geom_point(data = label_data, aes(x = Drug, median, label=paste(short_SOC, median_str, sep = ': ')))+
  geom_errorbar(data = label_data, aes(x = Drug, ymin=median, ymax=median, color = short_SOC),
                width = 0.5)+
  # geom_text(data = label_data, aes(x = Drug_num, y=median, label=paste(short_SOC, median_str,sep = '\n'), color = short_SOC),
  #           angle = 90, size = 3, hjust = 1)+
  scale_y_continuous(limits = c(min(label_data$median), 100), breaks = seq(14, 100, 14), labels = seq(14, 100, 14)/7, name = 'Time-to-Onset (weeks)')+
  coord_flip()+
  scale_color_manual(values = c("#ED0000", "#FFCD00", "#f79486", "#f9a828", "#FFDC91",'#fbfad3', "#a3de83","#d5def5", "#ADE2D0", "#00A087",
                                "#0387B1", "#925E9F", '#8594e4', '#6643b5', "#3C5488", "#6F99AD", "#A9A9A9",'#494953'),
                     breaks = levels(data_MedDRA$short_SOC)[1:18])+
  theme_bw()+
  # theme_classic2()+
  theme(axis.text = element_text(color = 'black'),
        # axis.text.x = element_text(angle=0, hjust = 1),
        panel.grid.minor = element_blank()
  )+
  labs(x = 'Time-to-Onset (weeks)', y = 'Density*(Outcome to Death (%))')
print(p)
p = ggplot()+
  # geom_point(data = label_data, aes(x = Drug, median, label=paste(short_SOC, median_str, sep = ': ')))+
  # geom_errorbar(data = label_data, aes(x = Drug, ymin=median, ymax=median, color = short_SOC), 
  #               width = 0.5)+
  geom_text(data = label_data, aes(x = Drug_num, y=median, label=median_str, color = short_SOC),
            angle = 45, size = 3, hjust = 1)+
  scale_y_continuous(limits = c(min(label_data$median), 100), breaks = seq(14, 100, 14), labels = seq(14, 100, 14)/7, name = 'Time-to-Onset (weeks)')+
  coord_flip()+
  scale_color_manual(values = c("#ED0000", "#FFCD00", "#f79486", "#f9a828", "#FFDC91",'#fbfad3', "#a3de83","#d5def5", "#ADE2D0", "#00A087",
                                "#0387B1", "#925E9F", '#8594e4', '#6643b5', "#3C5488", "#6F99AD", "#A9A9A9",'#494953'),
                     breaks = levels(data_MedDRA$short_SOC)[1:18])+
  theme_bw()+
  # theme_classic2()+
  theme(axis.text = element_text(color = 'black')
        # axis.text.x = element_text(angle=0, hjust = 1),
        # panel.grid.minor = element_blank()
  )+
  labs(x = 'Time-to-Onset (weeks)', y = 'Density*(Outcome to Death (%))')
print(p)
dev.off()


SourceData = createWorkbook()
SourceData_info(fig_data = Reac_term_out,
                colnames = colnames(Reac_term_out),
                cols = 1:ncol(Reac_term_out),
                sheet = 'sheet',
                widths = c(10, 20, 40))

saveWorkbook(SourceData, '毕业_Table_soc_pt_term.xlsx', overwrite = TRUE)


label_data$therapy = str_split_fixed(label_data$Drug, '_', 2)[,1]
label = label_data %>% group_by(therapy, short_SOC) %>% summarise(change = max(median) - min(median))

