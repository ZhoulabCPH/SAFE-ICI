# holdout test set

drug_reac = read.csv(paste(path, '/Original_data/drug_reac_holdout.csv', sep=''))

# length(unique(drug_reac[drug_reac$drugname_standard %in% c('PD-1/PD-L1', 'Combined', 'CTLA-4'), ]$primaryid))

ICI_holdout = drug_reac[drug_reac$drugname_standard %in% c('PD-1/PD-L1', 'Combined'), ]
ICI_holdout[ICI_holdout$cancer_type == '', 'cancer_type'] = 'Other'
ICI_holdout$event_dt_year = as.integer(ICI_holdout$event_dt_new/10000)


# ICI_holdout$event_dt_flag no flag
ICI_holdout[is.na(ICI_holdout$event_dt_flag) & !is.na(ICI_holdout$start_dt), 'duration'] = 
  as.numeric(difftime(as.Date(as.character(ICI_holdout[is.na(ICI_holdout$event_dt_flag) & !is.na(ICI_holdout$start_dt),]$event_dt_new), format='%Y%m%d'), 
                      as.Date(as.character(ICI_holdout[is.na(ICI_holdout$event_dt_flag) & !is.na(ICI_holdout$start_dt),]$start_dt), format='%Y%m%d'), units = 'days'))

ICI_holdout[ICI_holdout$duration <= 0 | is.na(ICI_holdout$duration), 'duration'] = NA

ICI_holdout$cancer_type = factor(ICI_holdout$cancer_type, levels = c('NSCLC', 'MEL', 'RCC', 'Other'))
ICI_holdout$outc_cod = factor(ICI_holdout$outc_cod, levels = c('DE'))
ICI_holdout$occr_country = factor(ICI_holdout$occr_country, levels = c('US', 'non_US'))
ICI_holdout$drugname_standard = factor(ICI_holdout$drugname_standard, levels = c('PD-1/PD-L1', 'Combined'))

ICI_holdout[ICI_holdout$outc_cod %in% c('DE'), 'outc_cod_num'] = 1
ICI_holdout[is.na(ICI_holdout$outc_cod_num), 'outc_cod_num'] = 0
ICI_holdout$Drug = paste(ICI_holdout$drugname_standard, ICI_holdout$cancer_type, sep='_')
ICI_holdout$Reac = ICI_holdout$PT_standard
ICI_holdout_Backup = ICI_holdout

ICI_holdout$age_num = round(ICI_holdout$age_num)
ICI_holdout_Backup$age_num = round(ICI_holdout_Backup$age_num)

####################
# table 1
ICI_holdout = ICI_holdout[!duplicated(ICI_holdout$primaryid) & 
                            ICI_holdout$duration <= 90 &
                            !is.na(ICI_holdout$duration),]

label(ICI_holdout$sex) <- "Sex"
label(ICI_holdout$age_num) <- "Age (years)"
label(ICI_holdout$age) <- "Age (class)"
label(ICI_holdout$cancer_type) <- "Cancer types"
label(ICI_holdout$duration) = 'Time-to-Onset (days)'
label(ICI_holdout$outc_cod) <- 'Outcome to Death in 90 days'
label(ICI_holdout$occr_country) = 'Region reporting'

# Paper version
strata <- c(list(`Total` = ICI_holdout[ICI_holdout$Drug %in% ICIs_item[1:8],]),
            split(ICI_holdout[ICI_holdout$Drug %in% ICIs_item[1:8],],
                  ICI_holdout[ICI_holdout$Drug %in% ICIs_item[1:8],]$drugname_standard),
            split(ICI_holdout[ICI_holdout$Drug %in% ICIs_item[1:8] & ICI_holdout$outc_cod_num == 1,],
                  ICI_holdout[ICI_holdout$Drug %in% ICIs_item[1:8] & ICI_holdout$outc_cod_num == 1,]$drugname_standard)
)

labels <- list(variables = list(age_num = render.varlabel(ICI_holdout$age_num),
                                age = render.varlabel(ICI_holdout$age),
                                sex = render.varlabel(ICI_holdout$sex), 
                                cancer_type = render.varlabel(ICI_holdout$cancer_type),
                                outc_cod = render.varlabel(ICI_holdout$outc_cod),
                                duration = render.varlabel(ICI_holdout$duration),
                                occr_country = render.varlabel(ICI_holdout$occr_country)
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

print(wilcox.test(ICI_holdout[!is.na(ICI_holdout$duration) & ICI_holdout$drugname_standard == 'Combined', ]$duration,
                  ICI_holdout[!is.na(ICI_holdout$duration) & ICI_holdout$drugname_standard == 'PD-1/PD-L1', ]$duration))

print(chisq.test(matrix(c(unname(table(ICI_holdout[ICI_holdout$drugname_standard == 'PD-1/PD-L1', ]$outc_cod_num)),
                          unname(table(ICI_holdout[ICI_holdout$drugname_standard == 'Combined', ]$outc_cod_num))), ncol = 2))$p.value)




# Chinese version
strata <- c(list(`Total` = ICI_holdout),
            split(ICI_clinical, ICI_holdout$drugname_standard))

labels <- list(variables = list(age_num = render.varlabel(ICI_holdout$age_num),
                                age = render.varlabel(ICI_holdout$age),
                                sex = render.varlabel(ICI_holdout$sex), 
                                cancer_type = render.varlabel(ICI_holdout$cancer_type),
                                outc_cod = render.varlabel(ICI_holdout$outc_cod),
                                duration = render.varlabel(ICI_holdout$duration),
                                occr_country = render.varlabel(ICI_holdout$occr_country)
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


# wilcox.test(ICI_holdout[ICI_holdout$durbation])


#########################################################################################################################
library('eulerr')
library('forestploter')
library('caret')

# table 3
Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration<=90 & !is.na(ICI_clinical_Backup$duration), ] %>% 
  group_by(Drug,Reac) %>% 
  summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))

table_data = data
table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]

table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
table_data$ROR = round(table_data$ROR,1)
table_data$EBGM = round(table_data$EBGM,1)

table_data$fatal = paste(table_data$Fatal_count, " (", round(table_data$Fatal_rate, 1), '%)',sep='')
table_data =  table_data[table_data$Fatal_rate >=20 , ]


w = table_data %>% group_by(Reac) %>% summarise(n = n())
table_data = table_data[table_data$Reac %in% w[w$n > 4, ]$Reac, ]

#######################################################################################################
# color_list = c('#ffffff', '#d9d9d9', '#add8e6', "#FFDC91",'#fffacd', "#d5def5", '#76eec6',"#f79486", '#b7e5b4', '#f0dbaf')

color_list = c("#fdface",'#d9d9d9', '#ffffff', '#b0d8e6','#ed8080', '#d5def5', '#fba07b','#ffeedd', '#67abed', '#cd95cd')


pdf('毕业_PT_VENN.pdf', height = 8, width = 8)

Reac_ref = table_data[,c('Reac', 'Fatal_rate')]
Reac_ref = unique(Reac_ref[order(-Reac_ref$Fatal_rate), ]$Reac)

for (i in 1:8){
  reac_list = table_data[table_data$Drug == ICIs_item[i],c('Reac', 'Fatal_rate')]
  reac_list = reac_list[order(-reac_list$Fatal_rate), ]
  
  temp = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                               ICI_clinical_Backup$outc_cod_num == 1 &
                               !is.na(ICI_clinical_Backup$duration) &
                               ICI_clinical_Backup$Drug == ICIs_item[i] & 
                               ICI_clinical_Backup$Reac %in% reac_list$Reac,]
  
  sum_m = length(unique(ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                                              ICI_clinical_Backup$outc_cod_num == 1 &
                                              !is.na(ICI_clinical_Backup$duration) &
                                              ICI_clinical_Backup$Drug == ICIs_item[i],]$primaryid))
  print(length(unique(temp$primaryid))/ sum_m)
  
  
  venn_list = list()
  for (item in reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac){
    venn_list = append(venn_list, list(temp[temp$Reac == item, ]$primaryid))
  }
  names(venn_list) = str_to_sentence(tolower(reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac))
  set.seed(1234)
  print(plot(euler(venn_list,
                   shape = "circle"),  
             quantities = list(type = c("counts"),cex=1), 
             labels=list(cex=1), 
             edges = list(col = "black", lex = 2),
             fills = list(fill = color_list[match(reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac, Reac_ref)],alpha=0.8), # 填充的颜色和透明度
             # legend = list(side = "right")      
  ))
}

dev.off()

##############################################################################################################################
##############################################################################################################################
# OR fatal forest adjusted by age sex

forest_data = data.frame()
for (i in 1:8){
  forest_data_temp = data.frame()
  for (Reac_term in unique(data[data$Drug == ICIs_item[i],]$Reac)){
    temp = ICI_clinical[ICI_clinical$Drug == ICIs_item[i] & ICI_clinical$duration <= 90 & !is.na(ICI_clinical$duration), ]
    temp[temp$primaryid %in% ICI_clinical_Backup[ICI_clinical_Backup$Reac == Reac_term,]$primaryid, 'Reac_ready'] = 1
    temp[is.na(temp$Reac_ready), 'Reac_ready'] = 0
    col_1 = sum(temp$Reac_ready==1)
    col_2 = sum(temp$Reac_ready==1 & temp$outc_cod_num == 1)
    col_3 = sum(temp$Reac_ready==0)
    col_4 = sum(temp$Reac_ready==0 & temp$outc_cod_num == 1)
    
    temp$Reac_ready = factor(temp$Reac_ready, levels = c(0, 1))
    temp$outc_cod_num = factor(temp$outc_cod_num, levels = c(0, 1))
    
    if (length(unique(temp$Reac_ready)) == 2 & length(unique(temp$outc_cod_num)) == 2){
      temp[temp$sex == 'Male', 'sex_num'] = 1
      temp[temp$sex != 'Male', 'sex_num'] = 0
      
      OR<-glm(formula=outc_cod_num ~Reac_ready+age_num+sex_num, data = temp, family = binomial)
      coef<-as.vector(coef(OR)[2])
      coef_CI<-as.vector(confint(OR)[2,])
      p_value = as.vector(summary(OR)$coefficients[,4][2])
      forest_data_temp = rbind(forest_data_temp, data.frame(Drug = ICIs_item[i], Reac = str_to_sentence(tolower(Reac_term)), 
                                                            OR = exp(coef), L_CI = exp(coef_CI[1]), U_CI = exp(coef_CI[2]), p_value = p_value,
                                                            col_1 = col_1, col_2 = col_2, col_3 = col_3, col_4 = col_4)
      )
    }
  }
  forest_data_temp$p_value_adjust = p.adjust(forest_data_temp$p_value, method = 'fdr')
  forest_data = rbind(forest_data, forest_data_temp)
}

# OR fatal forest adjusted by age sex duration

forest_data_duration = data.frame()
for (i in 1:8){
  forest_data_temp = data.frame()
  for (Reac_term in unique(data[data$Drug == ICIs_item[i],]$Reac)){
    temp = ICI_clinical[ICI_clinical$Drug == ICIs_item[i] & ICI_clinical$duration <= 90 & !is.na(ICI_clinical$duration), ]
    temp[temp$primaryid %in% ICI_clinical_Backup[ICI_clinical_Backup$Reac == Reac_term,]$primaryid, 'Reac_ready'] = 1
    temp[is.na(temp$Reac_ready), 'Reac_ready'] = 0
    col_1 = sum(temp$Reac_ready==1)
    col_2 = sum(temp$Reac_ready==1 & temp$outc_cod_num == 1)
    col_3 = sum(temp$Reac_ready==0)
    col_4 = sum(temp$Reac_ready==0 & temp$outc_cod_num == 1)
    
    temp$Reac_ready = factor(temp$Reac_ready, levels = c(0, 1))
    temp$outc_cod_num = factor(temp$outc_cod_num, levels = c(0, 1))
    
    if (length(unique(temp$Reac_ready)) == 2 & length(unique(temp$outc_cod_num)) == 2){
      temp[temp$sex == 'Male', 'sex_num'] = 1
      temp[temp$sex != 'Male', 'sex_num'] = 0
      
      OR<-glm(formula=outc_cod_num ~Reac_ready+age_num+sex_num+duration, data = temp, family = binomial)
      coef<-as.vector(coef(OR)[2])
      coef_CI<-as.vector(confint(OR)[2,])
      p_value = as.vector(summary(OR)$coefficients[,4][2])
      forest_data_temp = rbind(forest_data_temp, data.frame(Drug = ICIs_item[i], Reac = str_to_sentence(tolower(Reac_term)), 
                                                            OR = exp(coef), L_CI = exp(coef_CI[1]), U_CI = exp(coef_CI[2]), p_value = p_value,
                                                            col_1 = col_1, col_2 = col_2, col_3 = col_3, col_4 = col_4)
      )
    }
  }
  forest_data_temp$p_value_adjust = p.adjust(forest_data_temp$p_value, method = 'fdr')
  forest_data_duration = rbind(forest_data_duration, forest_data_temp)
}

###########

forest_data_out = merge(forest_data[, c(1:2, 7:10, 3:6, 11)], forest_data_duration[, c(1:2, 3:6, 11)], by = c('Drug', 'Reac'))
SourceData = createWorkbook()
SourceData_info(fig_data = forest_data_out,
                colnames = colnames(forest_data_out),
                cols = 1:ncol(forest_data_out),
                sheet = 'sheet',
                widths = c(10, 25, rep(12, ncol(forest_data_out)-2)))
saveWorkbook(SourceData, '毕业_Table_forest_data.xlsx', overwrite = TRUE)


Reac_drug_1 = forest_data_duration[!is.na(forest_data_duration$L_CI) & !is.na(forest_data_duration$U_CI) & forest_data_duration$p_value_adjust < 0.05, c('Drug', 'Reac')]
Reac_drug_2 = forest_data[!is.na(forest_data$L_CI) & !is.na(forest_data$U_CI) & forest_data$p_value_adjust < 0.05, c('Drug', 'Reac')]
Reac_drug = full_join(Reac_drug_1, Reac_drug_2, by = c('Drug', 'Reac'))

forest_data_outer_sig = merge(left_join(Reac_drug, forest_data, by = c('Drug', 'Reac')), left_join(Reac_drug, forest_data_duration, by = c('Drug', 'Reac')), c('Drug', 'Reac'))
pval_Q = c()
for (i in 1:dim(forest_data_outer_sig)[1]){
  ready = forest_data_outer_sig[i,]
  
  OR = as.vector(c(ready$OR.x, ready$OR.y))
  L_CI = as.vector(c(ready$L_CI.x, ready$L_CI.y))
  U_CI = as.vector(c(ready$U_CI.x, ready$U_CI.y))
  
  pval_Q = c(pval_Q, summary(metagen(OR, lower = L_CI, upper = U_CI, studlab = c('Simple', 'Full'), sm = "OR", transf = FALSE))$pval.Q)
}
forest_data_outer_sig$pval_Q = pval_Q
forest_data_outer_sig$OR_m = (forest_data_outer_sig$OR.x+forest_data_outer_sig$OR.y) / 2
forest_data_outer_sig = forest_data_outer_sig[order(-forest_data_outer_sig$OR_m), ]

forest_plot_self = function(forest_input, drug){
  tm <- forest_theme(
    base_size = 10, 
    ci_pch = 16,
    ci_col = "blue4",
    ci_fill = "blue4",
    ci_alpha = 0.8,
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = 0.2,
    
    refline_lwd = 1,
    refline_lty = "dashed",
    refline_col = "grey20",
    
    vertline_lwd = 1,
    vertline_lty = "dashed",
    vertline_col = "grey20",
    
    footnote_cex = 0.6,
    footnote_fontface = "italic",
    footnote_col = "red4"
  )
  
  forest_display = forest_input[forest_input$Drug == drug, c(2, 7:10, 17, 16, 15, 12)]
  colnames(forest_display) = c('related adverse events', 'WI\nTotal', 'WI\nDeath', 
                               'NWI\nTotal', 'NWI\nDeath', 
                               'log2(OR 95%CI)', 'OR 95%CI', 'FDR P', 'Hete P')
  
  p = forestploter::forest(
    forest_display,  
    est = log2(forest_input[forest_input$Drug == drug, 3]), 
    lower = log2(forest_input[forest_input$Drug == drug, 4]), 
    upper = log2(forest_input[forest_input$Drug == drug, 5]),
    sizes = forest_input[forest_input$Drug == drug, ]$size,
    ci_column = 6,
    ref_line = 0,
    arrow_lab = c("Low fatal risk", "High fatal risk"),
    xlim = c(-5.5, 5.5),
    ticks_at = seq(-4,4,2),
    theme = tm
  ) 
  
  p1 <- add_border(p, part = "header", where = "bottom")
  p1
}


pdf('毕业_PT_forest.pdf', height = 12, width = 8)
forest_input = forest_data_outer_sig[,c(1:11, 21)]
forest_input$OR.x = round(forest_input$OR.x, 2)
forest_input$L_CI_str = sprintf("%0.2f", round(forest_input$L_CI.x, 2))
forest_input$U_CI_str = sprintf("%0.2f", round(forest_input$U_CI.x, 2))
forest_input$p_value_adjust_str = as.character(signif(forest_input$p_value_adjust.x, 3))
forest_input[forest_input$p_value.x >=0.002, 'p_value_adjust_str'] = as.character(round(forest_input[forest_input$p_value.x >=0.002, ]$p_value_adjust.x, 3))
forest_input$OR_str = paste(forest_input$OR.x,' [', forest_input$L_CI_str, ', ', forest_input$U_CI_str, ']', sep='')
forest_input$forest_rcol = '                                       '
forest_input$size = ifelse(forest_input$p_value_adjust.x < 0.05, 0.8, 0.3)

print(forest_plot_self(forest_input, ICIs_item[1]))
print(forest_plot_self(forest_input, ICIs_item[2]))
print(forest_plot_self(forest_input, ICIs_item[3]))
print(forest_plot_self(forest_input, ICIs_item[4]))
print(forest_plot_self(forest_input, ICIs_item[5]))
print(forest_plot_self(forest_input, ICIs_item[6]))
print(forest_plot_self(forest_input, ICIs_item[7]))
print(forest_plot_self(forest_input, ICIs_item[8]))

forest_input = forest_data_outer_sig[,c(1:2, 12:20, 21)]
forest_input$OR.y = round(forest_input$OR.y, 2)
forest_input$L_CI_str = sprintf("%0.2f", round(forest_input$L_CI.y, 2))
forest_input$U_CI_str = sprintf("%0.2f", round(forest_input$U_CI.y, 2))
forest_input$p_value_adjust_str = as.character(signif(forest_input$p_value_adjust.y, 3))
forest_input[forest_input$p_value.y >=0.002, 'p_value_adjust_str'] = as.character(round(forest_input[forest_input$p_value.y >=0.002, ]$p_value_adjust.y, 3))
forest_input$OR_str = paste(forest_input$OR.y,' [', forest_input$L_CI_str, ', ', forest_input$U_CI_str, ']', sep='')
forest_input$forest_rcol = '                                       '
forest_input$size = ifelse(forest_input$p_value_adjust.y < 0.05, 0.8, 0.3)

print(forest_plot_self(forest_input, ICIs_item[1]))
print(forest_plot_self(forest_input, ICIs_item[2]))
print(forest_plot_self(forest_input, ICIs_item[3]))
print(forest_plot_self(forest_input, ICIs_item[4]))
print(forest_plot_self(forest_input, ICIs_item[5]))
print(forest_plot_self(forest_input, ICIs_item[6]))
print(forest_plot_self(forest_input, ICIs_item[7]))
print(forest_plot_self(forest_input, ICIs_item[8]))

dev.off()

#######################################################################################################################################

Fatal_rate_3 = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 & !is.na(ICI_clinical_Backup$duration), ] %>% 
  group_by(Drug,Reac) %>% 
  summarise(Fatal_rate = sum(outc_cod_num)/n()*100, Report_count = n(), Fatal_count = sum(outc_cod_num))

table_data = data
table_data$ICI = str_split_fixed(table_data$Drug, '_', 2)[,1]
table_data$Cancer = str_split_fixed(table_data$Drug, '_', 2)[,2]

table_data = merge(table_data, Fatal_rate_3, by = c('Drug', 'Reac'))
table_data$ROR = round(table_data$ROR,1)
table_data$EBGM = round(table_data$EBGM,1)

table_data$fatal = paste(table_data$Fatal_count, " (", round(table_data$Fatal_rate, 1), '%)',sep='')


pdf('毕业_PT_VENN_fatality.pdf', height = 8, width = 8)

for (i in c(2:3, 5:8)){
  color_list = c(colorRampPalette(c('#a3de83', '#ffffff'))(3)[1:2], colorRampPalette(c('#ffffff', '#ffa447'))(8))
  
  forest_input = forest_data[forest_data$Drug %in% ICIs_item[i] & 
                               !is.na(forest_data$L_CI) & 
                               !is.na(forest_data$U_CI) &
                               forest_data$p_value_adjust < 0.05, ]
  
  reac_list = table_data[table_data$Drug %in% ICIs_item[i] & table_data$Reac %in% toupper(forest_input$Reac),c('Reac', 'Fatal_rate')]
  reac_list = reac_list[order(-reac_list$Fatal_rate), ]
  
  temp = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                               # ICI_clinical_Backup$outc_cod_num == 1 &
                               !is.na(ICI_clinical_Backup$duration) &
                               ICI_clinical_Backup$Drug == ICIs_item[i],]
  
  venn_list = list()
  for (item in reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac){
    venn_list = append(venn_list, list(temp[temp$Reac == item, ]$primaryid))
  }
  names(venn_list) = str_to_sentence(tolower(reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac))
  set.seed(1234)
  print(length(venn_list))
  print(plot(euler(venn_list,
                   shape = "circle"),  
             quantities = list(type = c("counts"),cex=1), 
             labels=list(cex=1), 
             edges = list(col = "black", lex = 2),
             fills = list(fill = color_list[ as.integer(reac_list$Fatal_rate / 10) +1],alpha=0.8),
             legend = list(side = "right")
  ))
}
dev.off()

pdf('毕业_PT_VENN_fatality_0.pdf', height = 8, width = 8)

for (i in c(1,4)){
  color_list = c(colorRampPalette(c('#a3de83', '#ffffff'))(3)[1:2], colorRampPalette(c('#ffffff', '#ffa447'))(8))
  
  forest_input = forest_data[forest_data$Drug %in% ICIs_item[i] & 
                               !is.na(forest_data$L_CI) & 
                               !is.na(forest_data$U_CI) &
                               forest_data$p_value_adjust < 0.001, ]
  
  reac_list = table_data[table_data$Drug %in% ICIs_item[i] & table_data$Reac %in% toupper(forest_input$Reac),c('Reac', 'Fatal_rate')]
  reac_list = reac_list[order(-reac_list$Fatal_rate), ]
  
  temp = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                               # ICI_clinical_Backup$outc_cod_num == 1 &
                               !is.na(ICI_clinical_Backup$duration) &
                               ICI_clinical_Backup$Drug == ICIs_item[i],]
  
  venn_list = list()
  for (item in reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac){
    venn_list = append(venn_list, list(temp[temp$Reac == item, ]$primaryid))
  }
  names(venn_list) = str_to_sentence(tolower(reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac))
  set.seed(1234)
  print(length(venn_list))
  print(plot(euler(venn_list,
                   shape = "circle"),  
             quantities = list(type = c("counts"),cex=1), 
             labels=list(cex=1), 
             edges = list(col = "black", lex = 2),
             fills = list(fill = color_list[ as.integer(reac_list$Fatal_rate/ 10) +1],alpha=0.8),
             legend = list(side = "right")
  ))
}
dev.off()

pdf('毕业_PT_VENN_MG.pdf', height = 8, width = 8)

for (i in c(1:8)){
  color_list = c(colorRampPalette(c('#a3de83', '#ffffff'))(3)[1:2], colorRampPalette(c('#ffffff', '#ffa447'))(8))
  
  reac_list = table_data[table_data$Drug %in% ICIs_item[i] & 
                           table_data$Reac %in% c('RESPIRATORY FAILURE', 'MYASTHENIA GRAVIS', 'MYOSITIS', 'MYOCARDITIS'),
                         c('Reac', 'Fatal_rate')]
  reac_list = reac_list[order(-reac_list$Fatal_rate), ]
  
  temp = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                               # ICI_clinical_Backup$outc_cod_num == 1 &
                               !is.na(ICI_clinical_Backup$duration) &
                               ICI_clinical_Backup$Drug == ICIs_item[i],]
  
  venn_list = list()
  for (item in reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac){
    venn_list = append(venn_list, list(temp[temp$Reac == item, ]$primaryid))
  }
  names(venn_list) = str_to_sentence(tolower(reac_list[reac_list$Reac %in% toupper(temp$Reac), ]$Reac))
  set.seed(1234)
  print(length(venn_list))
  print(plot(euler(venn_list,
                   shape = "circle"),  
             quantities = list(type = c("counts"),cex=1), 
             labels=list(cex=1), 
             edges = list(col = "black", lex = 2),
             fills = list(fill = color_list[ as.integer(reac_list$Fatal_rate / 10) +1],alpha=0.8),
             legend = list(side = "right")
  ))
}
dev.off()

##############################################################################################################################
##############################################################################################################################
# hyper
# XGboost regression 
hyper_sub = data.frame()
for (i in 1:8){
  forest_input = forest_data[forest_data$Drug %in% ICIs_item[i] & 
                               !is.na(forest_data$L_CI) & 
                               !is.na(forest_data$U_CI) &
                               forest_data$p_value_adjust < 0.05,]
  
  logistic_data = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                                        # ICI_clinical_Backup$outc_cod_num == 1 &
                                        !is.na(ICI_clinical_Backup$duration) &
                                        ICI_clinical_Backup$Drug == ICIs_item[i] &
                                        ICI_clinical_Backup$Reac %in% toupper(forest_input$Reac),]
  logistic_data$value = 1
  logistic_data[logistic_data$sex == 'Male', 'sex_num'] = 1
  logistic_data[logistic_data$sex != 'Male', 'sex_num'] = 0
  
  logistic_data = reshape2::dcast(logistic_data, primaryid+age_num+sex_num+duration+outc_cod_num~Reac, value.var = 'value')
  logistic_data[,6:dim(logistic_data)[2]] = replace(logistic_data[,6:dim(logistic_data)[2]], is.na(logistic_data[,6:dim(logistic_data)[2]]), 0)
  
  set.seed(123)
  
  # xgboost fine
  fitControl <- trainControl(method = "cv",
                             number = 10,
                             verboseIter = FALSE,
                             search = "random" # 随机搜索，也可设定 "grid" 网格搜索
  )
  caret_xgb <- train(outc_cod_num ~. ,
                     data = logistic_data[,2:dim(logistic_data)[2]],
                     method="xgbTree",# xgbLinear
                     trainControl=fitControl
  )
  hyper_sub = rbind(hyper_sub, cbind(Drug = ICIs_item[i], caret_xgb$bestTune))
}

# XGboost regression 
hyper_sub_cv = data.frame()
out_sub = data.frame()
out_sub_train = data.frame()
holdout_sub = data.frame()
for (i in 1:8){
  print(i)
  forest_input = forest_data[forest_data$Drug %in% ICIs_item[i] & 
                               !is.na(forest_data$L_CI) & 
                               !is.na(forest_data$U_CI) &
                               forest_data$p_value_adjust < 0.05,]
  
  logistic_data = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                                        # ICI_clinical_Backup$outc_cod_num == 1 &
                                        !is.na(ICI_clinical_Backup$duration) &
                                        ICI_clinical_Backup$Drug == ICIs_item[i] &
                                        ICI_clinical_Backup$Reac %in% toupper(forest_input$Reac),]
  logistic_data$value = 1
  logistic_data[logistic_data$sex == 'Male', 'sex_num'] = 1
  logistic_data[logistic_data$sex != 'Male', 'sex_num'] = 0
  
  logistic_data = reshape2::dcast(logistic_data, primaryid+age_num+sex_num+duration+outc_cod_num~Reac, value.var = 'value')
  logistic_data[,6:dim(logistic_data)[2]] = replace(logistic_data[,6:dim(logistic_data)[2]], is.na(logistic_data[,6:dim(logistic_data)[2]]), 0)
  
  set.seed(123)
  folds <- createFolds(logistic_data$primaryid, k = 10, list = TRUE, returnTrain = FALSE)
  
  for (k in 1:length(folds)){
    idx = folds[[k]]
    
    set.seed(1234)
    fitControl <- trainControl(method = "cv",
                               number = 10,
                               verboseIter = FALSE,
                               search = "random" # 随机搜索，也可设定 "grid" 网格搜索
    )
    caret_xgb <- train(outc_cod_num ~. ,
                       data = logistic_data[,setdiff(2:dim(logistic_data)[2], 4)],
                       method="xgbTree",# xgbLinear
                       trainControl=fitControl
    )
    hyper_fold = caret_xgb$bestTune
    
    hyper_sub_cv = rbind(hyper_sub_cv, cbind(Drug = ICIs_item[i], fold = k, hyper_fold))
    
    dtrain <- xgb.DMatrix(data = as.matrix(logistic_data[-idx,c(2:3,6:dim(logistic_data)[2])],), label = logistic_data$outc_cod_num[-idx])
    params <- list(
      objective = "binary:logistic",
      eval_metric = "error",     
      max_depth = hyper_fold$max_depth,
      eta = hyper_fold$eta,
      gamma = hyper_fold$gamma,
      colsample_bytree = hyper_fold$colsample_bytree
    )
    xgb_model <- xgboost(data = dtrain, params = params, nrounds = hyper_fold$nrounds)
    pred_value = data.frame(Drug = ICIs_item[i],
                            out_de = logistic_data$outc_cod_num[idx],
                            duration = logistic_data$duration[idx],
                            primaryid = logistic_data$primaryid[idx],
                            pred = predict(xgb_model, xgb.DMatrix(data = as.matrix(logistic_data[idx,c(2:3,6:dim(logistic_data)[2])],))),
                            fold = k)
    out_sub = rbind(out_sub, pred_value)
    pred_value = data.frame(Drug = ICIs_item[i],
                            out_de = logistic_data$outc_cod_num[-idx],
                            duration = logistic_data$duration[-idx],
                            primaryid = logistic_data$primaryid[-idx],
                            pred = predict(xgb_model, xgb.DMatrix(data = as.matrix(logistic_data[-idx,c(2:3,6:dim(logistic_data)[2])],))),
                            fold = k)
    out_sub_train = rbind(out_sub_train, pred_value)
  }
  
  # all data model holdout test
  print('holdout')
  hyper = hyper_sub[hyper_sub$Drug == ICIs_item[i], ]
  
  dtrain <- xgb.DMatrix(data = as.matrix(logistic_data[,c(2:3,6:dim(logistic_data)[2])],), label = logistic_data$outc_cod_num)
  params <- list(
    objective = "binary:logistic",
    eval_metric = "error",     
    max_depth = hyper$max_depth,
    eta = hyper$eta,
    gamma = hyper$gamma,
    colsample_bytree = hyper$colsample_bytree
  )
  xgb_model <- xgboost(data = dtrain, params = params, nrounds = hyper$nrounds)
  
  logistic_holdout_data = ICI_holdout_Backup[ICI_holdout_Backup$duration <= 90 &
                                               # ICI_clinical_Backup$outc_cod_num == 1 &
                                               !is.na(ICI_holdout_Backup$duration) &
                                               ICI_holdout_Backup$Drug == ICIs_item[i] &
                                               ICI_holdout_Backup$Reac %in% toupper(forest_input$Reac),]
  if (dim(logistic_holdout_data)[1] > 0){
    logistic_holdout_data$value = 1
    logistic_holdout_data[logistic_holdout_data$sex == 'Male', 'sex_num'] = 1
    logistic_holdout_data[logistic_holdout_data$sex != 'Male', 'sex_num'] = 0
    
    logistic_holdout_data = reshape2::dcast(logistic_holdout_data, primaryid+age_num+sex_num+duration+outc_cod_num~Reac, value.var = 'value')
    logistic_holdout_data[,6:dim(logistic_holdout_data)[2]] = replace(logistic_holdout_data[,6:dim(logistic_holdout_data)[2]], is.na(logistic_holdout_data[,6:dim(logistic_holdout_data)[2]]), 0)
    
    for (name in colnames(logistic_data)){
      if (!name %in% colnames(logistic_holdout_data)){
        logistic_holdout_data[,name] = 0
      }
    }
    logistic_holdout_data = logistic_holdout_data[,colnames(logistic_data)]
    
    pred_hold_out_value = data.frame(Drug = ICIs_item[i],
                                     out_de = logistic_holdout_data$outc_cod_num,
                                     duration = logistic_holdout_data$duration,
                                     primaryid = logistic_holdout_data$primaryid,
                                     pred = predict(xgb_model, xgb.DMatrix(data = as.matrix(logistic_holdout_data[,c(2:3,6:dim(logistic_holdout_data)[2])])))
    )
    
    holdout_sub = rbind(holdout_sub, pred_hold_out_value)
  }
}

#######################################################################################
# XGboost data generated

hyper_sub$fold = 'holdout'
hyper = rbind(hyper_sub, hyper_sub_cv)
hyper$fold = factor(hyper$fold, levels = c(1:10, 'holdout'))
hyper = hyper[order(hyper$Drug, hyper$fold), ]

SourceData = createWorkbook()
SourceData_info(fig_data =hyper,
                colnames = colnames(hyper),
                cols = 1:ncol(hyper),
                sheet = 'sheet',
                widths = rep(12, ncol(hyper)))
saveWorkbook(SourceData, paste('毕业_Table_hyper_test.xlsx', sep = ''), overwrite = TRUE)


for (i in 1:8){
  print(i)
  forest_input = forest_data[forest_data$Drug %in% ICIs_item[i] & 
                               !is.na(forest_data$L_CI) & 
                               !is.na(forest_data$U_CI) &
                               forest_data$p_value_adjust < 0.05,]
  
  logistic_data = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                                        # ICI_clinical_Backup$outc_cod_num == 1 &
                                        !is.na(ICI_clinical_Backup$duration) &
                                        ICI_clinical_Backup$Drug == ICIs_item[i] &
                                        ICI_clinical_Backup$Reac %in% toupper(forest_input$Reac),]
  logistic_data$value = 1
  logistic_data[logistic_data$sex == 'Male', 'sex_num'] = 1
  logistic_data[logistic_data$sex != 'Male', 'sex_num'] = 0
  
  logistic_data = reshape2::dcast(logistic_data, primaryid+age_num+sex_num+duration+outc_cod_num~Reac, value.var = 'value')
  logistic_data[,6:dim(logistic_data)[2]] = replace(logistic_data[,6:dim(logistic_data)[2]], is.na(logistic_data[,6:dim(logistic_data)[2]]), 0)
  
  set.seed(123)
  folds <- createFolds(logistic_data$primaryid, k = 10, list = TRUE, returnTrain = FALSE)
  for (k in 1:10){
    logistic_data[folds[[k]], 'fold'] = k
  }
  logistic_data['Drug'] = ICIs_item[i] 
  logistic_data['out_de'] = logistic_data['outc_cod_num']
  logistic_data = logistic_data[, !colnames(logistic_data) %in% c('outc_cod_num')]
  
  logistic_data = merge(out_sub, logistic_data, by = c('Drug', 'primaryid', 'duration', 'out_de','fold'))
  logistic_data = logistic_data[order(logistic_data$fold, -logistic_data$pred), ]
  
  ########
  logistic_holdout_data = ICI_holdout_Backup[ICI_holdout_Backup$duration <= 90 &
                                               # ICI_clinical_Backup$outc_cod_num == 1 &
                                               !is.na(ICI_holdout_Backup$duration) &
                                               ICI_holdout_Backup$Drug == ICIs_item[i] &
                                               ICI_holdout_Backup$Reac %in% toupper(forest_input$Reac),]

  logistic_holdout_data$value = 1
  logistic_holdout_data[logistic_holdout_data$sex == 'Male', 'sex_num'] = 1
  logistic_holdout_data[logistic_holdout_data$sex != 'Male', 'sex_num'] = 0
  
  logistic_holdout_data = reshape2::dcast(logistic_holdout_data, primaryid+age_num+sex_num+duration+outc_cod_num~Reac, value.var = 'value')
  logistic_holdout_data[,6:dim(logistic_holdout_data)[2]] = replace(logistic_holdout_data[,6:dim(logistic_holdout_data)[2]], is.na(logistic_holdout_data[,6:dim(logistic_holdout_data)[2]]), 0)
  
  for (name in colnames(logistic_data)){
    if (!name %in% colnames(logistic_holdout_data)){
      logistic_holdout_data[,name] = 0
    }
  }
  logistic_holdout_data['Drug'] = ICIs_item[i] 
  logistic_holdout_data['fold'] = 'holdout'
  logistic_holdout_data = logistic_holdout_data[,colnames(logistic_holdout_data)]
  logistic_holdout_data['out_de'] = logistic_holdout_data['outc_cod_num']
  logistic_holdout_data = logistic_holdout_data[, !colnames(logistic_holdout_data) %in% c('outc_cod_num', 'pred')]
  
  logistic_holdout_data = merge(holdout_sub, logistic_holdout_data, by = c('Drug', 'primaryid', 'duration', 'out_de'))
  logistic_holdout_data = logistic_holdout_data[order(-logistic_holdout_data$pred), ]
  
  out_data = rbind(logistic_data, logistic_holdout_data)
  
  out_data$fold = factor(out_data$fold, levels = c(1:10, 'holdout'))
  out_data = out_data[order(out_data$Drug, out_data$fold), ]
  
  SourceData = createWorkbook()
  SourceData_info(fig_data =out_data,
                  colnames = colnames(out_data),
                  cols = 1:ncol(out_data),
                  sheet = 'sheet',
                  widths = c(10, 25, rep(12, ncol(out_data)-2)))
  saveWorkbook(SourceData, paste('毕业_Table_10_cv_test_',  gsub('/','_', ICIs_item[i]), '.xlsx', sep = ''), overwrite = TRUE)
}


####################################################################
# SHAP

SHAP_PLOT = function(i){
  forest_input = forest_data[forest_data$Drug %in% ICIs_item[i] & 
                               !is.na(forest_data$L_CI) & 
                               !is.na(forest_data$U_CI) &
                               forest_data$p_value_adjust < 0.05,]
  
  logistic_data = ICI_clinical_Backup[ICI_clinical_Backup$duration <= 90 &
                                        # ICI_clinical_Backup$outc_cod_num == 1 &
                                        !is.na(ICI_clinical_Backup$duration) &
                                        ICI_clinical_Backup$Drug == ICIs_item[i] &
                                        ICI_clinical_Backup$Reac %in% toupper(forest_input$Reac),]
  logistic_data$value = 1
  logistic_data[logistic_data$sex == 'Male', 'sex_num'] = 1
  logistic_data[logistic_data$sex != 'Male', 'sex_num'] = 0
  
  logistic_data = reshape2::dcast(logistic_data, primaryid+age_num+sex_num+duration+outc_cod_num~Reac, value.var = 'value')
  logistic_data[,6:dim(logistic_data)[2]] = replace(logistic_data[,6:dim(logistic_data)[2]], is.na(logistic_data[,6:dim(logistic_data)[2]]), 0)
  
  dtrain <- xgb.DMatrix(data = as.matrix(logistic_data[,c(2:3,6:dim(logistic_data)[2])],), label = logistic_data$outc_cod_num)
  hyper = hyper_sub[hyper_sub$Drug == ICIs_item[i], ]
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "error",     
    max_depth = hyper$max_depth,
    eta = hyper$eta,
    gamma = hyper$gamma,
    colsample_bytree = hyper$colsample_bytree
  )
  xgb_model <- xgboost(data = dtrain, params = params, nrounds = hyper$nrounds)
  
  # SHAP
  library(shapviz)
  shp <- shapviz(xgb_model, X_pred =as.matrix(logistic_data[,c(2:3,6:dim(logistic_data)[2])],))
  colnames(shp) = str_to_title(colnames(shp))
  sv_importance(shp, max_display = dim(logistic_data)[2]-3)+
    geom_vline(xintercept = 0.2, linetype ='dotdash', color = 'red')+
    theme_bw()+
    scale_x_continuous(breaks = seq(0, 1, 0.1))
}

pdf('毕业_Model_SHAP.pdf', height =8, width = 10)
for(i in 1:8){
  print(SHAP_PLOT(i))
}
dev.off()


#######################################################################################
data_generate_aupr = function(pred, group_name){
  perf = performance(pred, "prec", 'rec')
  prec = perf@"y.values"
  rec = perf@"x.values"
  perf = performance(pred, 'aucpr')
  aupr = perf@"y.values"
  
  test_temp = data.frame(prec, rec)
  colnames(test_temp) = c('prec', 'rec')
  idx = which.max(2*test_temp$prec*test_temp$rec / (test_temp$prec+test_temp$rec))
  
  part = data.frame(group_name, prec, rec, aupr, sprintf("%0.3f", aupr),
                    unlist(pred@cutoffs[1])[idx], test_temp$prec[idx], test_temp$rec[idx])
  colnames(part) = c('group', 'prec', 'rec', 'aupr', 'aupr_str', 'cut_off', 'prec_cut', 'rec_cut')
  part
}

aupr_plot = function(out, color_set){
  p = ggplot()+
    geom_line(data=f1_base, aes(x=recall, y=precision, group=f1), color='gray')+
    geom_text(data=f1_score_label,aes(x=0.95, y = y+0.03, label=f1), size=3, color='gray')+
    geom_line(data=out, aes(x=rec, y=prec, group=group, colour=group), linewidth=1, alpha=1)+
    geom_point(data=out[!duplicated(out[,c('group', 'cut_off', 'rec_cut', 'prec_cut')]),],
               aes(x=rec_cut, y=prec_cut, group=group, colour=group), size=2, shape=16)+
    # geom_text_repel(data=cut, aes(x=rec, y=prec, label = paste('Cut-off=', sprintf("%0.2f", cut_off), sep='')), size = 3)+
    scale_color_manual(name = '',
                       values = color_set,
                       breaks = unique(out$group))+
    scale_x_continuous(limits=c(0, 1), expand=c(0, 0), breaks=seq(0, 1, 0.2))+
    scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks=seq(0, 1, 0.2))+
    theme_bw()+
    theme(axis.ticks.length.x.bottom = unit(-0.1, "cm"),
          axis.ticks.length.y.left = unit(-0.1, "cm"),
          axis.text = element_text(color = 'black'),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          legend.justification = c(1,0)
    )+
    xlab("Recall")+
    ylab("Precision")
  
  strip = out[!duplicated(out$group),]
  p_strip<-ggplot(strip, aes(x=group, y=aupr, fill=group)) +
    geom_col(alpha=1)+
    scale_fill_manual(name = '',
                      values = color_set,
                      breaks = unique(out$group))+
    geom_text(aes(label=aupr_str), colour="black", angle=90, hjust = 1.5, size = 3)+
    geom_text(aes(label=group, y=0.2), colour="black", angle=90, size = 3)+
    scale_y_continuous(limits=c(0, 1.00), expand=c(0, 0), breaks=seq(0, 1, 0.2))+
    theme_bw()+
    theme_classic()+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(color = 'black'),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y.left = unit(0.1, "cm"),
          legend.position="none")+
    labs(x=NULL, y="AUPR")
  
  p_strip <- ggplotGrob(p_strip)
  return(p + annotation_custom(p_strip,xmin = 0.02,xmax = 0.7,ymin = 0.01,ymax = 0.7))
}

data_generate_auc = function(pred, group_name){
  perf = performance(pred, "tpr", 'fpr')
  tpr = perf@"y.values"
  fpr = perf@"x.values"
  perf = performance(pred, 'auc')
  auc = perf@"y.values"
  
  test_temp = data.frame(tpr, fpr)
  colnames(test_temp) = c('tpr', 'fpr')
  idx = which.max(test_temp$tpr - test_temp$fpr)
  
  part = data.frame(group_name, tpr, fpr, auc, sprintf("%0.3f", auc),
                    unlist(pred@cutoffs[1])[idx], test_temp$tpr[idx], test_temp$fpr[idx])
  colnames(part) = c('group', 'tpr', 'fpr', 'auc', 'auc_str', 'cut_off', 'tpr_cut', 'fpr_cut')
  part
}

auc_plot = function(out, color_set){
  p = ggplot()+
    geom_line(data=out, aes(x=fpr, y=tpr, group=group, colour=group), linewidth=1, alpha=1)+
    geom_point(data=out[!duplicated(out[,c('group', 'cut_off', 'tpr_cut', 'fpr_cut')]),],
               aes(x=fpr_cut, y=tpr_cut, group=group, colour=group), size=2, shape=16)+
    scale_color_manual(name = '',
                       values = color_set,
                       breaks = unique(out$group))+
    scale_x_continuous(limits=c(0, 1), expand=c(0, 0), breaks=seq(0, 1, 0.2))+
    scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks=seq(0, 1, 0.2))+
    geom_abline(intercept = 0, slope = 1, color = '#d5d5d5', linetype = 'dashed')+
    theme_bw()+
    theme(axis.ticks.length.x.bottom = unit(-0.1, "cm"),
          axis.ticks.length.y.left = unit(-0.1, "cm"),
          axis.text = element_text(color = 'black'),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          legend.justification = c(1,0)
    )+
    xlab("False positive rate")+
    ylab("True positive rate")
  
  strip = out[!duplicated(out$group),]
  p_strip<-ggplot(strip, aes(x=group, y=auc, fill=group)) +
    geom_col(alpha=1)+
    scale_fill_manual(name = '',
                      values = color_set,
                      breaks = unique(out$group))+
    geom_text(aes(label=auc_str), colour="black", angle=90, hjust = 1.5, size = 3)+
    geom_text(aes(label=group, y=0.2), colour="black", angle=90, size = 3)+
    scale_y_continuous(limits=c(0, 1.00), expand=c(0, 0), breaks=seq(0, 1, 0.2))+
    theme_bw()+
    theme_classic()+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(color = 'black'),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y.left = unit(0.1, "cm"),
          legend.position="none")+
    labs(x=NULL, y="AUROC")
  
  p_strip <- ggplotGrob(p_strip)
  return(p + annotation_custom(p_strip,xmin = 0.35,xmax = 1,ymin = 0.01,ymax = 0.7))
}

########################## 10-fold-cv

Data = data.frame(group='', prec=0, rec=0, aupr=0, aupr_str='', cut_off=0, prec_cut=0, rec_cut=0)[-1,]
for (i in 1:8){
  Data = rbind(Data, data_generate_aupr(pred = prediction(out_sub[out_sub$Drug == ICIs_item[i],]$pred, 
                                                          out_sub[out_sub$Drug == ICIs_item[i],]$out_de),
                                        group_name = ICIs_item[i]))
}
Data$group = factor(Data$group, levels = ICIs_item)
Data = Data[order(Data$group), ]

p_aupr_1 = aupr_plot(Data[Data$group %in% ICIs_item[1:4], ], color_set = c('#d96098', '#325288', '#24a19c',  '#d3d1d1'))
p_aupr_2 = aupr_plot(Data[Data$group %in% ICIs_item[5:8], ], color_set = c('#d96098', '#325288', '#24a19c',  '#d3d1d1'))

Data = data.frame(group='', tpr=0, fpr=0, auc=0, auc_str='', cut_off=0, tpr_cut=0, fpr_cut=0)[-1,]
for (i in 1:8){
  Data = rbind(Data, data_generate_auc(pred = prediction(out_sub[out_sub$Drug == ICIs_item[i],]$pred, 
                                                         out_sub[out_sub$Drug == ICIs_item[i],]$out_de),
                                       group_name = ICIs_item[i]))
}
Data$group = factor(Data$group, levels = ICIs_item)
Data = Data[order(Data$group), ]

p_auc_1 = auc_plot(Data[Data$group %in% ICIs_item[1:4], ], color_set = c('#d96098', '#325288', '#24a19c',  '#d3d1d1'))
p_auc_2 = auc_plot(Data[Data$group %in% ICIs_item[5:8], ], color_set = c('#d96098', '#325288', '#24a19c',  '#d3d1d1'))

dev.off()
pdf('毕业_Model_auc_aupr.pdf', height =8, width = 8)
print(p_auc_1 + p_aupr_1 + p_auc_2 + p_aupr_2)
dev.off()

##############################################################################################################################

pdf('毕业_Model_CV_KM.pdf', height = 5, width = 6)
for (i in 1:8){
  ready = out_sub[out_sub$Drug == ICIs_item[i], c('Drug', 'out_de', 'duration', 'primaryid', 'pred')]
  ready$pred_2 = (ready$pred > 0.5) + 1
  
  ready_baseline = ICI_clinical[ICI_clinical$Drug == ICIs_item[i] &
                                  ICI_clinical$duration <= 90 &
                                  # ICI_clinical_Backup$outc_cod_num == 1 &
                                  !is.na(ICI_clinical$duration) & 
                                  !ICI_clinical$primaryid %in% ready$primaryid,]
  ready_baseline$out_de = ready_baseline$outc_cod_num
  ready_baseline$pred = 0
  ready_baseline$pred_2 = 0
  
  ready = rbind(ready, ready_baseline[, c('Drug', 'out_de', 'duration', 'primaryid', 'pred', 'pred_2')])
  
  fit <- survfit(Surv(duration, out_de) ~ pred_2, data = ready)
  print(survdiff(Surv(duration, out_de) ~ pred_2, data = ready)$pvalue)
  print(survdiff(Surv(duration, out_de) ~ pred_2, data = ready[ready$pred_2 !=2, ])$pvalue)
  print(survdiff(Surv(duration, out_de) ~ pred_2, data = ready[ready$pred_2 !=0, ])$pvalue)
  
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
                 # palette = c('#302f32', 'darkgrey', '#fe7800'),
                 palette = c('#302f32', '#a2b5bb', '#fe7800'),
                 xlim = c(0, 90),
                 break.time.by = 14,
                 pval.coord = c(0, 0.2),
                 pval.prefix = "Log-rank test"
  )
  
  a = p$plot + 
    scale_x_continuous(breaks = seq(0, 90, 14), labels = seq(0, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    ggtitle(ICIs_item[i])+
    theme(axis.text = element_text(color = 'black'),
          panel.grid = element_blank(),
          legend.position = c(0.8, 1),
          legend.justification = c(0, 1),
          legend.box.just = "left"
    )+
    labs(y = '90 days survival probability')
  
  b = p$table + 
    scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 14), labels = seq(0, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    theme(axis.text = element_text(color = 'black'),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none'
    )
  
  print(a /b + plot_layout(heights = c(5,2),guides='collect') & theme(legend.position='right'))
}
dev.off()


####################################

Data = data.frame(group='', prec=0, rec=0, aupr=0, aupr_str='', cut_off=0, prec_cut=0, rec_cut=0)[-1,]
Data = rbind(Data, data_generate_aupr(pred = prediction(holdout_sub[holdout_sub$Drug %in%  ICIs_item[1:4],]$pred, 
                                                        holdout_sub[holdout_sub$Drug %in%  ICIs_item[1:4],]$out_de),
                                      group_name = 'PD-1/PD-L1'))
Data = rbind(Data, data_generate_aupr(pred = prediction(holdout_sub[holdout_sub$Drug %in%  ICIs_item[5:8],]$pred, 
                                                        holdout_sub[holdout_sub$Drug %in%  ICIs_item[5:8],]$out_de),
                                      group_name = 'Combined'))

Data$group = factor(Data$group, levels = c('PD-1/PD-L1', 'Combined'))
Data = Data[order(Data$group), ]
p_aupr = aupr_plot(Data, color_set = c('#d8b74c', '#2d9086'))


Data = data.frame(group='', tpr=0, fpr=0, auc=0, auc_str='', cut_off=0, tpr_cut=0, fpr_cut=0)[-1,]
Data = rbind(Data, data_generate_auc(pred = prediction(holdout_sub[holdout_sub$Drug %in%  ICIs_item[1:4],]$pred, 
                                                       holdout_sub[holdout_sub$Drug %in%  ICIs_item[1:4],]$out_de),
                                     group_name = 'PD-1/PD-L1'))
Data = rbind(Data, data_generate_auc(pred = prediction(holdout_sub[holdout_sub$Drug %in%  ICIs_item[5:8],]$pred, 
                                                       holdout_sub[holdout_sub$Drug %in%  ICIs_item[5:8],]$out_de),
                                     group_name = 'Combined'))

Data$group = factor(Data$group, levels = c('PD-1/PD-L1', 'Combined'))
Data = Data[order(Data$group), ]
p_auc = auc_plot(Data, color_set = c('#d8b74c', '#2d9086'))

dev.off()
pdf('毕业_Model_auc_aupr_holdout.pdf', height =4, width = 8)
print(p_auc + p_aupr)
dev.off()

pdf('毕业_Model_CV_KM_holdoout.pdf', height = 5, width = 6)
for (i in c(1,5)){
  ready = holdout_sub[holdout_sub$Drug %in% ICIs_item[i:(i+3)], c('Drug', 'out_de', 'duration', 'primaryid', 'pred')]
  ready$pred_2 = (ready$pred > 0.5) + 1
  
  ready_baseline = ICI_holdout[ICI_holdout$Drug %in% ICIs_item[i:(i+3)]&
                                 ICI_holdout$duration <= 90 &
                                 # ICI_clinical_Backup$outc_cod_num == 1 &
                                 !is.na(ICI_holdout$duration) & 
                                 !ICI_holdout$primaryid %in% ready$primaryid,]
  ready_baseline$out_de = ready_baseline$outc_cod_num
  ready_baseline$pred = 0
  ready_baseline$pred_2 = 0
  
  ready = rbind(ready, ready_baseline[, c('Drug', 'out_de', 'duration', 'primaryid', 'pred', 'pred_2')])
  
  fit <- survfit(Surv(duration, out_de) ~ pred_2, data = ready)
  print(survdiff(Surv(duration, out_de) ~ pred_2, data = ready)$pvalue)
  print(survdiff(Surv(duration, out_de) ~ pred_2, data = ready[ready$pred_2 !=2, ])$pvalue)
  print(survdiff(Surv(duration, out_de) ~ pred_2, data = ready[ready$pred_2 !=0, ])$pvalue)
  
  
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
                 palette = c('#525e75', '#a2b5bb', '#ff6363'),
                 xlim = c(0, 90),
                 break.time.by = 14,
                 pval.coord = c(0, 0.2),
                 pval.prefix = "Log-rank test"
  )
  
  a = p$plot + 
    scale_x_continuous(breaks = seq(0, 90, 14), labels = seq(0, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    # ggtitle(ICIs_item[i])+
    theme(axis.text = element_text(color = 'black'),
          panel.grid = element_blank(),
          legend.position = c(0.8, 1),
          legend.justification = c(0, 1),
          legend.box.just = "left"
    )+
    labs(y = '90 days survival probability')
  
  b = p$table + 
    scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 14), labels = seq(0, 90, 14)/7, name = 'Time-to-Onset (weeks)')+
    theme(axis.text = element_text(color = 'black'),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none'
    )
  
  print(a /b + plot_layout(heights = c(5,2),guides='collect') & theme(legend.position='right'))
}
dev.off()
