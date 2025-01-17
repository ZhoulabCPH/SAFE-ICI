# BiocManager::install("ggVennDiagram")
# devtools::install_github("junjunlab/jjPlot")
# devtools::install_github("junjunlab/jjAnno")
# devtools::install_github('erocoar/gghalves')
# BiocManager::install("graph",force = TRUE)
# BiocManager::install("RBGL",force = TRUE) 
# install.packages("Vennerable", repos="http://R-Forge.R-project.org")

library('openxlsx')
library('ggplot2')
# library('ggpolypath')
library('gghalves')
library('ggrepel')
library("scales")
library('ggsci')
library('ggforce')
library('patchwork')
library('data.table')
library('reshape2')
# library('tseries') #adf
# library('fUnitRoots') #adf
library('vcdExtra')
library('ggplotify')
library('ggpubr')
library('tidyverse')
library('pheatmap')
library('cowplot')
library('UpSetR')
# library('venn')
library('Vennerable')
library('VennDiagram')
# library('ggvenn')
library('stringr')
library('circlize')
library('ComplexHeatmap')
library('Hmisc')
library('corrplot')
library('ROCR')
library('ggiraphExtra')
library('jjPlot')
library('jjAnno')

library('table1')

library('survival')
library('survminer')

library('meta')
library('xgboost')

path = stringi::stri_reverse(str_split_fixed(stringi::stri_reverse(setwd(dirname(rstudioapi::getSourceEditorContext()$path))), "/", 2)[2])

min_N = 10

# prec = f1*rec/(2*rec - f1)

f1_base = data.frame(rep(seq(0.001, 1, 0.001), 5), rep(seq(1,9,2), each=1000)/10)
colnames(f1_base) = c('recall', 'f1')
f1_base$precision = f1_base$f1*f1_base$recall/(2*f1_base$recall - f1_base$f1)
f1_base$f1 = factor(f1_base$f1, levels = c(unique(f1_base$f1), 'F1='))
f1_base = f1_base[f1_base$precision<=1 & f1_base$precision>=0, ]

max_list = f1_base %>% group_by(f1) %>% summarise(recall = min(recall))
f1_base = rbind(f1_base, data.frame(max_list, precision = 1))

f1_score_label = f1_base %>% group_by(f1) %>% summarise(y = min(precision))
f1_score_label = rbind(f1_score_label, c('F1=', 0.9))
f1_score_label$y = as.numeric(f1_score_label$y)


d1 = data.frame(SOC = c(
  # 'Nivolumab',
  # 'Pembrolizumab',
  # 'Durvalumab',
  # 'Atezolizumab',
  # 'Avelumab',
  # 'Ipilimumab', 
  "Respiratory, thoracic and mediastinal disorders",
  "Gastrointestinal disorders",
  "Metabolism and nutrition disorders",
  "Hepatobiliary disorders",
  "Immune system disorders",
  "Endocrine disorders",
  "Nervous system disorders",
  "Cardiac disorders",
  "Neoplasms benign, malignant and unspecified (incl cysts and polyps)",
  "General disorders and administration site conditions",
  "Infections and infestations",
  "Renal and urinary disorders",
  "Injury, poisoning and procedural complications",
  "Blood and lymphatic system disorders",
  "Vascular disorders",
  "Skin and subcutaneous tissue disorders",
  "Musculoskeletal and connective tissue disorders",
  "Eye disorders",
  "Psychiatric disorders",
  "Ear and labyrinth disorders",
  "Reproductive system and breast disorders",
  "Pregnancy, puerperium and perinatal conditions",
  "Congenital, familial and genetic disorders",
  "Product issues"),
  short_SOC = c('Respiratory',
                'Gastrointestinal',
                'Metabolism & nutrition',
                'Hepatobiliary',
                'Immune',
                'Endocrine',
                'Nervous',
                'Cardiac',
                'Neoplasms',
                'General',
                'Infections & infestations',
                'Renal & urinary',
                'Injury & poisoning',
                'Blood',
                'Vascular',
                'Skin',
                'Musculoskeletal',
                'Eye',
                'Psychiatric',
                'Ear & labyrinth',
                'Reproductive & breast',
                'Pregnancy conditions',
                "Congenital & genetic",
                "Product issues"),
  col1 = c(
    # '#40A4D8', '#33BEB7', '#5e79c9', 
    # '#B2C224', '#FECC2F', '#be7aca', 
    '#fa75b0', '#ff8185', '#ffa255', '#fecc2f', '#FBA127', '#F66320', 
    '#c54585', '#f6a39b', '#DB3937', '#ffc2b0', '#e94363', '#A463D7', 
    '#0C5BCE', '#a2c86c', '#a8aabc', '#1d8c83', '#5e9795', '#00accc', 
    '#051937', '#004d7a', '#008793', '#00bf72', '#a8eb12', '#b5d5ff'))

ATC_MedDRA = read.csv(paste(path, '/Data/ATC_MedDRA.csv', sep=''), sep=',')
ATC_MedDRA$Reac = toupper(ATC_MedDRA$descendant_name)

# PT_review
AE_review = read.csv(paste(path, '/Data/AE_review.csv', sep=''), sep=',')
AE_review$Reac = toupper(AE_review$PT_exist)
AE_review = merge(AE_review, ATC_MedDRA, by='Reac')

AE_review = AE_review[AE_review$ancestor_class_id == 'PT', 
                      c("descendant_concept_code", "descendant_name", "descendant_class_id", 
                        "ancestor_concept_code", "ancestor_name")]
colnames(AE_review) = c('MedDRA ID (pre)', 'irAE terms', 'irAE terms class','MedDRA ID (after)', 'PT terms')
AE_review_PT = unique(AE_review)
AE_review_PT$Reac = toupper(AE_review_PT$`PT terms`)
AE_review_PT = merge(AE_review_PT, ATC_MedDRA_embed[, 1:2], by='Reac')
write.csv(unique(AE_review_PT[,2:7]), paste(path, '/Data/AE_review_PT.csv', sep=''), row.names = FALSE)
AE_review_PT = read.csv(paste(path, '/Data/AE_review_PT.csv', sep=''), sep=',')

ATC_RxNorm = read.csv(paste(path, '/Data/ATC_RxNorm.csv', sep=''), sep=',')

ATC_MedDRA = ATC_MedDRA[ATC_MedDRA$ancestor_class_id=='SOC' &
                          ATC_MedDRA$descendant_class_id =='PT', 
                        c('ancestor_name', 'descendant_name')]
colnames(ATC_MedDRA) = c('SOC', 'Reac')
ATC_MedDRA$Reac = toupper(ATC_MedDRA$Reac)
ATC_MedDRA_embed = merge(ATC_MedDRA, d1, by='SOC')

# drug_group
ICIs_item = c('PD-1/PD-L1_NSCLC', 'PD-1/PD-L1_MEL', 'PD-1/PD-L1_RCC', 'PD-1/PD-L1_Other', 
              'Combined_NSCLC', 'Combined_MEL', 'Combined_RCC', 'Combined_Other',
              'CTLA-4_NSCLC', 'CTLA-4_MEL', 'CTLA-4_RCC', 'CTLA-4')

Other_L01 = ATC_RxNorm[substr(ATC_RxNorm$concept_code_2, 1, 3) == 'L01', ]$concept_name_1
non_anti_tumor = setdiff(ATC_RxNorm$concept_name_1, ATC_RxNorm[substr(ATC_RxNorm$concept_code_2, 1, 3) == 'L01', ]$concept_name_1)
gluco = setdiff(ATC_RxNorm$concept_name_1, ATC_RxNorm[substr(ATC_RxNorm$concept_code_2, 1, 5) == 'H02AB', ]$concept_name_1)


# source data
SourceData_info = function(fig_data, colnames, sheet, cols, widths){
  
  colnames(fig_data) = colnames
  header_stytle = createStyle(fontSize = 10, fontName = "Times New Roman", textDecoration = "bold", halign = "left")
  body_header_stytle = createStyle(fontSize = 10, fontName = "Times New Roman", textDecoration = "bold", halign = "center", border='TopBottom')
  body_stytle = createStyle(fontSize = 10, fontName = "Times New Roman", halign = "right")
  
  addWorksheet(SourceData, sheet)
  
  writeDataTable(SourceData, sheet, fig_data, startRow = 1, tableStyle = "none")
  
  setColWidths(SourceData, sheet, cols = cols, widths = widths)
  
  addStyle(SourceData, sheet, style = header_stytle, rows = 1, cols = 1:ncol(fig_data), gridExpand = T)
  addStyle(SourceData, sheet, style = body_header_stytle, rows = 1, cols = 1:ncol(fig_data), gridExpand = T)
  addStyle(SourceData, sheet, style = body_stytle, rows = 2:(nrow(fig_data)+1), cols = 1:ncol(fig_data), gridExpand = T)
  addStyle(SourceData, sheet, style = createStyle(fontSize = 10, fontName = "Times New Roman", halign = "right", border='bottom'), 
           rows = nrow(fig_data)+1, cols = 1:ncol(fig_data), gridExpand = T)
}

