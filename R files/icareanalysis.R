# baseline icare paper - exploration
pacman::p_load(here, tidyverse, skimr, gtsummary, flextable, officer,
       scales, rstatix, GGally, janitor, DataExplorer, sjmisc,
       vtable, psych, car, lmtest, naniar, 
       GGally, ggpubr, margins, cowplot,
       MASS, pscl, AER, performance, apaTables, Hmisc, glmmTMB, margins, dplyr, SmartEDA, summarytools,
       ggstatsplot, rcompanion, dunn.test)

# set current directory
i_am("R files/icareanalysis.R")

# load data
df <- readRDS(here::here("R files/clean_self_assessment.RData"))

# filter by time 1 only. 
df <- df %>% filter(assessment_time == "Time 1")

# Functions --------------------------------------------------------------------------------------------------
# add assist interpretations to score columns
add_assist_interpretation <- function(data) {
  data %>%
    mutate(Interpretation = case_when(
      Drug == "Alcohol" & Score < 10 ~ "Low Risk",
      Drug == "Alcohol" & Score >= 10 & Score <= 26 ~ "Moderate Risk",
      Drug == "Alcohol" & Score > 26 ~ "High Risk",
      Drug != "Alcohol" & Score < 4 ~ "Low Risk",
      Drug != "Alcohol" & Score <= 26 ~ "Moderate Risk",
      Drug != "Alcohol" & Score > 26 ~ "High Risk",
      TRUE ~ "Unknown"  # In case of NA or unexpected values
    ))
}

# rename assist strings 
replace_drug_values <- function(data, column) {
  data %>%
    mutate(!!column := case_when(
      grepl("alcohol", !!sym(column)) ~ "Alcohol",
      grepl("smoking", !!sym(column)) ~ "Smoking",
      grepl("cannabis", !!sym(column)) ~ "Cannabis",
      grepl("cocaine", !!sym(column)) ~ "Cocaine",
      grepl("opioids", !!sym(column)) ~ "Opioids",
      grepl("hallucinogens", !!sym(column)) ~ "Hallucinogens",
      grepl("inhalants", !!sym(column)) ~ "Inhalants",
      grepl("amphetamines", !!sym(column)) ~ "Amphetamines",
      grepl("sedatives", !!sym(column)) ~ "Sedatives",
      grepl("other", !!sym(column)) ~ "Other",
      TRUE ~ !!sym(column)  # Keep the original value if none of the patterns match
    ))
}

# calculate min, median, and max columns for summaries
min_median_max <- function(data, variable, ...) {
  x <- data[[variable]]
  tibble(
    Min = min(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Max = max(x, na.rm = TRUE)
  )
}

add_n_level <- function(data, variable, ...) {
  
  x <- data[[variable]] 
  
  x %>% group_by(x) %>% summarise(n = n())
  
  
}

my_flextable <- function(ft) {
  ft %>%
    theme_booktabs() %>%
    flextable::padding(padding = 1, part = "all") %>%
    flextable::line_spacing(space = 0.9, part = "all") %>% 
    flextable::fontsize(size = 10, part = "all") %>% 
    flextable::font(fontname = "Times New Roman", 
         part = "all") %>%
    flextable::set_table_properties(
      layout = "autofit",   # enables AutoFit to window
      align = "left",       # aligns table to left margin
      width = 1             # 100% of page width
    )
}

add_bullets <- function(results, items, style = "List Paragraph", level = 1) {
  for (item in items) {
    results <- body_add_par(results, item, style = style)
  }
  results
}
add_paragraphs <- function(results, lines, style = "Normal") {
  for (line in lines) {
    results <- body_add_par(results, line, style = style)
  }
  results
}

# Sample Descriptives -----------------------------------------------------------------------------------------------
# Purpose: Describe sample characteristics and prevalence/severity of AOD use and psychological distress.

## Summary stats ------
key <- df %>%
  dplyr::select(k10_11_absenteeism, k10_12_presenteeism, alcohol_use_interpretation, alcohol_use_scores, max_drug_use_interpretation, 
                max_other_drug_score, k10_total, k10_total_interpretation, gender, age, agency_years,
                recent_alcohol_user, recent_drug_user, recent_any_user, poly_drug_user)
  

#summarytools::view(summarytools::dfSummary(key)) 
options(sci.pen = .99)
vars <- key %>% names()

#view(descr(key))
# specify the tested variable (coerce to numeric if needed)
key %>%
  dplyr::select(k10_11_absenteeism, k10_12_presenteeism, alcohol_use_scores,
                max_other_drug_score, k10_total, age, agency_years) %>%
  tidyr::pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::group_by(var) %>%
  rstatix::shapiro_test(value)

# SmartEDA::ExpData(key, type = 1)  # Basic overview
# SmartEDA::ExpNumStat(key, by = "A")  # Group-wise num stats


# Shapiro-Wilk for all relevant vars
options(scipen = 999)
key %>%
  dplyr::select(age, agency_years, alcohol_use_scores, k10_11_absenteeism, 
         k10_12_presenteeism, k10_total, max_other_drug_score) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>% 
  group_by(var) %>%
  shapiro_test(value) %>%
  print(n = 100)



## Characteristics of Sample ----
# Characteristics of sample table
 df %>% 
  dplyr::select(gender, age, agency_years, k10_total, k10_total_interpretation, absenteeism = k10_11_absenteeism, 
         presenteeism = k10_12_presenteeism, rtc_1_problem:rtc_3_readiness_category, 
         recent_alcohol_user, alcohol_use_scores, alcohol_use_interpretation, recent_drug_user, recent_any_user,
         max_other_drug_score, max_drug_use_interpretation, poly_drug_user,
         recent_drug_count) %>% 
  mutate(
    recent_alcohol_user = factor(case_when(recent_alcohol_user == "0" ~ "No",
                                           recent_alcohol_user == "1" ~ "Yes",
                                           TRUE ~ recent_alcohol_user)),
    recent_any_user = factor(case_when(recent_any_user == "0" ~ "No",
                                           recent_any_user == "1" ~ "Yes",
                                           TRUE ~ recent_any_user)),
    recent_drug_user = factor(case_when(recent_drug_user == "0" ~ "No",
                                        recent_drug_user == "1" ~ "Yes",
                                        TRUE ~ recent_drug_user)),
    poly_drug_user = factor(case_when(poly_drug_user == "0" ~ "No",
                                        poly_drug_user == "1" ~ "Yes",
                                        TRUE ~ poly_drug_user))) %>% 
  mutate_if(is.factor, fct_explicit_na) %>% 
   # rtc_1_problem = fct_explicit_na(rtc_1_problem)
    
  #) %>% 
  tbl_summary(
    statistic = list(all_continuous() ~ c("{mean} ({sd})"),
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ c(1,1),
                  all_categorical() ~ c(0, 1)),
    type = list(recent_drug_count = "continuous",
             recent_alcohol_user = "dichotomous",
             recent_drug_user = "dichotomous",
             poly_drug_user = "dichotomous",
             recent_any_user = "dichotomous"
             ),
    label = list(gender = "Gender",
                 age = "Age",
                 k10_total = "K10 Score",
                 k10_total_interpretation = "K10 Category",
                 absenteeism = "Absenteeism",
                 presenteeism = "Presenteeism",
                 agency_years = "Agency years",
                 rtc_1_problem = "Perceived problem with alcohol/other drugs",
                 rtc_2_difficulty = "Perceived difficulty changing/cutting down",
                 rtc_3_readiness_category = "Readiness to change category",
                 rtc_3_readiness = "Readiness to change score",
                 recent_alcohol_user = "Alcohol use (past 3 months)",
                 alcohol_use_scores = "Alcohol use score",
                 alcohol_use_interpretation = "Alcohol use risk category",
                 recent_drug_user = "Other drug use (past 3 months)",
                 mean_other_drug_score = "Individual average other drug score",
                 max_other_drug_score = "Individual maximum other drug score",
                 max_drug_use_interpretation = "Maximum other drug use risk category", 
                 poly_drug_user = "Polydrug use",
                 recent_drug_count = "Number of drugs used in 3 months",
                 recent_any_user = "Recent alcohol or other drug use"
    )
  ) %>% 
  bold_labels() %>% 
  modify_header(all_stat_cols() ~ "Mean (SD) / n (%)") %>% 
  add_n()  %>%
  add_stat(
    fns = all_continuous() ~ min_median_max,
    location = all_continuous() ~ "label"
  ) %>%
  modify_header(list(
    Min ~ "Min",
    Median ~ "Median",
    Max ~ "Max"
  )) %>% 
  modify_fmt_fun(
    list(
      Min ~ number_format(accuracy = 0.1, trim = FALSE),
      Median ~ number_format(accuracy = 0.1, trim = FALSE),
      Max ~ number_format(accuracy = 0.1, trim = FALSE)
    )
  ) %>% 
  bold_labels() %>% 
  modify_header(list(
    label ~ "**Characteristic**",
    stat_0 ~ "**Mean (SD) / n (%)**",
    Min ~ "**Min**",
    Median ~ "**Median**",
    Max ~ "**Max**"
  )) %>% 
  modify_footnote(everything() ~ NA) %>% 
  modify_table_styling(column = label, 
                       align = "left") %>% 
  gtsummary::as_flex_table() %>% 
  #align(j = 1, align = "left",part = "all") %>% 
  set_caption(caption ="Table X. Baseline characteristics of sample (n = 238)",
              align_with_table = FALSE) %>% 
  my_flextable() %>% 
  add_footer_lines('Note. Gender categories "Other" and "Prefer not to answer" were recoded as missing data in further analyses involving gender due to insufficient sample sizes') %>% 
  flextable::font(fontname = "Times New Roman", 
                  part = "all")

df %>% dplyr::select(rtc_1_problem
                     ) %>% pivot_longer(everything()) %>% 
  group_by(value) %>% 
  summarise(n = n()) %>% 
  mutate(p = n/sum(n)*100)
            

df <- df %>% mutate(gender = case_when(gender == "Prefer not to answer" ~ NA,
                                       gender == "Other" ~ NA,
                                       TRUE ~ gender),
                    gender = factor(gender))
  

# Substance use prevalance ----

(df %>% dplyr::select(starts_with("assist_total_"), -assist_total_scores_count) %>% 
  pivot_longer(
    cols = starts_with("assist_total_"),
    names_to = "substance",
    values_to = "score"
  ) %>%
  mutate(score = case_when(score < 1 ~ NA, TRUE ~ score)) %>% 
  filter(!is.na(score)) %>% 
  mutate(
    substance = str_remove(substance, "assist_total_") %>% str_to_title()
  ) %>% 
  mutate(
    `Risk` = case_when(
      substance == "Alcohol" & score < 10 ~ "Low risk",
      substance == "Alcohol" & score >= 10 & substance == "Alcohol" & score <= 26 ~ "Moderate risk",
      substance == "Alcohol" & score >26 ~ "High risk",
      substance != "Alcohol" & score < 4 ~ "Low risk",
      substance != "Alcohol" & score >= 4 & substance != "Alcohol" & score <= 26 ~ "Moderate risk",
      substance != "Alcohol" & score > 26 ~ "High risk",
      TRUE ~ "fix")) %>% group_by(substance, Risk) %>% 
  group_by(substance, Risk) %>% 
  summarise(
    n = n()) %>% 
  ungroup() %>% 
  group_by(substance) %>% 
  mutate(total = sum(n),
         percent = round(100*n/sum(n),1)) %>% 
  group_by(substance, Risk) %>% ungroup() %>% 
  dplyr::select(substance, total, Risk, n, percent) %>% 
  pivot_wider(names_from = Risk, values_from = c(n, percent) ) %>% 
  dplyr::select(Substance = substance, total, ends_with("Low risk"), ends_with("Moderate risk"), ends_with("High risk")) %>% 
  mutate(total = paste0(total, " (", number(total/220*100, accuracy = .1),"%)")) %>% 
  filter(Substance != 'Other') %>% 
  filter(Substance != "Smoking") %>% flextable() -> rates)
  
  
(typology <- data.frame(
  col_keys = c("Substance","total","n_Low risk","percent_Low risk","n_Moderate risk",  
  "percent_Moderate risk", "n_High risk","percent_High risk"),
  head = c("", "", rep("ASSIST scores", 6)),
  risk = c("", "", rep("Low risk", 2), rep("Moderate risk", 2), rep("High risk", 2)),
  names = c("Substance", "Total using", rep(c("n", "%"), 3))
)
)
(rates <- rates %>% set_header_df(mapping = typology, key = "col_keys")  %>% 
  merge_h(part = "header") %>% 
  set_caption("Table X. Rates of substance use in the past three months by risk category (n = 220)",
              align_with_table = FALSE) %>% 
  fontsize(size = 10, part = "all") %>%
  align(j = 1, align = "left", part = "all") %>% 
  align(j = -1, align = "center", part = "all") %>%
  bold(i = NULL, j = NULL, bold = TRUE, part = "header") %>% 
  padding(padding = 4, part = "all") %>% 
  bold(j = 1, part = "body", bold = TRUE) %>% 
  border_remove() %>% 
  
  #hline(border = fp_border(color = "black", width = 0.75), part = "body") %>%
  #vline(border = fp_border(color = "black", width = 0.5), part = "all") %>% 
  set_table_properties(
    layout = "autofit",   # enables AutoFit to window
    align = "left",       # aligns table to left margin
    width = 1             # 100% of page width
  ) %>% 
  theme_booktabs() %>% 
  hline_top(border = fp_border(color = "black", width = 1), part = "all") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "all") %>%
  hline(border = fp_border(color = "black", width = 1), i = c(1,2), part = "header") %>% 
  bold(part = "header", bold = TRUE) %>% 
  flextable::font(fontname = "Times New Roman", 
                  part = "all") %>% 
  align(align = "center", part = "all", j = -1) %>% 
    flextable::padding(padding = 1, part = "all") %>%
    flextable::line_spacing(space = 0.9, part = "all") %>% 
    flextable::fontsize(size = 10, part = "all") %>% 
    flextable::font(fontname = "Times New Roman", 
                    part = "all") %>%
    flextable::set_table_properties(
      layout = "autofit",   # enables AutoFit to window
      align = "left",       # aligns table to left margin
      width = 1             # 100% of page width
    ) %>% 
    align(j = -1, align = "center", part = "all") %>% 
    bold(bold = TRUE, part = "header"))


# td <- tempdir()
# tf <- tempfile("tables", td, ".docx")
# save_as_docx(rates, path = tf)  
# browseURL(tf)


# Rates and severity -----------------------------------------------------------------------------------------

## Alcohol use characteristics ----
  
# Summary based on alcohol risk scores
(alc <- df %>% 
    filter(recent_alcohol_user == 1) %>% 
    dplyr::select(gender, age, agency_years, 
                  alcohol_use_interpretation, 
                  k10_total, 
                  absenteeism = k10_11_absenteeism, 
                  presenteeism = k10_12_presenteeism
                  
                  ) %>% 
    tbl_summary(by = alcohol_use_interpretation,
                missing = "no",
                missing_text = "Missing",
                statistic = list(all_continuous() ~ c("{mean} ({sd})"),
                                 all_categorical() ~ "{n} ({p}%)"),
                percent = "row",
                digits = list(all_continuous() ~ c(1,1),
                              all_categorical() ~ c(0, 1)),
                type = list(recent_drug_count = "continuous",
                         #recent_alcohol_user = "categorical",
                         recent_drug_user = "categorical",
                         poly_drug_user = "categorical"),
                label = list(gender = "Gender",
                             age = "Age",
                             k10_total = "Psychological distress",
                             k10_total_interpretation = "K10 Category",
                             absenteeism = "Absenteeism",
                             presenteeism = "Presenteeism",
                             agency_years = "Agency years",
                             rtc_1_problem = "Perceived problem with alcohol/other drugs",
                             rtc_2_difficulty = "Perceived difficulty changing/cutting down",
                             rtc_3_readiness_category = "Readiness to change category",
                             rtc_3_readiness = "Readiness to change score",
                             recent_alcohol_user = "Alcohol use (past 3 months)",
                             alcohol_use_scores = "Alcohol use score",
                             alcohol_use_interpretation = "Alcohol use risk category",
                             recent_drug_user = "Other drug use (past 3 months)",
                             mean_other_drug_score = "Individual average other drug score",
                             max_other_drug_score = "Maximum other drug score",
                             max_drug_use_interpretation = "Maximum other drug use risk category", 
                             poly_drug_user = "Polydrug use",
                             recent_drug_count = "Number of drugs used in 3 months"
                )) %>% 
    bold_labels() %>% 
    add_n(statistic = "{N_nonmiss}",
          col_label = "n") %>% 
    modify_header(list(
      label ~ "Characteristic",
      stat_1 ~ "Low Risk\n(n = 77)",
      stat_2 ~ "Moderate Risk\n(n = 83)",
      stat_3 ~ "High Risk\n(n = 56)"
      )) %>% 
  add_p(test = list(all_continuous() ~ "kruskal.test",
      gender ~ "chisq.test"),
      
    
    #test.args = all_categorical() ~ list(simulate.p.value = TRUE),
    pvalue_fun = function(x) {
      formatted <- ifelse(x < 0.001, "<.001", sprintf("%.3f", x))
      sub("^0\\.", ".", formatted)
    }
    
    ) %>% 
    bold_p(t = .05) %>% 

    add_stat_label(label = list(all_continuous() ~ "Mean (SD)",
                                all_categorical() ~ "n (%)")) %>% 
    
    modify_footnote(everything() ~ NA) %>% 
  modify_header(p.value ~ "p") %>% gtsummary::as_flex_table() %>% 
  flextable::set_table_properties(
      layout = "autofit",   # enables AutoFit to window
      align = "left",       # aligns table to left margin
      width = 1             # 100% of page width
    ) %>% 
    align(j = 1, align = "left", part = "all") %>% 
    align(j = -1, align = "center", part = "all") %>% 
    set_caption(caption = "Table X. Sample characteristics by alcohol use risk levels (n = 216)", 
                align_with_table = "FALSE") %>% 
    #add_footer_lines(paste(c("Note. pvalues are from Kruskal-Wallis tests for continuous variables, and Chi-squared for categorical variables. Bolded values indicate statistically significant differences at p < .05"))) %>% 
    align(i = 1, align = "left", part = "footer") %>% 
    flextable::padding(padding = 1, part = "all") %>%
    flextable::line_spacing(space = 0.9, part = "all") %>% 
    flextable::fontsize(size = 10, part = "all") %>% 
    flextable::font(fontname = "Times New Roman", 
                    part = "all") %>% 
  bold(bold = T, i = 1, part = "header"))


(tests <- df %>% filter(recent_alcohol_user == 1))
(dunn.test(tests$age, tests$alcohol_use_interpretation,kw = TRUE, method = "bonferroni"))
(dunn.test(tests$agency_years, tests$alcohol_use_interpretation,kw = TRUE, method = "bonferroni"))
 (dunn.test(tests$k10_total, tests$alcohol_use_interpretation,kw = TRUE, method = "bonferroni"))
(dunn.test(tests$k10_11_absenteeism, tests$alcohol_use_interpretation,kw = TRUE, method = "bonferroni"))
 (dunn.test(tests$k10_12_presenteeism, tests$alcohol_use_interpretation,kw = TRUE, method = "bonferroni"))

 # overall test 
  report::report(stats::chisq.test(tests$gender, tests$alcohol_use_interpretation))
  
  gender <- table(tests$gender, tests$alcohol_use_interpretation)
  pairwiseNominalIndependence(
    gender,
    fisher = FALSE,
    gtest  = FALSE,
    chisq  = TRUE,
    method = "bonferroni"  # spelling: no capital B
  )
  
  prop.table(table(tests$gender, tests$alcohol_use_interpretation), margin = 1)

# # save table to word document
# save_as_docx(alc, path = tf)
# browseURL(tf)
## Other drug use characteristics ----

(drug <- df %>% 
  filter(recent_drug_user == 1) %>% 
  dplyr::select(gender, age, agency_years, 
                max_drug_use_interpretation, 
                k10_total, 
                absenteeism = k10_11_absenteeism, 
                presenteeism = k10_12_presenteeism
                
  ) %>% 
  tbl_summary(by = max_drug_use_interpretation,
              missing = "no",
              missing_text = "Missing",
              statistic = list(all_continuous() ~ c("{mean} ({sd})"),
                               all_categorical() ~ "{n} ({p}%)"),
              percent = "row",
              digits = list(all_continuous() ~ c(1,1),
                            all_categorical() ~ c(0, 1)),
              type = list(recent_drug_count = "continuous",
                       #recent_alcohol_user = "categorical",
                       recent_drug_user = "categorical",
                       poly_drug_user = "categorical"),
              label = list(gender = "Gender",
                           age = "Age",
                           k10_total = "Psychological distress",
                           k10_total_interpretation = "K10 Category",
                           absenteeism = "Absenteeism",
                           presenteeism = "Presenteeism",
                           agency_years = "Agency years",
                           rtc_1_problem = "Perceived problem with alcohol/other drugs",
                           rtc_2_difficulty = "Perceived difficulty changing/cutting down",
                           rtc_3_readiness_category = "Readiness to change category",
                           rtc_3_readiness = "Readiness to change score",
                           recent_alcohol_user = "Alcohol use (past 3 months)",
                           alcohol_use_scores = "Alcohol use score",
                           alcohol_use_interpretation = "Alcohol use risk category",
                           recent_drug_user = "Other drug use (past 3 months)",
                           mean_other_drug_score = "Individual average other drug score",
                           max_other_drug_score = "Maximum other drug score",
                           max_drug_use_interpretation = "Maximum other drug use risk category", 
                           poly_drug_user = "Polydrug use",
                           recent_drug_count = "Number of drugs used in 3 months"
              )) %>% 
  bold_labels() %>% 
  add_n(statistic = "{N_nonmiss}",
        col_label = "n") %>% 
  modify_header(list(
    label ~ "Characteristic",
    stat_1 ~ "Low Risk\n(n = 17)",
    stat_2 ~ "Moderate Risk\n(n = 46)",
    stat_3 ~ "High Risk\n(n = 5)"
  )) %>% 
  add_p(test = list(all_continuous() ~ "kruskal.test",
                    gender ~ "fisher.test"),
        
        
        #test.args = all_categorical() ~ list(simulate.p.value = TRUE),
        pvalue_fun = function(x) {
          formatted <- ifelse(x < 0.001, "<.001", sprintf("%.3f", x))
          sub("^0\\.", ".", formatted)
        }
        
  ) %>% 
  bold_p(t = .05) %>% 
  
  add_stat_label(label = list(all_continuous() ~ "Mean (SD)",
                              all_categorical() ~ "n (%)")) %>% 
  
  modify_footnote(everything() ~ NA) %>% 
  modify_header(p.value ~ "p") %>% gtsummary::as_flex_table() %>% 
  flextable::set_table_properties(
    layout = "autofit",   # enables AutoFit to window
    align = "left",       # aligns table to left margin
    width = 1             # 100% of page width
  ) %>% 
  align(j = 1, align = "left", part = "all") %>% 
  align(j = -1, align = "center", part = "all") %>% 
  set_caption(caption = "Table X. Sample characteristics by other drug use risk levels (n = 68)", 
              align_with_table = "FALSE") %>% 
  add_footer_lines(paste(c("Note. pvalues are from Kruskal-Wallis tests for continuous variables, and Chi-squared for categorical variables. Bolded values indicate statistically significant differences at p < .05"))) %>% 
  align(i = 1, align = "left", part = "footer") %>% 
  flextable::padding(padding = 1, part = "all") %>%
  flextable::line_spacing(space = 0.9, part = "all") %>% 
  flextable::fontsize(size = 10, part = "all") %>% 
  flextable::font(fontname = "Times New Roman", 
                  part = "all") %>% 
  bold(bold = T, i = 1, part = "header"))


tests <- df %>% filter(recent_drug_user == 1)
(fisher.test(table(tests$gender, tests$max_drug_use_interpretation)))
(dunn.test(tests$age, tests$max_drug_use_interpretation,kw = TRUE, method = "bonferroni"))
(dunn.test(tests$agency_years, tests$max_drug_use_interpretation,kw = TRUE, method = "bonferroni"))
(dunn.test(tests$k10_total, tests$max_drug_use_interpretation,kw = TRUE, method = "bonferroni"))
(dunn.test(tests$k10_11_absenteeism, tests$max_drug_use_interpretation,kw = TRUE, method = "bonferroni"))
dunn.test(tests$k10_12_presenteeism, tests$max_drug_use_interpretation,kw = TRUE, method = "bonferroni")


# save_as_docx(drug, path = tf)
# browseURL(tf)






# Associations -----------------------------------------------------------------------------------------------

# Psychological distress ----
df %>% filter(!(is.na(k10_total))) %>% 
  dplyr::select(gender, age, agency_years, k10_total_interpretation, 
                alcohol_use_scores, max_other_drug_score,
                absenteeism = k10_11_absenteeism, 
                presenteeism = k10_12_presenteeism ) %>% 
  
  tbl_summary(by = k10_total_interpretation,
              missing = "no",
              missing_text = "Missing",
              statistic = list(all_continuous() ~ c("{mean} ({sd})"),
                               all_categorical() ~ "{n} ({p}%)"),
              percent = "row",
              digits = list(all_continuous() ~ c(1,1),
                            all_categorical() ~ c(0, 1)),
              type = list(recent_drug_count = "continuous",
                       #recent_alcohol_user = "categorical",
                       recent_drug_user = "categorical",
                       poly_drug_user = "categorical"),
              label = list(gender = "Gender",
                           age = "Age",
                           k10_total = "Psychological distress",
                           k10_total_interpretation = "K10 Category",
                           absenteeism = "Absenteeism",
                           presenteeism = "Presenteeism",
                           agency_years = "Agency years",
                           rtc_1_problem = "Perceived problem with alcohol/other drugs",
                           rtc_2_difficulty = "Perceived difficulty changing/cutting down",
                           rtc_3_readiness_category = "Readiness to change category",
                           rtc_3_readiness = "Readiness to change score",
                           recent_alcohol_user = "Alcohol use (past 3 months)",
                           alcohol_use_scores = "Alcohol use score",
                           alcohol_use_interpretation = "Alcohol use risk category",
                           recent_drug_user = "Other drug use (past 3 months)",
                           mean_other_drug_score = "Individual average other drug score",
                           max_other_drug_score = "Maximum other drug score",
                           max_drug_use_interpretation = "Maximum other drug use risk category", 
                           poly_drug_user = "Polydrug use",
                           recent_drug_count = "Number of drugs used in 3 months"
              )) %>% 
  bold_labels() %>% 
  add_n(statistic = "{N_nonmiss}",
        col_label = "n") %>% 
 
  add_p(test = list(all_continuous() ~ "kruskal.test",
                    gender ~ "chisq.test"),
        
        
        #test.args = all_categorical() ~ list(simulate.p.value = TRUE),
        pvalue_fun = function(x) {
          formatted <- ifelse(x < 0.001, "<.001", sprintf("%.3f", x))
          sub("^0\\.", ".", formatted)
        }
        
  ) %>% 
  bold_p(t = .05) %>% 
  
  add_stat_label(label = list(all_continuous() ~ "Mean (SD)",
                              all_categorical() ~ "n (%)")) %>% 
  
  modify_footnote(everything() ~ NA)


  
tests <- df %>% filter(!is.na(k10_total))
dunn.test(tests$alcohol_use_scores, tests$k10_total_interpretation,kw = TRUE, method = "bonferroni")
dunn.test(tests$k10_11_absenteeism, tests$k10_total_interpretation,kw = TRUE, method = "bonferroni")
dunn.test(tests$k10_12_presenteeism, tests$k10_total_interpretation,kw = TRUE, method = "bonferroni")
report::report(stats::chisq.test(tests$gender, tests$alcohol_use_interpretation))

gender <- table(tests$gender, tests$alcohol_use_interpretation)
pairwiseNominalIndependence(
  gender,
  fisher = FALSE,
  gtest  = FALSE,
  chisq  = TRUE,
  method = "bonferroni"  # spelling: no capital B
)


tests <- df %>% filter(!is.na(k10_total))

dunn.test(tests$k10_11_absenteeism, tests$k10_total_interpretation,kw = TRUE, method = "bonferroni")
dunn.test(tests$k10_12_presenteeism, tests$k10_total_interpretation,kw = TRUE, method = "bonferroni")



## Correlations ----
# Pearson correlations (for normally distributed variables)
# spearman for those involving absenteeism and presenteeism

cor_vars <- df %>%
  dplyr::select(
    alcohol_use_scores, max_other_drug_score, 
    k10_total, absenteeism = k10_11_absenteeism, presenteeism = k10_12_presenteeism
  ) 
# define variable types 

continuous <- c("k10_total", "alcohol_use_scores", "max_other_drug_score")
count <- c('presenteeism', 'absenteeism', 'recent_drug_count')



glimpse(cor_vars)

# get all correlation pairs 
cor_pairs <- combn(names(cor_vars), 2, simplify = FALSE)
cor_pairs

# create function to get correct correlation type 
get_cors <- function(x_name, y_name, data) {
  x <- data[[x_name]]
  y <- data[[y_name]]
  
  # define variable types 
  is_x_count <- x_name %in% count
  is_y_count <- y_name %in% count
  is_x_continuous <- x_name %in% continuous
  is_y_continuous <- y_name %in% continuous
  
  # define correlation method 
  method <- case_when(
    is_x_count & is_y_count ~ "spearman", # use spearman
    is_x_count & is_y_continuous ~ "spearman", # use spearman
    is_x_continuous & is_y_count ~ "spearman", # use spearman
    is_x_continuous & is_y_continuous ~ "pearson",  # use pearson 
    TRUE ~ "spearman"  # if else, use spearman 
  )
  
  cor_type <- case_when(
    is_x_continuous & is_y_continuous ~ "pearson",
    TRUE ~ "spearman"
  )
  
  # calculate correlations 
  test <- suppressWarnings(cor.test(x, y, method = method))
  
  # return correlation coefficient, p-value, method used
  list(r = unname(test$estimate), p = test$p.value, method = method, type = cor_type)
}
  
  
  
# compute the correlations for all pairs
cor_results <- purrr::map_dfr(cor_pairs, function(vars) {
  res <- get_cors(vars[1], vars[2], cor_vars)
  tibble(
    var1 = vars[1],
    var2 = vars[2],
    cor = res$r,
    p = number(res$p, accuracy = .001),
    method = res$method,
    type = res$type#,
    #p2 = res$p
  )
})

print(cor_results, n = 80)
# Make symmetric by mirroring rows
cor_long <- cor_results %>%
  bind_rows(
    cor_results %>%
      rename(var1 = var2, var2 = var1)
  ) %>%
  mutate(
    cor = round(cor, 2),
    sig = p <= 0.05,
    stars = case_when(
      p < .001 ~ "***",
      p < .01 ~ "**", 
      p < .05 ~ "*",
      TRUE ~ ""),
    cor_ch = paste0(formatC(cor, format = "f", digits = 2), stars))
  
cor_long
# make wide for table
cor_matrix <- cor_long %>%
  dplyr::select(var1, var2, cor_ch) %>%
  pivot_wider(names_from = var2, values_from = cor_ch) %>% 
  dplyr::select(var1, alcohol_use_scores, max_other_drug_score, k10_total, absenteeism, presenteeism)
cor_matrix


# cor labels 
cor_labels <- c(
  k10_total = "Psychological Distress",
  presenteeism = "Presenteeism",
  absenteeism = "Absenteeism",
  alcohol_use_scores = "Alcohol Use Score",
  max_other_drug_score = "Max Other Drug Score")

cor_matrix

colnames(cor_matrix) <- c(" " = "var1", cor_labels[colnames(cor_matrix)[-1]])
cor_matrix[[1]] <- cor_labels[cor_matrix[[1]]]
cor_matrix[upper.tri(cor_matrix)] <- NA
cor_matrix 


# create ft 
cor_tbl <- flextable(cor_matrix)

cor_tbl %>% 
  fontsize(size = 10, part = "all") %>%
  align(j = 1, align = "left", part = "all") %>% 
  align(j = -1, align = "center", part = "all") %>%
  bold(i = NULL, j = NULL, bold = TRUE, part = "header") %>% 
  padding(padding = 4, part = "all") %>% 
  bold(j = 1, part = "body", bold = TRUE) %>% 
  border_remove() %>% 
  
  #hline(border = fp_border(color = "black", width = 0.75), part = "body") %>%
  #vline(border = fp_border(color = "black", width = 0.5), part = "all") %>% 
  theme_booktabs() %>% 
  hline_top(border = fp_border(color = "black", width = 1), part = "all") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "all") %>%
  #hline(border = fp_border(color = "black", width = 1), i = c(1,2), part = "header") %>% 
  bold(part = "header", bold = TRUE) %>% 
  flextable::font(fontname = "Times New Roman", 
                  part = "all") %>% 
  set_header_labels("var1" = " ") %>% 
  set_caption(caption = "Table X. Correlation matrix of key variables", align_with_table = "FALSE") %>% 
  add_footer_lines("Note. Bolded values indicate statistically significant correlations at * p < .05, ** p < .01, *** p < .001. Pearson correlations are reported for continuous-continuous pairs, and Spearman correlations were used for pairs involving count variables (absenteeism and presenteeism).") %>% 
  align(i = 1, align = "left", part = "footer") %>% 
  align(align = "center", part = "all", j = -1) %>% 
  flextable::padding(padding = 1, part = "all") %>%
  flextable::line_spacing(space = 0.9, part = "all") %>% 
  flextable::fontsize(size = 10, part = "all") %>% 
  flextable::font(fontname = "Times New Roman", 
                  part = "all") %>%
  flextable::set_table_properties(
    layout = "autofit",   # enables AutoFit to window
    align = "left",       # aligns table to left margin
    width = 1             # 100% of page width
  ) %>% 
  align(j = -1, align = "center", part = "all") %>% 
  bold(bold = TRUE, part = "header") -> cor_tbl


cor_tbl

 # tf <- tempfile("ex", here(), ".docx")
 # save_as_docx(cor_tbl, path = tf)
 # browseURL(tf)




# Workplace outcomes -----
glimpse(df)
df <- df %>% mutate(
  k10_11_absenteeism = round(k10_11_absenteeism, 1),
  k10_12_presenteeism = round(k10_12_presenteeism, 1)
) %>% rename(absenteeism = k10_11_absenteeism, 
             presenteeism = k10_12_presenteeism)

# Distribution plots

p1 <- df %>% dplyr::select(absenteeism, presenteeism) %>% 
  rename(Absenteeism = absenteeism, Presenteeism = presenteeism) %>% 
  pivot_longer(cols = c(Absenteeism, Presenteeism),
               names_to = "Variable",
               values_to = "Days") %>% 
  filter(!is.na(Days)) %>% 
gghistogram(
            x = "Days",
            y = "count",
            xlab = "Days",
            ylab = FALSE,
            #combine = TRUE,
            add = "median",
            facet.by = "Variable",
            #ggtheme = theme_pubr(),
            binwidth = 1
            
            ) + scale_x_continuous(breaks = seq(0, 30, by = 2))


p1

# very high zero count

## Check model type ----
df_sub <- df %>% dplyr::select(age, gender, agency_years, k10_total, rtc_3_readiness, 
                               recent_drug_user, recent_alcohol_user, recent_any_user, poly_drug_user,
                               alcohol_use_scores, max_other_drug_score, absenteeism, presenteeism)
                               

outcome <- "absenteeism"
predictors <- df_sub %>% dplyr::select(-absenteeism, -presenteeism) %>% names()

# run function to evaluate which count model handles outcome data best 
# poisson vs negative binomial vs zero inflated poisson vs zero inflated negative binomial, vs hurdle 
# choose best model type 

# function to run series of univariate models and check dispersion. 
compare_models <- function(var, data) {
  fmla <- as.formula(paste(outcome, "~", var))
  
  # print results:
  cat("\n============================\n")
  cat("Variable:", var, "\n")
  cat("============================\n\n")
  
  # fit models 
  # poisson
  pois <- glm(fmla, data = data, family = poisson())
  # negative binomial
  nb <- glm.nb(fmla, data = data)
  # zero inflated poisson
  zip <- zeroinfl(fmla, data = data, dist = "poisson")
  # zero inflated negative binomial
  zinb <- zeroinfl(fmla, data = data, dist = "negbin")
  # hurdle negative binomial 
  hurdle_nb <- hurdle(fmla, data = data, dist = "negbin", zero.dist = "binomial")
  
  # dispersion:
  disp_stat <- dispersiontest(pois)$statistic 
  cat('Dispersion (Poisson):', round(disp_stat, 2), '\n')
  
  # compare AIC / BIC
  aic_df <- AIC(pois, nb, zip, zinb, hurdle_nb)
  bic_df <- BIC(pois, nb, zip, zinb, hurdle_nb)
  
  cat("\nAIC Comparison:\n")
  print(aic_df)
  
  cat("\nBIC Comparison:\n")
  print(bic_df)
  
  # vuong tests
  cat("\nAIC Comparison:\n")
  print(aic_df)
  
  cat("\nBIC Comparison:\n")
  print(bic_df)
  
  # Vuong tests
  cat("\nVuong Test: ZIP vs Poisson\n")
  print(vuong(pois, zip))
  
  cat("\nVuong Test: ZINB vs NB\n")
  print(vuong(nb, zinb))
  
  cat("\nVuong Test: Hurdle vs NB\n")
  print(vuong(nb, hurdle_nb))


# get predicted vs observed zeros 
obs_zeros <- sum(data[[outcome]] == 0)
pred_zeros <- list(
  Poisson = sum(predict(pois, type = "response") < 1e-6),
  NB = sum(predict(nb, type = "response") < 1e-6),
  ZIP = sum(predict(zip, type = "response") < 1e-6),
  ZINB = sum(predict(zinb, type = "response") < 1e-6),
  Hurdle = sum(predict(hurdle_nb, type = "response") < 1e-6)
)

cat("\nObserved Zeros:", obs_zeros, "\n")
cat("Predicted Zeros:\n")
print(pred_zeros)

cat("\n--- Summary ---\n")
cat("- Use NB if overdispersion is present and excess zeros are not excessive.\n")
cat("- Use ZINB if you believe some people are 'always-zero'.\n")
cat("- Use Hurdle NB if everyone could have non-zero counts but a barrier (hurdle) prevents some.\n")
cat("- Use model with lowest AIC/BIC **and** conceptually appropriate assumptions.\n")
}

sink(here("R files/outputs/model_comparisons_absenteeism.txt"))
lapply(predictors, compare_models, data = df_sub %>% filter(!is.na(absenteeism)))
sink()

outcome <- "presenteeism"
predictors <- df_sub %>% dplyr::select(-presenteeism, -absenteeism) %>% names()

sink(here("R files/outputs/model_comparisons_presenteeism.txt"))
lapply(predictors, compare_models, data = df_sub %>% filter(!is.na(presenteeism)))
sink()

## Absenteeism ----


### univariate model ----
df_sub <- df %>% dplyr::select(age, gender, agency_years, 
                               alcohol_use_scores, max_other_drug_score, k10_total, absenteeism, presenteeism)
df_sub <- df_sub %>% mutate(
  alcohol_use_scores = case_when(is.na(alcohol_use_scores) ~ 0, TRUE ~ alcohol_use_scores),
  max_other_drug_score = case_when(is.na(max_other_drug_score) ~ 0, TRUE ~ max_other_drug_score)
)

#vtable(df_sub)
outcome <- "absenteeism"
predictors <- df_sub %>% dplyr::select(-absenteeism, -presenteeism) %>% names()
var_labels <- c(
  agency_years = "Years in agency",
  age = "Age",
  gender = "Gender",
  alcohol_use_scores = "Alcohol use risk",
  max_other_drug_score = "Maximum other drug use risk",
  k10_total = "Psychological distress",
  `max_other_drug_score:k10_total` = "Other drug use risk x Psychological distress"
)

run_uv_hurdle <- function(outcome, predictors, 
                          data, exponentiate = T) {
  
  # build both models seperately 
  model_tables <- map(predictors, ~ {
  fmla <- as.formula(paste0(outcome, " ~ ", .x, " | ", .x))
  model <- hurdle(fmla, data = data, dist = "negbin", zero.dist = "binomial")
    
  # zero part 
  zero_tbl <- tbl_regression(
    model,
    component = "zero_inflated",
    exponentiate = TRUE,
    estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
    pvalue_fun = function(x) {
      p <- style_pvalue(x, digits = 3)
      # Remove leading zero from both "0.03" and "<0.001"
      gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
    }, 
    tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    #modify_caption(paste("Zero-inflation model (logit) for", .x)) %>% 
    modify_header(label = "**Variable**", estimate = "**OR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )
    ) %>% 
    add_n() %>% 
    bold_p()
  
  # NB count part
  count_tbl <- tbl_regression(
    model,
    component = "conditional",
    exponentiate = T,
    estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
    pvalue_fun = function(x) {
      p <- style_pvalue(x, digits = 3)
      # Remove leading zero from both "0.03" and "<0.001"
      gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
    },
    
    tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    # modify_caption(paste("Count model (NB) for", .x)) %>% 
    modify_header(label = "**Variable**", estimate = "**IRR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )) %>% 

    bold_p()
  
  list(zero_tbl = zero_tbl, count_tbl = count_tbl)
  })

  # Separate and flatten count and zero tables
  zero_tbls  <- map(model_tables, "zero_tbl")
  count_tbls <- map(model_tables, "count_tbl")
  
  tbl_merge(
    tbls = list(
      tbl_stack(zero_tbls),
      tbl_stack(count_tbls) 
      ),
    tab_spanner = c("**Hurdle model (logit)**", "**Count model (NB)**")
  )
  
  
}

abs_hurdle_tbs <- run_uv_hurdle(outcome, predictors, data = df_sub,# %>% filter(!is.na(absenteeism)), 
                                exponentiate = T)

abs_hurdle_tbs

hurdle(absenteeism ~ gender, data = df_sub, dist = "negbin", zero.dist = "binomial") 
exp(0.2151)

abs_hurdle_tbs <- abs_hurdle_tbs %>% 
  modify_table_body(
    ~ .x %>%
      mutate(estimate_1 = number(estimate_1, accuracy = .01),
             estimate_1 = case_when(is.na(estimate_1) ~ "—", TRUE ~ estimate_1),
             estimate_2 = number(estimate_2, accuracy = .01),
             estimate_2 = case_when(is.na(estimate_2) ~ "—", TRUE ~ estimate_2),
             label = case_when(label == "0" ~ "No",
                               label == "1" ~ "Yes",
                               TRUE ~ label))) %>% 
  bold_labels() %>% 
  as_flex_table() %>% 
  delete_part("footer") 

abs_hurdle_tbs


# header_df_mv <- data.frame(
#   col_keys = names(abs_hurdle_tbs$body$dataset),
#   level_1 = c(" ", "", "Hurdle model (logit)", "Hurdle model (logit)", "Hurdle model (logit)",
#               "Count model (NB)", "Count model (NB)", "Count model (NB)"),
#   level_2 = c(" ", " ",  rep("Odds of Absenteeism", 3), rep("Frequency of Absenteeism", 3)),
#   level_3 = c("Variable", "N", "OR","95% CI", "p-value", "IRR", "95% CI", "p-value"),
#   stringsAsFactors = FALSE
# # )

# #abs_hurdle_tbs <- 
# abs_hurdle_tbs <- abs_hurdle_tbs %>% 
#   # set_header_df(
#   #   mapping = header_df_mv,
#   #   key = "col_keys") %>% 
#     merge_h(part = "header") %>% 
#     align(align = "center", part = "header") %>% 
#     set_caption(caption = "Table X. Univariate hurdle negative binomial regressions predicting absenteeism (n = 213)", 
#                 align_with_table = "FALSE") %>% 
#     bold(part = "header") %>% 
#     align(j = -1, align = "center", part = "all")  %>% 
#     flextable::font(fontname = "Times New Roman", 
#                     part = "all") %>% 
#     bold(part = "header") %>% 
#     hline(i = c(1:3), part = "header") %>%
#     hline_top(border = fp_border(color = "black"), part = "all")  %>% 
#   flextable::set_table_properties(
#     layout = "autofit",   # enables AutoFit to window
#     align = "left",       # aligns table to left margin
#     width = 1             # 100% of page width
#   ) %>% 
#   align(j = 1, align = "left", part = "all") %>% 
#   align(j = -1, align = "center", part = "all") %>% 
#   align(i = 1, align = "left", part = "footer") %>% 
#   flextable::padding(padding = 1, part = "all") %>%
#   flextable::line_spacing(space = 0.9, part = "all") %>% 
#   flextable::fontsize(size = 10, part = "all") %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") %>% 
#   bold(part = "header") 

# abs_hurdle_tbs 
# tf <- tempfile("ex", here(), ".docx")
# save_as_docx(abs_hurdle_tbs, path = tf)
# browseURL(tf)    

### multivariate model ----
abs_formula <- absenteeism ~ gender + alcohol_use_scores + max_other_drug_score + k10_total

abs_mv_mod <- hurdle(abs_formula, data = df_sub, #%>% filter(!is.na(absenteeism)), 
                     dist = "negbin", zero.dist = "binomial")

abs_mv_mod

hurdle(absenteeism ~ gender + alcohol_use_scores + max_other_drug_score + k10_total, data = df_sub, #%>% filter(!is.na(absenteeism)), 
       dist = "negbin", zero.dist = "binomial")

exp(0.03151)
# multicollinearity 
car::vif(lm(abs_formula, data = df_sub))

abs_mv_tbl <- list(
  zero = tbl_regression(abs_mv_mod, 
                        component = "zero_inflated",
                        exponentiate = TRUE,
                        estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                        pvalue_fun = function(x) {
                          p <- style_pvalue(x, digits = 3)
                          # Remove leading zero from both "0.03" and "<0.001"
                          gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                        }, 
                        tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    modify_header(label = "**Variable**", estimate = "**OR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )
    ) %>% 
    add_n() %>% 
    bold_p(t = .05) %>% 
    bold_labels(),
  count = tbl_regression(abs_mv_mod, 
                         exponentiate = TRUE,
                         component = "conditional",
                         estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                         pvalue_fun = function(x) {
                           p <- style_pvalue(x, digits = 3)
                           # Remove leading zero from both "0.03" and "<0.001"
                           gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                         },
                         
                         tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    # modify_caption(paste("Count model (NB) for", .x)) %>% 
    modify_header(label = "**Variable**", estimate = "**IRR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )) %>% 
    
    bold_p(t = .05) %>% 
    bold_labels()
) %>% 
  tbl_merge(
  tab_spanner = c("**Odds of Absenteeism**", 
                  "**Frequency of Absenteeism**")
)
abs_mv_tbl <- abs_mv_tbl %>% 
  modify_table_body(
    ~ .x %>% 
    mutate(
      estimate_1 = number(estimate_1, accuracy = .01),
      estimate_1 = case_when(is.na(estimate_1) ~ "—", TRUE ~ estimate_1),
      estimate_2 = number(estimate_2, accuracy = .01),
      estimate_2 = case_when(is.na(estimate_2) ~ "—", TRUE ~ estimate_2),
      label = case_when(label == "0" ~ "No",
                      label == "1" ~ "Yes",
                      TRUE ~ label))) %>% 
  as_flex_table() %>%  delete_part("footer") 
abs_mv_tbl
# header_df_mv <- data.frame(
#   col_keys = names(abs_mv_tbl$body$dataset),
#   level_1 = c(" ", "", "Hurdle model (logit)", "Hurdle model (logit)", "Hurdle model (logit)",
#               "Count model (NB)", "Count model (NB)", "Count model (NB)"),
#   level_2 = c(" ", " ",  rep("Odds of Absenteeism", 3), rep("Frequency of Absenteeism", 3)),
#   level_3 = c("Variable", "N", "OR","95% CI", "p-value", "IRR", "95% CI", "p-value"),
#   stringsAsFactors = FALSE)

# abs_mv_tbl %>% 
#   set_header_df(
#     mapping = header_df_mv,
#     key = "col_keys") %>% 
#   merge_h(part = "header") %>%
#   align(align = "center", part = "header") %>% 
#   set_caption(caption = "Table X. Multivariate hurdle negative binomial regression predicting absenteeism", 
#               align_with_table = "FALSE") %>% 
#   bold(bold = TRUE, part = "header") %>% 
#   hline_top(border = fp_border(color = "black"), part = "all") %>% 
#   hline(i = c(1:3), part = "header") %>% 
#   flextable::set_table_properties(
#     layout = "autofit",   # enables AutoFit to window
#     align = "left",       # aligns table to left margin
#     width = 1             # 100% of page width
#   ) %>% 
#   align(j = 1, align = "left", part = "all") %>% 
#   align(j = -1, align = "center", part = "all") %>% 
#   align(i = 1, align = "left", part = "footer") %>% 
#   flextable::padding(padding = 1, part = "all") %>%
#   flextable::line_spacing(space = 0.9, part = "all") %>% 
#   flextable::fontsize(size = 10, part = "all") %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") -> abs_mv_tbl
 
  

  
# abs_mv_tbl


# 
# save_as_docx(abs_mv_tbl, path = tf)
# browseURL(tf)



### interaction model ------------------------------------------------------------------------------------------
abs_alc_k10 <- absenteeism ~ gender + alcohol_use_scores*k10_total + max_other_drug_score
abs_drug_k10 <- absenteeism ~ gender + max_other_drug_score*k10_total + alcohol_use_scores

abs_alc_k10_mod <- hurdle(abs_alc_k10, data = df_sub, 
                          dist = "negbin", zero.dist = "binomial")
abs_drug_k10_mod <- hurdle(abs_drug_k10, data = df_sub, 
                           dist = "negbin", zero.dist = "binomial")
# no interaction 
summary(abs_mv_mod)
# with interaction 
summary(abs_alc_k10_mod)
summary(abs_drug_k10_mod)

# lrt 
lrtest(abs_mv_mod, abs_alc_k10_mod)
lrtest(abs_mv_mod, abs_drug_k10_mod)


abs_mv_mod <- hurdle(abs_formula, data = df_sub, #%>% filter(!is.na(absenteeism)), 
                     dist = "negbin", zero.dist = "binomial")

abs_mv_mod

# multicollinearity 
car::vif(lm(abs_formula, data = df_sub))

abs_mv_tbl <- list(
  zero = tbl_regression(abs_mv_mod, 
                        component = "zero_inflated",
                        exponentiate = TRUE,
                        estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                        pvalue_fun = function(x) {
                          p <- style_pvalue(x, digits = 3)
                          # Remove leading zero from both "0.03" and "<0.001"
                          gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                        }, 
                        tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    modify_header(label = "**Variable**", estimate = "**OR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )
    ) %>% 
    add_n() %>% 
    bold_p(t = .05) %>% 
    bold_labels(),
  count = tbl_regression(abs_mv_mod, 
                         exponentiate = TRUE,
                         component = "conditional",
                         estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                         pvalue_fun = function(x) {
                           p <- style_pvalue(x, digits = 3)
                           # Remove leading zero from both "0.03" and "<0.001"
                           gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                         },
                         
                         tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    # modify_caption(paste("Count model (NB) for", .x)) %>% 
    modify_header(label = "**Variable**", estimate = "**IRR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )) %>% 
    
    bold_p(t = .05) %>% 
    bold_labels()
) %>% 
  tbl_merge(
    tab_spanner = c("**Odds of Absenteeism**", 
                    "**Frequency of Absenteeism**")
  )
abs_mv_tbl <- abs_mv_tbl %>% 
  modify_table_body(
    ~ .x %>% 
      mutate(
        estimate_1 = number(estimate_1, accuracy = .01),
        estimate_1 = case_when(is.na(estimate_1) ~ "—", TRUE ~ estimate_1),
        estimate_2 = number(estimate_2, accuracy = .01),
        estimate_2 = case_when(is.na(estimate_2) ~ "—", TRUE ~ estimate_2),
        label = case_when(label == "0" ~ "No",
                          label == "1" ~ "Yes",
                          TRUE ~ label))) %>% 
  as_flex_table() %>%  delete_part("footer") 
# abs_mv_tbl
# header_df_mv <- data.frame(
#   col_keys = names(abs_mv_tbl$body$dataset),
#   level_1 = c(" ", "", "Hurdle model (logit)", "Hurdle model (logit)", "Hurdle model (logit)",
#               "Count model (NB)", "Count model (NB)", "Count model (NB)"),
#   level_2 = c(" ", " ",  rep("Odds of Absenteeism", 3), rep("Frequency of Absenteeism", 3)),
#   level_3 = c("Variable", "N", "OR","95% CI", "p-value", "IRR", "95% CI", "p-value"),
#   stringsAsFactors = FALSE)

# abs_mv_tbl %>% 
#   set_header_df(
#     mapping = header_df_mv,
#     key = "col_keys") %>% 
#   merge_h(part = "header") %>%
#   align(align = "center", part = "header") %>% 
#   set_caption(caption = "Table X. Multivariate hurdle negative binomial regression predicting absenteeism", 
#               align_with_table = "FALSE") %>% 
#   bold(bold = TRUE, part = "header") %>% 
#   hline_top(border = fp_border(color = "black"), part = "all") %>% 
#   hline(i = c(1:3), part = "header") %>% 
#   flextable::set_table_properties(
#     layout = "autofit",   # enables AutoFit to window
#     align = "left",       # aligns table to left margin
#     width = 1             # 100% of page width
#   ) %>% 
#   align(j = 1, align = "left", part = "all") %>% 
#   align(j = -1, align = "center", part = "all") %>% 
#   align(i = 1, align = "left", part = "footer") %>% 
#   flextable::padding(padding = 1, part = "all") %>%
#   flextable::line_spacing(space = 0.9, part = "all") %>% 
#   flextable::fontsize(size = 10, part = "all") %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") -> abs_mv_tbl




# abs_mv_tbl



# save_as_docx(abs_mv_tbl, path = tf)
# browseURL(tf)





## Presenteeism ----
### univariate model ----

#vtable(df_sub)
outcome <- "presenteeism"
predictors <- df_sub %>% dplyr::select(-absenteeism, -presenteeism) %>% names()
var_labels <- c(
  agency_years = "Years in agency",
  age = "Age",
  gender = "Gender",
  alcohol_use_scores = "Alcohol use risk",
  max_other_drug_score = "Maximum other drug use risk",
  k10_total = "Psychological distress"
)


pres_hurdle_tbs <- run_uv_hurdle(outcome, predictors, data = df_sub,
                                exponentiate = T)

pres_hurdle_tbs


pres_hurdle_tbs <- pres_hurdle_tbs %>% 
  modify_table_body(
    ~ .x %>%
      mutate(estimate_1 = number(estimate_1, accuracy = .01),
             estimate_1 = case_when(is.na(estimate_1) ~ "—", TRUE ~ estimate_1),
             estimate_2 = number(estimate_2, accuracy = .01),
             estimate_2 = case_when(is.na(estimate_2) ~ "—", TRUE ~ estimate_2),
             label = case_when(label == "0" ~ "No",
                               label == "1" ~ "Yes",
                               TRUE ~ label))) %>% 
  bold_labels() %>% 
  as_flex_table() %>% 
  delete_part("footer") 

pres_hurdle_tbs


# header_df_mv <- data.frame(
#   col_keys = names(abs_hurdle_tbs$body$dataset),
#   level_1 = c(" ", "", "Hurdle model (logit)", "Hurdle model (logit)", "Hurdle model (logit)",
#               "Count model (NB)", "Count model (NB)", "Count model (NB)"),
#   level_2 = c(" ", " ",  rep("Odds of Presenteeism", 3), rep("Frequency of Presenteeism", 3)),
#   level_3 = c("Variable", "N", "OR","95% CI", "p-value", "IRR", "95% CI", "p-value"),
#   stringsAsFactors = FALSE
# )

# #abs_hurdle_tbs <- 
# pres_hurdle_tbs <- pres_hurdle_tbs %>% 
#   set_header_df(
#     mapping = header_df_mv,
#     key = "col_keys") %>% 
#   merge_h(part = "header") %>% 
#   align(align = "center", part = "header") %>% 
#   set_caption(caption = "Table X. Univariate hurdle negative binomial regressions predicting presenteeism", 
#               align_with_table = "FALSE") %>% 
#   bold(part = "header") %>% 
#   align(j = -1, align = "center", part = "all")  %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") %>% 
#   bold(part = "header") %>% 
#   hline(i = c(1:3), part = "header") %>%
#   hline_top(border = fp_border(color = "black"), part = "all")  %>% 
#   flextable::set_table_properties(
#     layout = "autofit",   # enables AutoFit to window
#     align = "left",       # aligns table to left margin
#     width = 1             # 100% of page width
#   ) %>% 
#   align(j = 1, align = "left", part = "all") %>% 
#   align(j = -1, align = "center", part = "all") %>% 
#   align(i = 1, align = "left", part = "footer") %>% 
#   flextable::padding(padding = 1, part = "all") %>%
#   flextable::line_spacing(space = 0.9, part = "all") %>% 
#   flextable::fontsize(size = 10, part = "all") %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") %>% 
#   bold(part = "header") 

# pres_hurdle_tbs 
# # tf <- tempfile("ex", here(), ".docx")
# # save_as_docx(pres_hurdle_tbs, path = tf)
# browseURL(tf)    

### multivariate model ----
pres_formula <- presenteeism ~ age + agency_years + alcohol_use_scores + max_other_drug_score + k10_total

pres_mv_mod <- hurdle(pres_formula, data = df_sub,
                      zero.dist = "binomial", dist = "negbin")

pres_mv_mod
summary(pres_mv_mod)

# multicollinearity 
car::vif(lm(pres_formula, data = df_sub))

pres_mv_tbl <- list(
  zero = tbl_regression(pres_mv_mod, 
                        component = "zero_inflated",
                        exponentiate = TRUE,
                        estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                        pvalue_fun = function(x) {
                          p <- style_pvalue(x, digits = 3)
                          # Remove leading zero from both "0.03" and "<0.001"
                          gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                        }, 
                        tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    modify_header(label = "**Variable**", estimate = "**OR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )
    ) %>% 
    add_n() %>% 
    bold_p(t = .05) %>% 
    bold_labels(),
  count = tbl_regression(pres_mv_mod, 
                         exponentiate = TRUE,
                         component = "conditional",
                         estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                         pvalue_fun = function(x) {
                           p <- style_pvalue(x, digits = 3)
                           # Remove leading zero from both "0.03" and "<0.001"
                           gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                         },
                         
                         tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    # modify_caption(paste("Count model (NB) for", .x)) %>% 
    modify_header(label = "**Variable**", estimate = "**IRR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )) %>% 
    
    bold_p(t = .05) %>% 
    bold_labels()
) %>% 
  tbl_merge(
    tab_spanner = c("**Odds of Presenteeism**", 
                    "**Frequency of Presenteeism**")
  )
pres_mv_tbl
# pres_mv_tbl <- pres_mv_tbl %>% 
#   modify_table_body(
#     ~ .x %>% 
#       mutate(
#         estimate_1 = number(estimate_1, accuracy = .01),
#         estimate_1 = case_when(is.na(estimate_1) ~ "—", TRUE ~ estimate_1),
#         estimate_2 = number(estimate_2, accuracy = .01),
#         estimate_2 = case_when(is.na(estimate_2) ~ "—", TRUE ~ estimate_2),
#         label = case_when(label == "0" ~ "No",
#                           label == "1" ~ "Yes",
#                           TRUE ~ label))) %>% 
#   as_flex_table() %>%  delete_part("footer") 
# pres_mv_tbl
# header_df_mv <- data.frame(
#   col_keys = names(abs_mv_tbl$body$dataset),
#   level_1 = c(" ", "", "Hurdle model (logit)", "Hurdle model (logit)", "Hurdle model (logit)",
#               "Count model (NB)", "Count model (NB)", "Count model (NB)"),
#   level_2 = c(" ", " ",  rep("Odds of Presenteeism", 3), rep("Frequency of Presenteeism", 3)),
#   level_3 = c("Variable", "N", "OR","95% CI", "p-value", "IRR", "95% CI", "p-value"),
#   stringsAsFactors = FALSE)

# pres_mv_tbl %>% 
#   set_header_df(
#     mapping = header_df_mv,
#     key = "col_keys") %>% 
#   merge_h(part = "header") %>%
#   align(align = "center", part = "header") %>% 
#   set_caption(caption = "Table X. Multivariate hurdle negative binomial regression predicting presenteeism", 
#               align_with_table = "FALSE") %>% 
#   bold(bold = TRUE, part = "header") %>% 
#   hline_top(border = fp_border(color = "black"), part = "all") %>% 
#   hline(i = c(1:3), part = "header") %>% 
#   flextable::set_table_properties(
#     layout = "autofit",   # enables AutoFit to window
#     align = "left",       # aligns table to left margin
#     width = 1             # 100% of page width
#   ) %>% 
#   align(j = 1, align = "left", part = "all") %>% 
#   align(j = -1, align = "center", part = "all") %>% 
#   align(i = 1, align = "left", part = "footer") %>% 
#   flextable::padding(padding = 1, part = "all") %>%
#   flextable::line_spacing(space = 0.9, part = "all") %>% 
#   flextable::fontsize(size = 10, part = "all") %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") -> pres_mv_tbl




# pres_mv_tbl

# # save_as_docx(pres_mv_tbl, path = tf)
# # browseURL(tf)


### interaction model ------------------------------------------------------------------------------------------

pres_alc_k10 <- presenteeism ~ age + agency_years + alcohol_use_scores*k10_total + max_other_drug_score
pres_alc_k10_mod <- hurdle(pres_alc_k10, data = df_sub, 
                          dist = "negbin", zero.dist = "binomial")
pres_drug_k10 <- presenteeism ~ age + agency_years + max_other_drug_score*k10_total + alcohol_use_scores


pres_drug_k10_mod <- hurdle(pres_drug_k10, data = df_sub, 
                           dist = "negbin", zero.dist = "binomial")
# no interaction 
summary(pres_mv_mod)
# with interaction 
summary(pres_alc_k10_mod)
summary(pres_drug_k10_mod)

# lrt 
lrtest(pres_mv_mod, pres_alc_k10_mod)
lrtest(pres_mv_mod, pres_drug_k10_mod)

pres_mv_tbl2 <- list(
  zero = tbl_regression(pres_drug_k10_mod, 
                        component = "zero_inflated",
                        exponentiate = TRUE,
                        estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                        pvalue_fun = function(x) {
                          p <- style_pvalue(x, digits = 3)
                          # Remove leading zero from both "0.03" and "<0.001"
                          gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                        }, 
                        tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    modify_header(label = "**Variable**", estimate = "**OR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )
    ) %>% 
    add_n() %>% 
    bold_p(t = .05) %>% 
    bold_labels(),
  count = tbl_regression(pres_drug_k10_mod, 
                         exponentiate = TRUE,
                         component = "conditional",
                         estimate_fun = function(x) formatC(x, digits = 2, format = "f"),
                         pvalue_fun = function(x) {
                           p <- style_pvalue(x, digits = 3)
                           # Remove leading zero from both "0.03" and "<0.001"
                           gsub("(?<=<)0\\.|^0\\.", ".", p, perl = TRUE)
                         },
                         
                         tidy_fun = broom.helpers::tidy_zeroinfl
  ) %>% 
    # modify_caption(paste("Count model (NB) for", .x)) %>% 
    modify_header(label = "**Variable**", estimate = "**IRR**") %>% 
    modify_table_body(~ .x %>%
                        mutate(
                          label = ifelse(
                            row_type == "label" & variable %in% names(var_labels),
                            var_labels[variable],
                            label
                          ),
                          ci = ifelse(
                            !is.na(estimate) & (conf.high - conf.low > 50),
                            paste0(">50**"),
                            ci
                          )
                        )) %>% 
    
    bold_p(t = .05) %>% 
    bold_labels()
) %>% 
  tbl_merge(
    tab_spanner = c("**Odds of Presenteeism**", 
                    "**Frequency of Presenteeism**")
  )


pres_mv_tbl2 <- pres_mv_tbl2 %>% 
  modify_table_body(
    ~ .x %>% 
      mutate(
        estimate_1 = number(estimate_1, accuracy = .01),
        estimate_1 = case_when(is.na(estimate_1) ~ "—", TRUE ~ estimate_1),
        estimate_2 = number(estimate_2, accuracy = .01),
        estimate_2 = case_when(is.na(estimate_2) ~ "—", TRUE ~ estimate_2),
        label = case_when(label == "0" ~ "No",
                          label == "1" ~ "Yes",
                          TRUE ~ label))) %>% 
  as_flex_table() %>%  delete_part("footer") 
pres_mv_tbl2
# header_df_mv <- data.frame(
#   col_keys = names(abs_mv_tbl$body$dataset),
#   level_1 = c(" ", "", "Hurdle model (logit)", "Hurdle model (logit)", "Hurdle model (logit)",
#               "Count model (NB)", "Count model (NB)", "Count model (NB)"),
#   level_2 = c(" ", " ",  rep("Odds of Presenteeism", 3), rep("Frequency of Presenteeism", 3)),
#   level_3 = c("Variable", "N", "OR","95% CI", "p-value", "IRR", "95% CI", "p-value"),
#   stringsAsFactors = FALSE)

# pres_mv_tbl2 %>% 
#   set_header_df(
#     mapping = header_df_mv,
#     key = "col_keys") %>% 
#   merge_h(part = "header") %>%
#   align(align = "center", part = "header") %>% 
#   set_caption(caption = "Table X. Multivariate hurdle negative binomial regression predicting presenteeism", 
#               align_with_table = "FALSE") %>% 
#   bold(bold = TRUE, part = "header") %>% 
#   hline_top(border = fp_border(color = "black"), part = "all") %>% 
#   hline(i = c(1:3), part = "header") %>% 
#   flextable::set_table_properties(
#     layout = "autofit",   # enables AutoFit to window
#     align = "left",       # aligns table to left margin
#     width = 1             # 100% of page width
#   ) %>% 
#   align(j = 1, align = "left", part = "all") %>% 
#   align(j = -1, align = "center", part = "all") %>% 
#   align(i = 1, align = "left", part = "footer") %>% 
#   flextable::padding(padding = 1, part = "all") %>%
#   flextable::line_spacing(space = 0.9, part = "all") %>% 
#   flextable::fontsize(size = 10, part = "all") %>% 
#   flextable::font(fontname = "Times New Roman", 
#                   part = "all") -> pres_mv_tbl2




# pres_mv_tbl2

# # save_as_docx(pres_mv_tbl2, path = tf)
# # browseURL(tf)


# margins ----------------------------------------------------------------------------------------------------
library(marginaleffects)
library(emmeans)
df_sub <- df_sub %>% mutate(presenteeism1 = ifelse(presenteeism >= 1, 1, 0))
means <- df_sub %>%
  summarise(
    age_mean = mean(age, na.rm = TRUE),
    agency_years_mean = mean(agency_years, na.rm = TRUE),
    alcohol_mean = mean(alcohol_use_scores, na.rm = TRUE),
    k10_mean = mean(k10_total, na.rm = TRUE),
    k10_sd   = sd(k10_total, na.rm = TRUE))

k10_low  <- with(means, k10_mean - k10_sd)
k10_high <- with(means, k10_mean + k10_sd)
m <- glm(presenteeism1 ~ age + agency_years + max_other_drug_score*k10_total + alcohol_use_scores, data = df_sub, family = binomial)
summary(m)
tidy(m, exponentiate = TRUE, conf.int = .95)

preds <- avg_predictions(m, by = c("max_other_drug_score", "k10_total"), newdata = datagrid(max_other_drug_score = 0:40,
                                                                                   k10_total = c(k10_low, k10_high)))
plot_predictions(m, condition = list(max_other_drug_score = 0:40, k10_total = c(k10_low, k10_high)))

avg_comparisons(m, variables = "max_other_drug_score")
comparisons(m, variables = "max_other_drug_score", newdata = datagrid(k10_total = c(k10_low, k10_high)))
comparisons(m, hypothesis = "b2-b1=0", 
            variables = "max_other_drug_score",
            newdata = datagrid(k10_total = c(k10_low, k10_high)))


# on average:
#' moving from 0 to 1 in max other drug score is associated with a 1.9% increase in the probability of presenteeism 
#  moving from 0 to 1 in the max other drug score is associated with a change of 0.04 (4.0%) in presenteeism when psych distress is low 
#' but when psych distress is high, there is no significant change in the odds of presenteeism. 
#' comparisons are distinguishable. Psych distress has a moderating effect on the relationship between other drug use and presenteeism 
#' 
#' 


# plot -------------------------------------------------------------------------------------------------------
library(ggeffects)
library(marginaleffects)

# # smhr version 
# colourpal <- c("#C66951", "#d95f02", "#BF974D", "#e7298a", "#e6ab02")
# plot_predictions(m, condition = list(max_other_drug_score = 0:40, k10_total = c(k10_low, k10_high))) + 
#   #geom_line(size = 2) +
#   scale_y_continuous(labels = label_percent(accuracy = 1), breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
#   scale_colour_manual(values = c("#1b9e77","#7570b3"), labels = c("Low (mean + 1 SD)", "High (mean + 1 SD)"))+
#   scale_fill_manual(values = c("#1b9e77","#7570b3"), labels = NULL) +
#   guides(fill = "none")+
#   theme(legend.position = "right")+
#   labs(
#     x = "Maximum other drug use risk score (ASSIST)",
#     y = "Predicted probability of presenteeism",
#     color = "Psychological distress (K10+)",
#     fill = "Psychological distress (K10+)"
#   ) +
#   theme_pubr(base_size = 13) +
#   labs_pubr()
  


# preds <- preds %>% mutate(distress = if_else(k10_total < 13.9, "Low","High"),
#                           distress = factor(distress, levels = c("Low", "High"), labels = c("Low (mean - 1 SD)", "High (mean + 1 SD)"),
#                                             ordered = TRUE))
# preds
# preds %>% ggplot(aes(x = max_other_drug_score, y = estimate, color = distress, fill = distress)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = NA) 
#   scale_y_continuous(labels = label_percent(accuracy = 1), breaks = seq(0, 1, by = 0.1), limits = c(0, 1) )
# scale_x_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) + 
#   scale_colour_manual(values = c("#1b9e77","#7570b3")) +    
#   scale_fill_manual(values = c("#1b9e77","#7570b3")) +
  
#   labs(
#     x = "Maximum other drug use risk score (ASSIST)",
#     y = "Predicted probability of presenteeism",
#     color = "Psychological distress (K10+)",
#     fill = "Psychological distress (K10+)"
#   ) +
#   theme_pubr(base_size = 13) +
#   labs_pubr()


# existing version 

colourpal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#e6ab02")
preds <- preds %>% mutate(distress = if_else(k10_total < 13.9, "Low","High"),
                 distress = factor(distress, levels = c("Low", "High"), labels = c("Low (mean - 1 SD)", "High (mean + 1 SD)"),
                                   ordered = TRUE))
preds  
preds %>% ggplot(aes(x = max_other_drug_score, y = estimate, color = distress, fill = distress)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = NA) +
  scale_y_continuous(labels = label_percent(accuracy = 1), breaks = seq(0, 1, by = 0.1), limits = c(0, 1) ) +
  scale_x_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) + 
  scale_colour_manual(values = c("#1b9e77","#7570b3")) +    
  scale_fill_manual(values = c("#1b9e77","#7570b3")) +
    
    labs(
      x = "Maximum other drug use risk score (ASSIST)",
      y = "Predicted probability of presenteeism",
      color = "Psychological distress (K10+)",
      fill = "Psychological distress (K10+)"
    ) +
    theme_pubr(base_size = 13) +
    labs_pubr()



  

# old --------------------------------------------------------------------------------------------------------



# # marginal plot ----------------------------------------------------------------------------------------------
# drug_mean <- mean(df_sub$max_other_drug_score, na.rm = T)
# drug_sd <- sd(df_sub$max_other_drug_score, na.rm = T)
# drug_low <- drug_mean-drug_sd
# drug_high <- drug_mean+drug_sd

# mfx <- slopes(
#   pres_drug_k10_mod,
#   variables = "max_other_drug_score",
#   newdata = datagrid(
#     max_other_drug_score = c(drug_low, drug_high),  # Different drug use levels
#     k10_total = c(k10_low, k10_high),         # Low and high K10
#     # Hold others at means
#     age = mean(df_sub$age, na.rm = TRUE),
#     agency_years = mean(df_sub$agency_years, na.rm = TRUE),
#     alcohol_use_scores = mean(df_sub$alcohol_use_scores, na.rm = TRUE)), type = "zero")
  

# mfx$drug_level <- factor(mfx$max_other_drug_score,
#                          levels = c(drug_low, drug_high),
#                          labels = c("Low (mean - 1SD)", "High (mean + 1SD)"), ordered = TRUE)

# mfx$k10_level <- factor(mfx$k10_total,
#                         levels = c(k10_low, k10_high),
#                         labels = c("Low (mean - 1SD)", "High (mean + 1SD)"), ordered = TRUE)


# (p3 <- ggplot(mfx, aes(x = drug_level, y = estimate, 
#                 color = k10_level, group = k10_level)) +
#   geom_point(size = 3, position = position_dodge(width = 0.3)) +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
#                 width = 0.2, position = position_dodge(width = 0.3)) +
#     scale_y_continuous(breaks = seq(-0.01, 0.08, .01), limits = c(-0.01, 0.08)) +
#   geom_line(position = position_dodge(width = 0.3)) +
#   scale_color_manual(values = c("#1b9e77", "#7570b3")) +
#   labs(
#     x = "Maximum other drug use risk (ASSIST)",
#     y = "Marginal effect of other drug use risk on presenteeism\n(Change in probability per unit increase)",
#     color = "Psychological Distress (K10+)"
#   ) +
#   theme_pubr(base_size = 12) +
#   labs_pubr() +
#   #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5))

