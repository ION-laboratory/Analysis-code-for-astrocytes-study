cat("\014")
rm(list = ls())

# === Step 0: 设置统一路径变量 ===
base_path <- '/Users/baektony/Desktop/scRNA_seq/AD scRNA-seq/AD_data2_syn18485175'

# === Step 1: 加载包 ===
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

# === Step 2: 加载 Seurat 对象 ===
astro_obj <- readRDS(file.path(base_path, "astrocyte.Reactive_nonReactive.rds"))

# === Step 3: 小提琴图：HOPX in Reactive vs NonReactive ===
df <- FetchData(astro_obj, vars = c("HOPX", "cell_anno"))
df$cell_anno <- factor(df$cell_anno, levels = c("nonReactive", "Reactive"))

# 设置颜色
violin_colors1 <- c("nonReactive" = "grey40", "Reactive" = "#E64B35")
point_colors1 <- c("nonReactive" = "grey60", "Reactive" = "#F4A582")
point_borders1 <- c("nonReactive" = "black", "Reactive" = "#E64B35")

p1 <- ggplot(df, aes(x = cell_anno, y = HOPX)) +
  geom_violin(aes(fill = cell_anno), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = cell_anno, fill = cell_anno), 
              shape = 21, size = 1, stroke = 0.4, width = 0.2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = violin_colors1) +
  scale_color_manual(values = point_borders1) +
  labs(
    x = NULL,
    y = "HOPX expression (log-normalized)",
    title = "HOPX expression in Reactive vs Non-Reactive astrocytes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

print(p1)

ggsave(
  filename = file.path(base_path, "all_sample_HOPX_violin_Reactive_vs_NonReactive.pdf"),
  plot = p1,
  device = cairo_pdf,
  width = 5, height = 5, units = "in"
)


# === Step 4: Subset Reactive astrocytes and score ===
genes_pan <- toupper(c("Lcn2", "Steap4", "S1pr3", "Timp1", "Hspb1", "Cxcl10", 
                       "Cd44", "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap"))

genes_tox <- toupper(c("C3", "H2-T23", "Serping1", "H2-D1", "Ggta1", "Iigp1", "Gbp2", 
                       "Fbln5", "Ugt1a1", "Fkbp5", "Psmb8", "Srgn", "Amigo2", "Ifitm3", 
                       "Gsr", "Tapbp", "H2-T10"))

genes_pro <- toupper(c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", 
                       "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14", "Lgals3", "Flnc", 
                       "Tubb6", "Srxn1", "Tnfrsf12a", "Gcnt2"))

astro_reactive <- subset(astro_obj, subset = cell_anno == "Reactive")
astro_reactive <- AddModuleScore(astro_reactive, features = list(genes_pan), name = "pan_reactive")
astro_reactive <- AddModuleScore(astro_reactive, features = list(genes_tox), name = "neurotoxic")
astro_reactive <- AddModuleScore(astro_reactive, features = list(genes_pro), name = "neuroprotective")

scores <- FetchData(astro_reactive, vars = c("pan_reactive1", "neurotoxic1", "neuroprotective1"))

# === Step 5: 分类函数（含中间型） ===
assign_subgroup <- function(row, threshold = 0.05) {
  pan <- row["pan_reactive1"]
  tox <- row["neurotoxic1"]
  pro <- row["neuroprotective1"]
  
  if (pan < 0 & tox < 0 & pro < 0) {
    return("nonreactive")
  } else if (pan > 0 & (tox - pro) > threshold) {
    return("neurotoxic")
  } else if (pan > 0 & (pro - tox) > threshold) {
    return("neuroprotective")
  } else if (pan > 0) {
    return("intermediate")
  } else {
    return("unassigned")
  }
}

astro_reactive$reactive_subgroup <- apply(scores, 1, assign_subgroup)

# === Step 6: 提取表达数据用于绘图 ===
df_plot <- FetchData(astro_reactive, vars = c("HOPX", "reactive_subgroup")) %>%
  filter(reactive_subgroup %in% c("neurotoxic", "neuroprotective"))

df_plot$reactive_subgroup <- factor(df_plot$reactive_subgroup, 
                                    levels = c("neurotoxic", "neuroprotective"))

# 设置颜色
violin_colors2 <- c("neuroprotective" = "grey40", "neurotoxic" = "#E64B35")
point_colors2 <- c("neuroprotective" = "grey60", "neurotoxic" = "#F4A582")
point_borders2 <- c("neuroprotective" = "black", "neurotoxic" = "#E64B35")

p2 <- ggplot(df_plot, aes(x = reactive_subgroup, y = HOPX)) +
  geom_violin(aes(fill = reactive_subgroup), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = reactive_subgroup, fill = reactive_subgroup), 
              shape = 21, size = 1, stroke = 0.4, width = 0.2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = violin_colors2) +
  scale_color_manual(values = point_borders2) +
  labs(
    x = NULL,
    y = "HOPX expression (log-normalized)",
    title = "HOPX expression in Neurotoxic vs Neuroprotective Astrocytes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

print(p2)

ggsave(
  filename = file.path(base_path, "all_sample_HOPX_violin_toxic_vs_protective.pdf"),
  plot = p2,
  device = cairo_pdf,
  width = 5, height = 5, units = "in"
)

# 画三组数据
# Step 7: 合并 nonReactive、neurotoxic 和 neuroprotective 三组数据
df_nonreactive <- FetchData(astro_obj, vars = c("HOPX", "cell_anno")) %>%
  filter(cell_anno == "nonReactive") %>%
  mutate(subgroup = "nonReactive")

df_tox_pro <- FetchData(astro_reactive, vars = c("HOPX", "reactive_subgroup")) %>%
  filter(reactive_subgroup %in% c("neurotoxic", "neuroprotective")) %>%
  rename(subgroup = reactive_subgroup)

df_merge <- bind_rows(df_nonreactive, df_tox_pro)
df_merge$subgroup <- factor(df_merge$subgroup, levels = c("nonReactive", "neurotoxic", "neuroprotective"))

violin_colors3 <- c("nonReactive" = "grey40", 
                    "neurotoxic" = "#E64B35", 
                    "neuroprotective" = "#4DBBD5")

point_borders3 <- c("nonReactive" = "black", 
                    "neurotoxic" = "#E64B35", 
                    "neuroprotective" = "#4DBBD5")

p3 <- ggplot(df_merge, aes(x = subgroup, y = HOPX)) +
  geom_violin(aes(fill = subgroup), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = subgroup), 
              shape = 21, size = 1, stroke = 0.4, width = 0.2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     comparisons = list(c("nonReactive", "neurotoxic"),
                                        c("nonReactive", "neuroprotective"),
                                        c("neurotoxic", "neuroprotective"))) +
  scale_fill_manual(values = violin_colors3) +
  scale_color_manual(values = point_borders3) +
  labs(
    x = NULL,
    y = "HOPX expression (log-normalized)",
    title = "HOPX expression in nonReactive, neurotoxic and neuroprotective astrocytes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

print(p3)

ggsave(
  filename = file.path(base_path, "all_sample_HOPX_violin_nonReactive_toxic_protective.pdf"),
  plot = p3,
  device = cairo_pdf,
  width = 6, height = 5, units = "in"
)



group_stats_p3 <- df_merge %>%
  group_by(subgroup) %>%
  summarise(
    n = n(),
    mean = mean(HOPX),
    std = sd(HOPX),
    .groups = "drop"
  )

print(group_stats_p3)


# 分析AD病人数据中，上述三个值的差异 ------------------------------------------------------

# 提取所有 orig.ident
all_samples <- unique(astro_obj$orig.ident)

# 选择包含 "AD" 的样本名
ad_samples <- grep("AD", all_samples, value = TRUE)

# === Step B: 提取 AD 病人的 nonReactive astrocytes ===
df_nonreactive_ad <- FetchData(astro_obj, vars = c("HOPX", "cell_anno", "orig.ident")) %>%
  filter(orig.ident %in% ad_samples, cell_anno == "nonReactive") %>%
  mutate(subgroup = "nonReactive")

# === Step C: 提取 AD 病人的 neurotoxic 和 neuroprotective astrocytes ===
df_tox_pro_ad <- FetchData(astro_reactive, vars = c("HOPX", "reactive_subgroup", "orig.ident")) %>%
  filter(orig.ident %in% ad_samples, reactive_subgroup %in% c("neurotoxic", "neuroprotective")) %>%
  rename(subgroup = reactive_subgroup)

# === Step D: 合并数据并设置因子顺序 ===
df_ad_merge <- bind_rows(df_nonreactive_ad, df_tox_pro_ad)
df_ad_merge$subgroup <- factor(df_ad_merge$subgroup, levels = c("nonReactive", "neurotoxic", "neuroprotective"))

library(ggpubr)

violin_colors_ad <- c("nonReactive" = "grey40", 
                      "neurotoxic" = "#E64B35", 
                      "neuroprotective" = "#4DBBD5")

p_ad <- ggplot(df_ad_merge, aes(x = subgroup, y = HOPX)) +
  geom_violin(aes(fill = subgroup), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = subgroup), shape = 21, size = 1, stroke = 0.4,
              width = 0.2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("nonReactive", "neurotoxic"),
                                        c("nonReactive", "neuroprotective"),
                                        c("neurotoxic", "neuroprotective"))) +
  scale_fill_manual(values = violin_colors_ad) +
  scale_color_manual(values = violin_colors_ad) +
  labs(
    x = NULL,
    y = "HOPX expression (log-normalized)",
    title = "AD HOPX expr in Astrocytes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

print(p_ad)
ggsave(
  filename = file.path(base_path, "AD_HOPX_nonReactive_toxic_protective.pdf"),
  plot = p_ad,
  device = cairo_pdf,
  width = 6, height = 5, units = "in"
)



group_stats_ad <- df_ad_merge %>%
  group_by(subgroup) %>%
  summarise(
    n = n(),
    mean = mean(HOPX),
    std = sd(HOPX),
    .groups = "drop"
  )

print(group_stats_ad)




# 分析AD和WT中 HOPX表达差异 -------------------------------------------------------

# === Step 1: 提取所有 astrocyte 的 HOPX 表达和样本信息 ===
df_astro_all <- FetchData(astro_obj, vars = c("HOPX", "orig.ident"))

# === Step 2: 给每个样本打标签：AD or WT ===
df_astro_all$group <- ifelse(grepl("AD", df_astro_all$orig.ident), "AD", "WT")
df_astro_all$group <- factor(df_astro_all$group, levels = c("WT", "AD"))

# === Step 3: 绘图比较 ===
p_ad_vs_wt <- ggplot(df_astro_all, aes(x = group, y = HOPX)) +
  geom_violin(aes(fill = group), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = group), shape = 21, size = 1, stroke = 0.4,
              width = 0.2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("WT" = "grey40", "AD" = "#E64B35")) +
  scale_color_manual(values = c("WT" = "black", "AD" = "#E64B35")) +
  labs(
    x = NULL,
    y = "HOPX expression (log-normalized)",
    title = "HOPX expression in astrocytes (WT vs AD, cell-level)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

print(p_ad_vs_wt)

ggsave(
  filename = file.path(base_path, "AD_WT_HOPX_expr.pdf"),
  plot = p_ad_vs_wt,
  device = cairo_pdf,
  width = 6, height = 5, units = "in"
)

df_astro_all %>%
  group_by(group) %>%
  summarise(
    n_cells = n(),
    mean_HOPX = mean(HOPX),
    std_HOPX = sd(HOPX),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  print()



# 画出WT样本中各类星胶HOPX的表达qingk -------------------------------------------------
# 获取所有样本ID
all_samples <- unique(astro_obj$orig.ident)

# 选出不含 "AD" 的 WT 样本
wt_samples <- setdiff(all_samples, grep("AD", all_samples, value = TRUE))

df_nonreactive_wt <- FetchData(astro_obj, vars = c("HOPX", "cell_anno", "orig.ident")) %>%
  filter(orig.ident %in% wt_samples, cell_anno == "nonReactive") %>%
  mutate(subgroup = "nonReactive")

df_tox_pro_wt <- FetchData(astro_reactive, vars = c("HOPX", "reactive_subgroup", "orig.ident")) %>%
  filter(orig.ident %in% wt_samples, reactive_subgroup %in% c("neurotoxic", "neuroprotective")) %>%
  rename(subgroup = reactive_subgroup)

df_wt_merge <- bind_rows(df_nonreactive_wt, df_tox_pro_wt)
df_wt_merge$subgroup <- factor(df_wt_merge$subgroup, levels = c("nonReactive", "neurotoxic", "neuroprotective"))

p_wt <- ggplot(df_wt_merge, aes(x = subgroup, y = HOPX)) +
  geom_violin(aes(fill = subgroup), trim = FALSE, alpha = 0.4) +
  geom_jitter(aes(color = subgroup), shape = 21, size = 1, stroke = 0.4,
              width = 0.2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("nonReactive", "neurotoxic"),
                       c("nonReactive", "neuroprotective"),
                       c("neurotoxic", "neuroprotective")
                     ),
                     label = "p.format",
                     p.adjust.method = "BH") +
  scale_fill_manual(values = c("nonReactive" = "grey40", 
                               "neurotoxic" = "#E64B35", 
                               "neuroprotective" = "#4DBBD5")) +
  scale_color_manual(values = c("nonReactive" = "black", 
                                "neurotoxic" = "#E64B35", 
                                "neuroprotective" = "#4DBBD5")) +
  labs(
    x = NULL,
    y = "HOPX expression (log-normalized)",
    title = "HOPX expression in astrocyte subtypes (WT only)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

print(p_wt)

ggsave(
  filename = file.path(base_path, "WT_HOPX_expr in astro2.pdf"),
  plot = p_wt,
  device = cairo_pdf,
  width = 6, height = 5, units = "in"
)


# option，只有ealry AD和late AD分类的数据，才使用下面代码
# 
# # === Step 1: 提取 astrocyte 的 HOPX 表达和样本信息 ===
# df_astro_all <- FetchData(astro_obj, vars = c("HOPX", "orig.ident"))
# 
# # === Step 2: 给每个样本打标签：earlyAD or lateAD ===
# df_astro_all$group <- case_when(
#   grepl("^eAD", df_astro_all$orig.ident) ~ "earlyAD",
#   grepl("^L_AD", df_astro_all$orig.ident) ~ "lateAD",
#   TRUE ~ NA_character_  # 剔除 WT 等非 AD 样本
# )
# df_astro_all <- df_astro_all %>% filter(!is.na(group))  # 去掉NA
# df_astro_all$group <- factor(df_astro_all$group, levels = c("earlyAD", "lateAD"))
# 
# # === Step 3: 绘图比较 ===
# p_ad_stage <- ggplot(df_astro_all, aes(x = group, y = HOPX)) +
#   geom_violin(aes(fill = group), trim = FALSE, alpha = 0.4) +
#   geom_jitter(aes(color = group), shape = 21, size = 1, stroke = 0.4,
#               width = 0.2, alpha = 0.7) +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
#   stat_compare_means(method = "wilcox.test", label = "p.format") +
#   scale_fill_manual(values = c("earlyAD" = "#4DBBD5", "lateAD" = "#E64B35")) +
#   scale_color_manual(values = c("earlyAD" = "#4DBBD5", "lateAD" = "#E64B35")) +
#   labs(
#     x = NULL,
#     y = "HOPX expression (log-normalized)",
#     title = "HOPX expression in astrocytes (early vs late AD, cell-level)"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 14),
#     axis.text = element_text(size = 12),
#     axis.title.x = element_text(size = 13),
#     axis.title.y = element_text(size = 13),
#     legend.position = "none",
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black", size = 0.5)
#   )
# 
# print(p_ad_stage)
# 
# # === Step 4: 导出 PDF ===
# ggsave(
#   filename = file.path(base_path, "early_vs_late_AD_HOPX_expr.pdf"),
#   plot = p_ad_stage,
#   device = cairo_pdf,
#   width = 6, height = 5, units = "in"
# )
# 
# # === Step 5: 输出每组表达均值、标准差和细胞数 ===
# df_astro_all %>%
#   group_by(group) %>%
#   summarise(
#     n_cells = n(),
#     mean_HOPX = mean(HOPX),
#     std_HOPX = sd(HOPX),
#     .groups = "drop"
#   ) %>%
#   mutate(across(where(is.numeric), round, 3)) %>%
#   print()







# # 分析sample-level的neurotoxic比例 ---------------------------------------------
# 
# # === Step 1: 从完整 astro_obj 提取 subgroup 和 orig.ident ===
# astro_obj$group <- ifelse(grepl("^AD", astro_obj$orig.ident), "AD", "WT")
# 
# # 只在 Reactive 中的细胞计算 neurotoxic subgroup，其余标为 "other"
# astro_obj$reactive_subgroup <- NA
# reactive_cells <- colnames(astro_reactive)
# astro_obj$reactive_subgroup[reactive_cells] <- astro_reactive$reactive_subgroup
# 
# # === Step 2: 统计每个 sample 的 neurotoxic astrocytes 占比 ===
# df_neurotoxic_ratio <- astro_obj@meta.data %>%
#   group_by(orig.ident, group) %>%
#   summarise(
#     total = n(),
#     neurotoxic = sum(reactive_subgroup == "neurotoxic", na.rm = TRUE),
#     ratio = neurotoxic / total,
#     .groups = "drop"
#   )
# 
# # === Step 3: 统计每组的样本数、均值、标准差 ===
# summary_stats <- df_neurotoxic_ratio %>%
#   group_by(group) %>%
#   summarise(
#     n = n(),
#     mean_ratio = mean(ratio),
#     sd_ratio = sd(ratio),
#     .groups = "drop"
#   )
# 
# print(summary_stats)
# 
# # === Step 4: 绘制柱状图（neurotoxic占比） ===
# library(ggpubr)
# 
# p_ratio <- ggplot(df_neurotoxic_ratio, aes(x = group, y = ratio, fill = group)) +
#   geom_boxplot(width = 0.5, alpha = 0.4, outlier.shape = NA) +
#   geom_jitter(aes(color = group), shape = 21, size = 2, width = 0.2, stroke = 0.4, fill = "white") +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black") +
#   stat_compare_means(method = "wilcox.test", label = "p.format") +
#   scale_fill_manual(values = c("AD" = "#E64B35", "WT" = "grey40")) +
#   scale_color_manual(values = c("AD" = "#E64B35", "WT" = "black")) +
#   labs(
#     title = "Neurotoxic Astrocytes / Total Astrocytes (Sample-level)",
#     x = NULL,
#     y = "Proportion of Neurotoxic Astrocytes"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 14),
#     axis.text = element_text(size = 12),
#     legend.position = "none"
#   )
# 
# # 保存高清 PDF
# ggsave(
#   filename = file.path(base_path, "neurotoxic_ratio_of_total_Astrocytes_AD_vs_WT.pdf"),
#   plot = p_ratio,
#   device = cairo_pdf,
#   width = 5, height = 5, units = "in"
# )


# library(dplyr)
# 
# # 提取需要字段
# df_group3 <- FetchData(astro_obj, vars = c("HOPX", "orig.ident")) %>%
#   mutate(group = case_when(
#     grepl("^wt", orig.ident) ~ "WT",
#     grepl("^eAD", orig.ident) ~ "early_AD",
#     grepl("^L_AD", orig.ident) ~ "late_AD",
#     TRUE ~ NA_character_
#   )) %>%
#   filter(!is.na(group))
# 
# df_group3$group <- factor(df_group3$group, levels = c("WT", "early_AD", "late_AD"))
# 
# library(ggplot2)
# library(ggpubr)
# 
# # 定义成对比较组
# pairwise_comparisons <- list(
#   c("WT", "early_AD"),
#   c("WT", "late_AD"),
#   c("early_AD", "late_AD")
# )
# 
# p_group3_wilcox <- ggplot(df_group3, aes(x = group, y = HOPX)) +
#   geom_violin(aes(fill = group), trim = FALSE, alpha = 0.4) +
#   geom_jitter(aes(color = group), shape = 21, size = 1, stroke = 0.4,
#               width = 0.2, alpha = 0.7) +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
#   stat_compare_means(method = "wilcox.test",
#                      comparisons = pairwise_comparisons,
#                      label = "p.signif",
#                      p.adjust.method = "BH") +
#   scale_fill_manual(values = c("WT" = "grey40", "early_AD" = "#4DBBD5", "late_AD" = "#E64B35")) +
#   scale_color_manual(values = c("WT" = "black", "early_AD" = "#4DBBD5", "late_AD" = "#E64B35")) +
#   labs(
#     x = NULL,
#     y = "HOPX expression (log-normalized)",
#     title = "HOPX expression in WT, early AD, and late AD astrocytes (Wilcoxon test)"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 14),
#     axis.text = element_text(size = 12),
#     legend.position = "none"
#   )
# 
# print(p_group3_wilcox)







