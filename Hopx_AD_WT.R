cat("\014")
rm(list = ls())

# === Step 0: 设置统一路径变量 ===
base_path <- '/Users/baektony/Desktop/scRNA_seq/AD scRNA-seq/AD.data5.GSE185553.nature'

# === Step 1: 加载包 ===
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

# === Step 2: 加载 Seurat 对象 ===
astro_obj <- readRDS(file.path(base_path, "astrocyte.Reactive_nonReactive.rds"))

# === Step 3: 分 AD / WT 群体 ===
astro_AD <- subset(astro_obj, subset = grepl("^AD", orig.ident))
astro_WT <- subset(astro_obj, subset = grepl("^wt", orig.ident))

# === Step 4: 定义颜色样式 ===
violin_colors <- c("nonReactive" = "grey40", "Reactive" = "#E64B35")
point_colors <- c("nonReactive" = "grey60", "Reactive" = "#F4A582")
point_borders <- c("nonReactive" = "black", "Reactive" = "#E64B35")

# === Step 5: 函数：绘图和统计输出 ===
plot_violin_HOPX <- function(seurat_obj, group_label) {
  df <- FetchData(seurat_obj, vars = c("HOPX", "cell_anno")) %>%
    filter(cell_anno %in% c("nonReactive", "Reactive"))
  df$cell_anno <- factor(df$cell_anno, levels = c("nonReactive", "Reactive"))
  
  # 打印统计信息
  group_stats <- df %>%
    group_by(cell_anno) %>%
    summarise(
      n = n(),
      mean = mean(HOPX),
      sd = sd(HOPX),
      sem = sd / sqrt(n()),
      .groups = "drop"
    )
  print(paste0("====== ", group_label, " ======"))
  print(group_stats)
  
  # 画图
  p <- ggplot(df, aes(x = cell_anno, y = HOPX)) +
    geom_violin(aes(fill = cell_anno), trim = FALSE, alpha = 0.4) +
    geom_jitter(aes(color = cell_anno), shape = 21, size = 1, stroke = 0.4, width = 0.2, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = violin_colors) +
    scale_color_manual(values = point_borders) +
    labs(
      x = NULL,
      y = "HOPX expression (log-normalized)",
      title = paste0("HOPX expression in ", group_label, " astrocytes")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none"
    )
  
  return(p)
}

# === Step 6: 分别绘图 ===
p_ad <- plot_violin_HOPX(astro_AD, "AD")
p_wt <- plot_violin_HOPX(astro_WT, "WT")

# === Step 7: 拼图并保存 ===
p_combined <- p_ad + p_wt + plot_layout(ncol = 2)

ggsave(
  filename = file.path(base_path, "HOPX_violin_Reactive_vs_NonReactive_AD_WT.pdf"),
  plot = p_combined,
  device = cairo_pdf,
  width = 10, height = 5, units = "in"
)
