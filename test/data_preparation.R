# data_preparation.R
# 核心数据准备脚本：加载系统发育树、染色体计数和共线性信息
# Author: GitHub Copilot
# Date: 2025-04-17

library(ape)
library(dplyr)
library(readr)

# 加载系统发育树
dp_load_tree <- function(tree_file) {
  phy <- read.tree(tree_file)
  return(phy)
}

# 加载染色体计数数据，自动兼容常见列名
# 支持 species,count 或 species,chromosome_number 等
# 返回标准化为 species,chromosome_number 的数据框
dp_load_counts <- function(count_file) {
  counts <- read_csv(count_file, show_col_types = FALSE)
  # 自动兼容常见列名
  if (!"chromosome_number" %in% colnames(counts)) {
    possible_names <- c("count", "chr_num", "chrom_num", "chromosomeCount")
    found <- intersect(possible_names, colnames(counts))
    if (length(found) > 0) {
      colnames(counts)[colnames(counts) == found[1]] <- "chromosome_number"
    } else {
      stop("counts数据框必须包含chromosome_number列或常见别名（如count）")
    }
  }
  if (!"species" %in% colnames(counts)) {
    stop("counts数据框必须包含species列")
  }
  return(counts)
}

# 加载并处理共线性映射文件（all_bidirectional.tsv）
dp_load_synteny <- function(synteny_file) {
  synteny <- read_tsv(synteny_file, show_col_types = FALSE)
  # 计算每对物种的1:1映射比例
  synteny_summary <- synteny %>%
    group_by(species_A, species_B) %>%
    summarize(
      n_1to1 = sum(bidirectional_mapping_type == "1:1"),
      n_total = n(),
      ratio_1to1 = n_1to1 / n_total
    )
  return(list(raw = synteny, summary = synteny_summary))
}

# 数据一致性检查
dp_check_species_consistency <- function(phy, counts, synteny) {
  tree_species <- phy$tip.label
  count_species <- unique(counts$species)
  synteny_species <- unique(c(synteny$raw$species_A, synteny$raw$species_B))
  common_species <- Reduce(intersect, list(tree_species, count_species, synteny_species))
  if (length(common_species) == 0) {
    stop("No common species found across tree, counts, and synteny data!")
  }
  message("共有物种数: ", length(common_species))
  return(common_species)
}

# 主数据准备函数（支持外部参数传递）
data_prepare <- function(tree_file = NULL, count_file = NULL, synteny_file = NULL, args = NULL) {
  # 支持通过args列表传递参数
  if (!is.null(args)) {
    if (!is.null(args$tree_file)) tree_file <- args$tree_file
    if (!is.null(args$count_file)) count_file <- args$count_file
    if (!is.null(args$synteny_file)) synteny_file <- args$synteny_file
  }
  if (is.null(tree_file) || is.null(count_file) || is.null(synteny_file)) {
    stop("tree_file, count_file, synteny_file 必须全部指定！")
  }
  phy <- dp_load_tree(tree_file)
  counts <- dp_load_counts(count_file)
  synteny <- dp_load_synteny(synteny_file)
  common_species <- dp_check_species_consistency(phy, counts, synteny)
  # 可选：过滤数据，仅保留共有物种
  phy <- drop.tip(phy, setdiff(phy$tip.label, common_species))
  counts <- counts %>% filter(species %in% common_species)
  synteny$raw <- synteny$raw %>% filter(species_A %in% common_species, species_B %in% common_species)
  synteny$summary <- synteny$summary %>% filter(species_A %in% common_species, species_B %in% common_species)
  return(list(phy = phy, counts = counts, synteny = synteny))
}

# 用法示例：
# result <- data_prepare("tree_file.nwk", "chromosome_counts.csv", "all_bidirectional.tsv")
# str(result)

# ---- 命令行运行支持 ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat("用法: Rscript data_preparation.R tree_file.nwk chromosome_counts.csv all_bidirectional.tsv\n")
    quit(status = 1)
  }
  tree_file <- args[1]
  count_file <- args[2]
  synteny_file <- args[3]
  result <- data_prepare(tree_file = tree_file, count_file = count_file, synteny_file = synteny_file)
  # 可选：保存结果到RDS文件
  saveRDS(result, file = "data_preparation_result.rds")
  # 新增：保存关键数据为csv/tsv文件
  write.csv(result$counts, file = "prepared_counts.csv", row.names = FALSE)
  write.csv(result$synteny$summary, file = "prepared_synteny_summary.csv", row.names = FALSE)
  # phylo对象可选保存为newick格式
  write.tree(result$phy, file = "prepared_tree.nwk")
  cat("数据准备完成，结果已保存为 data_preparation_result.rds, prepared_counts.csv, prepared_synteny_summary.csv, prepared_tree.nwk\n")
}

