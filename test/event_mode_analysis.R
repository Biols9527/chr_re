# methods/event_mode_analysis.R
# 基于贝叶斯整合结果的事件识别与统计分析脚本
# Author: GitHub Copilot
# Date: 2025-04-17

library(ape)
library(dplyr)
library(ggplot2)

# 事件推断规则（优化：对状态取整，允许自定义阈值）
infer_event <- function(parent_state, child_state, threshold = 0) {
  if (is.na(parent_state) || is.na(child_state)) return(NA)
  p <- round(parent_state)
  c <- round(child_state)
  if (abs(p - c) <= threshold) return("none")
  if (p > c) return("fusion")
  if (p < c) return("fission")
  return("none")
}

# 通用事件分析函数（支持不同来源的节点状态）
event_mode_analysis_states <- function(states_df, phy, label, threshold = 0) {
  # 自动识别状态列
  state_col <- NULL
  for (col in c("state", "refined_state", "post_mean")) {
    if (col %in% colnames(states_df)) {
      state_col <- col
      break
    }
  }
  if (is.null(state_col)) stop("输入文件中未找到 state/refined_state/post_mean 等节点状态列！")
  node_states <- states_df[[state_col]]
  names(node_states) <- states_df$node
  edge_df <- data.frame(
    branch = paste(phy$edge[,1], phy$edge[,2], sep = "_"),
    parent = phy$edge[,1],
    child = phy$edge[,2],
    parent_label = NA,
    child_label = NA,
    parent_state = NA,
    child_state = NA,
    parent_state_rounded = NA,
    child_state_rounded = NA,
    event = NA
  )
  all_labels <- c(phy$tip.label, if (!is.null(phy$node.label)) phy$node.label else as.character((Ntip(phy)+1):(Ntip(phy)+Nnode(phy))))
  for (i in seq_len(nrow(edge_df))) {
    parent_idx <- edge_df$parent[i]
    child_idx <- edge_df$child[i]
    parent_label <- all_labels[parent_idx]
    child_label <- all_labels[child_idx]
    edge_df$parent_label[i] <- parent_label
    edge_df$child_label[i] <- child_label
    edge_df$parent_state[i] <- node_states[parent_label]
    edge_df$child_state[i] <- node_states[child_label]
    edge_df$parent_state_rounded[i] <- round(node_states[parent_label])
    edge_df$child_state_rounded[i] <- round(node_states[child_label])
    edge_df$event[i] <- infer_event(node_states[parent_label], node_states[child_label], threshold = threshold)
  }
  event_stats <- edge_df %>%
    group_by(event) %>%
    summarize(count = n()) %>%
    arrange(desc(count))
  write.csv(event_stats, file = paste0(label, "_event_type_stats.csv"), row.names = FALSE)
  write.csv(edge_df, file = paste0(label, "_event_branch_mapping.csv"), row.names = FALSE)
  p <- ggplot(event_stats, aes(x = reorder(event, -count), y = count, fill = event)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste0(label, " Event Type Frequency"), x = "Event Type", y = "Count") +
    theme(legend.position = "none")
  ggsave(paste0(label, "_event_type_frequency_barplot.pdf"), p, width = 6, height = 4)
}

# 主分析函数，分别对三种推断结果进行事件识别和统计
event_mode_analysis_all <- function(bayes_file, tree_file, iterative_file, preorder_file, threshold = 0) {
  phy <- read.tree(tree_file)
  # 贝叶斯后验均值
  bayes_df <- read.csv(bayes_file, stringsAsFactors = FALSE)
  bayes_states <- data.frame(node = bayes_df$node, state = bayes_df$post_mean)
  event_mode_analysis_states(bayes_states, phy, label = "bayesian", threshold = threshold)
  # 迭代式
  iterative_df <- read.csv(iterative_file, stringsAsFactors = FALSE)
  event_mode_analysis_states(iterative_df, phy, label = "iterative", threshold = threshold)
  # 前序修正
  preorder_df <- read.csv(preorder_file, stringsAsFactors = FALSE)
  event_mode_analysis_states(preorder_df, phy, label = "preorder", threshold = threshold)
  cat("三种推断结果的事件统计与可视化已全部完成。\n")
}

# ---- 命令行运行支持 ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) {
    cat("用法: Rscript methods/event_mode_analysis.R bayesian_synteny_states.csv prepared_tree.nwk synteny_states.csv synteny_states_preorder_refined.csv [threshold]\n")
    quit(status = 1)
  }
  bayes_file <- args[1]
  tree_file <- args[2]
  iterative_file <- args[3]
  preorder_file <- args[4]
  threshold <- ifelse(length(args) >= 5, as.numeric(args[5]), 0)
  event_mode_analysis_all(bayes_file, tree_file, iterative_file, preorder_file, threshold)
}

