# preorder_synteny_refinement.R
# 基于自上而下（前序遍历）的染色体数目状态修正脚本（修正版：叶节点保持观测值不变）
# Author: GitHub Copilot
# Date: 2025-04-17

library(ape)

# 自上而下递归修正函数（只修正内部节点，叶节点保持观测值）
do_preorder_refinement <- function(phy, node_states, refine_fun = NULL) {
  if (is.null(refine_fun)) {
    refine_fun <- function(parent_state, child_state) {
      if (abs(child_state - parent_state) > 5) {
        return(round((child_state + parent_state) / 2))
      } else {
        return(child_state)
      }
    }
  }
  preorder_rec <- function(node, node_states) {
    children <- phy$edge[phy$edge[,1] == node, 2]
    for (child in children) {
      # 只修正内部节点，叶节点保持观测值不变
      if (!(child %in% 1:Ntip(phy))) {
        node_states[child] <- refine_fun(node_states[node], node_states[child])
        node_states <- preorder_rec(child, node_states)
      }
    }
    return(node_states)
  }
  root_node <- Ntip(phy) + 1
  node_states <- preorder_rec(root_node, node_states)
  return(node_states)
}

# ---- 命令行运行支持 ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("用法: Rscript preorder_synteny_refinement.R prepared_tree.nwk synteny_states.csv\n")
    quit(status = 1)
  }
  phy <- read.tree(args[1])
  states_df <- read.csv(args[2], stringsAsFactors = FALSE)
  node_states <- states_df$state
  names(node_states) <- states_df$node
  refined_states <- do_preorder_refinement(phy, node_states)
  out_df <- data.frame(node = names(refined_states), refined_state = refined_states)
  write.csv(out_df, file = "synteny_states_preorder_refined.csv", row.names = FALSE)
  cat("自上而下修正完成，结果已保存为 synteny_states_preorder_refined.csv\n")
}

