# methods/iterative_synteny_reconstruction.R
# 迭代式共线性感知祖先状态与事件推断
# Author: GitHub Copilot
# Date: 2025-04-17

library(ape)
library(dplyr)

# 1. 定义事件推断规则（可根据实际需求扩展）
infer_event <- function(parent_state, child_state, synteny_info = NULL) {
  # 示例规则：可根据 parent_state, child_state, synteny_info 判断事件类型
  if (!is.null(synteny_info) && synteny_info == "WGD") {
    return("WGD")
  }
  if (parent_state > child_state) {
    return("fusion")
  }
  if (parent_state < child_state) {
    return("fission")
  }
  return("none")
}

# 2. 基线状态估计（可用共线性信息优化）
estimate_initial_states <- function(phy, counts, synteny_summary) {
  # 这里简单用众数或中位数作为初始估计，可结合synteny_summary优化
  node_states <- rep(median(counts$chromosome_number, na.rm = TRUE), Nnode(phy) + Ntip(phy))
  names(node_states) <- c(phy$tip.label, if (!is.null(phy$node.label)) phy$node.label else as.character((Ntip(phy)+1):(Ntip(phy)+Nnode(phy))))
  # 用已知叶节点数据覆盖初始状态
  for (tip in phy$tip.label) {
    val <- counts$chromosome_number[counts$species == tip]
    if (length(val) > 0 && !is.na(val[1])) node_states[tip] <- val[1]
  }
  return(node_states)
}

# 3. 迭代推断主函数（仅用postorder遍历）
iterative_synteny_reconstruction <- function(phy, counts, synteny, max_iter = 10) {
  if (!"chromosome_number" %in% colnames(counts)) {
    stop("counts数据框必须包含chromosome_number列")
  }
  node_states <- estimate_initial_states(phy, counts, synteny$summary)
  event_list <- list()
  for (iter in 1:max_iter) {
    prev_states <- node_states
    # 自下而上遍历
    postorder <- rev(ape::postorder(phy))
    for (node in postorder) {
      if (node > Ntip(phy)) {
        children <- phy$edge[phy$edge[,1] == node, 2]
        child_states <- node_states[children]
        # 可结合synteny$summary信息
        node_states[node] <- round(mean(child_states))
        # 推断事件
        for (child in children) {
          event <- infer_event(node_states[node], node_states[child])
          event_list[[paste(node, child, sep = "_")]] <- event
        }
      }
    }
    if (all(node_states == prev_states)) break
  }
  return(list(states = node_states, events = event_list))
}

# ---- 命令行运行支持 ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("用法: Rscript methods/iterative_synteny_reconstruction.R data_preparation_result.rds\n")
    quit(status = 1)
  }
  prep <- readRDS(args[1])
  result <- iterative_synteny_reconstruction(
    phy = prep$phy,
    counts = prep$counts,
    synteny = prep$synteny
  )
  saveRDS(result, file = "synteny_reconstruction_result.rds")
  # 也可输出为csv
  states_df <- data.frame(node = names(result$states), state = result$states)
  write.csv(states_df, file = "synteny_states.csv", row.names = FALSE)
  events_df <- data.frame(branch = names(result$events), event = unlist(result$events))
  write.csv(events_df, file = "synteny_events.csv", row.names = FALSE)
  cat("重建完成，结果已保存为 synteny_reconstruction_result.rds, synteny_states.csv, synteny_events.csv\n")
}
