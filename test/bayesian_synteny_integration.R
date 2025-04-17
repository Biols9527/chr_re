# bayesian_synteny_integration.R
# 综合多方法（迭代式+前序修正）结果的贝叶斯祖先状态重建示例脚本
# Author: GitHub Copilot
# Date: 2025-04-17

library(ape)
library(rjags)
library(coda)

# 读取数据
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("用法: Rscript bayesian_synteny_integration.R prepared_tree.nwk synteny_states.csv synteny_states_preorder_refined.csv\n")
  quit(status = 1)
}
tree_file <- args[1]
iterative_file <- args[2]
preorder_file <- args[3]

phy <- read.tree(tree_file)
iterative_df <- read.csv(iterative_file, stringsAsFactors = FALSE)
preorder_df <- read.csv(preorder_file, stringsAsFactors = FALSE)

# 整合两种方法的节点状态，计算均值和方差
all_nodes <- unique(c(iterative_df$node, preorder_df$node))
get_state <- function(df, node) {
  v <- df[df$node == node, 2]
  if (length(v) == 0) return(NA) else return(as.numeric(v[1]))
}
priors <- data.frame(
  node = all_nodes,
  mean = NA,
  sd = NA
)
for (i in seq_along(all_nodes)) {
  v1 <- get_state(iterative_df, all_nodes[i])
  v2 <- get_state(preorder_df, all_nodes[i])
  vals <- na.omit(c(v1, v2))
  priors$mean[i] <- mean(vals)
  priors$sd[i] <- ifelse(length(vals) > 1, sd(vals), 2) # 若只有一个值，给定较宽松的sd
}

# 判断节点类型
is_tip <- priors$node %in% phy$tip.label
node_type <- ifelse(is_tip, "tip", "internal")
species <- ifelse(is_tip, priors$node, "")

# 获取原始方法的状态
get_val <- function(df, node) {
  v <- df[df$node == node, 2]
  if (length(v) == 0) return(NA) else return(as.numeric(v[1]))
}
iterative_val <- sapply(priors$node, function(n) get_val(iterative_df, n))
preorder_val <- sapply(priors$node, function(n) get_val(preorder_df, n))

# 贝叶斯模型定义（每个节点状态为正态分布，先验为上述均值和方差）
# 这里只做简单示例：每个节点独立，实际可根据树结构扩展
model_string <- "model {
  for (i in 1:N) {
    state[i] ~ dnorm(prior_mean[i], prior_prec[i])
  }
}"

N <- nrow(priors)
data_jags <- list(
  N = N,
  prior_mean = priors$mean,
  prior_prec = 1/(priors$sd^2)
)

params <- c("state")
model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 2, quiet=TRUE)
update(model, 1000, progress.bar="none")
samples <- coda.samples(model, variable.names = params, n.iter = 5000, progress.bar="none")

# 提取后验均值和置信区间
posterior <- summary(samples)$statistics
posterior_ci <- summary(samples)$quantiles
result <- data.frame(
  node = priors$node,
  node_type = node_type,
  species = species,
  iterative_state = iterative_val,
  preorder_state = preorder_val,
  prior_mean = priors$mean,
  prior_sd = priors$sd,
  post_mean = posterior[,"Mean"],
  post_sd = posterior[,"SD"],
  ci2.5 = posterior_ci[,"2.5%"],
  ci97.5 = posterior_ci[,"97.5%"]
)
write.csv(result, file = "bayesian_synteny_states.csv", row.names = FALSE)
cat("贝叶斯整合分析完成，结果已保存为 bayesian_synteny_states.csv\n")

