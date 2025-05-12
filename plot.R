library(tidyverse)
library(readr)
library(stringr)
library(grid)
library(gridExtra)

file_list <- list.files(path = ".", pattern = "rolling_N200T1000.*\\.csv", full.names = TRUE)

extract_info <- function(filename) {
  dimension <- as.numeric(str_extract(filename, "(?<=_d)\\d+"))
  lag <- as.numeric(str_extract(filename, "(?<=_l)\\d+")) 
  list(d = dimension, l = lag)
}

df <- tibble()  # 初始化空 tibble


# 计算 RMSE 的函数
compute_rmse <- function(df) {
  rmse <- sqrt(colMeans((df - 1)^2))
  return(rmse)
}

# 批量读取、提取信息、计算 RMSE
results <- map_dfr(file_list, function(file) {
  df <- read_csv(file)[,2]
  info <- extract_info(file)
  rmse <- compute_rmse(df)
  tibble(
    dimension = info$d,
    lag = info$l,
    rmse = rmse
  )
})

for (file in file_list) {
  
  raw_vec <- read_csv(file, col_names = FALSE)[[2]]  # 每行是一个 replicate，读取单列
  
  info <- extract_info(file)
  dimension <- info$d
  lag <- info$l
  
  df_long <- tibble(
    replicate = seq_along(raw_vec),
    estimate = raw_vec,
    dimension = dimension,
    lag = lag,
    algorithm = "ORC-SMC",
    n_particles = 200
  )
  
  df <- bind_rows(df, df_long)
}



df_clean %>%
  group_by(dimension, lag) %>%
  summarise(rmse = sqrt(mean((estimate - 1)^2))) %>%
  ggplot(aes(x = factor(lag), y = rmse)) +
  geom_boxplot(aes(group = lag)) +
  facet_wrap(~ dimension, scales = "free_y") +
  labs(x = "Lag", y = "RMSE", title = "RMSE vs Lag across Dimensions") +
  theme_minimal()

ggplot(df_clean, aes(x = factor(lag), y = estimate)) +
  geom_boxplot() +
  facet_wrap(~ dimension, scales = "free_y") +
  labs(x = "Lag", y = "Estimate", title = "Estimate Distribution by Lag and Dimension") +
  theme_minimal()

library(ggplot2)

df_clean %>%
  ggplot(aes(x = factor(lag), y = estimate)) +
  geom_boxplot(aes(fill = factor(lag)), width = 0.6, outlier.size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "red") +  # 红叉表示平均值
  facet_wrap(~ dimension, scales = "free_y") +
  labs(x = "Lag", y = "Estimate", title = "Estimate Distributions by Dimension and Lag", fill = "Lag") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"  # 如果不需要颜色图例可以移除
  )



ggplot(results, aes(x = factor(lag), y = rmse, colour = dimension, group = dimension)) +
  geom_boxplot(outlier.shape = NA) +                  
  geom_jitter(width = 0.2, alpha = 0.4, color = "blue") +  
  #facet_wrap(~ dimension, scales = "free_y") +        
  labs(
    title = "RMSE over Lags by Dimension",
    x = "Lag",
    y = "RMSE"
  ) +
  theme_minimal(base_size = 14)

results$lag <- as.numeric(results$lag)
df_clean$lag <- as.numeric(df_clean$lag)

# 绘图
ggplot(results, aes(x = lag, y = rmse, color = factor(dimension), group = dimension)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  scale_x_continuous(breaks = c(4, 8, 16, 32, 64)) +
  labs(
    x = "Lag",
    y = expression("RMSE of " * hat(Z)[T]),
    color = "Dimension"
  ) +
  ggtitle("N200T1000d4")+
  theme_minimal()

df_clean <- df %>%
  filter(estimate != "x") %>%                  # 去掉非法数据
  mutate(
    estimate = as.numeric(estimate),           # 转换为 numeric
    dimension = factor(dimension),             # 因子化便于分面
    lag = factor(lag)                          # 因子化便于分组
  )

# 画图：大格子是 dimension，小格子里每个 lag 一组
ggplot(df_clean, aes(x = lag, y = estimate, fill = factor(lag))) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "red") +
  facet_wrap(~ dimension, scales = "free_y") +
  labs(
    x = "Lag", y = "Estimate",
    title = "Estimate distribution per Lag across Dimensions"
  ) +
  coord_cartesian(ylim = c(0, 3)) +  # 控制 Y 轴显示范围（根据你的数据调整）
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

df_clean %>%
  ggplot(aes(x = factor(lag), y = estimate, fill = factor(dimension))) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +  # 不显示离群值
  coord_cartesian(ylim = c(0, 2)) +  # 控制 y 轴显示范围，但不会删除数据
  labs(x = "Lag", y = "Estimate", fill = "Dimension") +
  theme_minimal() +
  ggtitle("Boxplot of Estimates for Different Dimensions and Lags")


ggplot(df_clean, aes(x = lag, y = rmse, color = factor(dimension), group = dimension)) +
  geom_line(linewidth = 1) +         # 连线
  geom_point(size = 2) +             # 点
  labs(
    x = "Lag", 
    y = "RMSE", 
    color = "Dimension", 
    title = "RMSE vs Lag for Different Dimensions"
  ) +
  theme_minimal()

#filtering 
library(ggplot2)
library(gridExtra)


j     <- 2
Times <- c(1, floor(Time/2), Time)

plots <- lapply(Times, function(t0){
  # 1) μ_true, σ²_true
  μ_true    <- smooth_means[t0, j]
  σ2_true   <- smooth_covs[j, j, t0]
  xs        <- seq(μ_true - 4*sqrt(σ2_true),
                   μ_true + 4*sqrt(σ2_true),
                   length.out = 200)
  
  df_true   <- data.frame(x = xs,
                          dens = dnorm(xs, mean = μ_true, sd = sqrt(σ2_true)),
                          type = "KF smoother")
  
  # 2) 
  particles_t0_j <- smooth_particles[, t0, j]
  
  # 2) 
  dens_est <- density(particles_t0_j)
  
  # 3) 
  df_part <- data.frame(
    x   = dens_est$x,
    dens= dens_est$y,
    type= "Orc-SMC"
  )
  
  # 3) 
  ggplot() +
   
    geom_line(data = df_part, aes(x = x, y = dens, color = type), size = 1) +
    
    geom_line(data = df_true, aes(x = x, y = dens, color = type), size = 1) +
    labs(
      title = paste0("Smoothing marginal, t=", t0, ", coor=", j),
      x     = bquote(x[.(t0)]), 
      y     = "density"
    ) +
    scale_color_manual(
      "", 
      values = c("Orc-SMC" = "steelblue", "KF smoother" = "red")
    ) +
    theme_minimal()
})


grid.arrange(
  grobs = plots,
  ncol = length(plots),
  top = textGrob(
    "N500T1000, d=2, lag=16",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

