library(dplyr)
df <- read.csv("C:/Users/14775/Downloads/df.csv")
data <- read.csv("C:/Users/14775/Downloads/N1000T100obs1_bpfiapf.csv")

final_df <- df %>%
  # The 'across' function applies the same operation to multiple columns
  # We select all columns except for 'date'
  mutate(across(-date, function(p) {
    # Calculate the difference of the logs: Δ log(p)
    diff_log_p <- c(NA, diff(log(p)))
    
    # Subtract the mean of these differences
    # na.rm = TRUE is important because the first element is NA
    y <- diff_log_p - mean(diff_log_p, na.rm = TRUE)
    
    return(y)
  }))

price_data <- df[, -1] # 选取除第一列外的所有列

# --- 步骤 1: 计算对数价格 ---
log_prices <- log(price_data)

# --- 步骤 2: 计算对数收益率 (一阶差分) ---
#    使用 lapply 将 diff 函数应用于每一列。结果的行数会比原始数据少1。
log_returns <- as.data.frame(lapply(log_prices, diff))

# --- 步骤 3: 减去均值 (De-meaning)，得到 y_t ---
#    scale() 函数非常适合这个任务：
#    center = TRUE: 从每列中减去该列的均值
#    scale = FALSE: 不对数据进行标准化（即不除以标准差）
y_t <- as.data.frame(scale(log_returns, center = TRUE, scale = FALSE))
#    给列重命名以便区分
colnames(y_t) <- paste0(colnames(df)[-1], "_y")


# --- 步骤 4 & 5: 平方后取对数，得到 w_t ---
#    这是最终可以输入到滤波器的格式
w_t <- log(y_t^2)
#    再次重命名
colnames(w_t) <- paste0(colnames(df)[-1], "_w")


# 3. 将处理好的数据与日期列合并到新的数据框中
#    由于计算差分后数据减少了一行，我们也需要去掉原始日期中的第一个日期
#    以确保数据框的行数能够对齐。

final_df <- cbind(date = df$date[-1], w_t)[,-1]

