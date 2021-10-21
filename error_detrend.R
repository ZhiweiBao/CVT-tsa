library(psych)
library(seasonal)

library(fpp2)

library(xts)

library(urca)
library(gridExtra)

library(seewave)

# data <- read.xlsx2("./data/CVT_data.xlsx", sheetIndex = 1, stringsAsFactors = FALSE)
data <- read.csv('./data/re_error.csv',sep=',',header=TRUE, row.names=1)
data <- as.data.frame(lapply(data,as.numeric))
time_index <- strptime("2018-01-31 23:45:00", "%Y-%m-%d %H:%M:%S") + 900 * 1:3552
# time_index <- strptime("2017-11-21 11:00:00", "%Y-%m-%d %H:%M:%S") + 86400 * 1:123
rownames(data) = time_index
data <- as.xts(data)


u_train <- data[,'A']


u_A <- data[,'A']
u_B <- data[,'B']
u_C <- data[,'C']

u1 <- autoplot(u_A) +
  xlab("") + ylab("UA") +
  labs(title="CVT三相电压相对误差") +
  theme(
    plot.title = element_text(hjust = 0.5))
u2 <- autoplot(u_B) +
  xlab("") + ylab("UB")
u3 <- autoplot(u_C) +
  xlab("时间") + ylab("UC")
grid.arrange(u1,u2,u3,nrow=3)

# u_A_train <- u_A[8257:11136]
# u_B_train <- u_B[8257:11136]
# u_C_train <- u_C[8257:11136]

u_A_train <- u_A[1001:2000]
u_B_train <- u_B[1001:2000]
u_C_train <- u_C[1001:2000]

u_A_test <- u_A[2001:3552]
u_B_test <- u_B[2001:3552]
u_C_test <- u_C[2001:3552]

time_index_forecast = time_index[2001:3552]
predict_step = length(time_index_forecast)

# lambda_A <- BoxCox.lambda(u_A_train, upper=100)
# lambda_B <- BoxCox.lambda(u_B_train, upper=10)
# lambda_C <- BoxCox.lambda(u_C_train, upper=10)


# lambda_A <- yj_trans(u_A_train)
# lambda_B <- BoxCox.lambda(u_B_train, upper=10)
# lambda_C <- BoxCox.lambda(u_C_train, upper=10)


# lambda_A <- log(u_A_train)
# lambda_B <- 0
# lambda_C <- 0
# 
# u_A_boxcox <- as.xts(data.frame(row.names = time_index[1:2000], values = BoxCox(u_A_train, lambda_A)))
# u_B_boxcox <- as.xts(data.frame(row.names = time_index[1:2000], values = BoxCox(u_B_train, lambda_B)))
# u_C_boxcox <- as.xts(data.frame(row.names = time_index[1:2000], values = BoxCox(u_C_train, lambda_C)))

# u1bc <- autoplot(u_A_boxcox) +
#   xlab("") + ylab("UA") +
#   labs(title="CVT三相电压幅值Box-Cox变换") +
#   theme(
#     plot.title = element_text(hjust = 0.5))
# u2bc <- autoplot(u_B_boxcox) +
#   xlab("") + ylab("UB")
# u3bc <- autoplot(u_C_boxcox) +
#   xlab("时间") + ylab("UC")
# grid.arrange(u1bc,u2bc,u3bc,nrow=3)

spec.pgram(ts_seasonal_A, plot = TRUE)
seewave::spec(ts_seasonal_A)

# ---------------model_A-------------------------------

u_A_mul <- msts(u_A_train, seasonal.periods=c(96))
model_A <- decompose(u_A_mul, type = "multiplicative")
model_A %>% autoplot() + ggtitle("经Box-Cox变换的A相电压时间序列MSTL分解图") + theme(plot.title = element_text(hjust = 0.5)) + xlab("周")
ts_seasonal_A <- seasonal(model_A)
ts_trend_A <- trendcycle(model_A)
ts_resid_A <- remainder(model_A)
ts_seasadj_A <- seasadj(model_A)

ts_resid_A_boxplot <- fivenum(ts_resid_A)
ts_resid_A_low <- ts_resid_A_boxplot[2] - 1.5 * (ts_resid_A_boxplot[4]- ts_resid_A_boxplot[2])
ts_resid_A_up <- ts_resid_A_boxplot[4] + 1.5 * (ts_resid_A_boxplot[4]- ts_resid_A_boxplot[2])

ts_resid_A_11 <- ts_resid_A[which(ts_resid_A > ts_resid_A_up),]

# checkresiduals(ts_resid_A)

# ---------------model_B-------------------------------

u_B_mul <- msts(u_B_train, seasonal.periods=c(96))
model_B <- decompose(u_B_mul, type = "multiplicative")
model_B %>% autoplot()
ts_seasonal_B <- seasonal(model_B)
ts_trend_B <- trendcycle(model_B)
ts_resid_B <- remainder(model_B)
ts_seasadj_B <- seasadj(model_B)

# ---------------model_C-------------------------------

u_C_mul <- msts(u_C_train, seasonal.periods=c(96))
model_C <- decompose(u_C_mul, type = "multiplicative")
model_C %>% autoplot()
ts_seasonal_C <- seasonal(model_C)
ts_trend_C <- trendcycle(model_C)
ts_resid_C <- remainder(model_C)
ts_seasadj_C <- seasadj(model_C)



ts_seasadj <- data.frame(ts_seasadj_A, ts_seasadj_B, ts_seasadj_C)
write.csv(ts_seasadj, './data/re_error_seasadj_1.csv',row.names = TRUE)







# ts_u_boxcox.mat <- data.frame(u_A_boxcox, u_B_boxcox, u_C_boxcox)
# VARselect(ts_u_boxcox.mat, lag.max = 60)
# model.var<-VAR(ts_u_boxcox.mat, p=2)

model_ts_test <- model_A
Ft <- max(0, 1-var(remainder(model_ts_test))/var(seasadj(model_ts_test)))
Fs1 <- max(0, 1-var(remainder(model_ts_test))/var(seasonal(model_ts_test)[,1] + remainder(model_ts_test)))
Fs2 <- max(0, 1-var(remainder(model_ts_test))/var(seasonal(model_ts_test)[,2] + remainder(model_ts_test)))
# Fs3 <- max(0, 1-var(remainder(model_ts_test))/var(seasonal(model_ts_test)[,3] + remainder(model_ts_test)))

urt.ts_seasadj_A <- ur.df(ts_seasadj_A)
urt.ts_seasadj_B <- ur.df(ts_seasadj_B)
urt.ts_seasadj_C <- ur.df(ts_seasadj_C)
summary(urt.ts_seasadj_A)

d_ts_seasadj_A <- diff(ts_seasadj_A)
d_ts_seasadj_B <- diff(ts_seasadj_B)
d_ts_seasadj_C <- diff(ts_seasadj_C)

urt.d_ts_seasadj_A <- ur.df(d_ts_seasadj_A)
urt.d_ts_seasadj_B <- ur.df(d_ts_seasadj_B)
urt.d_ts_seasadj_C <- ur.df(d_ts_seasadj_C)
summary(urt.d_ts_seasadj_A)


#-----------------------model-------------------------
ts_seasadj.mat <- data.frame(ts_seasadj_A, ts_seasadj_B, ts_seasadj_C)


lag_order = VARselect(ts_seasadj.mat, lag.max = 25)


model.vecm <- ca.jo(ts_seasadj.mat, ecdet="none", K=20)
jo.results <- summary(model.vecm)
# vecm.r2 <- cajorls(model.vecm, r=2)
# vecm.r2
# 
model.var <- vec2var(model.vecm,r=2)




# model.var <- VAR(ts_seasonal.mat, p=24, type = "both")

var.resid <- model.var$resid

# var.resid_A <- as.xts(data.frame(row.names = time_index[3:11136], values = var.resid[,1]))
# var.resid_B <- as.xts(data.frame(row.names = time_index[3:11136], values = var.resid[,2]))
# var.resid_C <- as.xts(data.frame(row.names = time_index[3:11136], values = var.resid[,3]))

checkresiduals(model.var$resid[,1])
# checkresiduals(model.var$resid[,2])
# checkresiduals(model.var$resid[,3])

#---------------------------predict----------------------

var.predict<-predict(model.var, n.ahead=predict_step, ci=0.95)

# fanchart(var.predict)

ts_seasadj_fcst <- var.predict$fcst
ts_seasadj_A_fcst <- ts_seasadj_fcst[1]
ts_seasadj_B_fcst <- ts_seasadj_fcst[2]
ts_seasadj_C_fcst <- ts_seasadj_fcst[3]

#----------------预测值-------------------

u_A_boxcox_fcst <- ts_seasadj_A_fcst[[1]][,1] + last(ts_seasonal_A[,1], predict_step/7) + last(ts_seasonal_A[,2], predict_step)
u_B_boxcox_fcst <- ts_seasadj_B_fcst[[1]][,1] + last(ts_seasonal_B[,1], predict_step/7) + last(ts_seasonal_B[,2], predict_step)
u_C_boxcox_fcst <- ts_seasadj_C_fcst[[1]][,1] + last(ts_seasonal_C[,1], predict_step/7) + last(ts_seasonal_C[,2], predict_step)

# u_A_fcst <- exp(u_A_boxcox_fcst)
# u_B_fcst <- exp(u_B_boxcox_fcst)
# u_C_fcst <- exp(u_C_boxcox_fcst)

u_A_fcst <- InvBoxCox(u_A_boxcox_fcst, lambda_A, biasadj=TRUE, fvar=ts_seasadj_A_fcst[[1]][,4]^2)
u_B_fcst <- InvBoxCox(u_B_boxcox_fcst, lambda_B, biasadj=TRUE, fvar=ts_seasadj_B_fcst[[1]][,4]^2)
u_C_fcst <- InvBoxCox(u_C_boxcox_fcst, lambda_C, biasadj=TRUE, fvar=ts_seasadj_C_fcst[[1]][,4]^2)


#-----------------上下限-------------------

u_A_boxcox_fcst_cpival <- ts_seasadj_A_fcst[[1]][,2:3] + last(ts_seasonal_A[,1], predict_step/7) + last(ts_seasonal_A[,2], predict_step)
u_B_boxcox_fcst_cpival <- ts_seasadj_B_fcst[[1]][,2:3] + last(ts_seasonal_B[,1], predict_step/7) + last(ts_seasonal_B[,2], predict_step)
u_C_boxcox_fcst_cpival <- ts_seasadj_C_fcst[[1]][,2:3] + last(ts_seasonal_C[,1], predict_step/7) + last(ts_seasonal_C[,2], predict_step)


u_A_fcst_cpival <- InvBoxCox(u_A_boxcox_fcst_cpival, lambda_A)
u_B_fcst_cpival <- InvBoxCox(u_B_boxcox_fcst_cpival, lambda_B)
u_C_fcst_cpival <- InvBoxCox(u_C_boxcox_fcst_cpival, lambda_C)


ts_u_A_fcst <- as.xts(data.frame(row.names = time_index_forecast, values = u_A_fcst, lower = u_A_fcst_cpival[,1], upper = u_A_fcst_cpival[,2]))
ts_u_B_fcst <- as.xts(data.frame(row.names = time_index_forecast, values = u_B_fcst, lower = u_B_fcst_cpival[,1], upper = u_B_fcst_cpival[,2]))
ts_u_C_fcst <- as.xts(data.frame(row.names = time_index_forecast, values = u_C_fcst, lower = u_C_fcst_cpival[,1], upper = u_C_fcst_cpival[,2]))


#------------------------原始值-------------------------------------------

ts_u_A <- u_A_train
ts_u_B <- u_B_train
ts_u_C <- u_C_train

#------------------------画图-------------------------------------

rb_u_A <- rbind.xts(ts_u_A, ts_u_A_fcst[,1])
rb_u_A_low <- rbind.xts(ts_u_A, ts_u_A_fcst[,2])
rb_u_A_up <- rbind.xts(ts_u_A, ts_u_A_fcst[,3])

rb_u_B <- rbind.xts(ts_u_B, ts_u_B_fcst[,1])
rb_u_B_low <- rbind.xts(ts_u_B, ts_u_B_fcst[,2])
rb_u_B_up <- rbind.xts(ts_u_B, ts_u_B_fcst[,3])

rb_u_C <- rbind.xts(ts_u_C, ts_u_C_fcst[,1])
rb_u_C_low <- rbind.xts(ts_u_C, ts_u_C_fcst[,2])
rb_u_C_up <- rbind.xts(ts_u_C, ts_u_C_fcst[,3])

p1 <- autoplot(rb_u_A) + 
  geom_ribbon(aes(ymin = rb_u_A_low, ymax = rb_u_A_up), fill = "slateblue") + 
  geom_line(aes(y = rb_u_A), color = "black") +
  xlab("") + ylab("UA") + 
  labs(title="CVT三相电压幅值") + 
  theme(
    # panel.background = element_rect(fill = "grey",colour = NA), 
    # panel.grid.minor = element_line(colour = "grey50"),
    # panel.grid.major = element_line(colour = "grey50"),
    # axis.line = element_line(colour = "black"),
    # panel.border = element_rect(fill = NA),
    # plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5)
  )

p2 <- autoplot(rb_u_B) + 
  geom_ribbon(aes(ymin = rb_u_A_low, ymax = rb_u_A_up), fill = "slateblue") + 
  geom_line(aes(y = rb_u_A), color = "black") +
  xlab("") + ylab("UB")

p3 <- autoplot(rb_u_C) + 
  geom_ribbon(aes(ymin = rb_u_A_low, ymax = rb_u_A_up), fill = "slateblue") + 
  geom_line(aes(y = rb_u_A), color = "black") +
  xlab("日期") + ylab("UC")

grid.arrange(p1,p2,p3,nrow=3)








