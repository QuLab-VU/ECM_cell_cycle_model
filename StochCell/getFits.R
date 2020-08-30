#load libraries
library(emg)
library(fracprolif)
library(xlsx)

setwd("C:\\Users\\Grim.m\\Documents\\pc9_cell_runs")
state=read.xlsx("tmp_state.xls", sheetName = "Sheet1")
total=read.xlsx("tmp_total.xls", sheetName = "Sheet1")

mitotic.lifespans <- subset(state, !End_of_Expt &
 !Death &
 !is.na(Lifespan))$Lifespan

censored.lifespans <- subset(state, End_of_Expt &
 !Death &
 !is.na(Lifespan))$Lifespan

est <- q.mle.emg.estimate(mitotic.lifespans, censored.lifespans)


r <- q.rates("emg", est)
m <- lm(AKTERK_max ~ Time_h, total)
untreated.division.rate <- coef(m)['Time_h']*log(2)
half.imt <- log(2)/untreated.division.rate/2
y <- total$AKTERK_treat - total$AKTERK_treat[1]
t <- total$Time_h
fit <- nls(y ~ log(split_tot(t, 1, 0,
 untreated.division.rate, r['d'],
 0.0, r['q'],
 0.0, 0.0, a, a,
 half.imt), 2),
 start=list(a=0.003),
 lower=list(a=0),
 algorithm="port")
a <- coef(fit)
dev.new(width=7, height=3)
par(mfrow=c(1,3))
d <- subset(state, !Death & !is.na(Lifespan))

plot(d$Birthtime,
 d$Lifespan,
 xlab="Birth time (h)",
 ylab="Lifespan (h)")

n <- length(d$Lifespan)
text(20,80,substitute("n" == n,list(n=n)))
p <- 100*length(subset(d, End_of_Expt)$Lifespan) / n
text(15,85,substitute(p * "% EoE", list(p=round(p))))

hist(mitotic.lifespans,
 xlab="Observed Intermitotic Times (h)",
 breaks=20,
 freq=FALSE)
curve(demg(x, coef(est)['mu'],
 coef(est)['sigma'],
 coef(est)['lambda']),
 add=TRUE, col='red', lwd=2)
text(35, 0.07, pos=4, labels=substitute("n"==n, list(n=length(mitotic.lifespans))))
text(35, 0.065, pos=4, labels=substitute(mu==m, list(m=round(coef(est)['mu'],2))))
text(35, 0.06, pos=4, labels=substitute(sigma==s, list(s=round(coef(est)['sigma'],2))))
text(35, 0.055, pos=4, labels=substitute(kappa==l, list(l=round(1/coef(est)['lambda'],2))))
ks <- ks.test(mitotic.lifespans, "pemg",
 mu=coef(est)['mu'],
 sigma=coef(est)['sigma'],
 lambda=coef(est)['lambda'])
text(35, 0.05, pos=4, labels=substitute("ks p-val"==p, list(p=round(ks$p.value,2))))
###################################
plot (t, y,
 xlab="Time (h)",
 ylab="Population Doublings")
model.total <- log(split_tot(t, 1, 0, untreated.division.rate,
 r['d'], 0.0, r['q'], 0, 0,
 a['a'], a['a'], half.imt), 2)
lines(t,model.total,col="green")
lines(t,model.total*split_fq(t, 1, 0, untreated.division.rate,
 r['d'], 0.0, r['q'], 0, 0,
 a['a'], a['a'], half.imt), col="red")
lines(t,model.total*split_fd(t, 1, 0, untreated.division.rate,
 r['d'], 0.0, r['q'], 0, 0,
 a['a'], a['a'], half.imt), col="blue")
dev.copy(pdf,'tmpplot.pdf')
dev.off()
save.image(file = "tmp_r.RData")
