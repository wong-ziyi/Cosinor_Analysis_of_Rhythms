#parameters setting
period=23.5
TimeTagsCol=1
ValueCol=c(2:4)
xtitle="Time(h)"
ytitle="Expression Level"

#input the raw data
raw<-read.csv("results.csv")
#processes raw data (mean SD calculation)
temp<-data.frame(Time=raw[, TimeTagsCol], Value=rowMeans(raw[, ValueCol]), SD=rowSds(as.matrix(raw[, ValueCol])))
#fit raw data to cosinor curve
fit_per_est<-cosinor.lm(Value~time(Time), period = period, data = temp)
#get fitted model
FitCurve<-data.frame(x=temp$Time, y=fit_per_est$coefficients[1]+fit_per_est$coefficients[2]*cos(2*pi*temp$Time/period+2*fit_per_est$coefficients[3]))
FitCurve<-data.frame(x=temp$Time, y=fit_per_est$fit$fitted.values)
#rhythm detection
res<-cosinor.detect(fit_per_est)
#Scatter plot
ForScatter<-melt(raw, id.vars=TimeTagsCol)
ggplot()+
  geom_point(aes(x=Time, y=value), ForScatter)+
  geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x=xtitle, y=ytitle)+
  expand_limits(y = 0)
#Mean SD plot
ggplot()+
  geom_point(aes(x=Time, y=Value), temp)+
  geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), temp)+
  geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x=xtitle, y=ytitle)+
  expand_limits(y = 0)

cor(temp$Value, FitCurve$y)^2 #R^2
cor.test(temp$Value, FitCurve$y) #P-value
