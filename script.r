library(latticeExtra)
library(dplyr)
library(tidyr)
library(lattice)
library(ggplot2)
library(bayesplot)
library(ammiBayes)
library(coda)
library(gridExtra)
library(ggpubr)

#
dados <- read.csv("data.csv", h=T)
head(dados)
dim(dados)
#
envy <- with(dados, paste(rep(1:8, each=18), rep(c("Y1", "Y2"), each=864/2), sep=""))
envy
geny <- with(dados, paste("G", Gen, sep=""))
geny

#
dados <- transform(dados, Env=factor(envy), Rep=factor(Rep), Gen=factor(geny))
#
dados <- transform(dados, Env=factor(Env), Rep=factor(Rep), Gen=factor(Gen))
#
#
model <- with(dados, ammiBayes(Y=GY, Gen=Gen, Env=Env, Rep=Rep, iter=4000,burn=500, jump=5, chains=3))
#
#
#save.image("model.RData")


#load("model.RData")



data.hist <- transform(dados, Gen=factor(Gen, levels=unique(geny)), 
											 Env=factor(Env, levels=unique(envy)))


# Geral

medias.cmc <- aggregate(GY ~ Gen, data=data.hist, FUN=mean)

medias.text <- round(medias.cmc$GY,2)


texto.hist <- lapply(medias.text, function(x){
	eval(parse(text = paste("expression(bar(x)==", x, ")", sep = "")))
})

texto.hist

textx <- round(do.breaks(range(data.hist$GY),8), 1)
textx 

histogram(~GY|Gen, data=data.hist,
			xlab="Yield", layout=c(3,6), main="Year 1",
			endpoints=range(data.hist$GY),			
			nint=8,
		  scales=list(alternating=1, x=list(at=textx)), ylim=c(0,70),
			medias.text=medias.text,
			panel=function(medias.text=medias.text,...){
				panel.histogram(...)
				panel.text(x=medias.text[panel.number()], y=56,
									 labels=texto.hist[[panel.number()]])
				panel.segments(x0=medias.text[panel.number()], y0=0,
											 x1=medias.text[panel.number()], y1=54, lty=2)
				
				})


# Year 1

medias.cmc <- aggregate(GY ~ Gen, data=data.hist, subset=Env%in%unique(envy)[1:8], FUN=mean)

medias.text <- round(medias.cmc$GY,2)


texto.hist <- lapply(medias.text, function(x){
	eval(parse(text = paste("expression(bar(x)==", x, ")", sep = "")))
})



textx <- round(do.breaks(range(data.hist$GY),8), 1)


fig1a <- histogram(~GY|Gen, data=data.hist, subset=Env%in%unique(envy)[1:8],
			xlab="Yield", layout=c(3,6), main="Year 1",
			endpoints=range(data.hist$GY),			
			nint=8,
		  scales=list(alternating=1, x=list(at=textx)), ylim=c(0,70),
			medias.text=medias.text,
			panel=function(medias.text=medias.text,...){
				panel.histogram(...)
				panel.text(x=medias.text[panel.number()], y=56,
									 labels=texto.hist[[panel.number()]])
				panel.segments(x0=medias.text[panel.number()], y0=0,
											 x1=medias.text[panel.number()], y1=54, lty=2)
				
				})


# Year 2


medias.cmc2 <- aggregate(GY ~ Gen, data=data.hist, subset=Env%in%unique(envy)[9:16], FUN=mean)

medias.text2 <- round(medias.cmc2$GY,2)


texto.hist2 <- lapply(medias.text2, function(x){
	eval(parse(text = paste("expression(bar(x)==", x, ")", sep = "")))
})



textx2 <- round(do.breaks(range(data.hist$GY),8), 1)



fig1b <- histogram(~GY|Gen, data=data.hist, subset=Env%in%unique(envy)[9:18],
			xlab="Yield", layout=c(3,6), main="Year 2",
			endpoints=range(data.hist$GY),			
nint=8,
		  scales=list(alternating=1, x=list(at=textx2)), ylim=c(0,70),
			medias.text=medias.text2,
			panel=function(medias.text=medias.text2,...){
				panel.histogram(...)
				panel.text(x=medias.text2[panel.number()], y=56,
									 labels=texto.hist2[[panel.number()]])
				panel.segments(x0=medias.text2[panel.number()], y0=0,
											 x1=medias.text2[panel.number()], y1=54, lty=2)
				
				})


pdf("fig1.pdf", width=20, height=11)
grid.arrange(fig1a, fig1b, ncol=2)
dev.off()

sink('summary.txt')
summary(model)
sink()

# Diagnosis

diag.gen <- diagnosis.ammiBayes(model, pars="Genotype")

diag.env1 <- diagnosis.ammiBayes(model, pars="Env.PC1")

diag.env2 <- diagnosis.ammiBayes(model, pars="Env.PC2")


gelman.diag(diag.gen)

raftery.diag(diag.gen)

# Bayesplot 

color_scheme_set("viridis")

pdf("tracegen.pdf", width=14)
mcmc_trace(diag.gen, facet_args=list(nrow=3)) + theme_minimal()
dev.off()

pdf("traceenv.pdf", width=14)
mcmc_trace(diag.env1, facet_args=list(nrow=4)) + theme_minimal()
dev.off()


mcmc_areas(diag.gen, pars=rownames(gen.ef)) + theme_minimal()

mcmc_hist(diag.gen)

# Prediction

pred <- predict(model)

dim(pred)
dim(dados)


head(dados)

dat.post <- data.frame(dados, pred)
head(dat.post)

post.mean <- aggregate(cbind(Mean, X2.5., X97.5.) ~ Env + Gen, data=dat.post, FUN=mean)
post.mean

length(unique(envy))
length(unique(geny))

16*18

write.csv(post.mean, "posterior_mean.csv")


post.mean.gen <- aggregate(cbind(Mean, X2.5., X97.5.) ~ Gen, data=dat.post, FUN=mean)
arrange(post.mean.gen, Gen)



# Genotype effect



gen.ef <- gen.effects(model)
gen.ef

x11()

pred1 <- ammiBayes.gen.plot(gen.ef)

pred2 <- mcmc_areas(diag.gen, pars=rownames(gen.ef), prob_outer=0.95, prob=0.5) + theme_minimal()


pdf("pred.pdf", width=15)
ggarrange(pred1, pred2, ncol=2, labels=c("A","B"))
#grid.arrange(pred1,pred2, ncol=2)
dev.off()


mcmc_intervals(diag.gen)

ammiBayes.mean.plot(model, pars.env=unique(envy)[1:8])

ammiBayes.mean.plot(model, col.text.env="transparent")

ammiBayes.mean.plot(model, col.text.gen="transparent")
x11()
ammiBayes.mean.plot(model)


inter <- detect.interaction(model, type="gen", conf=0.95)
inter

inter2 <- detect.interaction(model, type="env", conf=0.95)
inter2



ammiBayes.conf.plot(model, change.signal=TRUE, change.scale=0.1)

pdf("confplot.pdf", width=14, height=12)
ammiBayes.conf.plot(model, change.signal=TRUE, change.scale=0.1, 
										pars.gen=c("G.4", "G.3", "G.13","G.2", "G.12"),
										xlab="PC1 (81.14%)", ylab="PC2 (14.48%)")
dev.off()



y1 <- ammiBayes.conf.plot(model, pars.env=unique(envy)[1:8], change.signal=TRUE, change.scale=0.1)

y2 <- ammiBayes.conf.plot(model, pars.env=unique(envy)[9:16], change.signal=TRUE, change.scale=0.1)


pdf("confy.pdf", width=15)
ggarrange(y1, y2, ncol=2, labels=c("A","B"))
#grid.arrange(pred1,pred2, ncol=2)
dev.off()

citation("ammiBayes")

citation()


pdf("biplot_f.pdf", width=14, height=12)
ammiBayes.conf.plot(model, change.signal=TRUE, change.scale=0.1, xlab="PC1 (81.14%)", ylab="PC2 (14.48%)")
dev.off()


##### metan

library(metan)
head(dados)

summary(model)

pdf("boxplot.pdf", width=14, height=12)
bwplot(GY ~ Gen|Env, data=dados, scales=list(alternating=1), layout=c(4,4))

bwplot(GY ~ Gen|Rep, data=dados, scales=list(alternating=1), subset=Gen=="G.2")

bwplot(GY ~ Gen|Env, data=dados, scales=list(alternating=1), subset=Gen=="G.2")
dev.off()

grid.arrange(gr1, gr2, ncol=2)

ammi_model <- performs_ammi(dados, Env, Gen, Rep, resp = c(GY))
ammi_model

plot_scores(ammi_model, type = 2)

predicted <- predict(ammi_model)
predicted
predicted %>% 
   subset(TRAIT == "GY") %>% 
   make_mat(GEN, ENV, YpredAMMI) %>% 
  round_cols()


predg <- aggregate(YpredAMMI ~ GEN, data=predicted, FUN=mean)
head(predg)
predg[order(predg$V1),]
