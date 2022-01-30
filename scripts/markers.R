library(tidyverse)
library(spdep)
library(spatialreg)
library(pheatmap)
library(Hmisc)
DATA_PATH = 'data/nbhd_coord_schurch_2020/'
df = read_csv(paste0(DATA_PATH,'CRC_master.csv'))

d = df %>% select(X,Y,`CD44 - stroma`:`MMP12 - matrix metalloproteinase`,spots) %>%
  filter(spots == '15_A') %>%
  select(-spots)

MIN_DIST = 25

xy = d %>%
  dplyr::select(c(X,Y)) %>%
  rename(x = X,y=Y) %>%
  as.data.frame

set.seed(2022); xy = spatialRF::thinning(xy,minimum.distance = MIN_DIST)

d = d %>%
  semi_join(.,xy,by=c("X"="x","Y"="y")) %>%
  mutate(obs = 1:nrow(.)) %>%
  select(-c(X,Y)) %>%
  rename_with(~sub(" .*", "", .x)) %>%
  rename_with(~sub("-",".",.x))


col.rel.nb <- graph2nb(relativeneigh(as.matrix(xy)), sym=TRUE)

C = spdep::nb2mat(col.rel.nb,style="B",zero.policy=T)
listw = mat2listw(C)

pheatmap(cor(d))

m.1 <- errorsarlm(p53 ~ Cytokeratin,data = d,listw)

m.1$residuals

moran.test(m.1$residuals,listw)

moran.test(d$p53,listw)

m.2 <- lm(p53 ~ Cytokeratin, data=d)

m.3 <- errorsarlm(log(p53+1) ~ Cytokeratin,data = d,listw)
plot(m.3)

plot_resids <- function(m,listw) {
  resid <- m$residuals
  moran.test(resid,listw) %>% print
  hist(resid)
  resid = (resid - min(resid)) / (max(resid) - min(resid))
  ggplot(tibble(resid=resid,x=xy$x,y=xy$y)) +
    geom_point(aes(x=x,y=y,color=resid)) +
    scale_colour_gradient(low = "red", high = "blue")
}

plot_resids(m.3,listw)

out <- redun(~MMP12+CD44+CD45+T.bet+EGFR+CD30+MUC.1,data=d,nk=3)

out


m.4 <- errorsarlm(p53 ~ MMP12+CD44+CD45+T.bet+EGFR+CD30+MUC.1,data=d,listw)
start.time <- Sys.time()
m.5 <- errorsarlm(p53 ~ MMP12+CD44+CD45+T.bet+EGFR+CD30+MUC.1,data=d,listw=listw,Durbin=TRUE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
m.6 <-lm(p53 ~ MMP12+CD44+CD45+T.bet+EGFR+CD30+MUC.1,data=d)
plot_resids(m.4,listw)

plot_resids(m.5,listw)
plot_resids(m.6,listw)
summary(m.4)
summary(m.5)
summary(m.6)

m.7 <- lm(p53 ~ .,data=d)

out <- transace(d %>%select(-p53) %>% as.matrix)
d.new = cbind(out,p53=d$p53)

summary(m.7)

m.8 <- lm(p53 ~ .,data=d.new %>% as.data.frame)
anova(m.4,m.5)

summary(m.8)

m.9 <- errorsarlm(p53 ~ .,data=d,listw=listw,Durbin=TRUE)
m.10 <- errorsarlm(p53 ~ .,data=d.new %>% as.data.frame,listw=listw,Durbin=TRUE)

summary(m.9)

AIC(m.7)
AIC(m.8)
AIC(m.9)
AIC(m.10)

summary(m.8)

plot(m.8)
summary(m.10)


library(tidymodels)
