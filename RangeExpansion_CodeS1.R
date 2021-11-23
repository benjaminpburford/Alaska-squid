
# Supporting code for:

# Rapid range expansion of a marine ectotherm
# reveals the demographic and ecological consequences of
# short-term variability in seawater temperature and dissolved oxygen



library(tidyverse)
library(gsw)
library(rMR)
library(respirometry)
library(readxl)
library(mgcv)
library(lsmeans)
library(corrplot)
library(parameters)



  ######### Load data #########
# Set working directory to folder with data file only
setwd()
# Make list of the data file
file.list <- list.files()

# CODE MUST BE RUN SEQUENTIALLY



  ######### Observations of poleward range expansion and migratory capacity ##########

## Map of GOA and CCS with squid observations (Figure 3A)
#----

# Load map data
alaska <- map_data("world")

# Make squid observations (from Figure 2B) a dataframe
lat <- c(48.227,55.21373,57.04528,58.26137)	
lon <- c(-124.774,-133.6514,-135.3285,-154.048)
squid_obs <- data.frame(cbind(lat,lon))
names(squid_obs) = c("lat","lon")

# Define map
goaccsmap <- ggplot(data = alaska, aes(x=long, y=lat, group = group)) +
  geom_polygon(fill = "black", color="black") + 
  coord_fixed(xlim = c(-160, -108),  ylim = c(24, 61), ratio = 1.2) +
  theme_classic(base_size = 12) +
  xlab("Longitude (°W)")+
  ylab("Latitude (°N)")

# Plot and add points to map (takes a while) - Figure 3A
goaccsmap <- goaccsmap + geom_point(aes(x=squid_obs$lon[1],y=squid_obs$lat[1]),col="red")
goaccsmap <- goaccsmap + geom_point(aes(x=squid_obs$lon[2],y=squid_obs$lat[2]),col="red")
goaccsmap <- goaccsmap + geom_point(aes(x=squid_obs$lon[3],y=squid_obs$lat[3]),col="red")
goaccsmap <- goaccsmap + geom_point(aes(x=squid_obs$lon[4],y=squid_obs$lat[4]),col="red")
goaccsmap

ggsave("", height = 5, width = 5)

#----


## Calculate migratory capability from published data (required for Figure 3B)
#----

# Zeidberg 2004 = 0.21 m s-1
# Payne & O’Dor 2006 = 0.15 m s-1
# 50 % travel time during period of maturation, O'Dor 1988
# i.e., 12 h per day for last 3 months of life

# lower estimate
((0.15*3600)/1000)*12*90
# 583.2 km lifespan-1

# upper estimate
((0.21*3600)/1000)*12*90
# 816.48 km lifespan-1

# combination
((mean(c(0.21,0.15))*3600)/1000)*12*90
# 699.84 km lifespan-1

#----



  ######### Data organization and exploration for modeling GOA squid capture probability ########

## Assess correlation between regional factors (temperature) and the global environmenal factor (MEI)
#----

envt_data_goa <- read_excel(file.list[1],sheet = "envt_data_goa",col_types = rep("numeric", 13))

# Correlation matrix function
# Source: http://www.sthda.com/english/wiki/correlation-matrix-an-r-function-to-do-all-you-need
rquery.cormat<-function(x,
                        type=c('lower', 'upper', 'full', 'flatten'),
                        graph=TRUE,
                        graphType=c("correlogram", "heatmap"),
                        col=NULL, ...)
{
  library(corrplot)
  # Helper functions
  #+++++++++++++++++
  # Compute the matrix of correlation p-values
  cor.pmat <- function(x, ...) {
    mat <- as.matrix(x)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(method = "pearson",mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # Get lower triangle of the matrix
  getLower.tri<-function(mat){
    upper<-mat
    upper[upper.tri(mat)]<-""
    mat<-as.data.frame(upper)
    mat
  }
  # Get upper triangle of the matrix
  getUpper.tri<-function(mat){
    lt<-mat
    lt[lower.tri(mat)]<-""
    mat<-as.data.frame(lt)
    mat
  }
  # Get flatten matrix
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  # Define color
  if (is.null(col)) {
    col <- colorRampPalette(
      c("#67001F", "#B2182B", "#D6604D", "#F4A582",
        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
        "#4393C3", "#2166AC", "#053061"))(200)
    col<-rev(col)
  }
  
  # Correlation matrix
  cormat<-signif(cor(method = "spearman",x, use = "complete.obs", ...),2)
  pmat<-signif(cor.pmat(x, ...),2)
  # Reorder correlation matrix
  ord<-corrMatOrder(cormat, order="hclust")
  cormat<-cormat[ord, ord]
  pmat<-pmat[ord, ord]
  # Replace correlation coeff by symbols
  sym<-symnum(cormat, abbr.colnames=FALSE)
  # Correlogram
  if(graph & graphType[1]=="correlogram"){
    corrplot(cormat, type=ifelse(type[1]=="flatten", "lower", type[1]),
             tl.col="black", tl.srt=45,col=col,...)
  }
  else if(graphType[1]=="heatmap")
    heatmap(cormat, col=col, symm=TRUE)
  # Get lower/upper triangle
  if(type[1]=="lower"){
    cormat<-getLower.tri(cormat)
    pmat<-getLower.tri(pmat)
  }
  else if(type[1]=="upper"){
    cormat<-getUpper.tri(cormat)
    pmat<-getUpper.tri(pmat)
    sym=t(sym)
  }
  else if(type[1]=="flatten"){
    cormat<-flattenCorrMatrix(cormat, pmat)
    pmat=NULL
    sym=NULL
  }
  list(r=cormat, p=pmat, sym=sym)
}

# Correlation between MEI, Sitka temperature trends, Yakutat temperature trends, and GOA sst trends
envt_data_corr <- select(envt_data_goa,c(
  mei,sit_temp_c_trend,yak_temp_c_trend,goa_sst_c_trend))

envt_data_corr_results <- rquery.cormat(envt_data_corr, type="flatten", graph=FALSE)
envt_data_corr_results <- as.data.frame(envt_data_corr_results$r)
envt_data_corr_results <- filter(envt_data_corr_results,row == "sit_temp_c_trend" | column == "sit_temp_c_trend")
envt_data_corr_results
#               row           column  cor        p
#1              mei sit_temp_c_trend 0.54  4.3e-24
#2  goa_sst_c_trend sit_temp_c_trend 0.95 4.0e-139
#3 sit_temp_c_trend yak_temp_c_trend 0.95 4.9e-138

#----


## Plot environmental factors: squid catch frequency and temperature (Figure 4A)
#----

# Create raster data frame for filled gradient plots
goa_squid_temp <- envt_data_goa %>%
  select(year,month,sit_temp_c_trend) %>%
  filter(year>1998) %>%
  mutate(x = seq(1:length(year)))

goa_squid_temp1 <- goa_squid_temp %>%
  mutate(y = -.01,
         x = seq(1:length(year)))
goa_squid_temp2 <- goa_squid_temp %>%
  mutate(y = 0,
         x = seq(1:length(year)))
goa_squid_temp3 <- goa_squid_temp %>%
  mutate(y = .01,
         x = seq(1:length(year)))
goa_squid_temp4 <- goa_squid_temp %>%
  mutate(y = .02,
         x = seq(1:length(year)))
goa_squid_temp5 <- goa_squid_temp %>%
  mutate(y = .03,
         x = seq(1:length(year)))
goa_squid_temp6 <- goa_squid_temp %>%
  mutate(y = .04,
         x = seq(1:length(year)))
goa_squid_temp7 <- goa_squid_temp %>%
  mutate(y = .05,
         x = seq(1:length(year)))
goa_squid_temp8 <- goa_squid_temp %>%
  mutate(y = .06,
         x = seq(1:length(year)))

goa_squid_temp_raster <- bind_rows(goa_squid_temp1,
                                   goa_squid_temp2,
                                   goa_squid_temp3,
                                   goa_squid_temp4,
                                   goa_squid_temp5,
                                   goa_squid_temp6,
                                   goa_squid_temp7,
                                   goa_squid_temp8)

# create dataframe to use for labels
x_labels <- filter(goa_squid_temp1,month==1)

# Calculate catch frequency
abundance_data_goa <- read_excel(file.list[1],sheet = "abundance_data_goa")
abundance_data_goa$station <- as.factor(abundance_data_goa$station)
goa_squid_freq <- abundance_data_goa %>%
  group_by(year) %>%
  summarise(catch_freq = (sum(catch_bin)/length(catch_bin))) %>%
  ungroup()
# Join to dataframe of all years
year_seq <- data.frame(seq(1990,2017,1))
names(year_seq) <- c("year")
goa_squid_freq <- left_join(year_seq,goa_squid_freq,"year")
# Filter to relevant year and add column for plotting
goa_squid_freq <- goa_squid_freq %>%
  filter(year>1998) %>%
  mutate(x = seq(5,(length(year)*12),by=12)) %>%
  filter(catch_freq!="NA")

# Plot squid catch frequency and temperature trend (Figure 4A)
goa_squid_temp_raster %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = sit_temp_c_trend), interpolate=T)+
  scale_fill_gradientn(colours=c("blue","blanchedalmond","red"),
                       name=expression(atop("Temperature", paste("trend (°C)")))) +
  geom_hline(yintercept = 0,lty=2) +
  geom_line(data = goa_squid_freq,aes(x,catch_freq)) +
  geom_point(data = goa_squid_freq,aes(x,catch_freq)) +
  theme_classic(base_size = 15) +
  scale_y_continuous(limits = c(-.005,.045), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 245), expand = c(0,0),
                     breaks=x_labels$x,labels=x_labels$year) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size=11)) +
  ylab("Catch frequency (fraction of obs.)") +
  xlab("")

ggsave("", height = 3.5, width = 6.2)

#----


## Organize environmental factors and explore their corrlations
#----

## Calculate winter-spring (Jan-May) averages of climate indices (MEI, temperature) for each year
environmental_data <- envt_data_goa %>%
  filter(month<6) %>%
  group_by(year) %>%
  summarise(mei = mean(mei),
            sit_temp_c_trend = mean(sit_temp_c_trend)) %>%
  ungroup()


## Calculate least squares means of ccs squid data to add to dataframe
ccs_squid <- read_excel(file.list[1],sheet = "abundance_data_ccs")
#abundance_data_ccs <- read_excel(file.list[],sheet = "abundance_data_ccs")
ccs_squid <- filter(ccs_squid, sci_name=="Doryteuthis opalescens")
# Make station, the random effect, a factor
ccs_squid$station <- as.factor(ccs_squid$station)
# Add log-transformed catch column
ccs_squid$catch_ln <- log(ccs_squid$catch+1)
# lsmean calculation
squid.lm1 <- lm(catch_ln ~ as.factor(year) + as.numeric(station), data = ccs_squid)
anova(squid.lm1)
squid.rg1 <- ref.grid(squid.lm1)
ccs.year.effects <- as.data.frame(summary(squid.rg1, "year"))
# Add lsmean predictions to the dataframe
ccs_lsm <- ccs.year.effects %>%
  rename(ccs_lsm = prediction) %>%
  select(c(year,ccs_lsm))
# Join ccs squid with other environmental predictors
environmental_data <- left_join(environmental_data,ccs_lsm,by="year")


## Explore the association between regional and large-scale environmental factors

# Function for (absolute) correlations on the upper panels, with size proportional to the correlations
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use = "complete.obs",method = "spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Associations between environmental factors, all at no lag
pairs(select(environmental_data,-c(year)), lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)
# Correlation matrix function for p-values
# Gather relevant correlations
corr_results <- rquery.cormat(select(environmental_data,-c(year)), type="flatten")
corr_results <- as.data.frame(corr_results$r)
filter(corr_results,column=="sit_temp_c_trend"|row=="sit_temp_c_trend")
#      row           column  cor       p
#1 ccs_lsm sit_temp_c_trend 0.36 0.22000
#2     mei sit_temp_c_trend 0.66 0.00079
filter(corr_results,column=="ccs_lsm"|row=="ccs_lsm")
#      row           column  cor    p
#1 ccs_lsm              mei 0.28 0.25
#2 ccs_lsm sit_temp_c_trend 0.36 0.22


## Add 1 & 2 year lagged environmental predictors
environmental_data <- environmental_data %>%
  mutate(mei.1 = lag(mei),
         mei.2 = lag(mei.1),
         sit_temp_c_trend.1 = lag(sit_temp_c_trend),
         sit_temp_c_trend.2 = lag(sit_temp_c_trend.1),
         ccs_lsm.1 = lag(ccs_lsm),
         ccs_lsm.2 = lag(ccs_lsm.1)) %>%
  select(year,mei,mei.1,mei.2,sit_temp_c_trend,sit_temp_c_trend.1,
         sit_temp_c_trend.2,ccs_lsm,ccs_lsm.1,ccs_lsm.2)

#----


## Explore associations at different time lags, check for outliers, and assess collinearity
#----

## Add environmental predictors to goa catch data
goa_squid_pred <- left_join(abundance_data_goa,environmental_data,"year")
# Then filter to only include years where all envt variables were collected
goa_squid_pred <- goa_squid_pred %>%
  filter(year>1998 & year<2019) %>%
  select(-c(abundance_ref))

## Association between squid catch and environmental predictors (at 0-2 year lags)
# Gather relevant correlations
corr_results <- rquery.cormat(goa_squid_pred[5:14], type="flatten")
corr_results <- as.data.frame(corr_results$r)
corr_results <- filter(corr_results,row =="catch_bin" | column =="catch_bin")
# there is a weird ordering output - fix
corr_results <- arrange(corr_results,row,column)
row_fix <- c(as.character(corr_results[6:9,2]))
col_fix <- c(as.character(corr_results[6:9,1]))
cor_fix <- c(as.numeric(corr_results[6:9,3]))
p_fix <- c(as.numeric(corr_results[6:9,4]))
corr_fix <- as.data.frame(cbind(row_fix,col_fix,cor_fix,p_fix))
names(corr_fix) <- c("row","column","cor","p")
corr_fix$row <- as.factor(corr_fix$row)
corr_fix$column <- as.factor(corr_fix$column)
corr_fix$cor <- as.numeric(as.character(corr_fix$cor))
corr_fix$p <- as.numeric(as.character(corr_fix$p))
corr_results <- corr_results[-c(6:9),]
corr_results <- bind_rows(corr_results,corr_fix)
corr_results <- filter(corr_results,column!="catch"&column!="catch_ln"&column!="catch_bin")
corr_results <- arrange(corr_results,row,column)
# now plot the corelations and p-values to determine appropriate lags to include
corr_results %>%
  ggplot(aes(column,cor)) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point() +
  geom_text(aes(label=p),vjust=1.5,size=3) +
  facet_wrap(~row) +
  theme(axis.text.x = element_text(angle = 90))
# Lags with highest significant correlations = ccs_lsm.2, mei.1, sit_temp_c_trend.2
# In some cases the differene between lags is minor
# Including ccs_lsm.2, mei, and sit_temp_c_trend

## Outliers
# function for histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(goa_squid_pred[c(5,6,9,14)], panel = panel.smooth,
      cex = 1.5, pch = 24, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
# no outliers


## Collinearity
# first check correlations
pairs(goa_squid_pred[c(5,6,9,14)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)
# in addition, calculate variance inflation factors (GVIF)
car::vif(gam(catch_bin~s(ccs_lsm.2,k=3)+s(sit_temp_c_trend,k=3)+s(mei,k=3)+s(station,bs="re"),
             data=goa_squid_pred,family = binomial(link = "logit"),method="REML"))
car::vif(gam(catch_bin~s(ccs_lsm.2,k=3)+s(sit_temp_c_trend,k=3)+s(station,bs="re"),
             data=goa_squid_pred,family = binomial(link = "logit"),method="REML"))
# mei is strongly associated with nearshore temperature
# thus, only including temp and ccs squid in model (GVIF<3)

#----



  ######### Modeling GOA squid capture probability ######### 

## Model selection (Table S4A)
#----

# Distribution - nb with logit link based on exploratory data analysis
# Terms - ccs_lsm.2, sit_temp_c_trend, and station as random effect

## Smoothness of terms
length(unique(goa_squid_pred$year))/2
# k should be less than n/2, or 5 in this case
# keeping k to 3 should suffice, as it has minimal effect on edf and reml

## Run all possible model combinations
# models with temp
m1 <- gam(catch_bin~ccs_lsm.2+sit_temp_c_trend
          +s(station,bs="re"),
          data=goa_squid_pred,family = nb(link = "logit"),
          correlation = corAR1(form = ~ year),
          method = "REML")
m2 <- gam(catch_bin~s(ccs_lsm.2,k=3,bs="cr")+sit_temp_c_trend
          +s(station,bs="re"),
          data=goa_squid_pred,family = nb(link = "logit"),
          correlation = corAR1(form = ~ year),
          method = "REML")
m3 <- gam(catch_bin~ccs_lsm.2+s(sit_temp_c_trend,k=3,bs="cr")
          +s(station,bs="re"),
          data=goa_squid_pred,family = nb(link = "logit"),
          correlation = corAR1(form = ~ year),
          method = "REML")
m4 <- gam(catch_bin~s(ccs_lsm.2,k=3,bs="cr")+s(sit_temp_c_trend,k=3,bs="cr")
          +s(station,bs="re"),
          data=goa_squid_pred,family = nb(link = "logit"),
          correlation = corAR1(form = ~ year),
          method = "REML")


## Compare models using Bayesian Information Criterion (BIC)
mod_comp <- data.frame(AIC(m1,m2,m3,m4))
mod_comp <- mutate(mod_comp,model = c("m1","m2","m3","m4"))
mod_comp <- arrange(mod_comp,AIC)
mod_comp
#        df      AIC model
#1 22.82422 309.2299    m3
#2 22.82428 309.2300    m4
#3 21.45664 310.5527    m1
#4 21.45666 310.5527    m2

#----


## Model validation
#----

# Model summary and smooth fits
summary(m3)
plot(m3,pages=1,residuals = T)

# 1 plot residuals against fitted vaules to assess homogeneity
# 2 histogram of residuals to verify normailty, also QQ plot
par(mfrow=c(2,2))
gam.check(m3)
# 1 plot residuals against fitted vaules to assess homogeneity
plot(predict(m3,type="response"),residuals(m4))
# 2 histogram of residuals to verify normailty, also QQ plot
qq.gam(m3,rep=20,level=1)
plot(m3$linear.predictors,m3$residuals)
plot(predict(m3,type="response"),m3$y);abline(0,1,col=2)
# 3 plot residuals against explanatory variables used in model
plot(goa_squid_pred$ccs_lsm.2,m3$residuals)
plot(goa_squid_pred$sit_temp_c_trend,m3$residuals)
# 4 plot residuals against explanatory variables not used in model
plot(goa_squid_pred$mei,m3$residuals)
# 5 cook distance for influential observations
plot(cooks.distance(m3),m3$y)
dev.off()

# Over/underdispersion
E2 <- resid(m3, type = "pearson")
N  <- nrow(goa_squid_pred)
p  <- length(coef(m3)) +1 # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)
# underdispersed 0.4983829
# thus fit could mask truly significant relationships

#----


## Model interpretation (Figure 4C, Figure 4D, and Table S4B)
#----

# useful instructions https://stats.idre.ucla.edu/r/dae/logit-regression/


## Determine station with maximum probabiity of squid catch
# build dataframe to predict the probability of squid at each station with mean values of other predictors
newdata1 <- with(goa_squid_pred, data.frame(ccs_lsm.2 = mean(ccs_lsm.2),
                                            sit_temp_c_trend = mean(sit_temp_c_trend),
                                            station = unique(station)))
# predict prbability of squid at each station 
newdata1$stationP <- predict(m3, newdata = newdata1, type = "response")
# sort stations from min to max probabiities
newdata1 <- arrange(newdata1,stationP)
# probability at station with greatest likelihood of catch
maxP <- as.character(newdata1[59,3]) # "50", 0.04412924

max(newdata1$stationP)-min(newdata1$stationP)
# 0.03949071


## Plot catch prob ~ ccs squid abundance from two years prior (Figure 4C)

# build dataframe with range of ccs_lsm.2 and average temp @ station with maximum probability of catch
newdat1 <- data.frame(ccs_lsm.2 = rep(seq(from = min(goa_squid_pred$ccs_lsm.2),to = max(goa_squid_pred$ccs_lsm.2),length.out = 100),1),
                      sit_temp_c_trend = mean(goa_squid_pred$sit_temp_c_trend),
                      station = maxP)
# predict probability of squid catch
prdat1 <- stats::predict(m3,
                         newdata = newdat1,
                         type = "link",
                         se.fit = T)
prdat1 <- data.frame(within(prdat1, {
  PredictedProb <- plogis(fit)
  ll <- plogis(fit - (1.96 * se.fit))
  ul<- plogis(fit + (1.96 * se.fit))
}))
# add predicted probability 95% CI to dataframe for plotting 
newdat1$predicted <- prdat1$PredictedProb
newdat1$ul <- prdat1$ul
newdat1$ll <- prdat1$ll
# plot prediction
newdat1 %>%
  ggplot(aes(ccs_lsm.2, predicted)) +
  geom_hline(yintercept = 0,lty=2)+
  geom_point(aes(ccs_lsm.2,catch_bin),data = select(goa_squid_pred,c(ccs_lsm.2,catch_bin)),alpha=0.3)+
  geom_line(size=1)+
  geom_ribbon(aes(ymin = ll,ymax = ul),fill="blue",alpha=0.3)+
  theme_classic(base_size = 12)+
  xlab(expression(atop("Log-transformed CCS squid", paste("catch two years prior"))))+
  ylab("Probability of squid capture in GOA")+
  ylim(c(0,1))

ggsave("", height = 3.5, width = 3)

# range of ccs squid values
max(newdat1$ccs_lsm.2)-min(newdat1$ccs_lsm.2)
# 3.980521

# increase in catch probability from min to max ccs squid
max(newdat1$predicted)-min(newdat1$predicted)
# 0.1270977

# as a percent of the average catch prob at the station with the highest catch prob
((max(newdat1$predicted)-min(newdat1$predicted))/newdata1[59,4])*100
# 288.0123


## Plot catch prob ~ nearshore temperature (Figure 4D)

# build dataframe with range of temp and average ccs squid @ station with maximum probability of catch
newdat1 <- data.frame(sit_temp_c_trend = rep(seq(from = min(goa_squid_pred$sit_temp_c_trend),to = max(goa_squid_pred$sit_temp_c_trend),length.out = 100),1),
                      ccs_lsm.2 = mean(goa_squid_pred$ccs_lsm.2),
                      station = maxP)
# predict probability of squid catch
prdat1 <- stats::predict(m3,
                         newdata = newdat1,
                         type = "link",
                         se.fit = T)
prdat1 <- data.frame(within(prdat1, {
  PredictedProb <- plogis(fit)
  ll <- plogis(fit - (1.96 * se.fit))
  ul<- plogis(fit + (1.96 * se.fit))
}))
# add predicted probability 95% CI to dataframe for plotting 
newdat1$predicted <- prdat1$PredictedProb
newdat1$ul <- prdat1$ul
newdat1$ll <- prdat1$ll
# plot prediction
newdat1 %>%
  ggplot(aes(sit_temp_c_trend, predicted)) +
  geom_hline(yintercept = 0,lty=2)+
  geom_point(aes(sit_temp_c_trend,catch_bin),data = select(goa_squid_pred,c(sit_temp_c_trend,catch_bin)),alpha=0.3)+
  geom_line(size=1)+
  geom_ribbon(aes(ymin = ll,ymax = ul),fill="blue",alpha=0.3)+
  theme_classic(base_size = 12)+
  scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,.25)) +
  xlab(expression(atop("Average winter-spring", paste("temperature trend (°C)"))))+
  ylab("Probability of squid capture in GOA")

ggsave("", height = 3.5, width = 3)

# max vs lowest value
max(newdat1$predicted)-newdat1$predicted[1]
# 0.0293365

# max vs highest value
max(newdat1$predicted)-newdat1$predicted[100]
# 0.03154309


## Overall model goodness of fit - likelihood ratio test
# chi-square - difference in deviance between full and null model
with(m3, null.deviance - deviance) # 65.40523
# degrees of freedom
with(m3, df.null - df.residual) # 19.25505
# p-value
with(m3, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE)) # 6.208907e-07
# the deviance residual is -2*log likelihood - log likelihood is:
logLik(m3) # -131.7907 (df=22.82422)


## Model summary and smooth fits (Table S4B)
summary(m3)
# Parametric coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -6.5356     0.6275 -10.415  < 2e-16 ***
# ccs_lsm.2     0.5327     0.1810   2.942  0.00326 **
# Approximate significance of smooth terms:
#                        edf  Ref.df Chi.sq p-value  
# s(sit_temp_c_trend)  1.767  1.945  3.469 0.21211   
# s(station)          15.488 58.000 26.550 0.00321 **
parameters(m3)
# Fixed Effects component
#Parameter   | Coefficient |   SE |         95% CI |      z |      p
#(Intercept) |       -6.54 | 0.63 | [-7.77, -5.31] | -10.41 | < .001
#ccs_lsm.2   |        0.53 | 0.18 | [ 0.18,  0.89] |   2.94 | 0.003 
# Smooth Terms component
#Parameter                      | Coefficient |     z |     p
#Smooth term (sit_temp_c_trend) |        1.77 |  3.47 | 0.212
#Smooth term (station)          |       15.49 | 26.55 | 0.003

#----



  ######### Experimental components of metabolic index (phi) ##########

## Compare squid masses between temperature treatments (Table S7B) and proportion of females between temperature treatments
#----

# Mass comparison
pcrit_rmr_squid_metadata <- read_excel(file.list[1],sheet = "pcrit_rmr_squid_metadata")
pcrit_rmr_squid_metadata$target_temp_C <- as.factor(pcrit_rmr_squid_metadata$target_temp_C)

aov1 = aov(mass_g~target_temp_C, data = pcrit_rmr_squid_metadata)
summary(aov1)
print(TukeyHSD(aov1))

## Table S7B
#             diff        lwr        upr     p adj
# 10-7.5  8.056288   1.660851 14.4517250 0.0056123
# 13-7.5  1.426903  -4.968534  7.8223404 0.9730610
# 16-7.5  3.336237  -3.718986 10.3914586 0.6926895
# 19-7.5  4.493723  -2.003275 10.9907203 0.3204085
# 13-10  -6.629385 -12.948834 -0.3099355 0.0344644
# 16-10  -4.720051 -11.706465  2.2663626 0.3445218
# 19-10  -3.562565  -9.984776  2.8596462 0.5485803
# 16-13   1.909333  -5.077081  8.8957472 0.9443270
# 19-13   3.066820  -3.355391  9.4890308 0.6847555
# 19-16   1.157486  -5.922015  8.2369875 0.9915820

# Deutsch et al. 2015 found no effect of body mass on Pcrit for atlantic rock crab or amphipods
# For atlantic cod and sunbream, effects are only shown over large body mass range (Fig S1.)
# larger mass corresponds to marginally higher Pcrit
# squid used in the 10C treatment were only 1-8 g larger than other treatments on average
# this is an insignificant difference in body size

# Proportion of females comparison
pcrit_rmr_temp <- read_excel(file.list[1],sheet = "pcrit_rmr_temp")
#pcrit_rmr_temp <- read_excel(file.list[],sheet = "pcrit_rmr_temp")
pcrit_rmr_temp <- mutate(pcrit_rmr_temp, prop_female = total_female/total_squid)
pcrit_rmr_temp$target_temp_C <- as.factor(pcrit_rmr_temp$target_temp_C)

aov2 = aov(prop_female~target_temp_C, data = pcrit_rmr_temp)
summary(aov2)
print(TukeyHSD(aov2))
#              diff        lwr       upr     p adj
#10-7.5 -0.15052335 -0.8537602 0.5527135 0.9465861
#13-7.5 -0.02129630 -0.7245332 0.6819406 0.9999704
#16-7.5 -0.05681217 -0.8430549 0.7294306 0.9990723
#19-7.5  0.04717792 -0.6560590 0.7504148 0.9993070
#13-10   0.12922705 -0.5740098 0.8324639 0.9684621
#16-10   0.09371118 -0.6925316 0.8799539 0.9935764
#19-10   0.19770127 -0.5055356 0.9009382 0.8722015
#16-13  -0.03551587 -0.8217586 0.7507269 0.9998549
#19-13   0.06847422 -0.6347627 0.7717111 0.9970371
#19-16   0.10399009 -0.6822527 0.8902328 0.9904913

# No significant differences
# Moreover, Burford et al. 2019 found no effect of proportion of females on group routine metabolic rate (Fig. A1)

#----


## Compare Pcrit between temperature treatments (Table S5A)
#----
# Calculate average ± 1SD pcrits for each treatment
pcrit_rmr_temp %>%
  group_by(target_temp_C) %>%
  summarise(pcrit_mean = mean(pcrit_kPa_O2),
            pcrit_sd = sd(pcrit_kPa_O2)) %>%
  ungroup()
#target_temp_C pcrit_mean pcrit_sd
#<fct>              <dbl>    <dbl>
#1 7.5                 4.68    0.948
#2 10                  4.69    0.478
#3 13                  3.80    0.730
#4 16                  4.65    1.42 
#5 19                  5.77    1.31 

# Relevel with respect to the treatment with the lowest pcrit (13)
pcrit_rmr_temp <- within(pcrit_rmr_temp, target_temp_C <- relevel(target_temp_C, ref = '13'))

# Relate Pcrit to temperature treatment using linear regression
summary(lm(pcrit_kPa_O2~target_temp_C,pcrit_rmr_temp))

## Table S5A
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        3.7974     0.5689   6.675 9.11e-05 ***
#target_temp_C7.5   0.8868     0.8046   1.102   0.2990    
#target_temp_C10    0.8952     0.8046   1.113   0.2947    
#target_temp_C16    0.8505     0.8996   0.945   0.3691    
#target_temp_C19    1.9767     0.8046   2.457   0.0364 *
#----


## Use Arrhenius function to describe Pcrit as a function of temperature (Figure S1B)
#----
# See Deutsch et al. 2015 for derivation and examples 

# First, calcualte natural log of Pcrit and inverse temperature
# Inverse temperature requires Boltzmann's constant (8.617333262145e-5 eV K-1) and temperature in K
pcrit_rmr_temp <- pcrit_rmr_temp %>%
  mutate(ln_Y = log(pcrit_kPa_O2),
         mean_temp_K = mean_temp_C+273,
         inv_T = 1/(8.617333262145e-5*mean_temp_K))

# Then fit linear regression to calcualte temperature coefficient (E) and normalization constant (A)
lm1 <- lm(ln_Y~inv_T,data = pcrit_rmr_temp)
summary(lm1)
# Gather temperature coefficient (E)
E <- as.numeric(lm1$coefficients[2]) #-0.0998039
# Gather normalization constant (A)
A <- as.numeric(lm1$coefficients[1]) #5.582831

# Function to predict Pcrit given temperature (K)
# Uses E, A, and no mass scaling (^1)
squid_phi <- function(mass_g, mean_temp_K){
  ((A*(mean(mass_g)^1))/exp(-E/(8.617333262145e-5*mean_temp_K)))
}

# to determine the mass to use, compare model fit across a range of masses
test <- data.frame(mean_temp_K = pcrit_rmr_temp$mean_temp_C+273,
                   pcrit_kPa_O2 = pcrit_rmr_temp$pcrit_kPa_O2)
mu_pcrit = mean(pcrit_rmr_temp$pcrit_kPa_O2)
test <- test %>%
  mutate(g_38 = squid_phi(mass_g = 38,mean_temp_K = mean_temp_K),
         g_39 = squid_phi(mass_g = 39,mean_temp_K = mean_temp_K),
         g_40 = squid_phi(mass_g = 40,mean_temp_K = mean_temp_K),
         g_41 = squid_phi(mass_g = 41,mean_temp_K = mean_temp_K),
         g_42 = squid_phi(mass_g = 42,mean_temp_K = mean_temp_K),
         g_43 = squid_phi(mass_g = 43,mean_temp_K = mean_temp_K),
         g_44 = squid_phi(mass_g = 44,mean_temp_K = mean_temp_K),
         g_45 = squid_phi(mass_g = 45,mean_temp_K = mean_temp_K),
         g_46 = squid_phi(mass_g = 46,mean_temp_K = mean_temp_K),
         g_47 = squid_phi(mass_g = 47,mean_temp_K = mean_temp_K),
         g_48 = squid_phi(mass_g = 48,mean_temp_K = mean_temp_K),
         g_49 = squid_phi(mass_g = 49,mean_temp_K = mean_temp_K),
         g_50 = squid_phi(mass_g = 50,mean_temp_K = mean_temp_K),
         g_51 = squid_phi(mass_g = 51,mean_temp_K = mean_temp_K)) %>%
  gather(g_38:g_51,key="mass_g",value = "phi_pred") %>%
  mutate(phi = rep(pcrit_rmr_temp$pcrit_kPa_O2,14)) %>%
  group_by(mass_g) %>%
  summarise(tss = sum((phi - mu_pcrit)^2),
            rss = sum((phi - phi_pred)^2),
            R2 = 1-(rss/tss)) %>%
  ungroup()
# Plot fits to check
test %>%
  ggplot(aes(mass_g,R2)) +
  geom_point() +
  geom_hline(yintercept=0,lty=2)
# Gather best fit
test <- filter(test,R2 == max(R2))
test
#mass_g   tss   rss     R2
#<chr>  <dbl> <dbl>  <dbl>
# g_49    14.6  13.2 0.0997

# Arrhenius function's fit of Pcrit temperature dependency
# Need to convert target temperature back to being numeric
pcrit_rmr_temp <- pcrit_rmr_temp %>%
  mutate(target_temp_C = case_when(
    mean_temp_C>11&mean_temp_C<14 ~ 13,
    mean_temp_C<8 ~ 7.5,
    mean_temp_C>18 ~ 19,
    mean_temp_C>14&mean_temp_C<18 ~ 16,
    mean_temp_C>8&mean_temp_C<11 ~ 10))

# Redefine function to predict Pcrit given temperature (K) using 49 g mass
# As before, uses E, A, and no mass scaling (^1)
squid_phi <- function(mean_temp_K){
  ((A*(mean(49)^1))/exp(-E/(8.617333262145e-5*mean_temp_K)))
}

# Plot mean ± 1SD Pcrit over temperature with fit of Arrhenius function (Figure S1B)
pcrit_rmr_temp %>%
  group_by(target_temp_C) %>%
  summarise(pcrit_kPa_O2_mean = mean(pcrit_kPa_O2),
            pcrit_kPa_O2_sd = sd(pcrit_kPa_O2)) %>%
  ungroup() %>%
  ggplot(aes(target_temp_C,pcrit_kPa_O2_mean)) +
  geom_point(aes(target_temp_C,pcrit_kPa_O2_mean)) +
  geom_errorbar(aes(ymin = pcrit_kPa_O2_mean-pcrit_kPa_O2_sd,
                    ymax = pcrit_kPa_O2_mean+pcrit_kPa_O2_sd),width=0) +
  geom_line(data = pcrit_rmr_temp, aes(target_temp_C,squid_phi(mean_temp_K = target_temp_C+273)),lty=2) +
  scale_x_discrete(limits = seq(7,19,2)) +
  theme_classic(base_size = 15) +
  ylab(expression("Mean"~P[crit]~"± 1SD"~(kPa))) +
  xlab("Temperature (°C)")

ggsave("", height = 3.5, width = 5)

#----


## Determine factorial areobic scope at 14°C (Figure S1C) and compare with Pcrit
#----

# First calculate Q10 from RMR in Pcrit experiments
pcrit_rmr_temp %>% 
  ggplot(aes(mean_temp_C,routine_metabolic_rate_mass_normalized_mgO2gh)) + 
  geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x + I(x^2),col="black") +
  scale_x_continuous(limits = c(7, 19),breaks=seq(7,19,2)) +
  theme_classic(base_size = 15) +
  ylim(c(0.2,0.72))

# Fit quadratic regression to data in order to get a 10C window for Q10 (9-19°C)
lm3 <- lm(routine_metabolic_rate_mass_normalized_mgO2gh~mean_temp_C+I(mean_temp_C^2),pcrit_rmr_temp)
summary(lm3)
#                     Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)       0.4892386  0.1076355   4.545 0.000837 ***
#  mean_temp_C      -0.0469599  0.0174694  -2.688 0.021100 *  
#  I(mean_temp_C^2)  0.0030023  0.0006544   4.588 0.000780 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.03088 on 11 degrees of freedom
#Multiple R-squared:  0.9621,	Adjusted R-squared:  0.9552 
#F-statistic: 139.5 on 2 and 11 DF,  p-value: 1.529e-08

# Function for that fit
coeff_lm3 <-  as.numeric(lm3$coefficients)

RMR <- function(temp){
  coeff_lm3[3]*(temp^2)+(coeff_lm3[2]*temp)+coeff_lm3[1]
}

# Q10 determination
Q10_squid <- (RMR(temp=9)/RMR(temp=19))^(10/(9-19))
Q10_squid # 2.197749

# Now calculate maximum and minimum metabolic rate

# Load experimental aerobic scope data and data from literature
aerobic_scope <- read_excel(file.list[1],sheet = "aerobic_scope")
#aerobic_scope <- read_excel(file.list[],sheet = "aerobic_scope")

# Also include RMR from Pcrit trials at 13°C, normalize to 14°C
temp <- pcrit_rmr_temp %>%
  filter(target_temp_C==13) %>%
  mutate(rmr_or_mmr = "rmr",
         metabolic_rate_mass_temp_normalized_mgO2gh = routine_metabolic_rate_mass_normalized_mgO2gh/(Q10_squid^(1/(14/(mean_temp_C-14)))),
         normalized_temp_C = 14,
         q10 = Q10_squid)
aerobic_scope <- bind_rows(aerobic_scope,temp)  

# Gather MMR from experimental and literature data, typically occurred at the fastest speed
as_max <- aerobic_scope %>%
  filter(rmr_or_mmr == "mmr") %>%
  group_by(specimen_id) %>%
  summarize(MMR = max(metabolic_rate_mass_temp_normalized_mgO2gh,na.rm = T)) %>%
  ungroup()

# Gather RMR from experimental and literature data
as_min <- aerobic_scope %>%
  filter(rmr_or_mmr == "rmr") %>%
  rename(RMR = metabolic_rate_mass_temp_normalized_mgO2gh)

## Plot average MMR and RMR ± 1SD (Figure S1C)
ggplot()+
  geom_point(aes(2,mean(as_max$MMR)))+
  geom_errorbar(aes(2,mean(as_max$MMR)),
                ymin=(mean(as_max$MMR)-sd(as_max$MMR)),
                ymax=(mean(as_max$MMR)+sd(as_max$MMR)),width=0)+
  geom_point(aes(1,mean(as_min$RMR)))+
  geom_errorbar(aes(1,mean(as_min$RMR)),
                ymin=(mean(as_min$RMR)-sd(as_min$RMR)),
                ymax=(mean(as_min$RMR)+sd(as_min$RMR)),width=0)+
  scale_y_continuous(limits = c(0,2)) +
  scale_x_continuous(limits = c(.5,2.5),breaks = NULL) +
  geom_hline(yintercept = 0) +
  ylab(expression(atop("Mean metabolic rate ± 1SD", paste("(",mgO[2],g^-1,h^-1,")")))) +
  xlab("RMR                     MMR") +
  theme(axis.text.x = element_blank())+
  theme_classic(base_size = 15)

ggsave("", height = 4, width = 5)

# Calculate MMR and RMR
MMR_squid <- mean(as_max$MMR)
MMR_squid #1.397153
sd(as_max$MMR) #0.2713901
RMR_squid <- mean(as_min$RMR)
RMR_squid #0.3457932
sd(as_min$RMR) #0.05863407

# Factorial aerobic scope (FAS) calculated as quotient of average MMR and RMR at 14 °C
FAS_squid = MMR_squid/RMR_squid
FAS_squid
#4.040429

# Pcrit at 14°C determined from Arrhenius function
Pcrit_squid = squid_phi(14+273)
Pcrit_squid
#4.83586

# Pcrit for MMR, should be ~ 21
# However, data from multiple studies at 15°C was 18.53 (Seibel & Deutsch 2020)
Pcrit_squid*FAS_squid
#19.53895

#----



  ######### Environmental components of metabolic index (phi) #########

## Assess correlation between regional factors (temperature, oxygen) relevant for phi of CCS and the global environmenal factor (MEI)
#----

envt_data_ccs <- read_excel(file.list[1],sheet = "envt_data_ccs",col_types = rep("numeric", 18))

# Correlation between MEI, Monterey Bay temperature trends, and Point Reyes temperature trends
envt_data_corr <- envt_data_ccs %>%
  select(c(mei,mba_temp_c_trend,ptr_temp_c_trend))

rquery.cormat(envt_data_corr,type = "full")

envt_data_corr_results <- rquery.cormat(envt_data_corr, type="flatten", graph=FALSE)
envt_data_corr_results <- as.data.frame(envt_data_corr_results$r)
envt_data_corr_results <- filter(envt_data_corr_results,row == "mba_temp_c_trend" | column == "mba_temp_c_trend")
envt_data_corr_results
#               row           column  cor       p
#1              mei mba_temp_c_trend 0.53 8.7e-30
#2 mba_temp_c_trend ptr_temp_c_trend 0.54 3.6e-20

# Correlation between MEI, Monterey Bay mgLO2 trends, and Trinidad mgLO2 trends
envt_data_corr <- envt_data_ccs %>%
  select(c(mei,mba_mgL_O2_trend,trin_mgL_O2_trend))

rquery.cormat(envt_data_corr,type = "full")

envt_data_corr_results <- rquery.cormat(envt_data_corr, type="flatten", graph=FALSE)
envt_data_corr_results <- as.data.frame(envt_data_corr_results$r)
envt_data_corr_results <- filter(envt_data_corr_results,row == "mba_mgL_O2_trend" | column == "mba_mgL_O2_trend")
envt_data_corr_results
#                row           column  cor    p
#1 trin_mgL_O2_trend mba_mgL_O2_trend 0.12 0.14
#2               mei mba_mgL_O2_trend 0.33 0.23

#----


## Calculate CCS metabolic index (Figure S1A, and Table S5B)
#----

# Filtering and mutating for metabolic index calculation
squid_ecophys <- envt_data_ccs %>%
  select(year,month,mei,mba_temp_c_trend,mba_kPa_O2_trend) %>%
  mutate(phi_kPa_trend = mba_kPa_O2_trend/squid_phi(mean_temp_K = mba_temp_c_trend+273)) %>%
  filter(year>2002) %>%
  mutate(x = seq(1:length(year)))

# Create dataframe to use for plot labels
x_labels <- filter(squid_ecophys,month==1)


## Plot CCS seawater temperature and dissolved oxygen trends (Figure S1A)
squid_ecophys %>%
  ggplot(aes(x, mba_temp_c_trend)) +
  geom_line(col="red",size=1.5)+
  geom_line(aes(x, mba_kPa_O2_trend),col="blue",size=1.5) +
  scale_y_continuous(breaks = seq(10,20,2), limits = c(10, 20),
                     sec.axis = sec_axis(trans = ~.*1,name = expression(O[2]~"partial pressure trend (kPa)"))) +
  scale_x_continuous(limits = c(12, 197), expand = c(0,0),
                     breaks=x_labels$x,labels=x_labels$year) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(color = c("red")),
        axis.title.y.right = element_text(color = "blue", size=15)) +
  ylab("Temperature trend (°C)") +
  xlab("")

ggsave("", height = 3.8, width = 5.5)


## Summary statistics for metabolic index

# Calculate average phi for whole timeseries
mean(squid_ecophys$phi_kPa_trend,na.rm=T)
#3.51436
sd(squid_ecophys$phi_kPa_trend,na.rm=T)
#0.2607722

# Calculate average ± 1SD phi for each year
squid_ecophys %>%
  filter(year>2003) %>%
  filter(month<6) %>%
  group_by(year) %>%
  summarise(fas = mean(phi_kPa_trend,na.rm=T),
            sd_fas = sd(phi_kPa_trend,na.rm=T)) %>%
  ungroup()
#year    fas  sd_fas
#<dbl>  <dbl>   <dbl>
#  1  2004   2.78  0.0602
#2  2005   3.67  0.0800
#3  2006   3.83  0.0116
#4  2007   3.49  0.0413
#5  2008   3.40  0.121 
#6  2009   3.94  0.0221
#7  2010 NaN    NA     
#8  2011   3.63  0.0160
#9  2012   3.58  0.0294
#10  2013   3.48  0.0184
#11  2014   3.55  0.0190
#12  2015   3.17  0.0342
#13  2016   3.66  0.0231
#14  2017   3.34  0.0171
#15  2018   3.42  0.0297
#16  2019   3.58  0.0111

# Now compare 2015 metabolic index with all other years (Table S5B)
squid_ecophys_lm <- squid_ecophys %>%
  filter(year>2003) %>%
  filter(month<6)
squid_ecophys_lm$year <- as.factor(squid_ecophys_lm$year)
# Relevel with respect to 2015
squid_ecophys_lm <- within(squid_ecophys_lm, year <- relevel(year, ref = '2015'))
# Calculate significance of difference using a linear regression
summary(lm(phi_kPa_trend~year,squid_ecophys_lm))
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  3.16676    0.02063 153.523  < 2e-16 ***
#  year2004    -0.38546    0.02917 -13.214  < 2e-16 ***
#  year2005     0.50425    0.02917  17.286  < 2e-16 ***
#  year2006     0.66014    0.02917  22.630  < 2e-16 ***
#  year2007     0.32597    0.02917  11.174 2.70e-16 ***
#  year2008     0.22955    0.02917   7.869 8.19e-11 ***
#  year2009     0.77213    0.02917  26.469  < 2e-16 ***
#  year2011     0.46223    0.02917  15.845  < 2e-16 ***
#  year2012     0.40890    0.02917  14.017  < 2e-16 ***
#  year2013     0.31663    0.02917  10.854 8.76e-16 ***
#  year2014     0.38713    0.02917  13.271  < 2e-16 ***
#  year2016     0.48847    0.02917  16.745  < 2e-16 ***
#  year2017     0.16898    0.02917   5.793 2.72e-07 ***
#  year2018     0.24837    0.02917   8.514 6.53e-12 ***
#  year2019     0.41296    0.02917  14.156  < 2e-16 ***

#----


## Sensitivity analysis of CCS metabolic index (Figure S1D)
#----

# Assess how changes temp and O2 impact phi estimates
# Determine temp and O2 during winter-spring (Jan-May) 2015
envt_data_ccs %>%
  select(year,month,mba_temp_c_trend,mba_kPa_O2_trend) %>%
  filter(year==2015) %>%
  filter(month<6) %>%
  group_by(year) %>%
  summarise(temp_c = mean(mba_temp_c_trend),
            oxygen_kPa = mean(mba_kPa_O2_trend)) %>%
  ungroup() %>%
  mutate(phi = oxygen_kPa/squid_phi(mean_temp_K = temp_c+273)) # included to check
#year temp_c oxygen_kPa   phi
#<dbl>  <dbl>      <dbl> <dbl>
#1  2015   13.3       15.2  3.17

# Create dataframe with a range of each parameter (50 % inc or dec), centered around these values
sensitivity <- tibble(temp = rep(seq(13.3+(13.3/2),13.3-(13.3/2),length.out = 100),100),
                      oxygen = rep(seq(15.2+(15.2/2),15.2-(15.2/2),length.out = 100),each = 100))

# Determine phi
sensitivity <- mutate(sensitivity,phi = oxygen/squid_phi(mean_temp_K = temp+273))

# Plot heatmap of phi with respect to different temp and O2 combinations (Figure S1D)
sensitivity %>%
  ggplot(aes(oxygen, temp)) +
  geom_raster(aes(fill = phi), interpolate=F)+
  scale_fill_gradientn(colours=c("red","orange","yellow"),
                       name=expression(atop("Metabolic index",paste("(",phi,")")))) +
  geom_contour(aes(z= phi),col="black") +
  #geom_line(data=filter(sensitivity,round(phi,2)==2),aes(oxygen,temp),col="blue") + # for reference
  theme_classic(base_size = 13) +
  scale_x_continuous(limits = c(min(sensitivity$oxygen),max(sensitivity$oxygen)), expand = c(0,0),
                     breaks=seq(7,23,2)) +
  scale_y_continuous(limits = c(min(sensitivity$temp),max(sensitivity$temp)), expand = c(0,0),
                     breaks=seq(7,22,2)) +
  xlab("Oxygen partial pressure (kPa)") +
  ylab("Temperature (°C)")

ggsave("", height = 4, width = 5.7)


# How sensitive is the final estimate of phi to changes in the original timeseries?
# Increase temperature by 25% from the average value
# Then, decrease oxygen by 25% from the average value
# Finally, do both
temp_up <- mean(envt_data_ccs$mba_temp_c_trend,na.rm=T)/4
oxy_down <- mean(envt_data_ccs$mba_kPa_O2_trend,na.rm=T)/4

envt_data_ccs_test <- envt_data_ccs %>%
  mutate(temp_c = mba_temp_c_trend + temp_up,
         oxy_kPa = mba_kPa_O2_trend - oxy_down,
         phi = mba_kPa_O2_trend/squid_phi(mean_temp_K = mba_temp_c_trend+273),
         phi_temp_up = mba_kPa_O2_trend/squid_phi(mean_temp_K = temp_c+273),
         phi_oxy_down = oxy_kPa/squid_phi(mean_temp_K = mba_temp_c_trend+273),
         phi_both = oxy_kPa/squid_phi(mean_temp_K = temp_c+273))

mean(envt_data_ccs_test$phi,na.rm=T) #3.51436
mean(envt_data_ccs_test$phi_temp_up,na.rm=T) #3.367403
mean(envt_data_ccs_test$phi_oxy_down,na.rm=T) #2.635732
mean(envt_data_ccs_test$phi_both,na.rm=T) #2.525516

(3.51436-3.367403)/3.51436*100 #4.181615
(3.51436-2.635732)/3.51436*100 #25.00108
(3.51436-2.525516)/3.51436*100 #28.13724


# Now compare 2015 metabolic index with all other years (Table S5B)
envt_data_ccs_test_lm <- envt_data_ccs_test %>%
  filter(year>2014) %>%
  filter(month<6)
envt_data_ccs_test_lm$year <- as.factor(envt_data_ccs_test_lm$year)
# Relevel with respect to 2015
envt_data_ccs_test_lm <- within(envt_data_ccs_test_lm, year <- relevel(year, ref = '2015'))
# Calculate significance of difference using a linear regression
summary(lm(trin_phi_kPa_trend~year,envt_data_ccs_test_lm))
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  4.15090    0.01526 272.059  < 2e-16 ***
#year2016     0.15679    0.02158   7.267 4.98e-07 ***
#year2017     0.32364    0.02158  14.999 2.41e-12 ***
#year2018     0.04283    0.02158   1.985    0.061 .  
#year2019    -0.01824    0.02158  -0.845    0.408


#----


  ######### Data organization and exploration for modeling CCS squid abundance ######### 

## Plot environmental factors: squid abundance, competitor abundance, and metabolic index (Figure 4B)
#----

# Filter and mutate metabolic index data
squid_ecophys <- envt_data_ccs %>%
  select(year,month,mei,mba_temp_c_trend,mba_kPa_O2_trend) %>%
  mutate(phi_kPa_trend = mba_kPa_O2_trend/squid_phi(mean_temp_K = mba_temp_c_trend+273)) %>%
  filter(year>1998) %>%
  mutate(x = seq(1:length(year)))

# Create raster data frame for filled gradient plots
squid_ecophys1 <- squid_ecophys %>%
  mutate(y = -1,
         x = seq(1:length(year)))
squid_ecophys2 <- squid_ecophys %>%
  mutate(y = 0,
         x = seq(1:length(year)))
squid_ecophys3 <- squid_ecophys %>%
  mutate(y = 1,
         x = seq(1:length(year)))
squid_ecophys4 <- squid_ecophys %>%
  mutate(y = 2,
         x = seq(1:length(year)))
squid_ecophys5 <- squid_ecophys %>%
  mutate(y = 3,
         x = seq(1:length(year)))
squid_ecophys6 <- squid_ecophys %>%
  mutate(y = 4,
         x = seq(1:length(year)))
squid_ecophys7 <- squid_ecophys %>%
  mutate(y = 5,
         x = seq(1:length(year)))
squid_ecophys8 <- squid_ecophys %>%
  mutate(y = 6,
         x = seq(1:length(year)))
squid_ecophys9 <- squid_ecophys %>%
  mutate(y = 7,
         x = seq(1:length(year)))
squid_ecophys10 <- squid_ecophys %>%
  mutate(y = 8,
         x = seq(1:length(year)))

squid_ecophys_raster <- bind_rows(squid_ecophys1,
                                  squid_ecophys2,
                                  squid_ecophys3,
                                  squid_ecophys4,
                                  squid_ecophys5,
                                  squid_ecophys6,
                                  squid_ecophys7,
                                  squid_ecophys8,
                                  squid_ecophys9,
                                  squid_ecophys10,)

# Create dataframe to use for labels
x_labels <- filter(squid_ecophys1,month==1)

# Calculate least squares means of ccs squid data to add to dataframe
abundance_data_ccs <- read_excel(file.list[1],sheet = "abundance_data_ccs")

ccs_squid <- filter(abundance_data_ccs, sci_name=="Doryteuthis opalescens")
# Make station, the random effect, a factor
ccs_squid$station <- as.factor(ccs_squid$station)
# Add log-transformed catch column
ccs_squid$catch_ln <- log(ccs_squid$catch+1)
# lsmean calculation
squid.lm1 <- lm(catch_ln ~ as.factor(year) + as.numeric(station), data = ccs_squid)
anova(squid.lm1)
squid.rg1 <- ref.grid(squid.lm1)
ccs.year.effects <- as.data.frame(summary(squid.rg1, "year"))
# Add lsmean predictions to the dataframe
squid_lsm <- ccs.year.effects %>%
  select(c(year,prediction,SE)) %>%
  filter(year>1998) %>%
  mutate(x = seq(5,(length(year)*12),by=12))

# Calculate least squares means of ccs anchovy data to add to dataframe
ccs_fish <- filter(abundance_data_ccs, sci_name!="Doryteuthis opalescens")
ccs_fish <- ccs_fish %>%
  mutate(year.station = str_c(year,".",station)) %>%
  group_by(year.station) %>%
  summarise(catch = sum(catch),
            station = mean(station),
            year = mean(year)) %>%
  ungroup()
# Make station, the random effect, a factor
ccs_fish$station <- as.factor(ccs_fish$station)
# Add log-transformed catch column
ccs_fish$catch_ln <- log(ccs_fish$catch+1)
# lsmeans
fish.lm1 <- lm(catch_ln ~ as.factor(year) + as.numeric(station), data = ccs_fish)
anova(fish.lm1)
fish.rg1 <- ref.grid(fish.lm1)
ccs.year.effects <- as.data.frame(summary(fish.rg1, "year"))
# Add lsmean predictions to the dataframe
fish_lsm <- ccs.year.effects %>%
  select(c(year,prediction,SE)) %>%
  filter(year>1998) %>%
  mutate(x = seq(5,(length(year)*12),by=12))

# Plot squid abundance, competitor abundance, and metabolic index (Figure 4B)
squid_ecophys_raster %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = phi_kPa_trend), interpolate=T)+
  scale_fill_gradientn(colours=c("red","orange","yellow"),
                       name=expression(atop("Metabolic index",paste("(",phi,")")))) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point(data = squid_lsm,aes(x,prediction)) +
  geom_line(data = squid_lsm,aes(x,prediction)) +
  geom_errorbar(data = squid_lsm,aes(x,prediction,
                                     ymax=(prediction+SE),ymin=(prediction-SE)),width=0) +
  geom_point(data = fish_lsm,aes(x,prediction),col="blue3") +
  geom_line(data = fish_lsm,aes(x,prediction),col="blue3") +
  geom_errorbar(data = fish_lsm,aes(x,prediction,
                                    ymax=(prediction+SE),ymin=(prediction-SE)),width=0,col="blue3") +
  scale_y_continuous(limits = c(-.5, 7.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 245), expand = c(0,0),
                     breaks=x_labels$x,labels=x_labels$year) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size=11)) +
  ylab("Mean normalized catch ± SE") +
  xlab("")

ggsave("", height = 3.5, width = 6.4)

#----


## Organize environmental factors and explore their corrlations
#----

# Calculate winter-spring (Jan-May) averages of global predictor (MEI) for each year
environmental_data <- envt_data_ccs %>%
  filter(month<6) %>%
  group_by(year) %>%
  summarise(mei = mean(mei)) %>%
  ungroup()

# Calculate winter-spring averages for metabolic index year
squid_ecophysiology <- squid_ecophys %>%
  filter(month<6) %>%
  group_by(year) %>%
  summarise(phi_kPa_trend = mean(phi_kPa_trend)) %>%
  ungroup()
environmental_data <- left_join(environmental_data,squid_ecophysiology,by="year")

# Competitor (fish) abundance for each year
ccs.year.effects <- as.data.frame(summary(fish.rg1, "year"))
fish_lsm <- ccs.year.effects %>%
  rename(fish_lsm = prediction) %>%
  select(c(year,fish_lsm))
environmental_data <- left_join(environmental_data,fish_lsm,by="year")

# Function for (absolute) correlations on the upper panels, with size proportional to the correlations
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use = "complete.obs",method = "spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Associations between environmental factors, all at no lag
pairs(select(environmental_data,-c(year)), lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)
# Gather relevant correlations
corr_results <- rquery.cormat(select(environmental_data,-c(year)), type="flatten")
corr_results <- as.data.frame(corr_results$r)
filter(corr_results,column=="phi_kPa_trend"|row=="phi_kPa_trend")
#            row   column    cor    p
#1 phi_kPa_trend      mei -0.068 0.77
#2 phi_kPa_trend fish_lsm  0.200 0.89
filter(corr_results,column=="fish_lsm"|row=="fish_lsm")
#            row   column  cor    p
#1 phi_kPa_trend fish_lsm 0.20 0.89
#2           mei fish_lsm 0.23 0.41

# Add 1 & 2 year lagged environmental predictors
environmental_data <- environmental_data %>%
  mutate(mei.1 = lag(mei),
         mei.2 = lag(mei.1),
         phi_kPa_trend.1 = lag(phi_kPa_trend),
         phi_kPa_trend.2 = lag(phi_kPa_trend.1),
         fish_lsm.1 = lag(fish_lsm),
         fish_lsm.2 = lag(fish_lsm.1)) %>%
  select(year,mei,mei.1,mei.2,phi_kPa_trend,phi_kPa_trend.1,
         phi_kPa_trend.2,fish_lsm,fish_lsm.1,fish_lsm.2)

#----


## Explore associations at different time lags, check for outliers, and assess collinearity
#----

# Select relevant columns of squid abundance data
ccs_squid <- select(ccs_squid,year,station,catch_ln)
# Add environmental predictors to ccs catch data
ccs_squid_pred <- left_join(ccs_squid,environmental_data,"year")
# Then filter to exclude missing data
ccs_squid_pred <- filter(ccs_squid_pred,year>2004&year!=2011)

# Association between squid abundance and environmental predictors (at 0-2 year lags)
# Gather relevant correlations
corr_results <- rquery.cormat(ccs_squid_pred[3:12], type="flatten")
corr_results <- as.data.frame(corr_results$r)
corr_results <- filter(corr_results,row=="catch_ln"|column=="catch_ln")
corr_results <- arrange(corr_results,row,column)
# Now plot the corelations and p-values to determine appropriate lags to include
corr_results %>%
  ggplot(aes(column,cor)) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point() +
  geom_text(aes(label=p),vjust=1.5,size=3) +
  theme(axis.text.x = element_text(angle = 90))
# In some cases the differene between lags is minor
# Selected lags for model = phi_kPa_trend.1, fish_lsm.2, mei.1

# Outliers
# Function for histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(ccs_squid_pred[c(3,5,8,12)], panel = panel.smooth,
      cex = 1.5, pch = 24, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
# Outliers in catch - will need non-normal distribution and link, or transform

# Collinearity
# First check correlations
pairs(ccs_squid_pred[c(5,8,12)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)
# In addition, calculate variance inflation factors (GVIF)
car::vif(lm(catch_ln~phi_kPa_trend.1+fish_lsm.2+mei.1+station,
            data=ccs_squid_pred))
# Thus, including all factors in model (GVIF<3)

#----



  ######### Modeling CCS squid abundance ######### 

## Model selection (Table S4C)
#----

# Distribution is gaussian after log-transformation of response
# Terms - fish_lsm.2, fas_kPa_trend.1, mei.1, and station as random effect

# Smoothness of terms
length(unique(ccs_squid_pred$year))/2
# k should be less than n/2, or 6.5 in this case
# keeping k to 3 should suffice, as it has minimal effect on edf and reml
# run four possible combinations of smooth/no smooth terms
m1 <- gam(catch_ln~fish_lsm.2+phi_kPa_trend.1+mei.1
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m2 <- gam(catch_ln~s(fish_lsm.2,k=3,bs="cr")+phi_kPa_trend.1+mei.1
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m3 <- gam(catch_ln~fish_lsm.2+s(phi_kPa_trend.1,k=3,bs="cr")+mei.1
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m4 <- gam(catch_ln~fish_lsm.2+phi_kPa_trend.1+s(mei.1,k=3,bs="cr")
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m5 <- gam(catch_ln~s(fish_lsm.2,k=3,bs="cr")+s(phi_kPa_trend.1,k=3,bs="cr")+mei.1
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m6 <- gam(catch_ln~s(fish_lsm.2,k=3,bs="cr")+phi_kPa_trend.1+s(mei.1,k=3,bs="cr")
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m7 <- gam(catch_ln~fish_lsm.2+s(phi_kPa_trend.1,k=3,bs="cr")+s(mei.1,k=3,bs="cr")
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)
m8 <- gam(catch_ln~s(fish_lsm.2,k=3,bs="re")+s(phi_kPa_trend.1,k=3,bs="re")+s(mei.1,k=3,bs="re")
          +s(station,bs="re"),family = gaussian(),
          correlation = corAR1(form = ~ year),
          method = "REML",
          data=ccs_squid_pred)


## Compare models using Bayesian Information Criterion (AIC)
mod_comp <- data.frame(AIC(m1,m2,m3,m4,m5,m6,m7,m8))
mod_comp <- mutate(mod_comp,model = c("m1","m2","m3","m4","m5","m6","m7","m8"))
mod_comp <- arrange(mod_comp,AIC)
mod_comp
#        df      AIC model
#1 41.57057 1466.733    m5
#2 40.46314 1472.165    m3
#3 42.21402 1473.540    m7
#4 39.53979 1492.985    m2
#5 41.39092 1495.033    m6
#6 37.25717 1517.767    m8
#7 38.12006 1519.378    m1
#8 38.24667 1519.618    m4

#----


## Model validation
#----

# Model summary and smooth fits
summary(m5)
plot(m5,pages=1,residuals = T)

# 1 plot residuals against fitted vaules to assess homogeneity
# 2 histogram of residuals to verify normailty, also QQ plot
par(mfrow=c(2,2))
gam.check(m5)
# 1 plot residuals against fitted vaules to assess homogeneity
plot(predict(m5,type="response"),residuals(m5))
# 2 histogram of residuals to verify normailty, also QQ plot
qq.gam(m5,rep=20,level=1)
plot(m5$linear.predictors,m5$residuals)
plot(predict(m5,type="response"),m5$y);abline(0,1,col=2)
# 3 plot residuals against explanatory variables used in model
plot(ccs_squid_pred$fish_lsm.2,m5$residuals)
plot(ccs_squid_pred$phi_kPa_trend.1,m5$residuals)
# 4 plot residuals against explanatory variables not used in model
plot(ccs_squid_pred$mei.1,m5$residuals)
# 5 cook distance for influential observations
plot(cooks.distance(m5),m5$y)
dev.off()

#----


## Model interpretation (Figure 4E, Figure 4F, Figure 4G, and Table S4D)
#----

# Build dataframe to predict normalized squid catch at each station
newdata1 <- with(ccs_squid_pred, data.frame(fish_lsm.2 = mean(fish_lsm.2),
                                            phi_kPa_trend.1 = mean(phi_kPa_trend.1),
                                            mei.1 = mean(mei.1),
                                            station = unique(station)))
# Predict normalized squid catch at each station 
newdata1$stationP <- predict(m5, newdata = newdata1, type = "response")
# Gather the stations with the median catch over years
newdata1 <- arrange(newdata1,stationP)
# Catch at station with median catch
medP <- as.character(newdata1[(nrow(newdata1)/2+3),4])
medP #"165", 2.02


## Plot normalized ccs squid ~ fish range @ mean phi, mean mei, and at the station with median catch (Figure 4E)

# Build dataframe with range of fish_lsm.2 and average temp @ median station
newdat1 <- data.frame(fish_lsm.2 = seq(from = min(ccs_squid_pred$fish_lsm.2),to = max(ccs_squid_pred$fish_lsm.2),length.out = 100),
                      phi_kPa_trend.1 = mean(ccs_squid_pred$phi_kPa_trend.1),
                      mei.1 = mean(ccs_squid_pred$mei.1),
                      station = medP)
prdat1 <- stats::predict(m5,
                         newdata = newdat1,
                         type = "response",
                         se.fit = T)
newdat1$predicted <- prdat1$fit
newdat1$conf <- stats::qnorm(.95) * prdat1$se.fit

# Plot prediction
newdat1 %>% 
  ggplot(aes(fish_lsm.2, predicted)) +
  geom_hline(yintercept = 0,lty=2)+
  geom_point(aes(fish_lsm.2,catch_ln),data = select(ccs_squid_pred,c(fish_lsm.2,catch_ln)),alpha=0.3)+
  geom_line(size=1) +
  geom_ribbon(aes(ymin = (predicted-conf),ymax = (predicted+conf)), fill="red", alpha=0.3) +
  theme_classic(base_size = 12) +
  scale_y_continuous(limits = c(-1.5, 9),breaks = seq(0,9,3)) +
  xlab(expression(atop("Log-transformed CCS anchovy and", paste("sardine catch two years prior"))))+
  ylab("Log-transformed CCS squid catch")

ggsave("", height = 3.5, width = 3)

# range of fish values
max(newdat1$fish_lsm.2)-min(newdat1$fish_lsm.2)
# 3.691901

# change in squid over full range of fish
max(newdat1$predicted)-min(newdat1$predicted)
# 3.403859

# same, but as a % of mean squid value
((max(newdat1$predicted)-min(newdat1$predicted))/mean(newdat1$predicted))*100
# 201.8232


## Plot ccs squid ~ phi @ mean fish and mei (Figure 4F)

# Build dataframe with range of fish_lsm.2 and average temp @ median station
newdat1 <- data.frame(phi_kPa_trend.1 = seq(from = min(ccs_squid_pred$phi_kPa_trend.1),to = max(ccs_squid_pred$phi_kPa_trend.1),length.out = 100),
                      fish_lsm.2 = mean(ccs_squid_pred$fish_lsm.2),
                      mei.1 = mean(ccs_squid_pred$mei.1),
                      station = medP)
prdat1 <- stats::predict(m5,
                         newdata = newdat1,
                         type = "response",
                         se.fit = T)
newdat1$predicted <- prdat1$fit
newdat1$conf <- stats::qnorm(.95) * prdat1$se.fit

# Plot prediction
newdat1 %>% 
  ggplot(aes(phi_kPa_trend.1, predicted)) +
  geom_hline(yintercept = 0,lty=2)+
  geom_point(aes(phi_kPa_trend.1,catch_ln),data = select(ccs_squid_pred,c(phi_kPa_trend.1,catch_ln)),alpha=0.3)+
  geom_line(size=1) +
  geom_ribbon(aes(ymin = (predicted-conf),ymax = (predicted+conf)), fill="red", alpha=0.3) +
  theme_classic(base_size = 12) +
  scale_y_continuous(limits = c(-1.5, 9),breaks = seq(0,9,3)) +
  scale_x_continuous(limits = c(2.7, 4),breaks = seq(2.8,4,.4)) +
  xlab(expression(atop("Average winter-spring"~phi, paste("the year prior"))))+
  ylab("Log-transformed CCS squid catch")

#ggsave("~/Desktop/Denny Lab/Alaska market squid/Paper v2/American Naturalist/Revision/Figure 3F_rev.pdf", height = 3.5, width = 3)
#ggsave("", height = 3.5, width = 3)

# phi at maximum predicted value
filter(newdat1, predicted == max(newdat1$predicted))
#  phi_kPa_trend.1 fish_lsm.2      mei.1 station predicted      conf
#1        3.553023   1.160883 -0.3233262     165  2.039196 0.8988966

# change in squid over full range of phi
max(newdat1$predicted)-min(newdat1$predicted)
# 2.410164

# same, but as a % of mean squid value
((max(newdat1$predicted)-min(newdat1$predicted))/mean(newdat1$predicted))*100
# 178.0785


## Plot ccs squid ~ mei @ mean fish and phi (Figure 4G)

# Build dataframe with range of fish_lsm.2 and average temp @ median station
newdat1 <- data.frame(mei.1 = seq(from = min(ccs_squid_pred$mei.1),to = max(ccs_squid_pred$mei.1),length.out = 100),
                      fish_lsm.2 = mean(ccs_squid_pred$fish_lsm.2),
                      phi_kPa_trend.1 = mean(ccs_squid_pred$phi_kPa_trend.1),
                      station = medP)
prdat1 <- stats::predict(m5,
                         newdata = newdat1,
                         type = "response",
                         se.fit = T)
newdat1$predicted <- prdat1$fit
newdat1$conf <- stats::qnorm(.95) * prdat1$se.fit

# Plot prediction
newdat1 %>% 
  ggplot(aes(mei.1, predicted)) +
  geom_hline(yintercept = 0,lty=2)+
  geom_point(aes(mei.1,catch_ln),data = select(ccs_squid_pred,c(mei.1,catch_ln)),alpha=0.3)+
  geom_line(size=1) +
  geom_ribbon(aes(ymin = (predicted-conf),ymax = (predicted+conf)), fill="red", alpha=0.3) +
  theme_classic(base_size = 12) +
  scale_y_continuous(limits = c(-1.5, 9),breaks = seq(0,9,3)) +
  scale_x_continuous(limits = c(-1.75, 1.75),breaks = seq(-1.5,1.5,1)) +
  xlab(expression(atop("Average winter-spring ENSO", paste("the year prior"))))+
  ylab("Log-transformed CCS squid catch")

ggsave("", height = 3.5, width = 3)

# range of mei
max(newdat1$mei.1)-min(newdat1$mei.1)
# 3.186

# change in squid over full range of mei
max(newdat1$predicted)-min(newdat1$predicted)
# 0.5514878

# same, but as a % of mean squid value
((max(newdat1$predicted)-min(newdat1$predicted))/mean(newdat1$predicted))*100
# 27.89262


## Overall model summary statistics
m5
# Family: gaussian 
# Link function: identity 
# Estimated degrees of freedom:
# total = 39.83
summary.gam(m5)
#Parametric coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   2.0937     0.2772   7.553  4.1e-13 ***
#mei.1        -0.1731     0.1198  -1.445    0.149
#Approximate significance of smooth terms:
#                      Approximate significance of smooth terms:
#                      edf Ref.df       F  p-value    
#s(fish_lsm.2)       1.881  1.985 101.322  < 2e-16 ***
#s(phi_kPa_trend.1)  1.961  1.998  23.832 2.06e-10 ***
#s(station)         33.991 41.000   7.437  < 2e-16 ***

#R-sq.(adj) =  0.61

## Summary statistics (Table S4D) 
parameters(m5)
# Fixed Effects component
#Parameter   | Coefficient |   SE |        95% CI |     t |     df |      p
#(Intercept) |        2.09 | 0.28 | [ 1.54, 2.64] |  7.55 | 334.17 | < .001
#mei.1       |       -0.17 | 0.12 | [-0.41, 0.06] | -1.45 | 334.17 | 0.149 
# Smooth Terms component
#Parameter                     | Coefficient |      F |    df |      p
#Smooth term (fish_lsm.2)      |        1.88 | 101.75 |  1.98 | < .001
#Smooth term (phi_kPa_trend.1) |        1.96 |  23.83 |  2.00 | < .001
#Smooth term (station)         |       33.94 |   7.44 | 41.00 | < .001

#----



  ######### Comparative size and age at maturity ######### 
library(tidyverse)
library(SIBER)
library(viridis)
library(nlme)
library(lubridate)


## Plot size and age comparison between locations and years (Figure 5A)
#----

# Load data
age_size_maturity_goa_ccs <- read_excel(file.list[1],sheet = "age_size_maturity_goa_ccs")

age_size_maturity_goa_ccs <- age_size_maturity_goa_ccs %>%
  separate(capture_date,c("year","month","day")) %>%
  mutate(yearloc = str_c(year," ",capture_location))
  
# Prepare colors for plot
palette(c("red","black","#365C8D"))

# Use SIBER to make 95% bivariate ellipses 
dopa.siber.data<-age_size_maturity_goa_ccs[,c("age_at_maturity_days","size_at_maturity_length_mm", "yearloc")]
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="age_at_maturity_days"] <- "iso1"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="size_at_maturity_length_mm"] <- "iso2"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="yearloc"] <- "group"
dopa.siber.data$group <- as.numeric(factor(dopa.siber.data$group))
dopa.siber.data$community<-1
dopa.siber.data<-as.data.frame(dopa.siber.data)
dopa.siber.data<-subset(dopa.siber.data,iso1!="NA")
my.siber.data <- createSiberObject(dopa.siber.data)
communityMetricsML(my.siber.data)
names(my.siber.data)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, ci.mean = T, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# Plot age ~ size with 95% bivariate ellipses and colored by yearloc (Figure 5A)
par(mfrow=c(1,1))
setwd("")
pdf('',height = 4.5,width = 4)
plotSiberObject(my.siber.data,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(2,1),
                ylab="Age at maturity (days)",
                xlab="Size at maturity (dorsal mantle length, mm)",
                cex = 0.5
)
dev.off()

#----


## Compare age between locations and years (Figure S2, Table S6A)
#----

# Check normality of age at maturity in 2016 GOA
age_size_maturity_goa_ccs_1 <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2017 GOA" & yearloc != "2016 CCS")
shapiro.test(age_size_maturity_goa_ccs_1$age_at_maturity_days)
# W = 0.97831, p-value = 0.7491

# Check normality of age at maturity in 2016 CCS
age_size_maturity_goa_ccs_1 <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2017 GOA" & yearloc != "2016 GOA")
shapiro.test(age_size_maturity_goa_ccs_1$age_at_maturity_days)
# W = 0.93578, p-value = 0.1313

# Check normality of age at maturity in 2017 GOA
age_size_maturity_goa_ccs_1 <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2016 GOA" & yearloc != "2016 CCS")
shapiro.test(age_size_maturity_goa_ccs_1$age_at_maturity_days)
# W = 0.97152, p-value = 0.4016

# Make yearloc a factor, then organize with respect to 2016 Alaska squid
age_size_maturity_goa_ccs$yearloc <- as.factor(age_size_maturity_goa_ccs$yearloc)
age_size_maturity_goa_ccs <- within(age_size_maturity_goa_ccs, yearloc <- relevel(yearloc, ref = '2016 GOA'))

# ANCOVA of age ~ yearloc with size and sex as covariates, include all interactions
summary(aov(age_at_maturity_days~yearloc*size_at_maturity_length_mm*sex, data = age_size_maturity_goa_ccs))
#                                       Df Sum Sq Mean Sq F value   Pr(>F)    
#yearloc                                 2  64864   32432  31.522 6.12e-11 ***
#  size_at_maturity_length_mm            1   1796    1796   1.746   0.1900    
#sex                                     1    826     826   0.802   0.3729    
#yearloc:size_at_maturity_length_mm      2   7195    3598   3.497   0.0348 *  
#  yearloc:sex                           2   3483    1741   1.692   0.1903    
#size_at_maturity_length_mm:sex          1   1325    1325   1.288   0.2597    
#yearloc:size_at_maturity_length_mm:sex  2   5565    2783   2.705   0.0727 .  
#Residuals                              84  86426    1029   

# Plot covariation between size and age at maturity (Figure S2)
age_size_maturity_goa_ccs %>% 
  ggplot(aes(size_at_maturity_length_mm,age_at_maturity_days,color=yearloc)) + 
  scale_color_manual(values=c("black","red","#365C8D")) +
  geom_point() +
  geom_smooth(data = filter(age_size_maturity_goa_ccs,yearloc=="2016 GOA"),aes(size_at_maturity_length_mm,age_at_maturity_days),method = "lm") +
  geom_smooth(data = filter(age_size_maturity_goa_ccs,yearloc=="2016 CCS"),aes(size_at_maturity_length_mm,age_at_maturity_days),method = "lm") +
  geom_smooth(data = filter(age_size_maturity_goa_ccs,yearloc=="2017 GOA"),aes(size_at_maturity_length_mm,age_at_maturity_days),method = "lm") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank()) +
  ylab("Age at maturity (days)") +
  xlab("Size at maturity (dorsal mantle length, mm)")

ggsave("", height = 3, width = 4.5)

# Explore the significance of covariation between age and size at maturity
summary(lm(age_at_maturity_days~size_at_maturity_length_mm,filter(age_size_maturity_goa_ccs,yearloc=="2016 GOA")))
#            Estimate Std. Error t value Pr(>|t|)
#size_at_maturity_length_mm       -0.3437     0.3766  -0.913    0.369

summary(lm(age_at_maturity_days~size_at_maturity_length_mm,filter(age_size_maturity_goa_ccs,yearloc=="2016 CCS")))
#            Estimate Std. Error t value Pr(>|t|)
#size_at_maturity_length_mm       -0.2193     0.5128  -0.428  0.67312

summary(lm(age_at_maturity_days~size_at_maturity_length_mm,filter(age_size_maturity_goa_ccs,yearloc=="2017 GOA")))
#            Estimate Std. Error t value Pr(>|t|)
#size_at_maturity_length_mm        1.0025     0.3703   2.707   0.0101 *

# Cant do full ANOVA, violates non-homogeneity of regression slopes
# Will make pairwise comparisons between GOA 2016 and other locations


## First, compare GOA 2016 vs. CCS 2016, both have same age/size covariation so use ANOVA
# ANCOVA of age ~ yearloc with size and sex as covariates, include all interactions
summary(aov(age_at_maturity_days~yearloc*size_at_maturity_length_mm*sex, data = filter(age_size_maturity_goa_ccs, yearloc != "2017 GOA")))
#                    Df Sum Sq Mean Sq F value   Pr(>F)    
# yearloc             1  61908   61908  65.394 1.63e-10 ***
# size_at_maturity_length_mm              1    756     756   0.798   0.3761    
# sex                 1    660     660   0.697   0.4080    
# yearloc:size_at_maturity_length_mm      1     47      47   0.050   0.8239    
# yearloc:sex         1     69      69   0.073   0.7877    
# size_at_maturity_length_mm:sex          1   6039    6039   6.378   0.0149 *  
# yearloc:size_at_maturity_length_mm:sex  1    329     329   0.347   0.5585    
# Residuals          48  45442     947 

# No covariates interacted with independent variable, so eliminate interactions
summary(aov(age_at_maturity_days~yearloc+size_at_maturity_length_mm+sex, data = filter(age_size_maturity_goa_ccs, yearloc != "2017 GOA")))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# yearloc      1  61908   61908  61.997 2.02e-10 ***
# size_at_maturity_length_mm       1    756     756   0.757    0.388    
# sex          1    660     660   0.661    0.420    
# Residuals   52  51926     999 

# None of covariates significant, so can just do t-test
# Save test object for subsequent calcs
t1 <- t.test(age_at_maturity_days~yearloc, data = filter(age_size_maturity_goa_ccs, yearloc != "2017 GOA"))
t1
# t = 7.4445, df = 36.646, p-value = 7.754e-09
#  mean in group 2016 Alaska mean in group 2016 California 
#                   254.8125                      187.6250

# Calculate age difference
as.numeric(t1$estimate[1]-t1$estimate[2])
# 67.1875 days older in GOA than CCS on average
(as.numeric(t1$estimate[1]-t1$estimate[2])/t1$estimate[2])*100
# 35.80946 % older on average in GOA than CCS, as a percent of CCS average


## Next, compare 2016 GOA vs. 2017 GOA
# Multiple regression of age ~ yearloc with size and sex as covariates (allows for non-homogeneity of covariation between), include all interactions
summary(lm(age_at_maturity_days~yearloc*size_at_maturity_length_mm*sex, data = filter(age_size_maturity_goa_ccs, yearloc != "2016 CCS")))
#                                                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                                      445.8637   113.4493   3.930 0.000211 ***
#  yearloc2017 GOA                                 -398.6860   132.5087  -3.009 0.003748 ** 
#  size_at_maturity_length_mm                        -1.2473     0.7337  -1.700 0.093958 .  
#sexm                                            -216.2374   147.5882  -1.465 0.147778    
#yearloc2017 GOA:size_at_maturity_length_mm         2.5233     0.8826   2.859 0.005731 ** 
#  yearloc2017 GOA:sexm                             267.5981   173.2757   1.544 0.127435    
#size_at_maturity_length_mm:sexm                    1.4244     0.9432   1.510 0.135910    
#yearloc2017 GOA:size_at_maturity_length_mm:sexm   -1.9231     1.1462  -1.678 0.098246 .

# Covariate interacts with independent variable, so run multiple regression with significant covariates/interactions 
# Store model object for subsequent calculations
mr1 <- lm(age_at_maturity_days ~ yearloc + size_at_maturity_length_mm + yearloc:size_at_maturity_length_mm,
          data = filter(age_size_maturity_goa_ccs, yearloc != "2016 CCS"))
summary(mr1)
# (Table S6A)
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                308.3667    72.3747   4.261 6.42e-05 ***
# yearloc2017 Alaska        -232.5455    85.6711  -2.714  0.00841 ** 
# size_at_maturity_length_mm                      -0.3437     0.4632  -0.742  0.46064    
# yearloc2017 Alaska:size_at_maturity_length_mm    1.3462     0.5681   2.370  0.02065 *

# Calculate difference between age at average size of GOA 2016 and 2017 
# Average size size of GOA 2016
size <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2016 CCS" & yearloc != "2017 GOA") %>%
  group_by(yearloc) %>%
  summarise(size = mean(size_at_maturity_length_mm)) %>%
  ungroup()

# Function from model coefficients to calculate age at average size
age_1 <- function(size){
  mr1$coefficients[1]+(mr1$coefficients[3]*size)
}

# Replace with actual value
age_1 <- as.numeric(age_1(size$size))
# 254.8125 days

# Repeat for GOA 2017
size <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2016 CCS" & yearloc != "2016 GOA") %>%
  group_by(yearloc) %>%
  summarise(size = mean(size_at_maturity_length_mm)) %>%
  ungroup()
age_2 <- function(size){
  (mr1$coefficients[1]+mr1$coefficients[2])+((mr1$coefficients[3]+mr1$coefficients[4])*size)
}
age_2 <- as.numeric(age_2(size$size))
# 214.7625 days

# Calculate age difference
age_1-age_2
# 40.05 days older in 2016 than 2017 on average
((age_1-age_2)/age_2)*100
# 18.64851 % older in 2016 on average, as a percent of 2017 average

#----


## Compare size between locations and years (Table S6B)
#----

# Check normality of size at maturity in 2016 GOA
age_size_maturity_goa_ccs_1 <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2017 GOA" & yearloc != "2016 CCS")
shapiro.test(age_size_maturity_goa_ccs_1$size_at_maturity_length_mm)
# W = 0.92603, p-value = 0.0304

# Check normality of size at maturity in 2016 CCS
age_size_maturity_goa_ccs_1 <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2017 GOA" & yearloc != "2016 GOA")
shapiro.test(age_size_maturity_goa_ccs_1$size_at_maturity_length_mm)
# W = 0.95789, p-value = 0.3976

# Check normality of size at maturity in 2017 GOA
age_size_maturity_goa_ccs_1 <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2016 GOA" & yearloc != "2016 CCS")
shapiro.test(age_size_maturity_goa_ccs_1$size_at_maturity_length_mm)
# W = 0.95071, p-value = 0.08022


## First, compare 2016 GOA vs. 2016 CCS
# GOA 2016 non-normal, but use one-way anova because not particularly sensitive to non-normality
# Also, both have same age/size covariation (see previous section)
# ANCOVA of size ~ yearloc with age and sex as covariates, include all interactions
summary(aov(size_at_maturity_length_mm~yearloc*age_at_maturity_days*sex, data = filter(age_size_maturity_goa_ccs, yearloc != "2017 GOA")))
#                      Df Sum Sq Mean Sq F value   Pr(>F)    
# yearloc               1   3829    3829  21.094 3.18e-05 ***
# age_at_maturity_days              1    145     145   0.798   0.3761    
# sex                   1    377     377   2.074   0.1563    
# yearloc:age_at_maturity_days      1     24      24   0.134   0.7155    
# yearloc:sex           1     35      35   0.193   0.6626    
# age_at_maturity_days:sex          1    892     892   4.914   0.0314 *  
# yearloc:age_at_maturity_days:sex  1     41      41   0.226   0.6367    
# Residuals            48   8712     182 

# No covariates interacted with independent variable, so eliminate interactions
summary(aov(size_at_maturity_length_mm~yearloc+age_at_maturity_days+sex, data = filter(age_size_maturity_goa_ccs, yearloc != "2017 GOA")))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# yearloc      1   3829    3829  20.515 3.49e-05 ***
# age_at_maturity_days     1    145     145   0.776    0.382    
# sex          1    377     377   2.017    0.161    
# Residuals   52   9704     187

# None of covariates significant, so can just do t-test
t2 <- t.test(size_at_maturity_length_mm~yearloc, data = filter(age_size_maturity_goa_ccs, yearloc != "2017 GOA"))
t2
# t = 4.3064, df = 40.555, p-value = 0.0001025
# mean in group 2016 GOA mean in group 2016 CCS 
#                  155.8125                      139.1042
as.numeric(t2$estimate[1]-t2$estimate[2])
# 16.70833 mm larger in GOA
as.numeric((t2$estimate[1]-t2$estimate[2])/t2$estimate[2])*100
# 12.01138 % larger on average in GOA than CCS, as a percent of CCS average


## Next, compare 2016 GOA vs. 2017 GOA
# Need to do multiple regression to allow for non-homogeneity of slopes
summary(lm(size_at_maturity_length_mm~yearloc*age_at_maturity_days*sex, data = filter(age_size_maturity_goa_ccs, yearloc != "2016 CCS")))
#                                          Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                               193.0418    28.6272   6.743 5.24e-09 ***
#  yearloc2017 GOA                           -93.4422    33.8739  -2.759  0.00756 ** 
#  age_at_maturity_days                       -0.1526     0.1125  -1.357  0.17947    
#sexm                                      -56.7293    54.9339  -1.033  0.30564    
#yearloc2017 GOA:age_at_maturity_days        0.3276     0.1378   2.377  0.02043 *  
#  yearloc2017 GOA:sexm                       61.5395    60.5887   1.016  0.31360    
#age_at_maturity_days:sexm                   0.2389     0.2133   1.120  0.26697    
#yearloc2017 GOA:age_at_maturity_days:sexm  -0.2490     0.2436  -1.022  0.31047 

# Covariate interacts with independent variable, so run multiple regression with significant covariates/interactions 
# Store model object for subsequent calcs
mr2 <- lm(size_at_maturity_length_mm ~ yearloc + age_at_maturity_days + yearloc:age_at_maturity_days,
          data = filter(age_size_maturity_goa_ccs, yearloc != "2016 CCS"))
summary(mr2)
# (Table S6B)
#                              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 175.84531   24.10145   7.296 4.14e-10 ***
# yearloc2017 GOA          -71.87587   27.01961  -2.660  0.00973 ** 
# age_at_maturity_days                     -0.07862    0.09416  -0.835  0.40667    
# yearloc2017 GOA:age_at_maturity_days   0.23987    0.10959   2.189  0.03205 *

# Calculate difference between size at average age of GOA 2016 and 2017 
# Average age of GOA 2016
age <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2016 CCS" & yearloc != "2017 GOA") %>%
  group_by(yearloc) %>%
  summarise(age = mean(age_at_maturity_days)) %>%
  ungroup()

# Function from model coefficients to calculate age at average age
size_1 <- function(age){
  mr2$coefficients[1]+(mr2$coefficients[3]*age)
}

# Replace with actual value
size_1 <- as.numeric(size_1(age$age))
# 155.8125 mm

# Repeat for GOA 2017
age <- age_size_maturity_goa_ccs %>%
  filter(yearloc != "2016 CCS" & yearloc != "2016 GOA") %>%
  group_by(yearloc) %>%
  summarise(age = mean(age_at_maturity_days)) %>%
  ungroup()
size_2 <- function(age){
  (mr2$coefficients[1]+mr2$coefficients[2])+((mr2$coefficients[3]+mr2$coefficients[4])*age)
}
size_2 <- as.numeric(size_2(age$age))
# 138.6 mm

# Calculate age difference
size_1-size_2
# 17.2125 mm larger in 2016 than 2017 on average
((size_1-size_2)/size_2)*100
# 12.41883 % larger in 2016 on average, as a percent of 2017 average

#----


## Examine trends in size at maturity of CCS squid over time (Table S7A, Figure 5B)
#----

# Load data
size_maturity_ccs_time <- read_excel(file.list[1],sheet = "size_maturity_ccs_time")

size_maturity_ccs_time <- size_maturity_ccs_time %>%
  mutate(sex = replace(sex, sex=="f", c("female"))) %>%
  mutate(sex = replace(sex, sex=="m", c("male")))

# ANCOVA of size ~ time with sex as a covariate, include all interactions
summary(aov(size_at_maturity_length_mm~day_of_study_period*sex, data = size_maturity_ccs_time))
# Df Sum Sq Mean Sq F value   Pr(>F)    
# day_of_study_period           1  29659   29659 158.187  < 2e-16 ***
# sex            1   6654    6654  35.492 3.28e-09 ***
# day_of_study_period:sex       1    283     283   1.508     0.22    
# Residuals   1321 247676     187

# No covariates interacted with independent variable, so eliminate interactions
# Store model object for subsequent calcs
lm1 <- aov(size_at_maturity_length_mm~day_of_study_period+sex, data = size_maturity_ccs_time)
summary(lm1)
# (Table S7A)
#               Df Sum Sq Mean Sq F value  Pr(>F)    
# day_of_study_period           1  29659   29659  158.13 < 2e-16 ***
# sex            1   6654    6654   35.48 3.3e-09 ***
# Residuals   1322 247959     188

(lm1$coefficients[2])*365
# squid size decreasing -6.378559 mm/year
lm1$coefficients[3]
# males 5.354277 mm larger than females on average


## Plot size ~ time in CCS colored by sex (Figure 5B)
# Use model to make predictions and 95% CI
newdat1 <- data.frame(day_of_study_period = rep(seq(from = min(size_maturity_ccs_time$day_of_study_period),to = max(size_maturity_ccs_time$day_of_study_period),length.out = 100),2),
                      sex = c(rep("male",100),rep("female",100)))
prdat1 <- stats::predict(lm1,
                         newdata = newdat1,
                         type = "response",
                         se.fit = T)
newdat1$predicted <- prdat1$fit
newdat1$conf <- stats::qnorm(.95) * prdat1$se.fit

# Plot predictions from model, make sure sex is spelled out
newdat1 %>%
  ggplot(aes(day_of_study_period, predicted)) +
  geom_point(aes(day_of_study_period,size_at_maturity_length_mm,shape=sex),data = size_maturity_ccs_time,alpha=0.3)+
  geom_ribbon(aes(ymin = (predicted-conf),ymax = (predicted+conf), fill=sex), alpha=0.4) +
  geom_line(aes(col=sex),size=.7) +
  scale_color_manual(values = c("black","black")) +
  scale_shape_manual(values = c(19,2)) +
  scale_fill_manual(values = c("blue","red")) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = c(0.13, 0.13)) +
  scale_y_continuous(limits = c(65, 175), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 1490), expand = c(0,0),
                     breaks=c(0,365,730,1095,1460),labels=c(2016,2017,2018,2019,2020)) +
  ylab("Size at maturity (dorsal mantle length, mm)") +
  xlab("")

ggsave("", height = 3.7, width = 4)

#----


## Compare size between sexes in each location (Figure 5C)
#----

# Select relevant data from California dataset
squid_comp <- size_maturity_ccs_time %>%
  select(size_at_maturity_length_mm,sex) %>%
  mutate(year = factor("CCS 2016-2019"))

# Check normality
squid_comp_m <- filter(squid_comp,sex=="male")
squid_comp_f <- filter(squid_comp,sex=="female")
hist(squid_comp_m$size_at_maturity_length_mm)
shapiro.test(squid_comp_m$size_at_maturity_length_mm)
# W = 0.97785, p-value = 2.199e-11
hist(squid_comp_f$size_at_maturity_length_mm)
shapiro.test(squid_comp_f$size_at_maturity_length_mm)
# W = 0.98327, p-value = 0.001389

# Non-normal, but use one-way anova because not particularly sensitive to non-normality
# plus data are close to normal
# just be sure that significance is far below 0.05
summary(aov(size_at_maturity_length_mm~sex,squid_comp))
#               Df Sum Sq Mean Sq F value   Pr(>F)    
# sex            1   8012    8012   38.37 7.81e-10 ***
# Residuals   1323 276260     209

# To get means
t.test(squid_comp_m$size_at_maturity_length_mm,squid_comp_f$size_at_maturity_length_mm)
# 130.4806 (m)  124.6120 (f)

# Same with Alaska data
squid_comp1 <- age_size_maturity_goa_ccs %>%
  select(size_at_maturity_length_mm,sex) %>%
  mutate(year = factor("GOA 2016-2017"))

# Check normality
squid_comp1_m <- filter(squid_comp1,sex=="m")
squid_comp1_f <- filter(squid_comp1,sex=="f")
hist(squid_comp1_m$size_at_maturity_length_mm)
shapiro.test(squid_comp1_m$size_at_maturity_length_mm)
# W = 0.97615, p-value = 0.3903
hist(squid_comp1_f$size_at_maturity_length_mm)
shapiro.test(squid_comp1_f$size_at_maturity_length_mm)
# W = 0.92203, p-value = 0.004956

# Non-normal, but use one-way anova because not particularly sensitive to non-normality
# Plus data are close to normal
# Just be sure that significance is far below 0.05
summary(aov(size_at_maturity_length_mm~sex,squid_comp1))
#             Df Sum Sq Mean Sq F value Pr(>F)
# sex          1     80   79.67   0.301  0.585
# Residuals   94  24902  264.91
t.test(squid_comp1_m$size_at_maturity_length_mm,squid_comp1_f$size_at_maturity_length_mm) # to get means
# 143.6078 (m)  145.4333 (f)


# Plot comparison of size between sexes in both locations (Figure 5C)
# Bind data together
squid <- bind_rows(squid_comp,squid_comp1)

# Include full names for both sexes
squid %>%
  mutate(sex = replace(sex, sex=="f", c("female"))) %>%
  mutate(sex = replace(sex, sex=="m", c("male"))) %>%
  ggplot(aes(year,size_at_maturity_length_mm,fill=sex)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("blue","red")) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.15)) +
  xlab("") +
  ylab("Size at maturity (dorsal mantle length, mm)")

ggsave("",height = 3.7, width = 3.5)

#----


## Hatching season and size and age at maturity in GOA (Table S7C, Figure 5D, Table S7D, Figure S3)
#----

# Calculate birthdays and add size and age categories for plotting
birth <- age_size_maturity_goa_ccs %>%
  filter(capture_location=="GOA") %>%
  mutate(year = as.numeric(year),
         date = str_c(year,"-",month,"-",day),
         doy=yday(date),
         bday = ifelse(age_at_maturity_days<=doy, doy-age_at_maturity_days, 365-abs(doy-age_at_maturity_days)),
         byear = ifelse(age_at_maturity_days<=doy, year, year-1),
         bdate = as.Date(bday, origin = str_c(byear,"-",01,"-",01)),
         bmonth = as.numeric(format(as.Date(bdate), "%m")),
         size = ifelse(size_at_maturity_length_mm<=mean(size_at_maturity_length_mm), "< 144 mm", "> 144 mm"),
         age = ifelse(age_at_maturity_days<=mean(age_at_maturity_days), "< 233 days", "> 233 days"))

# Add season as a factor for plotting
birth <- birth %>%
  mutate(bseason = replace(bmonth, bmonth==12|bmonth<=2, c("winter"))) %>%
  mutate(bseason = replace(bseason, bmonth>2&bmonth<=5, c("spring"))) %>%
  mutate(bseason = replace(bseason, bmonth>5&bmonth<=8, c("summer"))) %>%
  mutate(bseason = replace(bseason, bmonth>8&bmonth<=11, c("fall")))


## Relate size at maturity to hatching season, starting with sex and year as covariates
summary(aov(size_at_maturity_length_mm~bseason*sex*factor(year),data = birth))
#                         Df Sum Sq Mean Sq F value   Pr(>F)    
#bseason                   3   7448  2482.7  15.933 9.87e-08 ***
#  sex                       1     10     9.7   0.062   0.8042    
#factor(year)              1   1088  1088.0   6.983   0.0105 *  
#  bseason:sex               3    106    35.3   0.227   0.8775    
#bseason:factor(year)      1      4     3.7   0.024   0.8780    
#sex:factor(year)          1    106   106.2   0.682   0.4123    
#bseason:sex:factor(year)  1     47    46.9   0.301   0.5851    
#Residuals                60   9349   155.8 

# No significant interactions so remove
summary(aov(size_at_maturity_length_mm~bseason+sex+factor(year),data = birth))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#bseason       3   7448  2482.7  17.047 2.64e-08 ***
#  sex           1     10     9.7   0.066  0.79754    
#factor(year)  1   1088  1088.0   7.471  0.00804 ** 
#  Residuals    66   9612   145.6 

# Sex not significant so remove for final model
lm2 <- aov(size_at_maturity_length_mm~bseason+factor(year),data = birth)
summary(lm2)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#bseason       3   7448  2482.7  17.288 2.03e-08 ***
#factor(year)  1   1088  1087.7   7.574  0.00761 ** 
#Residuals    67   9622   143.6
TukeyHSD(lm2)
# (Table S7C)
#                    diff        lwr        upr     p adj
#spring-fall    21.189258  11.090582  31.287935 0.0000034
#summer-fall    17.192157   7.607277  26.777037 0.0000708
#winter-fall   -22.441176 -46.043794   1.161441 0.0682302
#summer-spring  -3.997101 -12.747669   4.753466 0.6268131
#winter-spring -43.630435 -66.906733 -20.354136 0.0000322
#winter-summer -39.633333 -62.691363 -16.575304 0.0001446
#               diff       lwr       upr     p adj
#2017-2016 -5.500924 -11.17396 0.1721076 0.0571564

# Plot count of squid by birth month with respect to size (Figure 5D)
seasons <- data_frame(
  x = c(1.5,4,7,10,12.5),
  y = rep(17,5),
  lab = c("W","Sp","Su","F","W"),
  size = rep("< 144 mm",5)
)

birth %>%
  ggplot(aes(bmonth,fill=factor(size))) +
  facet_wrap(facets = ~year) +
  geom_histogram(binwidth=1,alpha = 0.5) +
  scale_x_continuous(limits = c(1,12.5),breaks = seq(1,12,1)) +
  scale_y_continuous(limits = c(0,17),breaks = seq(0,17,4)) +
  scale_fill_manual(values = c("red","blue"))+
  theme_classic(base_size = 12) +
  geom_vline(xintercept=c(2.5,5.5,8.5,11.5),size=0.2) +
  geom_text(data = seasons,aes(x,y,label=lab),size = 4) +
  labs(fill = expression(atop("Size at maturity", paste("(dorsal mantle length)")))) +
  theme(legend.title = element_text(size=9),
        legend.position = c(0.47, .3)) +
  ylab("Count") +
  xlab("Hatching month")

ggsave("",height = 3, width = 7)


## Relate age at maturity to hatching season, starting with sex and year as covariates
summary(aov(age_at_maturity_days~bseason*sex*year,data = birth))
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
#bseason           3  56037   18679  29.401 8.07e-12 ***
#  sex               1   2673    2673   4.207   0.0446 *  
#  year              1     14      14   0.023   0.8810    
#bseason:sex       3    933     311   0.489   0.6909    
#bseason:year      1    220     220   0.347   0.5581    
#sex:year          1   2266    2266   3.567   0.0638 .  
#bseason:sex:year  1    419     419   0.660   0.4197    
#Residuals        60  38119     635  

# No significant interactions so remove
summary(aov(age_at_maturity_days~bseason+sex+year,data = birth))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#bseason      3  56037   18679  29.383 3.49e-12 ***
#  sex          1   2673    2673   4.204   0.0443 *  
#  year         1     14      14   0.023   0.8810    
#Residuals   66  41957     636

# Year not significant so remove for final model
lm3 <- aov(age_at_maturity_days~bseason+sex,data = birth)
summary(lm3)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#bseason      3  56037   18679  29.818 2.32e-12 ***
#sex          1   2673    2673   4.266   0.0427 *  
#Residuals   67  41972     626
TukeyHSD(lm3)
# (Table S7D)
#                    diff        lwr        upr     p adj
#spring-fall     66.78730   45.69550  87.879095 0.0000000
#summer-fall     36.53464   16.51594  56.553338 0.0000523
#winter-fall    -42.17647  -91.47220   7.119261 0.1194106
#summer-spring  -30.25266  -48.52883 -11.976480 0.0002620
#winter-spring -108.96377 -157.57796 -60.349578 0.0000008
#winter-summer  -78.71111 -126.86943 -30.552790 0.0003177
#                 diff       lwr        upr     p adj
#male-female -12.22732 -24.07585 -0.3787932 0.0433038

# Plot count of squid by birth month with respect to age (Figure S3)
seasons <- data_frame(
  x = c(1.5,4,7,10,12.5),
  y = rep(13,5),
  lab = c("W","Sp","Su","F","W"),
  age = rep("< 233 days",5)
)

birth %>%
  ggplot(aes(bmonth,fill = factor(age))) +
  facet_wrap(facets = ~sex) +
  geom_histogram(binwidth=1,alpha = 0.5) +
  scale_x_continuous(limits = c(1,12.5),breaks = seq(1,12,1)) +
  scale_y_continuous(limits = c(0,13),breaks = seq(0,12,3)) +
  scale_fill_manual(values = c("red","blue"))+
  theme_classic(base_size = 12) +
  geom_vline(xintercept=c(2.5,5.5,8.5,11.5),size=0.2) +
  geom_text(data = seasons,aes(x,y,label=lab),size = 4) +
  labs(fill = c("Age at maturity")) +
  theme(legend.title = element_text(size=9)) +
  ylab("Count") +
  xlab("Hatching month")

ggsave("",height = 3, width = 9)

#----


## Size at maturity comparison with historical data (Figure 5E)
#----

# Convert to spelled out sexes
age_size_maturity_goa_ccs <- age_size_maturity_goa_ccs %>%
  mutate(sex = replace(sex, sex=="f", c("female"))) %>%
  mutate(sex = replace(sex, sex=="m", c("male")))

# Data from Fields (1965) of body size of spawning squid from 1948-1962
f_tot <- 3763
m_tot <- 3897
tot <- sum(f_tot,m_tot)
f_mean <- 143
m_mean <- 146.7
f_rand <- c(89,116,121,129,139,143,151,156,162,164)
m_rand <- c(107,111,120,125,137,142,152,161,169,173)
f_max <- 180
f_min <- 85
m_max <- 191
m_min <- 72

# Calculate upper and lower limits of 95% CI
# using t-didtribution b/c small sample size (df = 10-1)
f_low <- f_mean + (qt(0.025,9)*(sd(f_rand)/sqrt(f_tot)))
f_up <- f_mean + (qt(0.975,9)*(sd(f_rand)/sqrt(f_tot)))
m_low <- m_mean + (qt(0.025,9)*(sd(m_rand)/sqrt(m_tot)))
m_up <- m_mean + (qt(0.975,9)*(sd(m_rand)/sqrt(m_tot)))

# Calculate upper and lower limits of 95% CI for other data
squid_m <- filter(size_maturity_ccs_time,sex=="male"&size_at_maturity_length_mm!="NA")
squid_f <- filter(size_maturity_ccs_time,sex=="female"&size_at_maturity_length_mm!="NA")
asquid_m <- filter(age_size_maturity_goa_ccs,sex=="male"&size_at_maturity_length_mm!="NA")
asquid_f <- filter(age_size_maturity_goa_ccs,sex=="female"&size_at_maturity_length_mm!="NA")
af_low <- mean(asquid_f$size_at_maturity_length_mm) + (qnorm(0.025)*(sd(asquid_f$size_at_maturity_length_mm)/sqrt(nrow(asquid_f))))
af_up <- mean(asquid_f$size_at_maturity_length_mm) + (qnorm(0.975)*(sd(asquid_f$size_at_maturity_length_mm)/sqrt(nrow(asquid_f))))
am_low <- mean(asquid_m$size_at_maturity_length_mm) + (qnorm(0.025)*(sd(asquid_m$size_at_maturity_length_mm)/sqrt(nrow(asquid_m))))
am_up <- mean(asquid_m$size_at_maturity_length_mm) + (qnorm(0.975)*(sd(asquid_m$size_at_maturity_length_mm)/sqrt(nrow(asquid_m))))
cf_low <- mean(squid_f$size_at_maturity_length_mm) + (qnorm(0.025)*(sd(squid_f$size_at_maturity_length_mm)/sqrt(nrow(squid_f))))
cf_up <- mean(squid_f$size_at_maturity_length_mm) + (qnorm(0.975)*(sd(squid_f$size_at_maturity_length_mm)/sqrt(nrow(squid_f))))
cm_low <- mean(squid_m$size_at_maturity_length_mm) + (qnorm(0.025)*(sd(squid_m$size_at_maturity_length_mm)/sqrt(nrow(squid_m))))
cm_up <- mean(squid_m$size_at_maturity_length_mm) + (qnorm(0.975)*(sd(squid_m$size_at_maturity_length_mm)/sqrt(nrow(squid_m))))

# Use 95% CI to calculate SD
f_sd <- sqrt(f_tot)*((f_up-f_low)/3.92)
m_sd <- sqrt(m_tot)*((m_up-m_low)/3.92)

# Generate data with characteristics of Fields (1965) dataset
# Normal distribution with upper and lower limits
# Source for function:https://stackoverflow.com/questions/19343133/setting-upper-and-lower-limits-in-rnorm
mysamp <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}

# Generate data
set.seed(1234)
dat <- data.frame(year = factor(rep("CCS 1948-1962",tot)),
                  sex = factor(rep(c("male","female"),c(m_tot,f_tot))),
                  size_at_maturity_length_mm = c(mysamp(n=m_tot, m=m_mean, s=m_sd, lwr=m_min, upr=m_max, nnorm=5000),
                          c(mysamp(n=f_tot, m=f_mean, s=f_sd, lwr=f_min, upr=f_max, nnorm=5000))))

# Select relevant data from CCS dataset
squid_comp <- size_maturity_ccs_time %>%
  select(size_at_maturity_length_mm,sex) %>%
  mutate(year = factor("CCS 2016-2019"))

# Same with GOA data
squid_comp1 <- age_size_maturity_goa_ccs %>%
  select(size_at_maturity_length_mm,sex) %>%
  mutate(year = factor("GOA 2016-2017"))

# Bind them altogether
dat <- bind_rows(dat,squid_comp,squid_comp1)

# Calcualte means and add these and 95% CIs to dataframe
library(plyr)
cdat <- ddply(dat,.variables = c("year","sex"),summarise,dml.mean=mean(size_at_maturity_length_mm,na.rm=T))
cdat[1:2,3] <- c(f_mean,m_mean)
cdat <- cdat %>%
  mutate(dml.up = c(f_up,m_up,cf_up,cm_up,af_up,am_up),
         dml.low = c(f_low,m_low,cf_low,cm_low,af_low,am_low),
         dml = dml.mean,
         ci = dml.up-dml.low)
cdat
#           year    sex dml.mean   dml.up  dml.low      dml       ci
#1 CCS 1948-1962 female 143.0000 143.8706 142.1294 143.0000 1.741178
#2 CCS 1948-1962   male 146.7000 147.5611 145.8389 146.7000 1.722217
#3 CCS 2016-2019 female 124.6120 125.4909 123.7330 124.6120 1.757961
#4 CCS 2016-2019   male 130.4806 131.4534 129.5077 130.4806 1.945733
#5 GOA 2016-2017 female 145.4333 149.7665 141.1002 145.4333 8.666302
#6 GOA 2016-2017   male 143.6078 148.3968 138.8188 143.6078 9.577995

# Plot figure (Figure 5E)
dat%>%
  ggplot(aes(x=size_at_maturity_length_mm, fill=year))+
  facet_wrap(facets = ~sex) +
  geom_density(alpha=.3) +
  geom_vline(data = cdat,aes(xintercept=dml.mean,colour=year), size=1,alpha = 0.7) +
  geom_vline(data = cdat,aes(xintercept=dml.up,colour=year), size=0.4,linetype = "dotted") +
  geom_vline(data = cdat,aes(xintercept=dml.low,colour=year), size=0.4,linetype = "dotted") +
  scale_color_manual(values = c("burlywood","red","darkblue"))+
  scale_fill_manual(values = c("burlywood","red","darkblue"))+
  scale_y_continuous(limits = c(0,0.06)) +
  scale_x_continuous(limits = c(65,200),breaks = seq(70,200,30)) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = c(0.51, .79)) +
  ylab("Density") +
  xlab("Size at maturity (dorsal mantle length, mm)")

ggsave("", height = 3, width = 7)

#----


## Flexibility in size and age at maturity, comparison with D. gigas
#----

# Size of D. gigas, reported by Hoving et al. 2013
# 84.7 cm max (Chile05), 21.8 cm min (GB10)
# also 17.63 cm min for males and 18.21 cm min for females in 2015 (Hoving et al. 2019)
(84.7-21.8)/21.8*100
# 288.5321

# Age of D. gigas, reported by Hoving et al. 2013
# 474 day max (Chile0) 183 day min (GB10)
(474-183)/183*100
# 159.0164

# Size of D. opalescens
# 183 mm max GOA squid size, 68 mm min CCS size
(183-68)/68*100
# 169.1176

# Age of D. opalescens
# 320 days max GOA age, 137 days min CCS age
(320-137)/137*100
# 133.5766

#----



  ######### Trophic ecology and energy density ######### 
library(tidyverse)
library(SIBER)
library(siar)
library(spatstat.utils)
library(viridis)
library(MixSIAR)
library(rjags)
library(doBy)
library(scales)


## Explore trophic overlap between juvenile salmon (Figure S4)
#----

# Load data
isotopes <- read_excel(file.list[2],sheet = "isotopes")

# Establish colors for plot
palette(c("grey20", "grey40", "grey60"))

# Build dataframe for SIBER with iso1, iso2, group, and community columns 
dopa.siber.data<-isotopes[,c("d13C","d15N", "species")]

# Filter for species of interest
dopa.siber.data <- dopa.siber.data %>%
  filter(species=="Hatchery Chum" | species=="Hatchery Coho" | species=="Wild Juvenile Salmon")

# Rename columns
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="d13C"] <- "iso1"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="d15N"] <- "iso2"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="species"] <- "group"

# Format columns and add community
dopa.siber.data$group <- as.factor(dopa.siber.data$group)
dopa.siber.data$group <- as.integer(dopa.siber.data$group)
dopa.siber.data$community<-1
dopa.siber.data<-as.data.frame(dopa.siber.data)
dopa.siber.data$iso1<-as.numeric(dopa.siber.data$iso1)
dopa.siber.data$iso2<-as.numeric(dopa.siber.data$iso2)

# Now create Siber object
my.siber.data <- createSiberObject(dopa.siber.data)
communityMetricsML(my.siber.data)
names(my.siber.data)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2) 
group.hull.args      <- list(lty = 2, col = "grey20")

# Build isotope plot with 95% PI elipses (Figure S4)
par(mfrow=c(1,1))
setwd("")
pdf("",height = 4.5,width = 4)
plotSiberObject(my.siber.data,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
dev.off()


# Calculate niche-space overlap between juvenile salmon
spx <- split(my.siber.data$original.data$iso1, 
             my.siber.data$original.data$group)
spy <- split(my.siber.data$original.data$iso2,
             my.siber.data$original.data$group)

# Names are integers (5,6,11), corresponding to order in main dataframe 
levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"      
#[3] "Calanoid large"          "California market squid"
#[5] "Hatchery Chum"           "Hatchery Coho"          
#[7] "Juvenile sablefish"      "Krill"                  
#[9] "Pacific herring"         "Pyrosoma atlanticum"    
#[11] "Wild Juvenile Salmon"

# overlap between chum and coho
overlap <- overlap(as.numeric(spx[[1]]), spy[[1]], 
                   spx[[2]], spy[[2]],
                   steps = 1)
overlap$overlap / overlap$area1 # 0 % of chum
overlap$overlap / overlap$area2 # 0 % of coho
overlap$overlap / (overlap$area1 + overlap$area2) # 0 % of combined area

# overlap between chum and pink
overlap <- overlap(as.numeric(spx[[1]]), spy[[1]], 
                   spx[[3]], spy[[3]],
                   steps = 1)
overlap$overlap / overlap$area1 # 0.1912843 % of chum
overlap$overlap / overlap$area2 # 0.3754528 % of pink
overlap$overlap / (overlap$area1 + overlap$area2) # 0.1267223 % of combined area

# overlap between coho and pink
overlap <- overlap(as.numeric(spx[[2]]), spy[[2]], 
                   spx[[3]], spy[[3]],
                   steps = 1)
overlap$overlap / overlap$area1 # 0 % of chum
overlap$overlap / overlap$area2 # 0 % of pink
overlap$overlap / (overlap$area1 + overlap$area2) # 0 % of combined area

# There is enough overlap between two of the species to warrant the combination of all three
isotopes <- isotopes %>%
  mutate(species = replace(species,species=="Hatchery Chum",c("Juvenile salmon")),
         species = replace(species,species=="Hatchery Coho",c("Juvenile salmon")),
         species = replace(species,species=="Wild Juvenile Salmon",c("Juvenile salmon")))

levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"       "Calanoid large"         
#[4] "California market squid" "Juvenile sablefish"      "Juvenile salmon"        
#[7] "Krill"                   "Pacific herring"         "Pyrosoma atlanticum"    

#----


## Calculate trophic overlap between squid and other spp. in GOA (Figure 6A, Figure 6B)
#----

# Isotope niche space overlap for each species
# Calculated using the SIBER Package https://github.com/AndrewLJackson/SIBER/ 

# Establish colors for plot
palette(viridis(9))

# Build dataframe for SIBER with iso1, iso2, group, and community columns 
dopa.siber.data<-isotopes[,c("d13C","d15N", "species")]
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="d13C"] <- "iso1"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="d15N"] <- "iso2"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="species"] <- "group"
dopa.siber.data$group <- as.integer(as.factor(dopa.siber.data$group))
dopa.siber.data$community<-1
dopa.siber.data$iso1<-as.numeric(dopa.siber.data$iso1)
dopa.siber.data$iso2<-as.numeric(dopa.siber.data$iso2)
dopa.siber.data<-as.data.frame(dopa.siber.data)
dopa.siber.data<-filter(dopa.siber.data,!is.na(iso1))

# Create Siber object 
my.siber.data <- createSiberObject(dopa.siber.data)
communityMetricsML(my.siber.data)
names(my.siber.data)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2) 
group.hull.args      <- list(lty = 2, col = "grey20")


## Plot isotopic niche space of studied species in GOA (Figure 6A)

# Build isotope plot with 95% PI elipses 
par(mfrow=c(1,1))
setwd("")
pdf("",height = 4.5,width = 4)
plotSiberObject(my.siber.data,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
dev.off()

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(my.siber.data)
print(group.ML)
#       squid    coho_a   herring  sablefish calanoid pyrosome    krill
#TA   3.329239 1.2669503 2.1463116 3.5483500 4.458200 2.2675000 5.877581
#SEA  1.037651 0.3635120 0.6602749 0.6800329 1.110441 0.6687012 2.274349
#SEAc 1.069095 0.3837071 0.6969568 0.7000338 1.145142 0.6909912 2.373234
#     juv. salmon king_a
#TA   15.181300  3.456583
#SEA   2.116458  1.031860
#SEAc  2.130474  1.083453


## Calculate ellipses overlap between squid and other spp
spx <- split(my.siber.data$original.data$iso1, 
             my.siber.data$original.data$group)
spy <- split(my.siber.data$original.data$iso2,
             my.siber.data$original.data$group)

# Names are integers (1-9), corresponding to order in main dataframe 
levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"      
#[3] "Calanoid large"          "California market squid"
#[5] "Juvenile sablefish"      "Juvenile salmon"        
#[7] "Krill"                   "Pacific herring"        
#[9] "Pyrosoma atlanticum"

# Overlap between squid and adult coho
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[1]], spy[[1]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of adult coho

# Overlap between squid and adult king
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[2]], spy[[2]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of adult king

# Overlap between squid and copepods
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[3]], spy[[3]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of copepods

# Overlap between squid and sablefish
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[5]], spy[[5]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0.5612575 % of sablefish

# Overlap between squid and juv salmon
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[6]], spy[[6]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of juv salmon

# Overlap between squid and krill
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[7]], spy[[7]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of krill

# Overlap between squid and herring
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[8]], spy[[8]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of herring

# Overlap between squid and pyrosomes
overlap <- overlap(as.numeric(spx[[4]]), spy[[4]], 
                   spx[[9]], spy[[9]],
                   steps = 1)
overlap$overlap / overlap$area2 # 0 % of pyrosomes


## Plot overlapping niche space of squid and sablefish (Figure 6B)
palette(c("#2C728EFF","#21908CFF"))

# Build dataframe for SIBER with iso1, iso2, group, and community columns; 
dopa.siber.data<-isotopes[,c("d13C","d15N", "species")]
dopa.siber.data <- dopa.siber.data %>%
  filter(species=="California market squid" | species=="Juvenile sablefish")
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="d13C"] <- "iso1"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="d15N"] <- "iso2"
colnames(dopa.siber.data)[colnames(dopa.siber.data)=="species"] <- "group"
dopa.siber.data$group <- as.integer(as.factor(dopa.siber.data$group))
dopa.siber.data$community<-1
dopa.siber.data$iso1 <- as.numeric(dopa.siber.data$iso1)
dopa.siber.data$iso2 <- as.numeric(dopa.siber.data$iso2)
dopa.siber.data<-as.data.frame(dopa.siber.data)

# Create Siber object 
my.siber.data <- createSiberObject(dopa.siber.data)
communityMetricsML(my.siber.data)
names(my.siber.data)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2) 
group.hull.args      <- list(lty = 2, col = "grey20")

# Build isotope plot with 95% PI elipses 
par(mfrow=c(1,1))
setwd("")
pdf("",height = 4.5,width = 4)
plotSiberObject(my.siber.data,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
dev.off()

#----


## Determine the contribution of potential prey sources to GOA squid diets across years (Figure 6C)
#----

# Mixing models will give slightly different answer each time
setwd("")


## Prepare data
# Gather and save isotope data for predator (squid) mixture
Dopa.mixtures<-subset(isotopes, species=="California market squid")
Dopa.mixtures<-dplyr::select(Dopa.mixtures,c(d13C,d15N))
Dopa.mixtures$d13C <- as.numeric(Dopa.mixtures$d13C)
Dopa.mixtures$d15N <- as.numeric(Dopa.mixtures$d15N)
write.table(Dopa.mixtures, file="Dopa.mixtures.csv",sep=",")

# Mixtures (Predator) data 
mixDopa <- load_mix_data(filename="Dopa.mixtures.csv",
                         iso_names=c("d13C","d15N"),
                         factors=NULL, # looking at all years together
                         fac_random=NULL, # looking at all years together
                         fac_nested=NULL,
                         cont_effects=NULL)

# Gather and save isotope data for source (prey) species
# Must be in this format: "species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", and "n" (sample size each source)
Dopa.sources<-isotopes[,c("species","d13C","d15N")]
Dopa.sources$d13C <- as.numeric(Dopa.sources$d13C)
Dopa.sources$d15N <- as.numeric(Dopa.sources$d15N)
Dopa.sources<-Dopa.sources%>%
  filter(!is.na(d13C))%>%
  group_by(species) %>%
  summarise(Meand13C = mean(d13C),
            SDd13C = sd(d13C),
            Meand15N = mean(d15N),
            SDd15N = sd(d15N),
            n = length(d15N)) %>%
  ungroup() %>%
  filter(species!="California market squid"&species!="Pyrosoma atlanticum"&
           species!="Adult king salmon"&species!="Adult coho salmon"&
           species!="Juvenile sablefish")
Dopa.sources
#        species  Meand13C    SDd13C  Meand15N    SDd15N   n
#3  Calanoid large -20.73294 0.5690477  7.291471 0.6222984  34
#6 Juvenile salmon -18.99863 0.7696736 10.779804 0.8754055 153
#7           Krill -18.66356 1.2965069  9.662960 0.6068332  25
#8 Pacific herring -20.58193 0.5305993 11.834000 0.5252208  20
write.table(Dopa.sources, file="Dopa.sources.csv",sep=",")

# Source (prey) data
sourceDopa <- load_source_data(filename="Dopa.sources.csv",
                               source_factors=NULL,
                               conc_dep=FALSE,
                               data_type="means",
                               mixDopa)   

# Set trophic enrichment factors
# Needs to be in format: Species,Meand13C,SDd13C,Meand15N,SDd15N):
Species <- c("Calanoid large","Juvenile salmon","Krill","Pacific herring")
discr.Dopa<- data.frame(Species)
discr.Dopa$Meand13C <- 1
discr.Dopa$SDd13C <- 0.5
discr.Dopa$Meand15N <- 3.4
discr.Dopa$SDd15N <- 0.5
discr.Dopa
#          Species Meand13C SDd13C Meand15N SDd15N
#1  Calanoid large        1    0.5      3.4    0.5
#2 Juvenile salmon        1    0.5      3.4    0.5
#3           Krill        1    0.5      3.4    0.5
#4 Pacific herring        1    0.5      3.4    0.5
write.table(discr.Dopa, file="discr.Dopa.csv",sep=",")

# Trophic enrichment factor data
discrDopa <- load_discr_data(filename="discr.Dopa.csv", mixDopa)


## Define and run model
# Define model structure and write JAGS model file
resid_err <- TRUE
process_err <- TRUE
model_Dopa <- "MixSIAR_model_Dopa_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Dopa, resid_err, process_err, mixDopa, sourceDopa)

# Run the JAGS model ("normal" with chain length 100,000; burn in 50,000, and thin 50; with 3 chains): 
jags.uninf.Dopa <- run_model(run="normal",mixDopa,sourceDopa,discrDopa,model_Dopa,alpha.prior = 1, resid_err, process_err) #took 25 min to run.

# Process diagnostics, summary stats, and posterior plots
jags_output_dopa<- output_JAGS(jags.uninf.Dopa, mixDopa, sourceDopa) # DIC = 122.6205

# Gather data for plotting
DopaSummary<-jags.uninf.Dopa$BUGSoutput$summary
DopaSummary<-DopaSummary[c(37:40),]
DopaSummary<-as.matrix(DopaSummary[,c(1:7)])
colnames(DopaSummary)<- c("Mean","SD","ymin","Lower","Middle","Upper","ymax")
DopaSummary <- as.data.frame(DopaSummary)
DopaSummary <- mutate(DopaSummary,species = c("Copepods","Juvenile salmon","Krill","Pacific herring"))


## Plot prey contributions to squid diet across years (Figure 6C)
# Establish colors for plot
spp_colors = palette(viridis(9))
levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"      
#[3] "Calanoid large"          "California market squid"
#[5] "Juvenile sablefish"      "Juvenile salmon"        
#[7] "Krill"                   "Pacific herring"        
#[9] "Pyrosoma atlanticum"

DopaSummary <- mutate(DopaSummary,species = c("Large copepods","Juvenile salmon","Krill","Pacific herring"))
DopaSummary %>%
  ggplot(aes(species, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity",
               fill=c(spp_colors[3],spp_colors[6],spp_colors[7],spp_colors[8]),
               color="black",width=.6) +
  xlab("") +
  ylab("Proportion of squid diet") +
  scale_x_discrete(limits = c("Large copepods","Juvenile salmon","Krill","Pacific herring")) +
  ylim(c(0,1)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("", height = 3.9, width = 3)

#----


## Repeat but within years
#----


## Prepare data
# Gather and save isotope data for predator (squid) mixture
Dopa.mixtures<-subset(isotopes, species=="California market squid")
Dopa.mixtures <- separate(Dopa.mixtures,capture_date,c("year","month","day"),sep=c("-"))
Dopa.mixtures<-dplyr::select(Dopa.mixtures,c(d13C,d15N,year))
Dopa.mixtures$d13C <- as.numeric(Dopa.mixtures$d13C)
Dopa.mixtures$d15N <- as.numeric(Dopa.mixtures$d15N)
Dopa.mixtures$year <- as.numeric(Dopa.mixtures$year)
write.table(Dopa.mixtures, file="Dopa.mixtures.csv",sep=",")

# Mixtures (Predator) data
mixDopa <- load_mix_data(filename="Dopa.mixtures.csv",
                         iso_names=c("d13C","d15N"),
                         factors="year", # looking at years separately
                         fac_random=T, # looking at years separately
                         fac_nested=NULL,
                         cont_effects=NULL)
# Gathering and saving isotope data for source (prey) species is the same so no need to repeat


## Define and run model
# Define model structure and write JAGS model file
resid_err <- TRUE
process_err <- TRUE
model_Dopa <- "MixSIAR_model_Dopa_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Dopa, resid_err, process_err, mixDopa, sourceDopa)

# Run the JAGS model ("normal" with chain length 100,000; burn in 50,000, and thin 50; with 3 chains): 
jags.uninf.Dopa <- run_model(run="normal",mixDopa,sourceDopa,discrDopa,model_Dopa,alpha.prior = 1, resid_err, process_err) #took 25 min to run.

# Process diagnostics, summary stats, and posterior plots
jags_output_dopa<- output_JAGS(jags.uninf.Dopa, mixDopa, sourceDopa) # DIC = 92.91318

# Gather data for plotting
DopaSummary_year<-jags.uninf.Dopa$BUGSoutput$summary
DopaSummary_year<-DopaSummary_year[c(47:62),]
DopaSummary_year<-as.matrix(DopaSummary_year[,c(1:7)])
colnames(DopaSummary_year)<- c("Mean","SD","ymin","Lower","Middle","Upper","ymax")
DopaSummary_year <- as.data.frame(DopaSummary_year)
DopaSummary_year <- DopaSummary_year %>%
  mutate(species=c("Copepods","Copepods","Copepods",
                   "Juvenile salmon","Juvenile salmon","Juvenile salmon",
                   "Krill","Krill","Krill",
                   "Pacific herring","Pacific herring","Pacific herring",
                   "Copepods","Juvenile salmon","Krill","Pacific herring"),
         year=c(2016,2017,2018,
                2016,2017,2018,
                2016,2017,2018,
                2016,2017,2018,
                1,1,1,1))
DopaSummary_year$Middle<-as.numeric(DopaSummary_year$Middle)


## Plot prey contributions to squid diet within years
DopaSummary_year <- mutate(DopaSummary_year,
                           species=c("Large copepods","Large copepods","Large copepods",
                                     "Juvenile salmon","Juvenile salmon","Juvenile salmon",
                                     "Krill","Krill","Krill",
                                     "Pacific herring","Pacific herring","Pacific herring",
                                     "Large copepods","Juvenile salmon","Krill","Pacific herring"))
dodge <- position_dodge(width=0.9)
year_cols <- c("grey20", "grey40", "grey60")

DopaSummary_year %>%
  filter(year!=1) %>%
  ggplot(aes(species, Mean, alpha=factor(year), fill = species)) +
  geom_col(position = dodge) +
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD), position = dodge, width = 0) +
  scale_fill_manual(values = c(spp_colors[6],spp_colors[7],spp_colors[3],spp_colors[8])) +
  scale_color_manual(values = c(spp_colors[6],spp_colors[7],spp_colors[3],spp_colors[8])) +
  scale_alpha_manual(values = c(1,.7,.4)) +
  xlab("") +
  ylab("Proportion of squid diet") +
  scale_y_continuous(limits = c(0,1.15), breaks = seq(0,1,.25)) +
  theme_classic(base_size = 12) +
  guides(fill = FALSE) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, .6),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("", height = 3.9, width = 5)

#----


## Compare contribution of squid with other prey sources to adult Chinook salmon diets (Figure 6D)
#----


## Prepare data
# Gather and save isotope data for predator (squid) mixture
Dopa.mixtures<-subset(isotopes, species=="Adult king salmon")
Dopa.mixtures<-dplyr::select(Dopa.mixtures,c(d13C,d15N))
Dopa.mixtures$d13C <- as.numeric(Dopa.mixtures$d13C)
Dopa.mixtures$d15N <- as.numeric(Dopa.mixtures$d15N)
Dopa.mixtures <- filter(Dopa.mixtures, !is.na(d13C))
write.table(Dopa.mixtures, file="Dopa.mixtures.csv",sep=",")

# Mixtures (Predator) data set: 
mixDopa <- load_mix_data(filename="Dopa.mixtures.csv",
                         iso_names=c("d13C","d15N"),
                         factors=NULL,
                         fac_random=NULL,
                         fac_nested=NULL,
                         cont_effects=NULL)  

# Filter to only include squid collected 2017 and beyond
Dopa1 <- isotopes %>%
  separate(capture_date,c("year","month","day"),sep=c("-")) %>%
  mutate(year = as.numeric(year)) %>%
  filter(species!="California market squid" | year!=2016)

# Gather and save isotope data for source (prey) species
# Must be in this format: "species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", and "n" (sample size each source)
Dopa.sources<-Dopa1[,c("species","d13C","d15N")]
Dopa.sources$d13C <- as.numeric(Dopa.sources$d13C)
Dopa.sources$d15N <- as.numeric(Dopa.sources$d15N)
Dopa.sources<-Dopa.sources%>%
  filter(!is.na(d13C))%>%
  group_by(species) %>%
  summarise(Meand13C = mean(d13C),
            SDd13C = sd(d13C),
            Meand15N = mean(d15N),
            SDd15N = sd(d15N),
            n = length(d15N)) %>%
  ungroup() %>%
  filter(species!="Adult king salmon"&species!="Pyrosoma atlanticum"&
           species!="Adult coho salmon"&species!="Juvenile salmon"&
           species!="Juvenile sablefish")
Dopa.sources
#                species  Meand13C    SDd13C  Meand15N    SDd15N  n
#3          Calanoid large -20.73294 0.5690477  7.291471 0.6222984 34
#4 California market squid -18.87443 0.6268073 13.102400 0.5350956 25
#7                   Krill -18.66356 1.2965069  9.662960 0.6068332 25
#8         Pacific herring -20.58193 0.5305993 11.834000 0.5252208 20
write.table(Dopa.sources, file="Dopa.sources.csv",sep=",")

# Source (prey) data
sourceDopa <- load_source_data(filename="Dopa.sources.csv",
                               source_factors=NULL,
                               conc_dep=FALSE,
                               data_type="means",
                               mixDopa)   

# Set trophic enrichment factors
# Needs to be in format: Species,Meand13C,SDd13C,Meand15N,SDd15N
Species <- c("Calanoid large","California market squid","Krill","Pacific herring")
discr.Dopa<- data.frame(Species)
discr.Dopa$Meand13C <- 1
discr.Dopa$SDd13C <- 0.5
discr.Dopa$Meand15N <- 3.4
discr.Dopa$SDd15N <- 0.5
discr.Dopa
#                  Species Meand13C SDd13C Meand15N SDd15N
#1          Calanoid large        1    0.5      3.4    0.5
#2 California market squid        1    0.5      3.4    0.5
#3                   Krill        1    0.5      3.4    0.5
#4         Pacific herring        1    0.5      3.4    0.5
write.table(discr.Dopa, file="discr.Dopa.csv",sep=",")

# Trophic enrichment factor data
discrDopa <- load_discr_data(filename="discr.Dopa.csv", mixDopa)


## Define and run model
# Define model structure and write JAGS model file
resid_err <- TRUE
process_err <- TRUE
model_Dopa <- "MixSIAR_model_Dopa_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Dopa, resid_err, process_err, mixDopa, sourceDopa)

# Run the JAGS model ("normal" with chain length 100,000; burn in 50,000, and thin 50; with 3 chains): 
jags.uninf.Dopa <- run_model(run="normal",mixDopa,sourceDopa,discrDopa,model_Dopa,alpha.prior = 1, resid_err, process_err) #took 25 min to run.

# Process diagnostics, summary stats, and posterior plots
jags_output_dopa<- output_JAGS(jags.uninf.Dopa, mixDopa, sourceDopa) #DIC=123.1482

# Gather data for plotting
KingSummary<-jags.uninf.Dopa$BUGSoutput$summary
KingSummary<-KingSummary[c(24:27),]
KingSummary<-as.matrix(KingSummary[,c(1:7)])
colnames(KingSummary)<- c("mean","sd","ymin","Lower","Middle","Upper","ymax")
KingSummary <- as.data.frame(KingSummary)
KingSummary$Middle<-as.numeric(KingSummary$Middle)
KingSummary <- mutate(KingSummary, species = c("Copepods","California market squid","Krill","Pacific herring"))


## Plot prey contributions to chinook diet across years (Figure 6D)
# Align colors for plot 
levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"      
#[3] "Calanoid large"          "California market squid"
#[5] "Juvenile sablefish"      "Juvenile salmon"        
#[7] "Krill"                   "Pacific herring"        
#[9] "Pyrosoma atlanticum"

KingSummary <- mutate(KingSummary, species = c("Large copepods","California market squid","Krill","Pacific herring"))

KingSummary %>%
  ggplot(aes(species, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity",
               fill=c(spp_colors[3],spp_colors[4],spp_colors[7],spp_colors[8]),
               color="black",width=.6) +
  xlab("") +
  ylab("Proportion of Chinook salmon diet") +
  scale_x_discrete(limits = c("Large copepods","California market squid","Krill","Pacific herring")) +
  ylim(c(0,1)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("", height = 3.9, width = 3)

#----


## Repeat but for adult coho salmon (Figure 6E)
#----


## Prepare data
# Gather and save isotope data for predator (adult coho salmon) mixture
# Gather and save isotope data for predator (squid) mixture
Dopa.mixtures<-subset(isotopes, species=="Adult coho salmon")
Dopa.mixtures<-dplyr::select(Dopa.mixtures,c(d13C,d15N))
Dopa.mixtures$d13C <- as.numeric(Dopa.mixtures$d13C)
Dopa.mixtures$d15N <- as.numeric(Dopa.mixtures$d15N)
Dopa.mixtures <- filter(Dopa.mixtures, !is.na(d13C))
write.table(Dopa.mixtures, file="Dopa.mixtures.csv",sep=",")

# Mixtures (Predator) data set: 
mixDopa <- load_mix_data(filename="Dopa.mixtures.csv",
                         iso_names=c("d13C","d15N"),
                         factors=NULL,
                         fac_random=NULL,
                         fac_nested=NULL,
                         cont_effects=NULL)  
# Gathering and saving isotope data for source (prey) species is the same so no need to repeat


## Define and run model
# Define model structure and write JAGS model file
resid_err <- TRUE
process_err <- TRUE
model_Dopa <- "MixSIAR_model_Dopa_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Dopa, resid_err, process_err, mixDopa, sourceDopa)

# Run the JAGS model ("normal" with chain length 100,000; burn in 50,000, and thin 50; with 3 chains): 
jags.uninf.Dopa <- run_model(run="normal",mixDopa,sourceDopa,discrDopa,model_Dopa,alpha.prior = 1, resid_err, process_err) #took 25 min to run.

# Process diagnostics, summary stats, and posterior plots
jags_output_dopa<- output_JAGS(jags.uninf.Dopa, mixDopa, sourceDopa) #DIC=123.1482

# Gather data for plotting
CohoSummary<-jags.uninf.Dopa$BUGSoutput$summary
CohoSummary<-CohoSummary[c(22:25),]
CohoSummary<-as.matrix(CohoSummary[,c(1:7)])
colnames(CohoSummary)<- c("mean","sd","ymin","Lower","Middle","Upper","ymax")
CohoSummary <- as.data.frame(CohoSummary)
CohoSummary$Middle<-as.numeric(CohoSummary$Middle)
CohoSummary <- mutate(CohoSummary, species = c("Copepods","California market squid","Krill","Pacific herring"))


## Plot prey contributions to coho diet across years (Figure 6E)
# Align colors for plot 
levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"      
#[3] "Calanoid large"          "California market squid"
#[5] "Juvenile sablefish"      "Juvenile salmon"        
#[7] "Krill"                   "Pacific herring"        
#[9] "Pyrosoma atlanticum"

CohoSummary <- mutate(CohoSummary,species=c("Large copepods","California market squid","Krill","Pacific herring"))

CohoSummary %>%
  ggplot(aes(species, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity",
               fill=c(spp_colors[3],spp_colors[4],spp_colors[7],spp_colors[8]),
               color="black",width=.6) +
  xlab("") +
  ylab("Proportion of coho salmon diet") +
  scale_x_discrete(limits = c("Large copepods","California market squid","Krill","Pacific herring")) +
  ylim(c(0,1)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("", height = 3.9, width = 3)

#----


## Compare energy content of GOA squid with that of common forage fish (Figure 7A)
#----

# Load data
energy <- read_excel(file.list[1],sheet = "energy_density_proximal_comp")

# Compare energy density between juvenile salmon spp - can they be combined?
summary(aov(energy_density_kJ_g_dry_mass~species*length_mm, data = filter(energy,species=="Juvenile coho salmon"|species=="Juvenile chum salmon")))
#                Df Sum Sq Mean Sq F value   Pr(>F)    
# species         1 22.836  22.836  29.214 6.72e-06 ***
# length_mm          1  2.812   2.812   3.598   0.0672 .  
# species:length_mm  1  2.686   2.686   3.436   0.0733 .  
# Residuals      31 24.231   0.782

# Cant be combined so keep separate

# Compare energy density between the squid and the four forage species
# Check normality of energy density for all species
test <- filter(energy, species=="Juvenile coho salmon")
shapiro.test(test$energy_density_kJ_g_dry_mass)
#W = 0.94922, p-value = 0.4129
test <- filter(energy, species=="Juvenile chum salmon")
shapiro.test(test$energy_density_kJ_g_dry_mass)
#W = 0.91687, p-value = 0.1307
test <- filter(energy, species=="Juvenile sablefish")
shapiro.test(test$energy_density_kJ_g_dry_mass)
#W = 0.96674, p-value = 0.3434
test <- filter(energy, species=="California market squid")
shapiro.test(test$energy_density_kJ_g_dry_mass)
#W = 0.9047, p-value = 0.09562
test <- filter(energy, species=="Pacific herring")
shapiro.test(test$energy_density_kJ_g_dry_mass)
#W = 0.9452, p-value = 0.3547
rm(test)

# Make species a factor, then organize with respect to squid
energy$species <- as.factor(energy$species)
energy <- within(energy, species <- relevel(species, ref = 'California market squid'))

# ANCOVA of energy density ~ species with size as a covariate, include all interactions
summary(aov(energy_density_kJ_g_dry_mass~species*length_mm, data = energy))
#                Df Sum Sq Mean Sq F value   Pr(>F)    
# species         4  43.33  10.832  11.473 1.22e-07 ***
# length_mm          1   2.46   2.457   2.603  0.10999    
# species:length_mm  4  24.04   6.009   6.365  0.00014 ***
# Residuals      95  89.69   0.944

# Can't do ANCOVA b/c violates non-homogeneity of regression slopes
# Can't do mixed effects model because continuous covariate
# Need to do multiple regression to allow for non-homogeneity of slopes

# Explore the significance of covariation between energy density and size
summary(lm(energy_density_kJ_g_dry_mass~length_mm,filter(energy,species=="Juvenile coho salmon")))
# Estimate Std. Error t value Pr(>|t|)
# 0.00112    0.02870   0.039    0.969

summary(lm(energy_density_kJ_g_dry_mass~length_mm,filter(energy,species=="Juvenile chum salmon")))
# Estimate Std. Error t value Pr(>|t|)
# 0.06597    0.01950   3.382   0.0041 **

summary(lm(energy_density_kJ_g_dry_mass~length_mm,filter(energy,species=="Juvenile sablefish")))
# Estimate Std. Error t value Pr(>|t|)
# 0.025372   0.007098   3.574  0.00108 **

summary(lm(energy_density_kJ_g_dry_mass~length_mm,filter(energy,species=="California market squid")))
# Estimate Std. Error t value Pr(>|t|)
# -0.01245    0.00732  -1.701    0.111

summary(lm(energy_density_kJ_g_dry_mass~length_mm,filter(energy,species=="Pacific herring")))
# Estimate Std. Error t value Pr(>|t|)
# -0.02221    0.01110    -2.0   0.0628 .

# Plot non-homogeneity of covariation between size and energy density with only significant regression lines (Figure 7A)
levels(as.factor(isotopes$species))
#[1] "Adult coho salmon"       "Adult king salmon"       "Calanoid large"         
#[4] "California market squid" "Juvenile sablefish"      "Juvenile salmon"        
#[7] "Krill"                   "Pacific herring"         "Pyrosoma atlanticum"
levels(as.factor(energy$species))
#[1] "California market squid" "Juvenile chum salmon"    "Juvenile coho salmon"   
#[4] "Juvenile sablefish"      "Pacific herring" 
spp_colors <- palette(viridis(9))
energy %>% 
  ggplot(aes(length_mm,energy_density_kJ_g_dry_mass,col=species,shape=species)) + 
  geom_point() + 
  geom_smooth(aes(length_mm,energy_density_kJ_g_dry_mass,col=species),data = filter(energy,species=="Juvenile sablefish"),method="lm",se=F) +
  geom_smooth(aes(length_mm,energy_density_kJ_g_dry_mass,col=species),data = filter(energy,species=="Juvenile chum salmon"),method="lm",se=F) +
  scale_shape_manual(values = c(16,16,3,16,16)) +
  scale_color_manual(values = c(spp_colors[5],spp_colors[6],spp_colors[6],spp_colors[4],spp_colors[8])) +
  scale_y_continuous(limits = c(19,27), breaks = seq(19,27,2)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "") +
  ylab(expression("Energy density (kJ g dry"~mass^-1*")")) +
  xlab("Size (length, mm)")

ggsave("", height = 3, width = 4)

# Covariate interacts with independent variable, so run multiple regression with significant covariates/interactions 
# Store model object for subsequent calculations
mr1 <- lm(energy_density_kJ_g_dry_mass~species
          , data = energy)
summary(mr1)
#                                        Estimate Std. Error t value Pr(>|t|)    
#1(Intercept)                            23.847804   2.746401   8.683 1.07e-13 ***
#2speciesJuvenile chum salmon            -6.134130   3.169684  -1.935  0.05593 .  
#3speciesJuvenile coho salmon            -0.881465   4.346423  -0.203  0.83972    
#4speciesJuvenile sablefish             -10.008769   3.567020  -2.806  0.00609 ** 
#5speciesPacific herring                  4.212431   3.240700   1.300  0.19680    
#6length_mm                              -0.012452   0.021362  -0.583  0.56132    
#7speciesJuvenile chum salmon:length_mm   0.078421   0.034696   2.260  0.02609 *  
#8speciesJuvenile coho salmon:length_mm   0.013572   0.034454   0.394  0.69452    
#9speciesJuvenile sablefish:length_mm     0.037824   0.022352   1.692  0.09389 .  
#10speciesPacific herring:length_mm       -0.009756   0.022944  -0.425  0.67164

# Adjusted R-squared:  0.3845 
# F-statistic: 8.217 on 9 and 95 DF,  p-value: 6.036e-09

## Calculate difference between energy density at average size of squid in dataset (128 mm)

# Average size of squid in dataset (128 mm)
size <- energy %>%
  filter(species=="California market squid") %>%
  group_by(species) %>%
  summarise(size = mean(length_mm)) %>%
  ungroup()

# Function from model coefficients to calculate energy density at average size
energy_1 <- function(size){
  mr1$coefficients[1]+(mr1$coefficients[6]*size)
}

# Replace with actual value
energy_1 <- as.numeric(energy_1(size$size))
energy_1
# 22.25313

# Now repeat this process for each species, using energy-length_mm relationships

# Juvenile sablefish
energy_2 <- function(size){
  (mr1$coefficients[1]+mr1$coefficients[4])+((mr1$coefficients[9]+mr1$coefficients[6])*size)
}
energy_2 <- as.numeric(energy_2(size$size))
energy_2
# 17.08818
energy_1-energy_2
# 5.164945

# Pacific herring
energy_3 <- function(size){
  (mr1$coefficients[1]+mr1$coefficients[5])+((mr1$coefficients[10]+mr1$coefficients[6])*size)
}
energy_3 <- as.numeric(energy_3(size$size))
energy_3
# 25.21616
energy_1-energy_3
# -2.963031

# Juvenile coho
energy_4 <- function(size){
  (mr1$coefficients[1]+mr1$coefficients[3])+((mr1$coefficients[8]+mr1$coefficients[6])*size)
}
energy_4 <- as.numeric(energy_4(size$size))
energy_4
# 23.10973
energy_1-energy_4
# -0.8566061

# Juvenile chum
energy_5 <- function(size){
  (mr1$coefficients[1]+mr1$coefficients[2])+((mr1$coefficients[7]+mr1$coefficients[6])*size)
}
energy_5 <- as.numeric(energy_5(size$size))
energy_5
# 26.16179
energy_1-energy_5
# -3.908666

#----


## Compare lipid and protein content between the squid, sablefish, and herring (Figure 7B)
#----

# Protein and lipid content of squid
squid_sum <- energy %>%
  filter(species=="California market squid") %>%
  mutate(protein_perc_wet_mass = as.numeric(protein_perc_wet_mass),
         lipid_perc_wet_mass = as.numeric(lipid_perc_wet_mass)) %>%
  filter(!is.na(protein_perc_wet_mass) | !is.na(lipid_perc_wet_mass)) %>%
  mutate(n = 1) %>%
  group_by(species) %>%
  summarise(protein_m = mean(protein_perc_wet_mass,na.rm=T),
            protein_sd = sd(protein_perc_wet_mass,na.rm=T),
            squid_count = sum(n),
            squid_length_m = mean(length_mm),
            squid_length_sd = sd(length_mm),
            lipid_sd = sd(lipid_perc_wet_mass,na.rm=T),
            lipid_m = mean(lipid_perc_wet_mass,na.rm=T)) %>%
  ungroup()

# Squid lipid
dl_mean <- squid_sum$lipid_m
dl_sd <- squid_sum$lipid_sd
dl_tot <- squid_sum$squid_count
# Squid protein
dp_mean <- squid_sum$protein_m
dp_sd <- squid_sum$protein_sd
dp_tot <- squid_sum$squid_count

# Data from Vollenweider et al. (2011) of protein and lipid content

# Herring lipid
hl_mean <- 7.8 # Appendix A - % Lipid May
hl_sd <- 0.6 # Appendix A - % Lipid May
hl_tot <- 17+21+36 # Table 2 - May 2001, 2002, 2003
# Herring protein
hp_mean <- 16.3 # Appendix A - % Protein May
hp_sd <- 0.7 # Appendix A - % Protein May
hp_tot <- 17+21+36 # Table 2 - May 2001, 2002, 2003

# Sablefish lipid
sl_mean <- 15.5 # Appendix A - % Lipid May
sl_sd <- 6.6 # Appendix A - % Lipid May
sl_tot <- 17 # Table 2 - May 2004
# Sablefish protein
sp_mean <- 13.3 # Appendix A - % Protein May
sp_sd <- 0.6 # Appendix A - % Protein May
sp_tot <- 17 # Table 2 - May 2004

# Calculate upper and lower limits of 95% CI
# using t-didtribution b/c small sample size (df = n-1)

# Squid
dl_low <- dl_mean + (qt(0.025,(dl_tot-1))*(dl_sd/sqrt(dl_tot)))
dl_up <- dl_mean + (qt(0.975,(dl_tot-1))*(dl_sd/sqrt(dl_tot)))
dp_low <- dp_mean + (qt(0.025,(dp_tot-1))*(dp_sd/sqrt(dp_tot)))
dp_up <- dp_mean + (qt(0.975,(dp_tot-1))*(dp_sd/sqrt(dp_tot)))
# Herring
hl_low <- hl_mean + (qt(0.025,(hl_tot-1))*(hl_sd/sqrt(hl_tot)))
hl_up <- hl_mean + (qt(0.975,(hl_tot-1))*(hl_sd/sqrt(hl_tot)))
hp_low <- hp_mean + (qt(0.025,(hp_tot-1))*(hp_sd/sqrt(hp_tot)))
hp_up <- hp_mean + (qt(0.975,(hp_tot-1))*(hp_sd/sqrt(hp_tot)))
# Sablefish
sl_low <- sl_mean + (qt(0.025,(sl_tot-1))*(sl_sd/sqrt(sl_tot)))
sl_up <- sl_mean + (qt(0.975,(sl_tot-1))*(sl_sd/sqrt(sl_tot)))
sp_low <- sp_mean + (qt(0.025,(sp_tot-1))*(sp_sd/sqrt(sp_tot)))
sp_up <- sp_mean + (qt(0.975,(sp_tot-1))*(sp_sd/sqrt(sp_tot)))

# Join together in dataframe
proximal <- data.frame(species = c(rep("California market squid",2),
                                   rep("Juvenile sablefish",2),
                                   rep("Pacific herring",2)),
                       type = rep(c("lipid","protein"),3),
                       mean = c(dl_mean,dp_mean,sl_mean,sp_mean,hl_mean,hp_mean),
                       up.ci = c(dl_up,dp_up,sl_up,sp_up,hl_up,hp_up),
                       low.ci = c(dl_low,dp_low,sl_low,sp_low,hl_low,hp_low)) 

# Plot (Figure 7B)
proximal %>%
  ggplot(aes(species,mean,color = species)) +
  facet_wrap(facets = ~type) +
  geom_point() +
  geom_errorbar(aes(ymin=low.ci, ymax=up.ci), width=.2) +
  scale_color_manual(values = c(spp_colors[5],spp_colors[4],spp_colors[8])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0,20), breaks = seq(0,20,5)) +
  ylab(c("Mean % of wet mass ± 95 % CI")) +
  xlab(c(""))

ggsave("", height = 3.5, width = 6)

#----


## All done!


