library(mgcv)
library(tidyverse)
library(lubridate)
library(mgcViz)
library(viridis)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(gratia)
library(gganimate)
library(gifski)
library(zoo)
library(ggsci)
library(ggbreak)
library(rgdal)
library(scales)
library(sf)
library(itsadug)

Silalatshani_df <- readRDS("GAMM_data\\Silalatshani_df.rds")

# Columns 
Silalatshani_df$sensor <- as.factor(Silalatshani_df$sensor)
Silalatshani_df <- Silalatshani_df %>% 
  mutate_at(vars(date), dplyr::funs(year = lubridate::year, month = lubridate::month, day =lubridate::day))


# Clean data
Silalatshani_df_nona <- Silalatshani_df[!is.na(Silalatshani_df$GPP) & (Silalatshani_df$GPP > 0) & 
                                !is.na(Silalatshani_df$ETa) & (Silalatshani_df$ETa > 0) &
                                !is.na(Silalatshani_df$Ea) & (Silalatshani_df$Ea > 0) & 
                                !is.na(Silalatshani_df$Ta) & (Silalatshani_df$Ta > 0) &
                                (Silalatshani_df$year > 2012),]

#Silalatshani_df_nona <- Silalatshani_df_nona %>% mutate(GPP = replace(GPP, GPP>15, NA))


# List of formulae with different time lags for rainfall (30, 60, 90, 120 day)
GPP_f_30 <- log(GPP) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y, k=50, bs='gp') +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("gp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_30, bs="cr", k=30) 

GPP_f_60 <- log(GPP) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y, k=50, bs='gp') +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("gp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_60, bs="cr", k=30) 

GPP_f_90 <- log(GPP) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y, k=50, bs='gp') +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("gp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_90, bs="cr", k=30) 

GPP_f_120 <- log(GPP) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y, k=50, bs='gp') +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("gp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_120, bs="cr", k=30)

GPP_fs <- list(GPP_f_30, GPP_f_60, GPP_f_90, GPP_f_120)

#Loop through different lags and test with anova
GPP_Silalatshani_GAM1 <- bam(GPP_f_30, discrete=TRUE, nthreads=8, data=Silalatshani_df_nona)
GPP_Silalatshani_GAM2 <- bam(GPP_f_60, discrete=TRUE, nthreads=8, data=Silalatshani_df_nona)
GPP_Silalatshani_GAM3 <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Silalatshani_df_nona)
GPP_Silalatshani_GAM4 <- bam(GPP_f_120, discrete=TRUE, nthreads=8, data=Silalatshani_df_nona)

anova(GPP_Silalatshani_GAM1, GPP_Silalatshani_GAM2, GPP_Silalatshani_GAM3, GPP_Silalatshani_GAM4, test='F')

# Run GAMs for biomass (GPP)
GPP_Silalatshani_GAM <- GPP_Silalatshani_GAM3

resids_Silalatshani_GPP <- residuals.gam(GPP_Silalatshani_GAM)
valRho_Silalatshani_GPP <- acf(resids_Silalatshani_GPP, plot=FALSE)$acf[2]
GPP_Silalatshani_GAM <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Silalatshani_df_nona, rho = valRho_Silalatshani_GPP)
GPP_Silalatshani_GAM1_vis <- getViz(GPP_Silalatshani_GAM)

check(GPP_Silalatshani_GAM1_vis)

GPP_Silalatshani_year_plot <- plot_smooth(GPP_Silalatshani_GAM, view = 'year', transform = exp)

##
GPP_Silalatshani_year_plot_gg <- ggplot(data=GPP_Silalatshani_year_plot$fv, aes(x=year, y=fit))+
  geom_line(aes(col="GPP"))+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y=expression(paste("gC/m"^2,paste("/day"))))+theme_minimal() +
  scale_color_manual(values=c("black")) + theme(legend.title=element_blank())

GPP_Silalatshani_space_plot <- plot(sm(GPP_Silalatshani_GAM1_vis, 3), trans = function(x){ 
  GPP_Silalatshani_GAM$family$linkinv(coef(GPP_Silalatshani_GAM)[1] + x) 
})
GPP_Silalatshani_space_plot_gg <- ggplot(data=GPP_Silalatshani_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name=expression(paste("GPP gC/m"^2)))+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()


##

GPP_Silalatshani_yearxspace <- plotSlice(sm(GPP_Silalatshani_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

GPP_Silalatshani_yearxspace_plot <- ggplot()+geom_raster(data=GPP_Silalatshani_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=GPP_Silalatshani_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="GPP (gC/m2)")+coord_equal()+guides(alpha="none", colour="none")

GPP_Silalatshani_yearxmonth <- plot(sm(GPP_Silalatshani_GAM1_vis, 4), trans = function(x){ 
  GPP_Silalatshani_GAM$family$linkinv(coef(GPP_Silalatshani_GAM)[1] + x) 
}) 

GPP_Silalatshani_yearxmonth_plot <- ggplot(data=GPP_Silalatshani_yearxmonth$data$fit, aes(x=y, y=x, fill=tz)) +
  geom_raster()+
  geom_contour(aes(x=y, y=x, z=tz, col="red"), show.legend = FALSE)+
  scale_fill_viridis_c(labels=label_wrap(5))+
  scale_y_continuous(breaks=seq(2013, 2021, by=1))+
  scale_x_continuous(breaks=seq(1, 12, by=1))+
  labs(y="Year", x="Month", fill="GPP (gC/m2)") + theme_bw()

# Run GAMs for WUE, evapotranspiration (ETa), evaporation (Ea), transpiration (Ta) 
# Stick with 90 day lag

Silalatshani_df_nona$WUE <- Silalatshani_df_nona$GPP / Silalatshani_df_nona$ETa

WUE_f_90 <- log(WUE) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y) +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("tp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_90, bs="cr", k=30)

ETa_f_90 <- log(ETa) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y) +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("tp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_90, bs="cr", k=30)

Ta_f_90 <- log(Ta) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y) +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("tp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_90, bs="cr", k=30)

Ea_f_90 <- log(Ea) ~ 
  # smooth term for year
  s(year, bs="cr", k=9) + 
  # cyclic term for season
  s(month, bs="cc", k=12) +  
  # smooth term for spatial interaction
  s(x,y) +                  
  # seasonal within year
  ti(month, year, bs = c("cc", "cr"), k = c(12,9)) +  
  # space x time
  ti(x, y, year, d = c(2, 1), bs = c("tp", "cr"), 
     k = c(50,9)) +  
  # random effect for ls
  #s(sensor, bs = 're') + 
  # rainfall lag
  s(chirps_90, bs="cr", k=30)

WUE_Silalatshani_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona)
resids_Silalatshani_WUE <- residuals.gam(WUE_Silalatshani_GAM)
valRho_Silalatshani_WUE <- acf(resids_Silalatshani_WUE , plot=FALSE)$acf[2]
WUE_Silalatshani_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona, rho = valRho_Silalatshani_WUE)

WUE_Silalatshani_GAM1_vis <- getViz(WUE_Silalatshani_GAM)

ETa_Silalatshani_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona)
resids_Silalatshani_ETa <- residuals.gam(ETa_Silalatshani_GAM)
valRho_Silalatshani_ETa <- acf(resids_Silalatshani_ETa , plot=FALSE)$acf[2]
ETa_Silalatshani_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona, rho = valRho_Silalatshani_ETa)

ETa_Silalatshani_GAM1_vis <- getViz(ETa_Silalatshani_GAM)

Ta_Silalatshani_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona)
resids_Silalatshani_Ta <- residuals.gam(Ta_Silalatshani_GAM)
valRho_Silalatshani_Ta <- acf(resids_Silalatshani_Ta , plot=FALSE)$acf[2]
Ta_Silalatshani_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona, rho = valRho_Silalatshani_Ta)

Ta_Silalatshani_GAM1_vis <- getViz(Ta_Silalatshani_GAM)

Ea_Silalatshani_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona)
resids_Silalatshani_Ea <- residuals.gam(Ea_Silalatshani_GAM)
valRho_Silalatshani_Ea <- acf(resids_Silalatshani_Ea , plot=FALSE)$acf[2]
Ea_Silalatshani_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Silalatshani_df_nona, rho = valRho_Silalatshani_Ea)

Ea_Silalatshani_GAM1_vis <- getViz(Ea_Silalatshani_GAM)

WUE_Silalatshani_year_plot <- plot_smooth(WUE_Silalatshani_GAM, view = 'year', transform = exp)

WUE_Silalatshani_year_plot_gg <- ggplot(data=WUE_Silalatshani_year_plot$fv, aes(x=year, y=fit))+
  geom_line()+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y="WUE gC/mm")+theme_minimal()

##
WUE_Silalatshani_space_plot <- plot(sm(WUE_Silalatshani_GAM1_vis, 3), trans = function(x){ 
  WUE_Silalatshani_GAM$family$linkinv(coef(WUE_Silalatshani_GAM)[1] + x) 
})
WUE_Silalatshani_space_plot_gg <- ggplot(data=WUE_Silalatshani_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name="WUE gC/mm")+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()
##

ETa_Silalatshani_year_plot <- plot_smooth(ETa_Silalatshani_GAM, view = 'year', transform = exp)

##
ETa_Silalatshani_space_plot <- plot(sm(ETa_Silalatshani_GAM1_vis, 3), trans = function(x){ 
  ETa_Silalatshani_GAM$family$linkinv(coef(ETa_Silalatshani_GAM)[1] + x) 
})
ETa_Silalatshani_space_plot_gg <- ggplot(data=ETa_Silalatshani_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="ETa (mm)")+
  theme(element_blank()) + coord_equal()

ETa_Silalatshani_yearxspace <- plotSlice(sm(ETa_Silalatshani_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

ETa_Silalatshani_yearxspace_plot <- ggplot()+geom_raster(data=ETa_Silalatshani_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=ETa_Silalatshani_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="ETa mm/day")+coord_equal()+guides(alpha="none", colour="none")

##

Ea_Silalatshani_year_plot <- plot_smooth(Ea_Silalatshani_GAM, view = 'year', transform = exp)

##
Ea_Silalatshani_space_plot <- plot(sm(Ea_Silalatshani_GAM1_vis, 3), trans = function(x){ 
  Ea_Silalatshani_GAM$family$linkinv(coef(Ea_Silalatshani_GAM)[1] + x) 
})
Ea_Silalatshani_space_plot_gg <- ggplot(data=Ea_Silalatshani_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ea (mm)")+
  theme(element_blank()) + coord_equal()
##


Ta_Silalatshani_year_plot <- plot_smooth(Ta_Silalatshani_GAM, view = 'year', transform = exp)


##
Ta_Silalatshani_space_plot <- plot(sm(Ta_Silalatshani_GAM1_vis, 3), trans = function(x){ 
  Ea_Silalatshani_GAM$family$linkinv(coef(Ea_Silalatshani_GAM)[1] + x) 
})
Ta_Silalatshani_space_plot_gg <- ggplot(data=Ta_Silalatshani_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ta (mm)")+
  theme(element_blank()) + coord_equal()
##

GPP_Silalatshani_plot_dat <- data.frame(var="GPP", value=GPP_Silalatshani_year_plot$fv$fit,
                                ll=GPP_Silalatshani_year_plot$fv$ll,
                                ul=GPP_Silalatshani_year_plot$fv$ul,
                                year=GPP_Silalatshani_year_plot$fv$year)

WUE_Silalatshani_plot_dat <- data.frame(var="WUE", value=WUE_Silalatshani_year_plot$fv$fit,
                                  ll=WUE_Silalatshani_year_plot$fv$ll,
                                  ul=WUE_Silalatshani_year_plot$fv$ul,
                                  year=WUE_Silalatshani_year_plot$fv$year)

E_Silalatshani_plot_dat <- data.frame(var="E", value=Ea_Silalatshani_year_plot$fv$fit,
                                 ll=Ea_Silalatshani_year_plot$fv$ll,
                                 ul=Ea_Silalatshani_year_plot$fv$ul,
                                 year=Ea_Silalatshani_year_plot$fv$year)
T_Silalatshani_plot_dat <- data.frame(var="T", value=Ta_Silalatshani_year_plot$fv$fit,
                                ll=Ta_Silalatshani_year_plot$fv$ll,
                                ul=Ta_Silalatshani_year_plot$fv$ul,
                                year=Ta_Silalatshani_year_plot$fv$year)
ET_Silalatshani_plot_dat <- data.frame(var="ET", value=ETa_Silalatshani_year_plot$fv$fit,
                                 ll=ETa_Silalatshani_year_plot$fv$ll,
                                 ul=ETa_Silalatshani_year_plot$fv$ul,
                                 year=ETa_Silalatshani_year_plot$fv$year)
All_Silalatshani_plot_dat <- rbind(GPP_Silalatshani_plot_dat, WUE_Silalatshani_plot_dat,
                             E_Silalatshani_plot_dat, T_Silalatshani_plot_dat, ET_Silalatshani_plot_dat)

All_Silalatshani_plot_dat$var <- as.factor(All_Silalatshani_plot_dat$var)
All_Silalatshani_plot_dat$var <- factor(All_Silalatshani_plot_dat$var, levels=(c("WUE", "GPP", "ET", "E", "T")))
All_Silalatshani_plot_dat$scheme <- "Silalatshani"
saveRDS(All_Silalatshani_plot_dat, "All_Silalatshani_year.rds")

ET_Silalatshani_plot_dat <- rbind(E_Silalatshani_plot_dat, T_Silalatshani_plot_dat, ET_Silalatshani_plot_dat)

ET_Silalatshani_plot_dat$var <- as.factor(ET_Silalatshani_plot_dat$var)
ET_Silalatshani_plot_dat$var <- factor(ET_Silalatshani_plot_dat$var, levels=(c("ET", "E", "T")))

ET_Silalatshani_year_plot <- ggplot(data= ET_Silalatshani_plot_dat, aes(x=year)) + 
  geom_line(aes(y=value, col=var))+
  geom_ribbon(aes(ymin=ll, ymax=ul, fill=var), alpha=0.2)+
  scale_x_continuous(breaks=seq(2013, 2021, by=1))+
  scale_color_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  theme(legend.title = element_blank())+
  labs(x="Year", y="mm/day")+
  theme_minimal() +theme(legend.title=element_blank())
ET_Silalatshani_year_plot

##
laya <- rbind(c(1,2,3))
gs <- list(ETa_Silalatshani_space_plot_gg, 
           Ta_Silalatshani_space_plot_gg, Ea_Silalatshani_space_plot_gg)
grid.arrange(grobs=gs, layout=laya, ncol=3)

ETa_Silalatshani_oneplot_dat <- data.frame(x=ETa_Silalatshani_space_plot$data$fit$x,
                                     y=ETa_Silalatshani_space_plot$data$fit$y,
                                     z=ETa_Silalatshani_space_plot$data$fit$z,
                                     var=c("Evapotranspiration"))

Ta_Silalatshani_oneplot_dat <- data.frame(x=Ta_Silalatshani_space_plot$data$fit$x,
                                      y=Ta_Silalatshani_space_plot$data$fit$y,
                                      z=Ta_Silalatshani_space_plot$data$fit$z,
                                      var=c("Transpiration"))

Ea_Silalatshani_oneplot_dat <- data.frame(x=Ea_Silalatshani_space_plot$data$fit$x,
                                     y=Ea_Silalatshani_space_plot$data$fit$y,
                                     z=Ea_Silalatshani_space_plot$data$fit$z,
                                     var=c("Evaporation"))

ET_Silalatshani_oneplot_dat <- rbind(ETa_Silalatshani_oneplot_dat, Ta_Silalatshani_oneplot_dat,
                                Ea_Silalatshani_oneplot_dat)

ET_Silalatshani_oneplot_dat$var <- factor(ET_Silalatshani_oneplot_dat$var,
                                     levels=c("Evapotranspiration","Evaporation","Transpiration"))
ET_Silalatshani_oneplot <- ggplot(data=ET_Silalatshani_oneplot_dat, aes(x=x, y=y))+geom_raster(aes(fill=z))+
  geom_contour(aes(z=z), colour="black") + facet_grid(rows=vars(var)) + coord_equal() +
  scale_fill_viridis_c(name="mm")+theme(element_blank(), axis.text = element_blank())


##

p <- ggdraw() +
  #draw_plot(WUE_Silalatshani_year_plot_gg, x=0, y=0.8, width=0.9, height=0.2) +
  draw_plot(GPP_Silalatshani_year_plot_gg, x=0, y=0.6, width=1, height=0.4) +
  #draw_plot(GPP_Silalatshani_space_plot_gg, x=0.6, y=0.6,width = 0.4, height=0.2)+
  #draw_plot(WUE_Silalatshani_space_plot_gg, x=0.6, y=0.8,width = 0.4, height=0.2)+
  draw_plot(ET_Silalatshani_year_plot, x=0, y=0,width = 1, height=0.6)#+
  #draw_plot(ET_Silalatshani_oneplot, x=0.6, y=0,width = 0.4, height=0.6) 

title <- ggdraw() + draw_label("Silalatshani", fontface='bold')
Silalatshani_GAMplot <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
Silalatshani_GAMplot
ggsave('Silalatshani_GAMplot.jpeg', dpi=900, height=10, width=10, device='jpeg')
save(Silalatshani_GAMplot, file='Plots/Silalatshani_GAMplot.rdata')