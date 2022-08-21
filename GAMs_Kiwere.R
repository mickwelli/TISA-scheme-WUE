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

Kiwere_df <- readRDS("GAMM_data\\Kiwere_df.rds")

# Columns 
Kiwere_df$sensor <- as.factor(Kiwere_df$sensor)
Kiwere_df <- Kiwere_df %>% 
  mutate_at(vars(date), dplyr::funs(year = lubridate::year, month = lubridate::month, day =lubridate::day))


# Clean data
Kiwere_df_nona <- Kiwere_df[!is.na(Kiwere_df$GPP) & (Kiwere_df$GPP > 0) & 
                                !is.na(Kiwere_df$ETa) & (Kiwere_df$ETa > 0) &
                                !is.na(Kiwere_df$Ea) & (Kiwere_df$Ea > 0) & 
                                !is.na(Kiwere_df$Ta) & (Kiwere_df$Ta > 0) &
                                (Kiwere_df$year > 2012),]

#Kiwere_df_nona <- Kiwere_df_nona %>% mutate(GPP = replace(GPP, GPP>15, NA))


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
GPP_Kiwere_GAM1 <- bam(GPP_f_30, discrete=TRUE, nthreads=8, data=Kiwere_df_nona)
GPP_Kiwere_GAM2 <- bam(GPP_f_60, discrete=TRUE, nthreads=8, data=Kiwere_df_nona)
GPP_Kiwere_GAM3 <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Kiwere_df_nona)
GPP_Kiwere_GAM4 <- bam(GPP_f_120, discrete=TRUE, nthreads=8, data=Kiwere_df_nona)

anova(GPP_Kiwere_GAM1, GPP_Kiwere_GAM2, GPP_Kiwere_GAM3, GPP_Kiwere_GAM4, test='F')

# Run GAMs for biomass (GPP)
GPP_Kiwere_GAM <- GPP_Kiwere_GAM3

resids_Kiwere_GPP <- residuals.gam(GPP_Kiwere_GAM)
valRho_Kiwere_GPP <- acf(resids_Kiwere_GPP, plot=FALSE)$acf[2]
GPP_Kiwere_GAM <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Kiwere_df_nona, rho = valRho_Kiwere_GPP)
GPP_Kiwere_GAM1_vis <- getViz(GPP_Kiwere_GAM)

check(GPP_Kiwere_GAM1_vis)

GPP_Kiwere_year_plot <- plot_smooth(GPP_Kiwere_GAM, view = 'year', transform = exp)

##
GPP_Kiwere_year_plot_gg <- ggplot(data=GPP_Kiwere_year_plot$fv, aes(x=year, y=fit))+
  geom_line(aes(col="GPP"))+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y=expression(paste("gC/m"^2,paste("/day"))))+theme_minimal() +
  scale_color_manual(values=c("black")) + theme(legend.title=element_blank())

GPP_Kiwere_space_plot <- plot(sm(GPP_Kiwere_GAM1_vis, 3), trans = function(x){ 
  GPP_Kiwere_GAM$family$linkinv(coef(GPP_Kiwere_GAM)[1] + x) 
})
GPP_Kiwere_space_plot_gg <- ggplot(data=GPP_Kiwere_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name=expression(paste("GPP gC/m"^2)))+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()


##

GPP_Kiwere_yearxspace <- plotSlice(sm(GPP_Kiwere_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

GPP_Kiwere_yearxspace_plot <- ggplot()+geom_raster(data=GPP_Kiwere_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=GPP_Kiwere_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="GPP (gC/m2)")+coord_equal()+guides(alpha="none", colour="none")

GPP_Kiwere_yearxmonth <- plot(sm(GPP_Kiwere_GAM1_vis, 4), trans = function(x){ 
  GPP_Kiwere_GAM$family$linkinv(coef(GPP_Kiwere_GAM)[1] + x) 
}) 

GPP_Kiwere_yearxmonth_plot <- ggplot(data=GPP_Kiwere_yearxmonth$data$fit, aes(x=y, y=x, fill=tz)) +
  geom_raster()+
  geom_contour(aes(x=y, y=x, z=tz, col="red"), show.legend = FALSE)+
  scale_fill_viridis_c(labels=label_wrap(5))+
  scale_y_continuous(breaks=seq(2013, 2021, by=1))+
  scale_x_continuous(breaks=seq(1, 12, by=1))+
  labs(y="Year", x="Month", fill="GPP (gC/m2)") + theme_bw()

# Run GAMs for WUE, evapotranspiration (ETa), evaporation (Ea), transpiration (Ta) 
# Stick with 90 day lag

Kiwere_df_nona$WUE <- Kiwere_df_nona$GPP / Kiwere_df_nona$ETa

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

WUE_Kiwere_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona)
resids_Kiwere_WUE <- residuals.gam(WUE_Kiwere_GAM)
valRho_Kiwere_WUE <- acf(resids_Kiwere_WUE , plot=FALSE)$acf[2]
WUE_Kiwere_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona, rho = valRho_Kiwere_WUE)

WUE_Kiwere_GAM1_vis <- getViz(WUE_Kiwere_GAM)

ETa_Kiwere_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona)
resids_Kiwere_ETa <- residuals.gam(ETa_Kiwere_GAM)
valRho_Kiwere_ETa <- acf(resids_Kiwere_ETa , plot=FALSE)$acf[2]
ETa_Kiwere_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona, rho = valRho_Kiwere_ETa)

ETa_Kiwere_GAM1_vis <- getViz(ETa_Kiwere_GAM)

Ta_Kiwere_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona)
resids_Kiwere_Ta <- residuals.gam(Ta_Kiwere_GAM)
valRho_Kiwere_Ta <- acf(resids_Kiwere_Ta , plot=FALSE)$acf[2]
Ta_Kiwere_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona, rho = valRho_Kiwere_Ta)

Ta_Kiwere_GAM1_vis <- getViz(Ta_Kiwere_GAM)

Ea_Kiwere_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona)
resids_Kiwere_Ea <- residuals.gam(Ea_Kiwere_GAM)
valRho_Kiwere_Ea <- acf(resids_Kiwere_Ea , plot=FALSE)$acf[2]
Ea_Kiwere_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Kiwere_df_nona, rho = valRho_Kiwere_Ea)

Ea_Kiwere_GAM1_vis <- getViz(Ea_Kiwere_GAM)

WUE_Kiwere_year_plot <- plot_smooth(WUE_Kiwere_GAM, view = 'year', transform = exp)

WUE_Kiwere_year_plot_gg <- ggplot(data=WUE_Kiwere_year_plot$fv, aes(x=year, y=fit))+
  geom_line()+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y="WUE gC/mm")+theme_minimal()

##
WUE_Kiwere_space_plot <- plot(sm(WUE_Kiwere_GAM1_vis, 3), trans = function(x){ 
  WUE_Kiwere_GAM$family$linkinv(coef(WUE_Kiwere_GAM)[1] + x) 
})
WUE_Kiwere_space_plot_gg <- ggplot(data=WUE_Kiwere_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name="WUE gC/mm")+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()
##

ETa_Kiwere_year_plot <- plot_smooth(ETa_Kiwere_GAM, view = 'year', transform = exp)

##
ETa_Kiwere_space_plot <- plot(sm(ETa_Kiwere_GAM1_vis, 3), trans = function(x){ 
  ETa_Kiwere_GAM$family$linkinv(coef(ETa_Kiwere_GAM)[1] + x) 
})
ETa_Kiwere_space_plot_gg <- ggplot(data=ETa_Kiwere_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="ETa (mm)")+
  theme(element_blank()) + coord_equal()

ETa_Kiwere_yearxspace <- plotSlice(sm(ETa_Kiwere_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

ETa_Kiwere_yearxspace_plot <- ggplot()+geom_raster(data=ETa_Kiwere_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=ETa_Kiwere_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="ETa mm/day")+coord_equal()+guides(alpha="none", colour="none")

##

Ea_Kiwere_year_plot <- plot_smooth(Ea_Kiwere_GAM, view = 'year', transform = exp)

##
Ea_Kiwere_space_plot <- plot(sm(Ea_Kiwere_GAM1_vis, 3), trans = function(x){ 
  Ea_Kiwere_GAM$family$linkinv(coef(Ea_Kiwere_GAM)[1] + x) 
})
Ea_Kiwere_space_plot_gg <- ggplot(data=Ea_Kiwere_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ea (mm)")+
  theme(element_blank()) + coord_equal()
##


Ta_Kiwere_year_plot <- plot_smooth(Ta_Kiwere_GAM, view = 'year', transform = exp)


##
Ta_Kiwere_space_plot <- plot(sm(Ta_Kiwere_GAM1_vis, 3), trans = function(x){ 
  Ea_Kiwere_GAM$family$linkinv(coef(Ea_Kiwere_GAM)[1] + x) 
})
Ta_Kiwere_space_plot_gg <- ggplot(data=Ta_Kiwere_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ta (mm)")+
  theme(element_blank()) + coord_equal()
##

GPP_Kiwere_plot_dat <- data.frame(var="GPP", value=GPP_Kiwere_year_plot$fv$fit,
                                ll=GPP_Kiwere_year_plot$fv$ll,
                                ul=GPP_Kiwere_year_plot$fv$ul,
                                year=GPP_Kiwere_year_plot$fv$year)

WUE_Kiwere_plot_dat <- data.frame(var="WUE", value=WUE_Kiwere_year_plot$fv$fit,
                                  ll=WUE_Kiwere_year_plot$fv$ll,
                                  ul=WUE_Kiwere_year_plot$fv$ul,
                                  year=WUE_Kiwere_year_plot$fv$year)

E_Kiwere_plot_dat <- data.frame(var="E", value=Ea_Kiwere_year_plot$fv$fit,
                                 ll=Ea_Kiwere_year_plot$fv$ll,
                                 ul=Ea_Kiwere_year_plot$fv$ul,
                                 year=Ea_Kiwere_year_plot$fv$year)
T_Kiwere_plot_dat <- data.frame(var="T", value=Ta_Kiwere_year_plot$fv$fit,
                                ll=Ta_Kiwere_year_plot$fv$ll,
                                ul=Ta_Kiwere_year_plot$fv$ul,
                                year=Ta_Kiwere_year_plot$fv$year)
ET_Kiwere_plot_dat <- data.frame(var="ET", value=ETa_Kiwere_year_plot$fv$fit,
                                 ll=ETa_Kiwere_year_plot$fv$ll,
                                 ul=ETa_Kiwere_year_plot$fv$ul,
                                 year=ETa_Kiwere_year_plot$fv$year)
All_Kiwere_plot_dat <- rbind(GPP_Kiwere_plot_dat, WUE_Kiwere_plot_dat,
                             E_Kiwere_plot_dat, T_Kiwere_plot_dat, ET_Kiwere_plot_dat)

All_Kiwere_plot_dat$var <- as.factor(All_Kiwere_plot_dat$var)
All_Kiwere_plot_dat$var <- factor(All_Kiwere_plot_dat$var, levels=(c("WUE", "GPP", "ET", "E", "T")))
All_Kiwere_plot_dat$scheme <- "Kiwere"
saveRDS(All_Kiwere_plot_dat, "All_Kiwere_year.rds")

ET_Kiwere_plot_dat <- rbind(E_Kiwere_plot_dat, T_Kiwere_plot_dat, ET_Kiwere_plot_dat)

ET_Kiwere_plot_dat$var <- as.factor(ET_Kiwere_plot_dat$var)
ET_Kiwere_plot_dat$var <- factor(ET_Kiwere_plot_dat$var, levels=(c("ET", "E", "T")))

ET_Kiwere_year_plot <- ggplot(data= ET_Kiwere_plot_dat, aes(x=year)) + 
  geom_line(aes(y=value, col=var))+
  geom_ribbon(aes(ymin=ll, ymax=ul, fill=var), alpha=0.2)+
  scale_x_continuous(breaks=seq(2013, 2021, by=1))+
  scale_color_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  theme(legend.title = element_blank())+
  labs(x="Year", y="mm/day")+
  theme_minimal() +theme(legend.title=element_blank())
ET_Kiwere_year_plot

##
laya <- rbind(c(1,2,3))
gs <- list(ETa_Kiwere_space_plot_gg, 
           Ta_Kiwere_space_plot_gg, Ea_Kiwere_space_plot_gg)
grid.arrange(grobs=gs, layout=laya, ncol=3)

ETa_Kiwere_oneplot_dat <- data.frame(x=ETa_Kiwere_space_plot$data$fit$x,
                                     y=ETa_Kiwere_space_plot$data$fit$y,
                                     z=ETa_Kiwere_space_plot$data$fit$z,
                                     var=c("Evapotranspiration"))

Ta_Kiwere_oneplot_dat <- data.frame(x=Ta_Kiwere_space_plot$data$fit$x,
                                      y=Ta_Kiwere_space_plot$data$fit$y,
                                      z=Ta_Kiwere_space_plot$data$fit$z,
                                      var=c("Transpiration"))

Ea_Kiwere_oneplot_dat <- data.frame(x=Ea_Kiwere_space_plot$data$fit$x,
                                     y=Ea_Kiwere_space_plot$data$fit$y,
                                     z=Ea_Kiwere_space_plot$data$fit$z,
                                     var=c("Evaporation"))

ET_Kiwere_oneplot_dat <- rbind(ETa_Kiwere_oneplot_dat, Ta_Kiwere_oneplot_dat,
                                Ea_Kiwere_oneplot_dat)

ET_Kiwere_oneplot_dat$var <- factor(ET_Kiwere_oneplot_dat$var,
                                     levels=c("Evapotranspiration","Evaporation","Transpiration"))
ET_Kiwere_oneplot <- ggplot(data=ET_Kiwere_oneplot_dat, aes(x=x, y=y))+geom_raster(aes(fill=z))+
  geom_contour(aes(z=z), colour="black") + facet_grid(rows=vars(var)) + coord_equal() +
  scale_fill_viridis_c(name="mm")+theme(element_blank(), axis.text = element_blank())


##

p <- ggdraw() +
  #draw_plot(WUE_Kiwere_year_plot_gg, x=0, y=0.8, width=0.9, height=0.2) +
  draw_plot(GPP_Kiwere_year_plot_gg, x=0, y=0.6, width=1, height=0.4) +
  #draw_plot(GPP_Kiwere_space_plot_gg, x=0.6, y=0.6,width = 0.4, height=0.2)+
  #draw_plot(WUE_Kiwere_space_plot_gg, x=0.6, y=0.8,width = 0.4, height=0.2)+
  draw_plot(ET_Kiwere_year_plot, x=0, y=0,width = 1, height=0.6)#+
  #draw_plot(ET_Kiwere_oneplot, x=0.6, y=0,width = 0.4, height=0.6) 

title <- ggdraw() + draw_label("Kiwere", fontface='bold')
Kiwere_GAMplot <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
Kiwere_GAMplot
ggsave('Kiwere_GAMplot.jpeg', dpi=900, height=10, width=10, device='jpeg')
save(Kiwere_GAMplot, file='Plots/Kiwere_GAMplot.rdata')