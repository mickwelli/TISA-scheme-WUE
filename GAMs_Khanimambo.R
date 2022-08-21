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

Khanimambo_df <- readRDS("GAMM_data\\Khanimambo_df.rds")

# Columns 
Khanimambo_df$sensor <- as.factor(Khanimambo_df$sensor)
Khanimambo_df <- Khanimambo_df %>% 
  mutate_at(vars(date), dplyr::funs(year = lubridate::year, month = lubridate::month, day =lubridate::day))


# Clean data
Khanimambo_df_nona <- Khanimambo_df[!is.na(Khanimambo_df$GPP) & (Khanimambo_df$GPP > 0) & 
                                !is.na(Khanimambo_df$ETa) & (Khanimambo_df$ETa > 0) &
                                !is.na(Khanimambo_df$Ea) & (Khanimambo_df$Ea > 0) & 
                                !is.na(Khanimambo_df$Ta) & (Khanimambo_df$Ta > 0) &
                                (Khanimambo_df$year > 2012),]

#Khanimambo_df_nona <- Khanimambo_df_nona %>% mutate(GPP = replace(GPP, GPP>15, NA))


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
GPP_Khanimambo_GAM1 <- bam(GPP_f_30, discrete=TRUE, nthreads=8, data=Khanimambo_df_nona)
GPP_Khanimambo_GAM2 <- bam(GPP_f_60, discrete=TRUE, nthreads=8, data=Khanimambo_df_nona)
GPP_Khanimambo_GAM3 <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Khanimambo_df_nona)
GPP_Khanimambo_GAM4 <- bam(GPP_f_120, discrete=TRUE, nthreads=8, data=Khanimambo_df_nona)

anova(GPP_Khanimambo_GAM1, GPP_Khanimambo_GAM2, GPP_Khanimambo_GAM3, GPP_Khanimambo_GAM4, test='F')

# Run GAMs for biomass (GPP)
GPP_Khanimambo_GAM <- GPP_Khanimambo_GAM3

resids_Khanimambo_GPP <- residuals.gam(GPP_Khanimambo_GAM)
valRho_Khanimambo_GPP <- acf(resids_Khanimambo_GPP, plot=FALSE)$acf[2]
GPP_Khanimambo_GAM <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Khanimambo_df_nona, rho = valRho_Khanimambo_GPP)
GPP_Khanimambo_GAM1_vis <- getViz(GPP_Khanimambo_GAM)

check(GPP_Khanimambo_GAM1_vis)

GPP_Khanimambo_year_plot <- plot_smooth(GPP_Khanimambo_GAM, view = 'year', transform = exp)

##
GPP_Khanimambo_year_plot_gg <- ggplot(data=GPP_Khanimambo_year_plot$fv, aes(x=year, y=fit))+
  geom_line(aes(col="GPP"))+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y=expression(paste("gC/m"^2,paste("/day"))))+theme_minimal() +
  scale_color_manual(values=c("black")) + theme(legend.title=element_blank())

GPP_Khanimambo_space_plot <- plot(sm(GPP_Khanimambo_GAM1_vis, 3), trans = function(x){ 
  GPP_Khanimambo_GAM$family$linkinv(coef(GPP_Khanimambo_GAM)[1] + x) 
})
GPP_Khanimambo_space_plot_gg <- ggplot(data=GPP_Khanimambo_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name=expression(paste("GPP gC/m"^2)))+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()


##

GPP_Khanimambo_yearxspace <- plotSlice(sm(GPP_Khanimambo_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

GPP_Khanimambo_yearxspace_plot <- ggplot()+geom_raster(data=GPP_Khanimambo_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=GPP_Khanimambo_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="GPP (gC/m2)")+coord_equal()+guides(alpha="none", colour="none")

GPP_Khanimambo_yearxmonth <- plot(sm(GPP_Khanimambo_GAM1_vis, 4), trans = function(x){ 
  GPP_Khanimambo_GAM$family$linkinv(coef(GPP_Khanimambo_GAM)[1] + x) 
}) 

GPP_Khanimambo_yearxmonth_plot <- ggplot(data=GPP_Khanimambo_yearxmonth$data$fit, aes(x=y, y=x, fill=tz)) +
  geom_raster()+
  geom_contour(aes(x=y, y=x, z=tz, col="red"), show.legend = FALSE)+
  scale_fill_viridis_c(labels=label_wrap(5))+
  scale_y_continuous(breaks=seq(2013, 2021, by=1))+
  scale_x_continuous(breaks=seq(1, 12, by=1))+
  labs(y="Year", x="Month", fill="GPP (gC/m2)") + theme_bw()

# Run GAMs for WUE, evapotranspiration (ETa), evaporation (Ea), transpiration (Ta) 
# Stick with 90 day lag

Khanimambo_df_nona$WUE <- Khanimambo_df_nona$GPP / Khanimambo_df_nona$ETa

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

WUE_Khanimambo_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona)
resids_Khanimambo_WUE <- residuals.gam(WUE_Khanimambo_GAM)
valRho_Khanimambo_WUE <- acf(resids_Khanimambo_WUE , plot=FALSE)$acf[2]
WUE_Khanimambo_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona, rho = valRho_Khanimambo_WUE)

WUE_Khanimambo_GAM1_vis <- getViz(WUE_Khanimambo_GAM)

ETa_Khanimambo_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona)
resids_Khanimambo_ETa <- residuals.gam(ETa_Khanimambo_GAM)
valRho_Khanimambo_ETa <- acf(resids_Khanimambo_ETa , plot=FALSE)$acf[2]
ETa_Khanimambo_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona, rho = valRho_Khanimambo_ETa)

ETa_Khanimambo_GAM1_vis <- getViz(ETa_Khanimambo_GAM)

Ta_Khanimambo_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona)
resids_Khanimambo_Ta <- residuals.gam(Ta_Khanimambo_GAM)
valRho_Khanimambo_Ta <- acf(resids_Khanimambo_Ta , plot=FALSE)$acf[2]
Ta_Khanimambo_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona, rho = valRho_Khanimambo_Ta)

Ta_Khanimambo_GAM1_vis <- getViz(Ta_Khanimambo_GAM)

Ea_Khanimambo_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona)
resids_Khanimambo_Ea <- residuals.gam(Ea_Khanimambo_GAM)
valRho_Khanimambo_Ea <- acf(resids_Khanimambo_Ea , plot=FALSE)$acf[2]
Ea_Khanimambo_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Khanimambo_df_nona, rho = valRho_Khanimambo_Ea)

Ea_Khanimambo_GAM1_vis <- getViz(Ea_Khanimambo_GAM)

WUE_Khanimambo_year_plot <- plot_smooth(WUE_Khanimambo_GAM, view = 'year', transform = exp)

WUE_Khanimambo_year_plot_gg <- ggplot(data=WUE_Khanimambo_year_plot$fv, aes(x=year, y=fit))+
  geom_line()+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y="WUE gC/mm")+theme_minimal()

##
WUE_Khanimambo_space_plot <- plot(sm(WUE_Khanimambo_GAM1_vis, 3), trans = function(x){ 
  WUE_Khanimambo_GAM$family$linkinv(coef(WUE_Khanimambo_GAM)[1] + x) 
})
WUE_Khanimambo_space_plot_gg <- ggplot(data=WUE_Khanimambo_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name="WUE gC/mm")+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()
##

ETa_Khanimambo_year_plot <- plot_smooth(ETa_Khanimambo_GAM, view = 'year', transform = exp)

##
ETa_Khanimambo_space_plot <- plot(sm(ETa_Khanimambo_GAM1_vis, 3), trans = function(x){ 
  ETa_Khanimambo_GAM$family$linkinv(coef(ETa_Khanimambo_GAM)[1] + x) 
})
ETa_Khanimambo_space_plot_gg <- ggplot(data=ETa_Khanimambo_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="ETa (mm)")+
  theme(element_blank()) + coord_equal()

ETa_Khanimambo_yearxspace <- plotSlice(sm(ETa_Khanimambo_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

ETa_Khanimambo_yearxspace_plot <- ggplot()+geom_raster(data=ETa_Khanimambo_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=ETa_Khanimambo_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="ETa mm/day")+coord_equal()+guides(alpha="none", colour="none")

##

Ea_Khanimambo_year_plot <- plot_smooth(Ea_Khanimambo_GAM, view = 'year', transform = exp)

##
Ea_Khanimambo_space_plot <- plot(sm(Ea_Khanimambo_GAM1_vis, 3), trans = function(x){ 
  Ea_Khanimambo_GAM$family$linkinv(coef(Ea_Khanimambo_GAM)[1] + x) 
})
Ea_Khanimambo_space_plot_gg <- ggplot(data=Ea_Khanimambo_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ea (mm)")+
  theme(element_blank()) + coord_equal()
##


Ta_Khanimambo_year_plot <- plot_smooth(Ta_Khanimambo_GAM, view = 'year', transform = exp)


##
Ta_Khanimambo_space_plot <- plot(sm(Ta_Khanimambo_GAM1_vis, 3), trans = function(x){ 
  Ea_Khanimambo_GAM$family$linkinv(coef(Ea_Khanimambo_GAM)[1] + x) 
})
Ta_Khanimambo_space_plot_gg <- ggplot(data=Ta_Khanimambo_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ta (mm)")+
  theme(element_blank()) + coord_equal()
##

GPP_Khanimambo_plot_dat <- data.frame(var="GPP", value=GPP_Khanimambo_year_plot$fv$fit,
                                ll=GPP_Khanimambo_year_plot$fv$ll,
                                ul=GPP_Khanimambo_year_plot$fv$ul,
                                year=GPP_Khanimambo_year_plot$fv$year)

WUE_Khanimambo_plot_dat <- data.frame(var="WUE", value=WUE_Khanimambo_year_plot$fv$fit,
                                  ll=WUE_Khanimambo_year_plot$fv$ll,
                                  ul=WUE_Khanimambo_year_plot$fv$ul,
                                  year=WUE_Khanimambo_year_plot$fv$year)

E_Khanimambo_plot_dat <- data.frame(var="E", value=Ea_Khanimambo_year_plot$fv$fit,
                                 ll=Ea_Khanimambo_year_plot$fv$ll,
                                 ul=Ea_Khanimambo_year_plot$fv$ul,
                                 year=Ea_Khanimambo_year_plot$fv$year)
T_Khanimambo_plot_dat <- data.frame(var="T", value=Ta_Khanimambo_year_plot$fv$fit,
                                ll=Ta_Khanimambo_year_plot$fv$ll,
                                ul=Ta_Khanimambo_year_plot$fv$ul,
                                year=Ta_Khanimambo_year_plot$fv$year)
ET_Khanimambo_plot_dat <- data.frame(var="ET", value=ETa_Khanimambo_year_plot$fv$fit,
                                 ll=ETa_Khanimambo_year_plot$fv$ll,
                                 ul=ETa_Khanimambo_year_plot$fv$ul,
                                 year=ETa_Khanimambo_year_plot$fv$year)
All_Khanimambo_plot_dat <- rbind(GPP_Khanimambo_plot_dat, WUE_Khanimambo_plot_dat,
                             E_Khanimambo_plot_dat, T_Khanimambo_plot_dat, ET_Khanimambo_plot_dat)

All_Khanimambo_plot_dat$var <- as.factor(All_Khanimambo_plot_dat$var)
All_Khanimambo_plot_dat$var <- factor(All_Khanimambo_plot_dat$var, levels=(c("WUE", "GPP", "ET", "E", "T")))
All_Khanimambo_plot_dat$scheme <- "Khanimambo"
saveRDS(All_Khanimambo_plot_dat, "All_Khanimambo_year.rds")

ET_Khanimambo_plot_dat <- rbind(E_Khanimambo_plot_dat, T_Khanimambo_plot_dat, ET_Khanimambo_plot_dat)

ET_Khanimambo_plot_dat$var <- as.factor(ET_Khanimambo_plot_dat$var)
ET_Khanimambo_plot_dat$var <- factor(ET_Khanimambo_plot_dat$var, levels=(c("ET", "E", "T")))

ET_Khanimambo_year_plot <- ggplot(data= ET_Khanimambo_plot_dat, aes(x=year)) + 
  geom_line(aes(y=value, col=var))+
  geom_ribbon(aes(ymin=ll, ymax=ul, fill=var), alpha=0.2)+
  scale_x_continuous(breaks=seq(2013, 2021, by=1))+
  scale_color_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  theme(legend.title = element_blank())+
  labs(x="Year", y="mm/day")+
  theme_minimal() +theme(legend.title=element_blank())
ET_Khanimambo_year_plot

##
laya <- rbind(c(1,2,3))
gs <- list(ETa_Khanimambo_space_plot_gg, 
           Ta_Khanimambo_space_plot_gg, Ea_Khanimambo_space_plot_gg)
grid.arrange(grobs=gs, layout=laya, ncol=3)

ETa_Khanimambo_oneplot_dat <- data.frame(x=ETa_Khanimambo_space_plot$data$fit$x,
                                     y=ETa_Khanimambo_space_plot$data$fit$y,
                                     z=ETa_Khanimambo_space_plot$data$fit$z,
                                     var=c("Evapotranspiration"))

Ta_Khanimambo_oneplot_dat <- data.frame(x=Ta_Khanimambo_space_plot$data$fit$x,
                                      y=Ta_Khanimambo_space_plot$data$fit$y,
                                      z=Ta_Khanimambo_space_plot$data$fit$z,
                                      var=c("Transpiration"))

Ea_Khanimambo_oneplot_dat <- data.frame(x=Ea_Khanimambo_space_plot$data$fit$x,
                                     y=Ea_Khanimambo_space_plot$data$fit$y,
                                     z=Ea_Khanimambo_space_plot$data$fit$z,
                                     var=c("Evaporation"))

ET_Khanimambo_oneplot_dat <- rbind(ETa_Khanimambo_oneplot_dat, Ta_Khanimambo_oneplot_dat,
                                Ea_Khanimambo_oneplot_dat)

ET_Khanimambo_oneplot_dat$var <- factor(ET_Khanimambo_oneplot_dat$var,
                                     levels=c("Evapotranspiration","Evaporation","Transpiration"))
ET_Khanimambo_oneplot <- ggplot(data=ET_Khanimambo_oneplot_dat, aes(x=x, y=y))+geom_raster(aes(fill=z))+
  geom_contour(aes(z=z), colour="black") + facet_grid(rows=vars(var)) + coord_equal() +
  scale_fill_viridis_c(name="mm")+theme(element_blank(), axis.text = element_blank())


##

p <- ggdraw() +
  #draw_plot(WUE_Khanimambo_year_plot_gg, x=0, y=0.8, width=0.9, height=0.2) +
  draw_plot(GPP_Khanimambo_year_plot_gg, x=0, y=0.6, width=1, height=0.4) +
  #draw_plot(GPP_Khanimambo_space_plot_gg, x=0.6, y=0.6,width = 0.4, height=0.2)+
  #draw_plot(WUE_Khanimambo_space_plot_gg, x=0.6, y=0.8,width = 0.4, height=0.2)+
  draw_plot(ET_Khanimambo_year_plot, x=0, y=0,width = 1, height=0.6)#+
  #draw_plot(ET_Khanimambo_oneplot, x=0.6, y=0,width = 0.4, height=0.6) 

title <- ggdraw() + draw_label("Khanimambo", fontface='bold')
Khanimambo_GAMplot <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
Khanimambo_GAMplot
ggsave('Khanimambo_GAMplot.jpeg', dpi=900, height=10, width=10, device='jpeg')
save(Khanimambo_GAMplot, file='Plots/Khanimambo_GAMplot.rdata')