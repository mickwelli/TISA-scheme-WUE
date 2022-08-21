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

Magozi_df <- readRDS("GAMM_data\\Magozi_df.rds")

# Columns 
Magozi_df$sensor <- as.factor(Magozi_df$sensor)
Magozi_df <- Magozi_df %>% 
  mutate_at(vars(date), dplyr::funs(year = lubridate::year, month = lubridate::month, day =lubridate::day))


# Clean data
Magozi_df_nona <- Magozi_df[!is.na(Magozi_df$GPP) & (Magozi_df$GPP > 0) & 
                                !is.na(Magozi_df$ETa) & (Magozi_df$ETa > 0) &
                                !is.na(Magozi_df$Ea) & (Magozi_df$Ea > 0) & 
                                !is.na(Magozi_df$Ta) & (Magozi_df$Ta > 0) &
                                (Magozi_df$year > 2012),]

#Magozi_df_nona <- Magozi_df_nona %>% mutate(GPP = replace(GPP, GPP>15, NA))


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
GPP_Magozi_GAM1 <- bam(GPP_f_30, discrete=TRUE, nthreads=8, data=Magozi_df_nona)
GPP_Magozi_GAM2 <- bam(GPP_f_60, discrete=TRUE, nthreads=8, data=Magozi_df_nona)
GPP_Magozi_GAM3 <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Magozi_df_nona)
GPP_Magozi_GAM4 <- bam(GPP_f_120, discrete=TRUE, nthreads=8, data=Magozi_df_nona)

anova(GPP_Magozi_GAM1, GPP_Magozi_GAM2, GPP_Magozi_GAM3, GPP_Magozi_GAM4, test='F')

# Run GAMs for biomass (GPP)
GPP_Magozi_GAM <- GPP_Magozi_GAM3

resids_Magozi_GPP <- residuals.gam(GPP_Magozi_GAM)
valRho_Magozi_GPP <- acf(resids_Magozi_GPP, plot=FALSE)$acf[2]
GPP_Magozi_GAM <- bam(GPP_f_90, discrete=TRUE, nthreads=8, data=Magozi_df_nona, rho = valRho_Magozi_GPP)
GPP_Magozi_GAM1_vis <- getViz(GPP_Magozi_GAM)

check(GPP_Magozi_GAM1_vis)

GPP_Magozi_year_plot <- plot_smooth(GPP_Magozi_GAM, view = 'year', transform = exp)

##
GPP_Magozi_year_plot_gg <- ggplot(data=GPP_Magozi_year_plot$fv, aes(x=year, y=fit))+
  geom_line(aes(col="GPP"))+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y=expression(paste("gC/m"^2,paste("/day"))))+theme_minimal() +
  scale_color_manual(values=c("black")) + theme(legend.title=element_blank())

GPP_Magozi_space_plot <- plot(sm(GPP_Magozi_GAM1_vis, 3), trans = function(x){ 
  GPP_Magozi_GAM$family$linkinv(coef(GPP_Magozi_GAM)[1] + x) 
})
GPP_Magozi_space_plot_gg <- ggplot(data=GPP_Magozi_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name=expression(paste("GPP gC/m"^2)))+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()


##

GPP_Magozi_yearxspace <- plotSlice(sm(GPP_Magozi_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

GPP_Magozi_yearxspace_plot <- ggplot()+geom_raster(data=GPP_Magozi_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=GPP_Magozi_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="GPP (gC/m2)")+coord_equal()+guides(alpha="none", colour="none")

GPP_Magozi_yearxmonth <- plot(sm(GPP_Magozi_GAM1_vis, 4), trans = function(x){ 
  GPP_Magozi_GAM$family$linkinv(coef(GPP_Magozi_GAM)[1] + x) 
}) 

GPP_Magozi_yearxmonth_plot <- ggplot(data=GPP_Magozi_yearxmonth$data$fit, aes(x=y, y=x, fill=tz)) +
  geom_raster()+
  geom_contour(aes(x=y, y=x, z=tz, col="red"), show.legend = FALSE)+
  scale_fill_viridis_c(labels=label_wrap(5))+
  scale_y_continuous(breaks=seq(2013, 2021, by=1))+
  scale_x_continuous(breaks=seq(1, 12, by=1))+
  labs(y="Year", x="Month", fill="GPP (gC/m2)") + theme_bw()

# Run GAMs for WUE, evapotranspiration (ETa), evaporation (Ea), transpiration (Ta) 
# Stick with 90 day lag

Magozi_df_nona$WUE <- Magozi_df_nona$GPP / Magozi_df_nona$ETa

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

WUE_Magozi_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona)
resids_Magozi_WUE <- residuals.gam(WUE_Magozi_GAM)
valRho_Magozi_WUE <- acf(resids_Magozi_WUE , plot=FALSE)$acf[2]
WUE_Magozi_GAM <- bam(WUE_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona, rho = valRho_Magozi_WUE)

WUE_Magozi_GAM1_vis <- getViz(WUE_Magozi_GAM)

ETa_Magozi_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona)
resids_Magozi_ETa <- residuals.gam(ETa_Magozi_GAM)
valRho_Magozi_ETa <- acf(resids_Magozi_ETa , plot=FALSE)$acf[2]
ETa_Magozi_GAM <- bam(ETa_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona, rho = valRho_Magozi_ETa)

ETa_Magozi_GAM1_vis <- getViz(ETa_Magozi_GAM)

Ta_Magozi_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona)
resids_Magozi_Ta <- residuals.gam(Ta_Magozi_GAM)
valRho_Magozi_Ta <- acf(resids_Magozi_Ta , plot=FALSE)$acf[2]
Ta_Magozi_GAM <- bam(Ta_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona, rho = valRho_Magozi_Ta)

Ta_Magozi_GAM1_vis <- getViz(Ta_Magozi_GAM)

Ea_Magozi_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona)
resids_Magozi_Ea <- residuals.gam(Ea_Magozi_GAM)
valRho_Magozi_Ea <- acf(resids_Magozi_Ea , plot=FALSE)$acf[2]
Ea_Magozi_GAM <- bam(Ea_f_90, discrete=TRUE, nthreads=8, data = Magozi_df_nona, rho = valRho_Magozi_Ea)

Ea_Magozi_GAM1_vis <- getViz(Ea_Magozi_GAM)

WUE_Magozi_year_plot <- plot_smooth(WUE_Magozi_GAM, view = 'year', transform = exp)

WUE_Magozi_year_plot_gg <- ggplot(data=WUE_Magozi_year_plot$fv, aes(x=year, y=fit))+
  geom_line()+geom_ribbon(aes(ymin=ll, ymax=ul), alpha=0.2)+scale_x_continuous(breaks=seq(2013,2021,1))+
  labs(x="Year", y="WUE gC/mm")+theme_minimal()

##
WUE_Magozi_space_plot <- plot(sm(WUE_Magozi_GAM1_vis, 3), trans = function(x){ 
  WUE_Magozi_GAM$family$linkinv(coef(WUE_Magozi_GAM)[1] + x) 
})
WUE_Magozi_space_plot_gg <- ggplot(data=WUE_Magozi_space_plot$data$fit, aes(x=x, y=y, z=tz))+
  geom_raster(aes(fill=tz)) + geom_contour(colour="black") + scale_fill_viridis_c(name="WUE gC/mm")+
  theme(element_blank(), axis.text = element_blank()) +coord_equal()
##

ETa_Magozi_year_plot <- plot_smooth(ETa_Magozi_GAM, view = 'year', transform = exp)

##
ETa_Magozi_space_plot <- plot(sm(ETa_Magozi_GAM1_vis, 3), trans = function(x){ 
  ETa_Magozi_GAM$family$linkinv(coef(ETa_Magozi_GAM)[1] + x) 
})
ETa_Magozi_space_plot_gg <- ggplot(data=ETa_Magozi_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="ETa (mm)")+
  theme(element_blank()) + coord_equal()

ETa_Magozi_yearxspace <- plotSlice(sm(ETa_Magozi_GAM1_vis,5), fix=list("year"=seq(2013,2021)), trans = exp)+labs(title = "GPP") +
  guides(fill=guide_colourbar("GPP"))

ETa_Magozi_yearxspace_plot <- ggplot()+geom_raster(data=ETa_Magozi_yearxspace$data$fit, aes(x=x, y=y, fill=tz), alpha=0.8)+
  scale_fill_viridis_c(labels=label_wrap(5))+theme_void()+
  facet_wrap(~.fx.year)+geom_contour(data=ETa_Magozi_yearxspace$data$fit, aes(x=x, y=y, z=tz, col='red'))+
  labs(fill="ETa mm/day")+coord_equal()+guides(alpha="none", colour="none")

##

Ea_Magozi_year_plot <- plot_smooth(Ea_Magozi_GAM, view = 'year', transform = exp)

##
Ea_Magozi_space_plot <- plot(sm(Ea_Magozi_GAM1_vis, 3), trans = function(x){ 
  Ea_Magozi_GAM$family$linkinv(coef(Ea_Magozi_GAM)[1] + x) 
})
Ea_Magozi_space_plot_gg <- ggplot(data=Ea_Magozi_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ea (mm)")+
  theme(element_blank()) + coord_equal()
##


Ta_Magozi_year_plot <- plot_smooth(Ta_Magozi_GAM, view = 'year', transform = exp)


##
Ta_Magozi_space_plot <- plot(sm(Ta_Magozi_GAM1_vis, 3), trans = function(x){ 
  Ea_Magozi_GAM$family$linkinv(coef(Ea_Magozi_GAM)[1] + x) 
})
Ta_Magozi_space_plot_gg <- ggplot(data=Ta_Magozi_space_plot$data$fit, aes(x=x, y=y, z=z))+
  geom_raster(aes(fill=z)) + geom_contour(colour="black") + scale_fill_viridis_c(name="Ta (mm)")+
  theme(element_blank()) + coord_equal()
##

GPP_Magozi_plot_dat <- data.frame(var="GPP", value=GPP_Magozi_year_plot$fv$fit,
                                ll=GPP_Magozi_year_plot$fv$ll,
                                ul=GPP_Magozi_year_plot$fv$ul,
                                year=GPP_Magozi_year_plot$fv$year)

WUE_Magozi_plot_dat <- data.frame(var="WUE", value=WUE_Magozi_year_plot$fv$fit,
                                  ll=WUE_Magozi_year_plot$fv$ll,
                                  ul=WUE_Magozi_year_plot$fv$ul,
                                  year=WUE_Magozi_year_plot$fv$year)

E_Magozi_plot_dat <- data.frame(var="E", value=Ea_Magozi_year_plot$fv$fit,
                                 ll=Ea_Magozi_year_plot$fv$ll,
                                 ul=Ea_Magozi_year_plot$fv$ul,
                                 year=Ea_Magozi_year_plot$fv$year)
T_Magozi_plot_dat <- data.frame(var="T", value=Ta_Magozi_year_plot$fv$fit,
                                ll=Ta_Magozi_year_plot$fv$ll,
                                ul=Ta_Magozi_year_plot$fv$ul,
                                year=Ta_Magozi_year_plot$fv$year)
ET_Magozi_plot_dat <- data.frame(var="ET", value=ETa_Magozi_year_plot$fv$fit,
                                 ll=ETa_Magozi_year_plot$fv$ll,
                                 ul=ETa_Magozi_year_plot$fv$ul,
                                 year=ETa_Magozi_year_plot$fv$year)
All_Magozi_plot_dat <- rbind(GPP_Magozi_plot_dat, WUE_Magozi_plot_dat,
                             E_Magozi_plot_dat, T_Magozi_plot_dat, ET_Magozi_plot_dat)

All_Magozi_plot_dat$var <- as.factor(All_Magozi_plot_dat$var)
All_Magozi_plot_dat$var <- factor(All_Magozi_plot_dat$var, levels=(c("WUE", "GPP", "ET", "E", "T")))
All_Magozi_plot_dat$scheme <- "Magozi"
saveRDS(All_Magozi_plot_dat, "All_Magozi_year.rds")

ET_Magozi_plot_dat <- rbind(E_Magozi_plot_dat, T_Magozi_plot_dat, ET_Magozi_plot_dat)

ET_Magozi_plot_dat$var <- as.factor(ET_Magozi_plot_dat$var)
ET_Magozi_plot_dat$var <- factor(ET_Magozi_plot_dat$var, levels=(c("ET", "E", "T")))

ET_Magozi_year_plot <- ggplot(data= ET_Magozi_plot_dat, aes(x=year)) + 
  geom_line(aes(y=value, col=var))+
  geom_ribbon(aes(ymin=ll, ymax=ul, fill=var), alpha=0.2)+
  scale_x_continuous(breaks=seq(2013, 2021, by=1))+
  scale_color_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00", "#009E73"))+
  theme(legend.title = element_blank())+
  labs(x="Year", y="mm/day")+
  theme_minimal() +theme(legend.title=element_blank())
ET_Magozi_year_plot

##
laya <- rbind(c(1,2,3))
gs <- list(ETa_Magozi_space_plot_gg, 
           Ta_Magozi_space_plot_gg, Ea_Magozi_space_plot_gg)
grid.arrange(grobs=gs, layout=laya, ncol=3)

ETa_Magozi_oneplot_dat <- data.frame(x=ETa_Magozi_space_plot$data$fit$x,
                                     y=ETa_Magozi_space_plot$data$fit$y,
                                     z=ETa_Magozi_space_plot$data$fit$z,
                                     var=c("Evapotranspiration"))

Ta_Magozi_oneplot_dat <- data.frame(x=Ta_Magozi_space_plot$data$fit$x,
                                      y=Ta_Magozi_space_plot$data$fit$y,
                                      z=Ta_Magozi_space_plot$data$fit$z,
                                      var=c("Transpiration"))

Ea_Magozi_oneplot_dat <- data.frame(x=Ea_Magozi_space_plot$data$fit$x,
                                     y=Ea_Magozi_space_plot$data$fit$y,
                                     z=Ea_Magozi_space_plot$data$fit$z,
                                     var=c("Evaporation"))

ET_Magozi_oneplot_dat <- rbind(ETa_Magozi_oneplot_dat, Ta_Magozi_oneplot_dat,
                                Ea_Magozi_oneplot_dat)

ET_Magozi_oneplot_dat$var <- factor(ET_Magozi_oneplot_dat$var,
                                     levels=c("Evapotranspiration","Evaporation","Transpiration"))
ET_Magozi_oneplot <- ggplot(data=ET_Magozi_oneplot_dat, aes(x=x, y=y))+geom_raster(aes(fill=z))+
  geom_contour(aes(z=z), colour="black") + facet_grid(rows=vars(var)) + coord_equal() +
  scale_fill_viridis_c(name="mm")+theme(element_blank(), axis.text = element_blank())


##

p <- ggdraw() +
  #draw_plot(WUE_Magozi_year_plot_gg, x=0, y=0.8, width=0.9, height=0.2) +
  draw_plot(GPP_Magozi_year_plot_gg, x=0, y=0.6, width=1, height=0.4) +
  #draw_plot(GPP_Magozi_space_plot_gg, x=0.6, y=0.6,width = 0.4, height=0.2)+
  #draw_plot(WUE_Magozi_space_plot_gg, x=0.6, y=0.8,width = 0.4, height=0.2)+
  draw_plot(ET_Magozi_year_plot, x=0, y=0,width = 1, height=0.6)#+
  #draw_plot(ET_Magozi_oneplot, x=0.6, y=0,width = 0.4, height=0.6) 

title <- ggdraw() + draw_label("Magozi", fontface='bold')
Magozi_GAMplot <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
Magozi_GAMplot
ggsave('Magozi_GAMplot.jpeg', dpi=900, height=10, width=10, device='jpeg')
save(Magozi_GAMplot, file='Plots/Magozi_GAMplot.rdata')