#' Magozi_df
#' 
#' @description ET, T, E, and GPP data for the Magozi irrigation scheme.
#' 
#' @format A data frame comprising 1648248 rows and 35 columns. 
#' \describe{
#'    \item{x}{Longitude}
#'    \item{y}{Latitude}
#'    \item{time}{Time of satellite image acquisition}
#'    \item{blue}{Value for Landsat blue band}
#'    \item{green}{Value for Landsat green band}
#'    \item{red}{Value for Landsat red band}
#'    \item{nir}{Value for Landsat near infrared band}
#'    \item{swir1}{Value for Landsat shortwave infrared 1 band}
#'    \item{swir2}{Value for Landsat shortwave infrared 2 band}
#'    \item{sensor}{Landsat sensor- 7 or 8}
#'    \item{Rh}{Relative humidity}
#'    \item{daily_Ta}{Daily air temperature}
#'    \item{NDVI_vals}{NDVI values}
#'    \item{LAI_vals}{Leaf Area Index values}
#'    \item{SAVI_vals}{Soil Adjusted Vegetation Index values}
#'    \item{GCVI}{Green Chlorophyll Vegetation Index values}
#'    \item{DOY}{Day of Year}
#'    \item{Dlength}{Daylength (hours)}
#'    \item{Lat}{Latitude}
#'    \item{DAS}{Days After Sowing}
#'    \item{Crop}{Crop type (estimated but not used)}
#'    \item{LUE}{Predicted Light Use Efficiency}
#'    \item{SWdown}{Downwelling shortwave radiation}
#'    \item{PAR}{Photosynthetically Active Radiation}
#'    \item{fAPAR}{fraction of Absorbed Photosynthetically Active Radiation}
#'    \item{GPP}{Gross Primary Productivity prediction}
#'    \item{ETa}{Actual evapotranspiration}
#'    \item{Ea}{Actual evaporation}
#'    \item{Ta}{Actual transpiration}
#'    \item{date}{Date of image acquisition}
#'    \item{chirps}{CHIRPS rainfall (mm) on day of acquisition}
#'    \item{chirps_30}{CHIRPS rainfall (mm) with 30 day lag}
#'    \item{chirps_60}{CHIRPS rainfall (mm) with 60 day lag}
#'    \item{chirps_90}{CHIRPS rainfall (mm) with 90 day lag}
#'    \item{chirps_120}{CHIRPS rainfall (mm) with 120 day lag}
#'    }
"Magozi_df"