source('~/Current Projects/SBG/LPJ/Comparisons/scripts/MODIS_comparison.R')


lpj <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc'
month <- 12
year <- 2015
out <- '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/LPJ_version021_Dec_2015_MODIS_comparisons.pdf'


pdf(out, height = 11, width = 8.5, paper = 'letter')
MODIS_comparison(lpj.nc = lpj,
                 month = month,
                 year = year)
dev.off()

