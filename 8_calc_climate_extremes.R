library(pracma) # for detrend function
library(plyr)
library(lubridate)
library(ggplot2)

###############################################################################
# Function to calculate if data values exceed (or are lower than) a certain 
# percentile.
is_extreme <- function(data_vec, prob, greater=TRUE, data_subset=NULL, thresholds_only=FALSE) {
    if (is.null(data_subset)) {data_subset<- rep(TRUE, length(data_vec))}
    if (thresholds_only == TRUE) {
        return(quantile(data_vec[data_subset], prob=prob/100, na.rm=TRUE))
    }
    if (greater) {
        return(data_vec > quantile(data_vec[data_subset], prob=prob/100, na.rm=TRUE))
    } else {
        return(data_vec < quantile(data_vec[data_subset], prob=prob/100, na.rm=TRUE))
    }
}

###############################################################################
# Function for formatting text p values
format_p <- function(p_val) {
    if (p_val < .001) {
        return('<.001')
    } else {
        return(round(p_val, digits=3))
    }
}

###############################################################################
# Function for adding equation to ggplot facets, plots modified from:  
# http://bit.ly/11SrhE2, and from http://bit.ly/10zlf8f
eqnfunc_slope <- function(d, model_formula, two_lines=FALSE) {
    m <- lm(formula(model_formula), data=d)
    l <- list(b=format(abs(coef(m)[2]), digits=2),
              pslope=format_p(summary(m)$coefficients[2, 4]))
    if (coef(m)[2] >= 0)  {
        eq <- substitute("slope" == b*","*~~"p"[slope]==~pslope, l)
    } else {
        eq <- substitute("slope" == - b*","*~~"p"[slope]==~pslope, l)
    }
    c(eqn=as.character(as.expression(eq)))
}

###############################################################################
# Function to test for significant linear trend
signif_slope <- function(d, model_formula) {
    m <- lm(formula(model_formula), data=d)
    p <- summary(m)$coefficients[2, 4]
    b <- coef(m)[2]
    if (p < .05) {
        return(data.frame(signif_slope=TRUE, sign=sign(b)))
    } else {
        return(data.frame(signif_slope=FALSE, sign=sign(b)))
    }
}

###############################################################################
# Start data processing
load('pentad_vg_plot_ppt.RData')

vg_plot_ppt <- vg_plot_ppt[vg_plot_ppt$date < as.Date('2014/1/1'), ]

# Add in continent names and countries
sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
vg_plot_ppt <- merge(vg_plot_ppt, sites, by.x='sitecode', by.y='sitecode')
vg_plot_ppt$SiteNamePretty <- ordered(vg_plot_ppt$SiteNamePretty, levels=sites$SiteNamePretty)

#TODO: MUST DETREND

ppt_one <- vg_plot_ppt[vg_plot_ppt$plot_ID == 'VGBBS1', ]

ppt <- ddply(vg_plot_ppt, .(sitecode, plot_ID, plot_num), transform,
             date=date,
             ppt_dt=detrend(ppt),
             pentad=pentad)

ppt <- ddply(ppt, .(sitecode, plot_ID), transform,
             ppt_dt_gt_90=is_extreme(ppt_dt, 90, data_subset=(ppt_dt > 0)),
             ppt_dt_gt_95=is_extreme(ppt_dt, 95, data_subset=(ppt_dt > 0)),
             ppt_dt_gt_99=is_extreme(ppt_dt, 99, data_subset=(ppt_dt > 0)))
ppt$Year <- year(ppt$date)

###############################################################################
# 90th percentiles
gt_90 <- ddply(ppt, .(SiteNamePretty, plot_num, Year), summarize,
               num_days=sum(ppt_dt_gt_90, na.rm=TRUE),
               Continent=continent[1])
gt_90 <- ddply(gt_90, .(SiteNamePretty, Year), summarize,
               num_days=mean(num_days, na.rm=TRUE),
               Continent=Continent[1])
labeldata <- ddply(gt_90, .(SiteNamePretty), eqnfunc_slope, 'num_days ~ order(Year)')
# labeldata$eqn <- gsub('order[(]Year[)]', 'Year', labeldata$eqn)
# labeldata$eqn <- gsub('num_days', 'Num. Days', labeldata$eqn)
pct_90_plot <- ggplot(gt_90, aes(Year, num_days, colour=Continent)) +
    geom_line() + xlab('Date') + facet_wrap(~SiteNamePretty) +
    ylab('Days above 90th percentile') + 
    geom_smooth(method="lm", se=TRUE) +
    scale_y_continuous(breaks=c(0, 3, 6, 9, 12)) +
    theme_grey(base_size=18) + guides(colour=FALSE) +
    transparent_opts
    #geom_text(data=labeldata, aes(x=1970, y=13, label=eqn), parse=TRUE, 
    #          colour='black', hjust=0, size=8) + ylim(c(0,14))
ggsave('precip_90th_above.png', width=width, height=height, dpi=dpi,
       bg='transparent')
# What is the 90th percentile per site?
ddply(ppt, .(SiteNamePretty), summarize, 
      ppt_dt_gt_95=is_extreme(ppt_dt, 95, data_subset=(ppt_dt > 0), thresholds_only=TRUE))

###############################################################################
# 95th percentiles
gt_95 <- ddply(ppt, .(SiteNamePretty, plot_num, Year), summarize, num_days=sum(ppt_dt_gt_95, na.rm=TRUE))
# labeldata <- ddply(gt_95, .(SiteNamePretty, plot_num), eqnfunc_slope, 'num_days ~ order(Year)')
# labeldata$eqn <- gsub('order[(]Year[)]', 'Year', labeldata$eqn)
# labeldata$eqn <- gsub('num_days', 'Num. Days', labeldata$eqn)
pct_95_plot <- ggplot(gt_95, aes(Year, num_days, colour=plot_num)) +
    geom_line() + xlab('Date') + facet_wrap(~SiteNamePretty) +
    ylab('Days above 95th percentile') + 
    geom_smooth(method="lm", se=TRUE) +
    scale_y_continuous(breaks=c(0, 3, 6, 9, 12)) +
    theme_grey(base_size=18) + guides(colour=FALSE) +
    transparent_opts
    #geom_text(data=labeldata, aes(x=1970, y=13, label=eqn), parse=TRUE, 
    #          colour='black', hjust=0, size=8) + ylim(c(0,14))
ggsave('precip_95th_above.png', width=width, height=height, dpi=dpi,
       bg='transparent')
# What is the 95th percentile per site?
ddply(ppt, .(SiteNamePretty), summarize, 
      ppt_dt_gt_95=is_extreme(ppt_dt, 95, data_subset=(ppt_dt > 0), thresholds_only=TRUE))

###############################################################################
# 99th percentiles
gt_99 <- ddply(ppt, .(SiteNamePretty, plot_num, Year), summarize, num_days=sum(ppt_dt_gt_99, na.rm=TRUE))
# labeldata <- ddply(gt_99, .(SiteNamePretty, plot_num), eqnfunc_slope, 'num_days ~ order(Year)')
# labeldata$eqn <- gsub('order[(]Year[)]', 'Year', labeldata$eqn)
# labeldata$eqn <- gsub('num_days', 'Num. Days', labeldata$eqn)
pct_99_plot <- ggplot(gt_99, aes(Year, num_days, colour=plot_num)) +
    geom_line() + xlab('Date') + facet_wrap(~SiteNamePretty) +
    ylab('Days above 99th percentile') + 
    geom_smooth(method="lm", se=TRUE) +
    scale_y_continuous(breaks=c(0, 3, 6, 9, 12)) +
    theme_grey(base_size=18) + guides(colour=FALSE) +
    transparent_opts
    #geom_text(data=labeldata, aes(x=1970, y=13, label=eqn), parse=TRUE, 
    #          colour='black', hjust=0, size=8) + ylim(c(0,14))
ggsave('precip_99th_above.png', width=width, height=height, dpi=dpi,
       bg='transparent')
# What is the 99th percentile per site?
ddply(ppt, .(SiteNamePretty), summarize, 
      ppt_dt_gt_99=is_extreme(ppt_dt, 99, data_subset=(ppt_dt > 0), thresholds_only=TRUE))

###############################################################################
# Num dry days
num_dry <- ddply(ppt, .(SiteNamePretty, plot_num, Year), summarize,
                 num_days=sum(ppt < 5, na.rm=TRUE),
                 Continent=continent[1])
num_dry <- ddply(num_dry, .(SiteNamePretty, Year), summarize,
                 num_days=mean(num_days, na.rm=TRUE),
                 Continent=Continent[1])
num_dry_plot <- ggplot(num_dry, aes(Year, num_days, colour=Continent)) +
    geom_line() + xlab('Date') + facet_wrap(~SiteNamePretty) +
    ylab('Number of dry days') + 
    geom_smooth(method="lm", se=TRUE) +
    theme_grey(base_size=18) + guides(colour=FALSE) +
    transparent_opts
    #geom_text(data=labeldata, aes(x=1970, y=13, label=eqn), parse=TRUE, 
    #          colour='black', hjust=0, size=8) + ylim(c(0,14))
ggsave('precip_num_dry_days.png', width=width, height=height, dpi=dpi,
       bg='transparent')

# How many dry days per site?
ddply(num_dry, .(SiteNamePretty), summarize, 
      mean_num_dry=mean(num_dry, na.rm=TRUE))

# Significant trends?
num_dry_signif <- ddply(num_dry, .(SiteNamePretty), signif_slope,
                       'num_days ~ order(Year)')
num_dry_signif <- ddply(num_dry_signif, .(SiteNamePretty), summarize,
                        num_signif=sum(signif_slope),
                        num_neg=sum(sign == -1),
                        num_zer=sum(sign == 0),
                        num_pos=sum(sign == 1))
num_dry_signif

###############################################################################
# Total annual precip
ppt_ann <- ddply(ppt, .(SiteNamePretty, plot_num, Year), summarize, total=sum(ppt, na.rm=TRUE))
ppt_ann_plot <- ggplot(ppt_ann, aes(Year, total, colour=plot_num)) +
    geom_line() + xlab('Date') + facet_wrap(~SiteNamePretty) +
    ylab('Total annual precipitation (mm)') + 
    geom_smooth(method="lm", se=TRUE) +
    theme_grey(base_size=18) + guides(colour=FALSE) +
    transparent_opts
    #geom_text(data=labeldata, aes(x=1970, y=13, label=eqn), parse=TRUE, 
    #          colour='black', hjust=0, size=8) + ylim(c(0,14))
ggsave('precip_annual_total.png', width=width, height=height, dpi=dpi,
       bg='transparent')

###############################################################################
# Test if slopes are significant at 95%
gt_90_signif <- ddply(gt_90, .(SiteNamePretty), signif_slope,
                      'num_days ~ order(Year)')
gt_90_signif <- ddply(gt_90_signif, .(SiteNamePretty), summarize,
                      num_signif=sum(signif_slope),
                      num_neg=sum(sign == -1),
                      num_zer=sum(sign == 0),
                      num_pos=sum(sign == 1))
gt_90_signif

gt_95_signif <- ddply(gt_95, .(SiteNamePretty, plot_num), signif_slope,
                      'num_days ~ order(Year)')
gt_95_signif <- ddply(gt_95_signif, .(SiteNamePretty), summarize,
                      num_signif=sum(signif_slope),
                      num_neg=sum(sign == -1),
                      num_zer=sum(sign == 0),
                      num_pos=sum(sign == 1))
gt_95_signif

gt_99_signif <- ddply(gt_99, .(SiteNamePretty, plot_num), signif_slope,
                      'num_days ~ order(Year)')
gt_99_signif <- ddply(gt_99_signif, .(SiteNamePretty), summarize,
                      num_signif=sum(signif_slope),
                      num_neg=sum(sign == -1),
                      num_zer=sum(sign == 0),
                      num_pos=sum(sign == 1))
gt_99_signif
