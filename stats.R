#cols: consensus_compound_id,num_peaks,num_times_fragmented,num_times_matched    
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname='mss', user='bmartin')

# What degree of fragmentation (num peaks) is required for identification ?


# Get the statistics about attempts at identification of compounds.
rs = dbSendQuery(con, 'SELECT * FROM cpd_ident_stats')
df = fetch(rs, n=-1)

#---------------------------------------------------------------------
#  With the distribution of the proportion of peaks fragmented,
#  find 20th percentile.
#---------------------------------------------------------------------
frag_20_pctile = quantile( df$num_times_fragmented/df$num_peaks, 0.2)

#---------------------------------------------------------------------
#  With the distribution of the proportion of fragmented peaks that are
#  identified, find 20th percentile.
#---------------------------------------------------------------------
match_20_pctile = quantile( df$num_times_matched/df$num_times_fragmented, 0.2)

#---------------------------------------------------------------------
#  Since there are so many compounds that were identified only once,
#  perhaps the statistics are biased a little bit from some false identified.
#---------------------------------------------------------------------
df2 = subset(df, num_times_matched > 1)

mod_frag_20_pctile = quantile( df2$num_times_fragmented/df2$num_peaks, 0.2)
mod_match_20_pctile = 
   quantile( df2$num_times_matched/df2$num_times_fragmented, 0.2)

frag_20_pctile
match_20_pctile
mod_frag_20_pctile
mod_match_20_pctile

