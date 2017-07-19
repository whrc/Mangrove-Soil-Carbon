# Kylen Solvik, ksolvik@whrc.org
# Short script to model all soil profiles using linear, log-linear, log-log, and quadratic models.
# Finds pvalues for each model in order to compare model types against each other. 
# You may encounter "perfect fit" warnings for sites with only ~3 horizons. 

# Read soil profile CSV
soil_csv = read.csv("./soil_profiles.csv")

# Set path to output
out_csv = "./depth_model_outputs.csv"

# Set up empty dataframe for saving outputs
df_cols = c("SiteName","SiteNum","LinearP","LinearDirection","LogLinearP","LogLinearDirection",
          "LogLogP","LogLogDirection","QuadP")
out_df=data.frame(matrix(ncol=length(df_cols),nrow=0))
colnames(out_df) = df_cols

# Simple function to get the direction of the trend line slope
get_slope_dir = function(slope) {
  if(slope>0) {
    return(1)
  }
  if(slope<0) {
    return(-1)
  }
  if(slope==0) {
    return(0)
  }
}

# Loop over all sites
for (sitenum in unique(soil_csv$SiteNum)) {

    # Get horizons
    horizons=soil_csv[soil_csv$SiteNum==sitenum,]

    # Only proceed if there is more than 2 horizons
    if(nrow(horizons)>2) {
      out_row=c(levels(horizons$Site.name[1])[horizons$Site.name[1]],sitenum)
      
      # Get depths and carbon densities 
      depth=horizons[,"mid_point"]
      cd=horizons[,"CD_calc"]
      # For log-log model, exclude any horizons where carbon density == 0
      nonzero_cd=horizons[horizons$CD_calc>0,"CD_calc"]
      nonzero_cd_depth=horizons[horizons$CD_calc>0,"mid_point"]
      
      # Linear model
      m=lm(cd ~ depth)
      direction=get_slope_dir(m$coefficients[2])
      pval=summary(m)$coefficients[2,4]
      out_row=c(out_row,pval,direction)
      
      # Log-Linear model
      m=lm(cd ~ log(depth))
      direction=get_slope_dir(m$coefficients[2])
      pval=summary(m)$coefficients[2,4]
      out_row=c(out_row,pval,direction)
      
      # Log-Log model
      m=lm(log(nonzero_cd) ~ log(nonzero_cd_depth))
      direction=get_slope_dir(m$coefficients[2])
      pval=summary(m)$coefficients[2,4]
      out_row=c(out_row,pval,direction)
      
      # Quadratic model
      m=lm(cd ~ poly(depth,2))
      fstat=summary(m)$fstatistic
      pval=pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE) 
      out_row=c(out_row,pval)

      # Append model results to output dataframe
      out_df[nrow(out_df)+1,]=out_row
    }    
}

# Write output dataframe to csv
write.csv(out_df,out_csv)

print("All Done")
