# Mangrove-Soil-Carbon

R code and spatial predictions for WHRC-TNC project **mapping soil carbon stocks under mangrove forests of the world**. Output maps at 30 m resolution are available for download from Google Earth Engine following the links below. To learn more about soil organic carbon mapping visit **[this tutorial](http://gsif.isric.org/doku.php/wiki:soil_organic_carbon)**.

![Global map of SOCS under mangrove forests of the world](https://github.com/whrc/Mangrove-Soil-Carbon/blob/master/img/mSOC_combinedLayout_sm.jpg "Output predictions of soil organic carbon stock under mangrove forests of the world.")

Files in this repository include:

* `rmatrix_OCD.csv` = regression matrix used to fit a spatial prediction model to map soil organic carbon density,
* `mangroves_SOC_points.gpkg` = soil carbon mangrove forests data base in GeoPackage format,

Input soil data for mangove forests, predictions at 100 m resolution (for 0-1 m and 0-2 m depths) and tiled 30 m resolution predictions can be downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/OCYUIT without restrictions.

Soil organic carbon under mangrove forests grids at 30 m **can also be downloaded from** https://code.earthengine.google.com/?asset=users/gfiske/SOCS/Mangroves_SOC_0_100cm_30m

Please cite as:

* Sanderman, J., Hengl, T., Fiske, G., Solvik, K., Adame, M., Benson, L., et al., (2017?) **"A global map of mangrove forest soil carbon at 30 m spatial resolution"**. submitted to Environmental Research Letters (Special Issue).


