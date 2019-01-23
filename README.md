# Mangrove-Soil-Carbon

R code and spatial predictions for WHRC-TNC project [**mapping soil carbon stocks under mangrove forests of the world**](https://doi.org/10.1088/1748-9326/aabe1c). Output maps at 30 m resolution are available for download from Google Earth Engine following the links below. To learn more about soil organic carbon mapping visit **[this tutorial](https://envirometrix.github.io/PredictiveSoilMapping/SOC-chapter.html)**.

![Global map of SOCS under mangrove forests of the world](https://github.com/whrc/Mangrove-Soil-Carbon/blob/master/img/mSOC_combinedLayout_sm.jpg "Output predictions of soil organic carbon stock under mangrove forests of the world.")

Data can be visualized as a web-based map here: https://storage.googleapis.com/gfiske1/global_mangrove/index_w_slider.html. Updated version (October, 2018) of this map can be browsed [here](https://www.arcgis.com/apps/MapSeries/index.html?appid=fe214a492f114bde8b3aa1d54ef23224). Spatial predictions can be downloaded from https://doi.org/10.5281/zenodo.1469347

Files in this repository include:

* `rmatrix_OCD.csv` = regression matrix used to fit a spatial prediction model to map soil organic carbon density,
* `mangroves_SOC_points.gpkg` = soil carbon mangrove forests data base in GeoPackage format,

Input soil data for mangove forests, predictions at 100 m resolution (for 0-1 m and 0-2 m depths) and tiled 30 m resolution predictions can be downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/OCYUIT without restrictions. In the October 2018 update we also use the point data from [Rovai et al. (2018)](https://doi.org/10.1038/s41558-018-0162-5) (see `41558_2018_162_MOESM2_ESM.csv`).

Soil organic carbon under mangrove forests grids at 30 m **can also be downloaded from** https://code.earthengine.google.com/?asset=users/gfiske/SOCS/Mangroves_SOC_0_100cm_30m_Dec_15_2017

Please cite as:

* Sanderman, J., Hengl, T., Fiske, G., Solvik, K., Adame, M., Benson, L., et al., (2018) [**"A global map of mangrove forest soil carbon at 30 m spatial resolution"**](https://doi.org/10.1088/1748-9326/aabe1c). Environmental Research Letters, 13(5), 055002.


