# Masking tools
Command line Python tools to generate gridded country and region masks at any resolution. This includes (1) binary masks, (2) percent weighted by area of gridcell in country/region of interest, (3) percent weighted by population of gridcell in country/region of interest.

Here’s how the code works:
* Install a Python conda environment using the geopandas.yml file. 
* You can then run python at the command line to generate pretty much any mask you can imagine. A few examples below:

This run will generate a mask at 2x2.5 degrees that has values of 1 for all grid cells in the World Bank Middle East and North Africa region, plus Mongolia (ISO3 code):

python make_country_landmask.py -grid "2.0x2.5" -custom "REGION_WB:Middle East & North Africa" -country "MNG" -o mideast_plus_mongolia.npy

This will generate a mask of Europe, cut off a bit west of Russia:

python make_country_landmask.py -grid "2.0x2.5" -continent Europe -lat '33,72' -lon 'm26,37' -o Europe_no_russia.npy

If your outfile ends in .csv, it will save as csv. If .npy, an npy file. If .nc or .nc4, a netcdf file

## Making masks for regions, capturing coastline percentages, and using curvilinear grids

The code is documented, so hopefully you can figure out how to express the masks you want. You can select UN subregions using commands of the type: "SUBREGION:Caribbean,Central America” for the -custom command line argument. This argument would make a mask including both the Caribbean and Central America definitions. 

If your grids are non-standard, you can supply any netcdf file which includes latitude and longitude dimensions for the -grid command line argument.

If you want to make a mask that captures the percent of a gridcell occupied by the country/region of interest, here's an example:

python make_country_landmask.py -grid "1x1" -grid2Agg "4.0x5.0" -country "CAN,USA:50" -o /hpc/group/shindell/ap851/masks/canada_and_alaska_4x5_percent_mask.npy

The idea here is that it first makes a 1x1 mask, then sums that mask up into the 4x5 grid. This provides an estimate of the percent of the bigger grid occupied by your country/region You can pick a finer grid for a more accurate 4x5 percent mask, at the cost of population time.

## Population-weighted masking

If you want to make a mask weighted by population, here's an example:

python make_country_landmask.py -grid '0.5x0.5' -aggByPop True -continent "South America" -country "FRA" -lon "m180,m30" -custom "SUBREGION:Caribbean,Central America" -path2pop /hpc/group/shindell/ap851/masks/gpwv4/gpw_v4_une_atotpopbt_cntm_2pt5_min.nc -o /hpc/group/shindell/ap851/masks/south_and_central_america_0p5x0p5_pop_weighted_mask.nc4

The idea here is that you capture what percent of a given gridcell belongs to your country/countries or region of interest. If a population lays over a border, that will reduce the percentage of the gridcell. Since no one lives off of a coastline, coastal gridcells evaluate to 1.

This feature requires external data. Note that if you want to use the population-weighted feature, you'll need to download the 2.5 arcminute Gridded Population of the World v4 dataset and supply a path to it. You'll also need a land-sea mask at the 2.5 arcminute grid. While the repository does provide you the code to generate this, you can also download a pre-generated mask from the release tab of this repository.

