# Compile phenology manuscript

# Define variables
TEXFILE  = phenology
IMGDIR   = ./img
OUTDIR   = ./out
DATDIR   = ./dat

# Include .pdf here to ensure it is always built
# latexmk always run, make cannot easily track dependencies in .aux, .bib etc.
.PHONY : $(TEXFILE).pdf all clean get ts

# Depends on final PDF, which starts dependency chain
all : $(TEXFILE).pdf 

# R scripts

# Plot data prep
$(DATDIR)/plots.rds $(DATDIR)/trees.rds $(OUTDIR)/prep_vars.tex : prep.R $(DATDIR)/plots_v2.12.csv $(DATDIR)/stems_iluaii_v2.12.csv tex_func.R
	@echo Plot data preparation
	Rscript $<

# Create abundance matrices
$(DATDIR)/ba_clust_mat.rds $(DATDIR)/ab_plot_mat.rds : abund.R $(DATDIR)/trees.rds
	@echo Build abundance matrices
	Rscript $<

# Land cover classification extract
$(DATDIR)/zambia_landcover/lcc.tif : lcc.R 
	@echo Prepare Zambia land cover classification
	Rscript $<

# BioClim extract
$(DATDIR)/bioclim.rds $(DATDIR)/bioclim_zambia.rds : bioclim.R $(DATDIR)/plots.rds $(DATDIR)/africa_countries/africa.shp $(DATDIR)/wc2.1_30s_bio/*.tif
	@echo BioClim data extraction
	Rscript $<

# TRMM extract
$(DATDIR)/trmm.rds $(OUTDIR)/trmm_vars.tex : trmm.R $(DATDIR)/plots.rds $(DATDIR)/trmm_ts.rds tex_func.R
	@echo TRMM data extraction
	Rscript $<

# MODIS extract
$(DATDIR)/modis.rds $(OUTDIR)/modis_vars.tex : modis.R $(DATDIR)/plots.rds $(DATDIR)/trmm.rds $(DATDIR)/africa_countries/africa.shp tex_func.R
	@echo MODIS data extraction
	Rscript $<

# Diversity
$(IMGDIR)/plot_dist_hist.pdf $(IMGDIR)/basal_area_dom_hist.pdf $(IMGDIR)/ward_sil_mean.pdf $(IMGDIR)/ward_sil.pdf $(DATDIR)/taxon.rds $(DATDIR)/div.rds $(OUTDIR)/clust_summ.tex $(DATDIR)/indval.rds $(OUTDIR)/diversity_vars.tex : div.R $(DATDIR)/plots.rds $(DATDIR)/trees.rds $(DATDIR)/ba_clust_mat.rds $(DATDIR)/ab_plot_mat.rds tex_func.R plot_func.R
	@echo Diversity metrics
	Rscript $<

# Create visualisations
$(IMGDIR)/site_map.pdf $(IMGDIR)/site_clim.pdf $(IMGDIR)/site_loc.pdf $(IMGDIR)/phen_bivar.pdf $(OUTDIR)/clust_summ.tex $(IMGDIR)/dens_lag.pdf $(IMGDIR)/phen_dens_clust.pdf : vis.R $(DATDIR)/plots.rds $(DATDIR)/div.rds $(DATDIR)/modis.rds $(DATDIR)/bioclim.rds $(DATDIR)/indval.rds $(DATDIR)/africa_countries/africa.shp $(DATDIR)/bioclim_zambia.rds $(DATDIR)/zambia_landcover/lcc.tif plot_func.R
	@echo Visualisation and descriptive tables
	Rscript $<

# Models
$(IMGDIR)/mod_slopes.pdf $(IMGDIR)/mod_marg.pdf $(OUTDIR)/mod_stat.tex $(OUTDIR)/tukey_terms.tex $(OUTDIR)/models_vars.tex : models.R $(DATDIR)/plots.rds $(DATDIR)/div.rds $(DATDIR)/modis.rds $(DATDIR)/bioclim.rds plot_func.R tex_func.R
	@echo Models
	Rscript $<

# Compile all latex variables to one file
$(OUTDIR)/vars.tex : $(OUTDIR)/prep_vars.tex\
	$(OUTDIR)/trmm_vars.tex\
	$(OUTDIR)/modis_vars.tex\
	$(OUTDIR)/diversity_vars.tex\
	$(OUTDIR)/models_vars.tex
	@echo Compile LaTeX variables
	cat $^ > $@

# Convert .drawio to .pdf
$(IMGDIR)/schematic.pdf : drawio/schematic.drawio
	@echo Compile drawio images
	./drawio_export.sh $< $@

# Compile main tex file and show errors
$(TEXFILE).pdf : $(TEXFILE).tex\
	$(OUTDIR)/vars.tex\
	$(IMGDIR)/plot_loc.pdf\
	$(IMGDIR)/ts_example.pdf\
	$(IMGDIR)/mod_slopes.pdf\
	$(IMGDIR)/mod_marg.pdf\
	$(IMGDIR)/schematic.pdf\
	$(IMGDIR)/phen_dens_clust.pdf\
	$(OUTDIR)/all_mod_sel.tex\
	$(OUTDIR)/lsq_terms.tex\
	$(OUTDIR)/clust_summ.tex
	@echo Compile manuscript
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make -bibtex $<

# Clean up stray intermediary files
clean :
	@echo Clean LaTeX files
	latexmk -C

# Re-create time series and trait values, needs external data, takes a long time
get : 
	@echo Extract raw time series data
	Rscript trmm_get.R
	Rscript modis_get.R

