# Compile phenology manuscript

# Define variables
TEXFILE = phenology
IMGDIR = ./img
OUTDIR = ./out
DATDIR = ./dat

# Include .pdf here to ensure it is always built
# latexmk always run, make cannot easily track dependencies in .aux, .bib etc.
.PHONY : $(TEXFILE).pdf all clean get ts

# Depends on final PDF, which starts dependency chain
all : $(TEXFILE).pdf 

# Plot data prep
$(DATDIR)/plots.rds $(DATDIR)/trees.rds $(OUTDIR)/prep_vars.tex : prep.R $(DATDIR)/plots_v2.12.csv $(DATDIR)/stems_iluaii_v2.12.csv tex_func.R
	@echo Plot data preparation
	Rscript $<

# Create abundance matrices
$(DATDIR)/ba_clust_mat.rds $(DATDIR)/ab_plot_mat.rds : abund.R $(DATDIR)/trees.rds
	@echo Build abundance matrices
	Rscript $<

# BioClim extract
$(DATDIR)/bioclim.rds $(DATDIR)/bioclim_zambia.rds : bioclim.R $(DATDIR)/plots.rds $(DATDIR)/africa_countries/africa.shp $(DATDIR)/wc2.1_30s_bio/*.tif
	@echo BioClim data extraction
	Rscript $<

# Season length raster extract
$(DATDIR)/slen_zambia.tif : slen.R $(DATDIR)/africa_countries/africa.shp $(DATDIR)/modis_zambia.tif $(DATDIR)/zambia_landcover/lcc.tif
	@echo Zambia season length raster extraction
	Rscript $<

# Time series statistics 
$(DATDIR)/stat_avg.rds $(DATDIR)/stat_all.rds $(OUTDIR)/stat_vars.tex : stat.R $(DATDIR)/plots.rds $(DATDIR)/trmm_ts.rds $(DATDIR)/era_ts.rds $(DATDIR)/modis_ts.rds tex_func.R
	@echo Time series statistics extraction
	Rscript $<

# Diversity
$(IMGDIR)/plot_dist_hist.pdf $(IMGDIR)/basal_area_dom_hist.pdf $(IMGDIR)/ward_sil_mean.pdf $(IMGDIR)/ward_sil.pdf $(DATDIR)/taxon.rds $(DATDIR)/div.rds $(OUTDIR)/clust_summ.tex $(DATDIR)/indval.rds $(OUTDIR)/diversity_vars.tex : div.R $(DATDIR)/plots.rds $(DATDIR)/trees.rds $(DATDIR)/ba_clust_mat.rds $(DATDIR)/ab_plot_mat.rds tex_func.R plot_func.R
	@echo Diversity metrics
	Rscript $<

# Create visualisations
$(IMGDIR)/box_facet_map_mat.pdf $(IMGDIR)/site_map.pdf $(IMGDIR)/phen_bivar.pdf : vis.R $(DATDIR)/plots.rds $(DATDIR)/div.rds $(DATDIR)/stat_all.rds $(DATDIR)/stat_avg.rds $(DATDIR)/bioclim.rds $(DATDIR)/indval.rds $(DATDIR)/africa_countries/africa.shp $(DATDIR)/bioclim_zambia.rds plot_func.R $(DATDIR)/trmm_ts.rds $(DATDIR)/bioclim_zambia.rds $(DATDIR)/slen_zambia.tif
	@echo Visualisation and descriptive tables
	Rscript $<

# Models
$(IMGDIR)/mod_slopes.pdf $(IMGDIR)/mod_marg.pdf $(OUTDIR)/dredge.tex $(OUTDIR)/models_vars.tex : models.R $(DATDIR)/plots.rds $(DATDIR)/div.rds $(DATDIR)/stat_all.rds plot_func.R tex_func.R
	@echo Models
	Rscript $<

# Compile all latex variables to one file
$(OUTDIR)/vars.tex : $(OUTDIR)/prep_vars.tex\
	$(OUTDIR)/diversity_vars.tex\
	$(OUTDIR)/models_vars.tex\
	$(OUTDIR)/stat_vars.tex
	@echo Compile LaTeX variables
	cat $^ > $@

# Convert .drawio to .pdf
$(IMGDIR)/schematic.pdf : drawio/schematic.drawio
	@echo Compile drawio images
	./drawio_export.sh $< $@

# Convert time series illustation from .svg to .png
$(IMGDIR)/ts_illus.png : drawio/ts_illus.svg
	@echo Compile svg images
	rsvg-convert -w 2586 -h 1709 $@ $<

# Compile main tex file and show errors
$(TEXFILE).pdf : $(TEXFILE).tex\
	$(OUTDIR)/vars.tex\
	$(OUTDIR)/clust_summ.tex\
	$(OUTDIR)/dredge_adj.tex\
	$(IMGDIR)/schematic_rev_2023-09-01.pdf\
	$(IMGDIR)/site_map.pdf\
	$(IMGDIR)/box_facet_map_mat.pdf\
	$(IMGDIR)/ts_illus.png\
	$(IMGDIR)/mod_slopes.pdf\
	$(IMGDIR)/mod_marg.pdf\
	$(IMGDIR)/boxplots.pdf\
	$(IMGDIR)/phen_bivar.pdf
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
	Rscript era_get.R
