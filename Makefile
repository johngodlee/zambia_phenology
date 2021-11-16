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

# Data prep
$(DATDIR)/plots.rds $(DATDIR)/sites_loc.rds $(DATDIR)/plot_id_lookup.rds $(OUTDIR)/data_prep_vars.tex $(DATDIR)/ba_clust_mat.rds : data_prep.R functions.R $(DATDIR)/plots_v2.12.csv
	@echo Data preparation
	Rscript $<
# Feeds into modis_get.R, trmm_get.R, try_get.R, vipphen.R
# but too big to run in main Makefile

# MODIS extract
$(DATDIR)/plots_phen.rds $(IMGDIR)/ts_example.pdf $(OUTDIR)/modis_extract_vars.tex : modis_extract.R functions.R $(DATDIR)/plots.rds $(DATDIR)/vipphen.rds $(DATDIR)/evi.rds modis_get.R
	@echo MODIS data extraction
	Rscript $<

# TRMM extract
$(DATDIR)/plots_trmm.rds $(OUTDIR)/trmm_extract_vars.tex : trmm_extract.R functions.R $(DATDIR)/plots_phen.rds $(DATDIR)/trmm.rds trmm_get.R
	@echo TRMM data extraction
	Rscript $<

# TRY extract
$(DATDIR)/plots_try.rds : try_extract.R $(DATDIR)/plots_trmm.rds $(DATDIR)/ba_clust_mat.rds try_get.R $(DATDIR)/traits_sp.rds
	@echo TRY data extraction
	Rscript $<

# Diversity
$(DATDIR)/plots_div.rds $(OUTDIR)/clust_summ.tex $(OUTDIR)/diversity_vars.tex : diversity.R functions.R $(DATDIR)/plots_try.rds $(DATDIR)/plot_id_lookup.rds $(DATDIR)/ba_clust_mat.rds
	@echo Diversity data
	Rscript $<

# Analysis prep
$(DATDIR)/plots_anal.rds $(IMGDIR)/phen_dens_clust.pdf $(IMGDIR)/plot_loc.pdf $(OUTDIR)/analysis_vars.tex : analysis.R functions.R $(DATDIR)/plots_div.rds $(DATDIR)/vipphen_stack.rds 
	@echo Analysis preparation
	Rscript $<

# Models
$(IMGDIR)/mod_slopes.pdf $(IMGDIR)/mod_marg.pdf $(OUTDIR)/mod_stat.tex $(OUTDIR)/all_mod_sel.tex $(OUTDIR)/models_vars.tex : models.R functions.R $(DATDIR)/plots_anal.rds 
	@echo Models
	Rscript $<

# Compile all latex variables to one file
$(OUTDIR)/vars.tex : $(OUTDIR)/data_prep_vars.tex\
	$(OUTDIR)/modis_extract_vars.tex\
	$(OUTDIR)/trmm_extract_vars.tex\
	$(OUTDIR)/diversity_vars.tex\
	$(OUTDIR)/analysis_vars.tex\
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

# Re-create time series and trait values, needs external data
get : 
	@echo Extract raw time series data
	Rscript zambia_clim_get.R
	Rscript modis_get.R
	Rscript trmm_get.R
	Rscript try_get.R
	Rscript vipphen.R

