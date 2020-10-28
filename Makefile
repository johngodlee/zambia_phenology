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
$(DATDIR)/plots.rds $(DATDIR)/plot_id_lookup.rds $(OUTDIR)/data_prep_vars.tex : data_prep.R $(DATDIR)/plots_v2.7.csv
	Rscript $<

# MODIS extract
$(IMGDIR)/ts_example.pdf $(DATDIR)/plots_phen.rds $(OUTDIR)/modis_extract_vars.tex : modis_extract.R functions.R $(DATDIR)/plots.rds $(DATDIR)/vipphen.rds $(DATDIR)/evi.rds
	Rscript $<

# TRMM extract
$(DATDIR)/plots_trmm.rds $(OUTDIR)/trmm_extract_vars.tex : trmm_extract.R functions.R $(DATDIR)/plots_phen.rds $(DATDIR)/trmm.rds
	Rscript $<

# Diversity
$(DATDIR)/plots_div.rds $(OUTDIR)/indval.tex $(OUTDIR)/diversity_vars.tex : diversity.R functions.R $(DATDIR)/plots_trmm.rds $(DATDIR)/plot_id_lookup.rds $(DATDIR)/stems_latest_v2.7.csv
	Rscript $<

# Analysis prep
$(DATDIR)/plots_anal.rds $(OUTDIR)/analysis_vars.tex $(IMGDIR)/phen_dens_clust.pdf $(IMGDIR)/plot_loc.pdf : analysis.R functions.R $(DATDIR)/plots_div.rds $(DATDIR)/vipphen_stack.rds 
	Rscript $<

# Models
$(IMGDIR)/mod_slopes.pdf $(IMGDIR)/mod_marg.pdf $(OUTDIR)/mod_stat.tex $(OUTDIR)/all_mod_sel.tex : models.R $(DATDIR)/plots_anal.rds 
	Rscript $<

# Compile all latex variables to one file
$(OUTDIR)/vars.tex : $(OUTDIR)/data_prep_vars.tex $(OUTDIR)/modis_extract_vars.tex $(OUTDIR)/trmm_extract_vars.tex $(OUTDIR)/diversity_vars.tex $(OUTDIR)/analysis_vars.tex
	cat $^ > $@

# Convert .drawio to .pdf
$(IMGDIR)/schematic.pdf : drawio/schematic.drawio
	./drawio_export.sh $< $@

# Format some tables
$(OUTDIR)/indval_fmt.tex $(OUTDIR)/all_mod_sel_fmt.tex : $(OUTDIR)/indval.tex $(OUTDIR)/all_mod_sel.tex 
	./table_fmt.sh $< $@

# Compile main tex file and show errors
$(TEXFILE).pdf : $(TEXFILE).tex\
	$(OUTDIR)/vars.tex\
	$(IMGDIR)/plot_loc.pdf\
	$(IMGDIR)/ts_example.pdf\
	$(IMGDIR)/mod_slopes.pdf\
	$(IMGDIR)/mod_marg.pdf\
	$(IMGDIR)/nsca.pdf\
	$(IMGDIR)/schematic.pdf\
	$(IMGDIR)/phen_dens_clust.pdf\
	$(OUTDIR)/all_mod_sel.tex\
	$(OUTDIR)/indval_fmt.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make -bibtex $<

# Only run time series data getting
get : 
	Rscript modis_get.R
	Rscript trmm_get.R

# Clean up stray intermediary files
clean :
	latexmk -C

# Re-create time series
ts : 
	Rscript modis_get.R
	Rscript trmm_get.R
	Rscript vipphen.R

