#!/usr/bin/env sh

# Multirow response variables in Tukey's table
sed '9s/^Cumulative\sEVI\s\&/\\multirow{3}{*}{Cumulative EVI} \&/g' out/lsq_terms.tex |\
	sed '10s/^\s\+Cumulative\sEVI\s\&/\&/g' |\
	sed '11s/^\s\+Cumulative\sEVI\s\&/\&/g' |\
	sed '13s/^Season\slength\s\&/\\multirow{3}{*}{Season length} \&/g' |\
	sed '14s/^\s\+Season\slength\s\&/\&/g' |\
	sed '15s/^\s\+Season\slength\s\&/\&/g' |\
	sed '17s/^Green-up\srate\s\&/\\multirow{3}{*}{Green-up rate} \&/g' |\
	sed '18s/^\s\+Green-up\srate\s\&/\&/g' |\
	sed '19s/^\s\+Green-up\srate\s\&/\&/g' |\
	sed '21s/^Senescence\srate\s\&/\\multirow{3}{*}{Senescence rate} \&/g' |\
	sed '22s/^\s\+Senescence\srate\s\&/\&/g' |\
	sed '23s/^\s\+Senescence\srate\s\&/\&/g' |\
	sed '25s/^Green-up\slag\s\&/\\multirow{3}{*}{Green-up lag} \&/g' |\
	sed '26s/^\s\+Green-up\slag\s\&/\&/g' |\
	sed '27s/^\s\+Green-up\slag\s\&/\&/g' |\
	sed '29s/^Senescence\slag\s\&/\\multirow{3}{*}{Senescence lag} \&/g' |\
	sed '30s/^\s\+Senescence\slag\s\&/\&/g' |\
	sed '31s/^\s\+Senescence\slag\s\&/\&/g' > out/lsq_terms_fmt.tex
