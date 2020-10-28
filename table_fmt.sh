#!/usr/bin/env sh

# Multirow indicator species table clusters
sed '9s/^1\s\&/\\multirow{3}{*}{1} \&/g' out/indval.tex |\
	sed '10s/^\s\+1\s\&/\&/g' |\
	sed '11s/^\s\+1\s\&/\&/g' |\
	sed '13s/^2\s\&/\\multirow{3}{*}{2} \&/g' |\
	sed '14s/^\s\+2\s\&/\&/g' |\
	sed '15s/^\s\+2\s\&/\&/g' |\
	sed '17s/^3\s\&/\\multirow{3}{*}{3} \&/g' |\
	sed '18s/^\s\+3\s\&/\&/g' |\
	sed '19s/^\s\+3\s\&/\&/g' > out/indval_fmt.tex

# 
