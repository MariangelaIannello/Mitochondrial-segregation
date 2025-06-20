#!/bin/bash
if [ -d back-translated.fas ]
	then rm -r back-translated.fas
	fi
mkdir back-translated.fas

if [ -d back-translated_masked.fas ]
	then rm -r back-translated_masked.fas
	fi
mkdir back-translated_masked.fas

if [ ! -d R ]
	then mkdir R
	fi

cd R
Rscript ../back-translated_masked.R
cd ..
