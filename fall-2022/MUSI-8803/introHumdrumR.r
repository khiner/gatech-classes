# Intro to humdrumR - Fall 2022
###############################
## Installations: uncomment install lines as necessary 
# install.packages('devtools') 
## Install humdrumR from github:
# devtools::install_github("Computational-Cognitive-Musicology-Lab/humdrumR", build_vignettes = TRUE)

## Load libraries
library(humdrumR)

## setwd
setwd(humdrumRroot)

## read-in a single humdrum file

readHumdrum('HumdrumData/BachChorales/chor001.krn') -> chor1

## Equivalent of 'extract'
ex1 <- chor1[[,4]] # index by column position -- double index required
subset(chor1, Spine==4) # same as prior 

# table (count) all rhythms in the melody
with(ex1, recip(Token)) |> table()

# same as:
with(ex1, recip(Token), table(.))

# table (count) all scale degrees (absolute - no direction/8ve) in the melody 
## equivalent of 'deg'
with(ex1, deg(Token,simple=TRUE) |> table())
# table (count) all melodic intervals in the melody
## equivalent of 'mint'
with(ex1, mint(Token)) |> table()


# load local corpus (#change to your own local file path with the bhchorales!!)
mychors <- readHumdrum('~/Desktop/Music_Files/KernFiles/ClassicalMusic/bach/bhchorale/chor')
## equivalent of 'census'
census(mychors)

# what chords most commonly end a phrase in the bhchorale corpus?     
with(mychors[[,1]], subset=Token%~%";", table(Token)) #like 'grepl' in R

# what are the top 5 most common chord bigrams?
with(mychors[[,1]], paste(Token, Token[lag=1])) |> table() |> sort() |> tail(n=5) 

# get rid of semicolon & ask same question
with(mychors[[,1]], gsub(";", "", Token), paste(., .[lag=1])) |> table() |> sort() |> tail(n=5) 
# (notice by default it will not cross file boundaries)

# graph all pitches in 'semits' in the melodic part
with(mychors[[,5]], semits(Token), hist(.))#(., main="Pitch Distribution in Semits"))

# graph all pitches in 'semits' by each part
with(mychors[[,3:6]], semits(Token) |> hist(), by=Spine)

# table the interval distribution of the melody part in semitones
with(mychors[[,5]], mint(Token, deparser=semits), table(.))

# graph the interval distribution of the melody part in semitones
with(mychors[[,5]], mint(Token, deparser=semits, initial=NULL)) |> table() |> barplot()

# examine melody over chord (simple scale degree over chords)
# requires creating new 'harm' field so that harm and kern 'views' can overlap
# requires overwriting (save new output - this altered dimensions of columns)
mychors <- foldExclusive(mychors, "harm", "kern")
with(mychors[[,5]], deg(Token, simple=TRUE), paste(Harm, .), table(.)) |> sort()
