# "group by" functions
# for those new to R (or pandas) I suggest reading "implicit loops" in the "Learning Stats with R" document (p.267 of the pdf)

# For many of you, you'll want to use loops. If you're using R you should rarely need loops.
# Here are some helpful examples.

# sapply - use when you want to apply a function to each item inside an entire array/list/column/dataframe
## typically the equivalent of simply passing the function INSIDE at the end of "with"
readHumdrum('HumdrumData/BachChorales/chor001.krn') -> chor1
## example syntax: within(data, format(Token), function(.))
## "as" humdrum data:
within(chor1, kern(Token), toupper(.))
## "as" R data
with(chor1[,1], kern(Token), nchar(.))

# often we want to apply some function over groups or subgroups of our data (e.g., average test scores by Grade level)
# tapply() and by()
## Average rhythmic duration (rhythmic value) by part. Use argument "by" to apply a grouping
with(chor1, duration(Token), by = Spine, mean(.))

## Average duration value (rhythm value) per measure (across all spines)
with(chor1, duration(Token), by = Bar, mean(.))

## find any notes higher than "cc" in the Soprano part using boolean indexing
with(chor1[,4], semits(Token) > 12) # outputs boolean mask

### to see the 'humdrum' view you could use subset:
subset(chor1, semits(Token) > 12 & Spine == 4)

### to see only the tokens themselves and not the boolean mask:
with(chor1, subset = semits(Token) >= 12 & Spine == 4, semits(Token))

# Example read-in MCFlow dataset (local corpus downloaded from github)
mcflow <- readHumdrum('~/Documents/MyGits/Public/MCFlow/Humdrum/*.rap')

## see which Fields are available by default:
fields(mcflow)

## examine stress (2nd col) in Verses compared to Chorus:
### first need to "line up" stress and recip. Use "foldExclusive"
### (note that rests are coded as "." in stress (they are neither 0 nor 1))
mcstress <- foldExclusive(mcflow, 'stress', 'recip')
with(mcstress, table(Stress, Token))

### now let's group by Formal label (section)
with(mcstress, table(Stress, Token), by = Formal)

### verses were numbered individually... to collapse them and count again, do:
gsub('[0-9]', '', Formal)
with(mcstress, table(Stress, Token), by = sub('[0-9]', '', Formal))
