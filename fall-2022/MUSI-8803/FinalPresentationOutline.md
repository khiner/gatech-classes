## Chord-Ratio Methodology 

### Step 1: Understanding our data

We used several methods to get a broad understanding of our chord data before formally testing our hypothesis.

(Show some charts here with things like total chords, distribution of most popular chords, & other high-level data exploration stuff to convey a sense of what we're working with.)

### Chord-Ratio Methodology I (ratio of songs with _at least one_ progression usage)

Initially, we simply graphed how many times the I-V-vi-IV progression occurred _at least once_ for each year (1991-2022).

More specifically, for each year $y$, we calculated the ratio

$R_y = \dfrac{\text{Number of songs containing at least one I-V-vi-IV in year } y}{\text{Number of songs first appearing in the Billboard Hot 100 in year } y}$

(Chart showing this)

### (Slide 2) Chord-Ratio Methodology I (ratio of songs with _at least one_ progression usage)

(Same chart)

Upon first glance, it’s clear that I-V-vi-IV usage peaked around 2008 and has steadily declined since.

However, we weren’t convinced this was the best way to formally test our hypothesis of declining popularity.
After all, a progression is arguably more "popular" if songs repeat it many times vs. only using it once.

### Chord-Ratio Methodology II (ratio of quadgram occurrences)

The method we use in our final analysis is to calculate, for each year, the ratio of the total number of I-V-vi-IV quadram occurrences (with rotations), over the total number of quadgrams:

$R_y = \dfrac{\text{Number of I-V-vi-IV occurrences in year } y}{\text{Total number of quadgrams across all songs in year } y}$

(Chart showing this)

### (Slide 2) Chord-Ratio Methodology II (ratio of quadgram occurrences)

(Side-by-side of Method 1 chart and Method 2 chart)

There are some clear differences between the graphs using these two methods, but we can see the overall trend it very similar.


## Other interesting findings

  * The quadram data shows a clear declining trend in 4-chord progression utilization since 1991 (both in the number of unique 4-chord progressions, and in total 4-chord progression usage).
    - This may be partly because of the proliferation of genres such as hip-hop and rap that are primarily loop-based.
    (Show chart)
  * ...
