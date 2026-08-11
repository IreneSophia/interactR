# `interactR`

**An easy-to-use R package for extracting behavioural features from social interactions**

<img src="https://github.com/IreneSophia/interactR/blob/main/logo/InteractRLogo_2026-08-01.svg" width="200" alt="logo">

`interactR` provides an accessible workflow for researchers interested in studying social interactions. The package extracts a range of behavioural features from data derived from video and audio recordings or collected in the virtual reality environment `VERSE`.  The pipeline produces analysis-ready measures with minimal coding.

Whether you are interested in dwell times to areas of interest, emotional facial expressions, kinematic analysis or interpersonal synchrony, `interactR` provides a unified pipeline that streamlines preprocessing and feature extraction.

This R package is currently under active development and is being maintained by **Irene Sophia Plank**, with new features and improvements being added on an ongoing basis.

## Features

### Gaze Patterns and Joint Attention

Measure where participants direct their attention and identify moments of shared attention between interaction partners.

Outputs include:

- Dwell times to specified areas of interest
- Joint attention between two interaction partners

### Emotional Facial Expressions (EFE)

Detect and classify emotional facial expressions, either with your own categorisation or using the emotion categorisation framework proposed by Aldenhoven et al. (2026).

Outputs include:

- Emotion facial expressions over time
- Frequency of emotional facial expressions

### Head Gesture Detection

Automatically detect head movements using a zero-crossing approach, inspired by [Hale, Ward et al. (2020)](https://doi.org/10.1007/s10919-019-00320-3).

Outputs include:

- Nodding events
- Head-shaking events
- Frequency measures
- Plots showing extracted Head Gestures across time

### Speech Patterns and Turn Taking

If the interactions was also captured with audio recordings, then speech and turn-taking features can be extracted as well. This data can either be captured within the virtual reality environment `VERSE` or as audio recordings which are processed with the attached `praat` script: `praatScript/featSpeech.praat`. `praat` is a free, open-source computer software package widely used for speech analysis and synthesis in phonetics.

The `praat` script uses the [`uhm-o-meter`](https://sites.google.com/view/uhm-o-meter/home) created by [De Jong, N.H., Pacilly, J., & Heeren, W. (2021)](https://www.google.com/url?q=https%3A%2F%2Fdoi.org%2F10.1080%2F0969594X.2021.1951162&sa=D&sntz=1&usg=AOvVaw0lcSmN4aCMwjZ4iRQRR6Ab). This algorithm detects syllables, speaking and sounding instances. The `praat` script was extended to also extract pitch and intensity from the audio recordings.

Regardless of whether you use data collected in `VERSE` or recordings processed with `praat`, `interactR` offers functions to further process this data and extract relevant speech and turn-taking features. This pipeline was described in [Plank et al. (2023)](https://doi.org/10.3389/fpsyt.2023.1257569).

Outputs include: 

- Turn-taking gaps
- Number of turns
- Pitch and intensity variance
- Articulation rate
- Silence-to-turn ratio
- Phonation balance

### Interpersonal Synchrony (IPS)

Measure the behavioural adaptation between two interaction partners using two different methods: 

- **Wavelet Coherence (WTC)**, localised coherence in both time and frequency space. WTC decomposes the signals into distinct periodicities (e.g., distinguishing fast from slow nods) and evaluates how strongly the two time series co-vary at each specific frequency band across the interaction.
- **Windowed Lagged Cross-Correlation (WLCC)**, localised linear correlation in the time domain only, evaluated over moving temporal windows. It treats the signal as a whole aggregate waveform without separating it into underlying frequency components, focusing instead on how the overall shape of the signal matches up when shifted in time.

Outputs include:

- Peak synchrony values for WLCC
- Lead-lag relationships for WLCC
- Aggregated phase and coherence values in specified frequency bands for WTC
- Binned phase values in specified frequency bands for WTC
- Summary statistics across the interaction
- Rose plots showing phase distribution for one dyad

For instance, WLCC of smiling may indicate how much participants feel or communicate Joy with each other.  For head gestures, WLCC can inform the best lag with the highest adaptation values, while WTC can show difference in rhythmic adaptation of nodding. 

### Pseudosynchrony Comparison

Observed synchrony can arise either because interaction partners are genuinely coordinating with one another or simply by chance.

To help distinguish between these possibilities, all synchrony measures can be compared with **pseudo IPS** values. Pseudo IPS distributions are created in two ways: 

- Shuffling interaction partners to create pseudo dyads of people who did not interact with each other
- Shuffling segments of the time course data of one interaction partner

This allows researchers to determine whether observed synchrony exceeds chance levels.

Outputs include:

- Pseudosynchrony distributions
- Statistical comparisons (frequentist, Bayesian or based on permutation tests)
- Visual comparison revealing the best lag to distinguish observed from pseudo IPS

## Who Is This Package For?

`interactR` is designed for researchers in:

- Psychology
- Communication science
- Behavioural science
- Human-computer interaction
- Clinical research
- Social interaction research

No advanced signal-processing or programming expertise is required.

## Getting Started

In the `vignettes` folder, you will soon find examples for the extraction of different features, paired with videos allowing you to compare what was extracted with the behaviour of the interaction partners. 

### Data Structure

Each data collection session is uniquely identified by the combination of two variables: 

- `Time`
- `Dyad`

In the case of dyadic interactions, `Dyad` consists of the two interaction partners separated by a hyphen, e.g., `sub01-sub02`. Preprocessing dataframes contain the data of the interaction partners in long format, specified by `Identifier`. 

Solo data also has a `Dyad` column, but only consisting of the `Identifier` of the one person combined with the suffix `-solo`. Upcoming versions will also offer options for interactions between more than two interaction partners. 

### Created Files

Throughout the preprocessing, large data is saved as `arrow` files while aggregated data (one row per interaction partner) is saved as `csv` files. Columns in the `csv` file start with the feature group name, except for dyadic variables which are the same for both interaction partners of the same dyad and have the additional prefix `Dyad`. 


## Installation

```r
# install.packages("pak")

pak::pak("IreneSophia/interactR")
```

## Citation

If you use `interactR` in your research, please cite:

```text
Citation information will be added upon publication.
```

