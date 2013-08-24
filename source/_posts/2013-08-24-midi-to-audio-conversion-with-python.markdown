---
layout: post
title: "MIDI To Audio Conversion With Python"
date: 2013-08-24 14:54
comments: true
categories: MIDI, Audio, MIR, Python, Machine Learning, Data
---

Last year I was working on a Music Information Retrieval (MIR) problem, trying to identify a song's key from audio recordings.  With many machine learning problems, one of the most difficult tasks is finding good data sets to train and test with.  While there are some good [annotated music datasets](http://isophonics.net/datasets) out there, sometimes you need to generate your own test data.  This post includes a simple Python script we used to generate test data based on MIDI files.

<!-- more -->

For those interested in the actual song key identification research, the source code for that is available [here](http://github.com/devonbryant/uccs-mir-key).  The algorithm we used was based on Hidden Markov model classifiers and Beat-Synchronous Chromagram audio features.

## Why MIDI?

Extracting annotations from a MIDI file is often much simpler than other formats because the timing information, instruments, and notes are all encoded in the file itself.  Unfortunately the default MIDI soundbanks on most computers are pretty terrible and don't sound like real instruments or musicians, so you don't want to use them for evaluating or training MIR systems.  However, there are many realistic soundfonts available online that can make a MIDI recording sound close to real instruments.  To avoid overfitting our system, we wanted to use real audio training sets with diverse instrument sounds (rock, jazz, classical, metal, etc.) so that the key identification was mostly genre-agnostic.

## Conversion Script

The following Python script takes a directory of MIDI files and a directory of SF2 soundfont files and generates corresponding audio files (wav, aiff, mp3, etc.).  The selection of which soundfont is used for each MIDI file is random.  To use this script, you need to have [FluidSynth](http://sourceforge.net/apps/trac/fluidsynth/) installed.

{% gist 1810984 %}