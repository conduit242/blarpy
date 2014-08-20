#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#
# Copyright 2014 Robert Bird
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# blar.py
# Version 0.1
# A simple program for converting log files to a streaming UTF genome and analyzing it.
# Syntax is: python blar.py -i <file_name> or python blar.py -h for help & options

import sys
import xxhash
import numpy
import fileinput
import locale
import codecs
import unicodedata
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import operator
import time
import argparse


# Is odd...roughly twice as fast as mod test

def is_odd(num):
    return num & 0x1


# Moving average

def moving_average(values, window):
    weights = numpy.repeat(1.0, window) / window
    sma = numpy.convolve(values, weights, 'valid')
    return sma


# Calculate vector magnitude

def magnitude(V):
    return numpy.sqrt(sum([x * x for x in numpy.nditer(V)]))


# Vector normalizer

def normalize(V):
    v_m = magnitude(V)
    return [numpy.divide(vi, v_m) for vi in numpy.nditer(V)]


# Feature hash input. Has the effect of quantizing inputs with graceful random collisions

def feature_hash_string(s, window, dim):

    start = time.clock()

    # Generate window-char Markov chains & create feature hash vector

    v = {}
    for x in range(0, dim):
        v[x] = 0
    length = len(s)
    max_num = 2.0 ** 64
    for x in range(0, length - window):
        key = xxhash.xxh64(s[x:x + window]) % dim
        v[key] += 0x1

    return numpy.asarray(v.values())


# Use random projection for LSH and output a UTF char for the hash

def locality_hash_vector(v, width):

    start = time.clock()
    hash = numpy.zeros(width, dtype=int)
    for x in range(0, width):
        projection = numpy.dot(PROJECTION_VECTORS[x], v)
        if projection < 0:
            hash[x] = 0
        else:
            hash[x] = 0x1

    # Return unicode char equal to the LSH

    return unichr(int(''.join(map(str, hash)), 2))


# MAIN
#

# Handle encodings for STDOUT vs other
# if sys.stdout.isatty():
#    default_encoding = sys.stdout.encoding
# else:
#    default_encoding = locale.getpreferredencoding()

# Process command line arguments

parser = argparse.ArgumentParser()

parser.add_argument(
    '-a',
    action='store',
    dest='a_width',
    default=0,
    type=int,
    help='Set the nucleotide alphabet width in bits',
    )

parser.add_argument(
    '-f',
    action='store',
    dest='f_width',
    default=0,
    type=int,
    help='Set the feature hash vector width in array length',
    )

parser.add_argument(
    '-v',
    action='store',
    dest='v_width',
    default=0,
    type=int,
    help='Set the n-gram vectorization length for line processing',
    )

parser.add_argument(
    '-c',
    action='store',
    dest='c_width',
    default=0,
    type=int,
    help='Set the width of codons in nucleotides',
    )

parser.add_argument(
    '-C',
    action='store_true',
    dest='C',
    default=False,
    help='Display charts; default off'
    )

parser.add_argument(
    '-m',
    action='store',
    dest='ma_window',
    default=100,
    type=int,
    help='Set the window size for moving average display',
    )

parser.add_argument(
    '-g',
    action='store',
    dest='granularity',
    default=2, 
    type=int,
    help='Sets the granularity factor for auto-tuning, from 1-N'
    )

parser.add_argument(
    '-t',
    action='store',
    dest='i_threshold',
    default=0.,
    type=float,
    help='Sets the threshold for printing interesting items in std. deviations, default 2.0'
    )

parser.add_argument(
    '-p',
    action='store_false',
    dest='hard',
    required=False,
    default=True,
    help='Do not use hard projection vectors',
    )

parser.add_argument(
    '-i',
    action='store',
    dest='f_name',
    required=True,
    help='Input file name'
    )

args = vars(parser.parse_args())

# File handling

f = open(args['f_name'], 'r')

# Read in the file once and build a list of line offsets for display use later

line_offset = []
offset = 0
file_length = 0
for line in f:
    line_offset.append(offset)
    offset += len(line)
    file_length += 0x1
print '# of lines: ', file_length
f.seek(0)

# Assign argument values; autotune feature_width & alphabet_width
# based on file length re: Shannon, et al.

granularity = args['granularity']
ma_window = args['ma_window']

if args['c_width'] > 0:
    codon_width = args['c_width']
else:
    codon_width = 4

if args['v_width'] > 0:
    vect_width = args['v_width']
else:
    vect_width = 4

if args['f_width'] > 0:
    feature_width = args['f_width']
else:
    print 'Autotuning feature hash vector width'
    feature_width = granularity * int(math.ceil(math.log(file_length,
            2)))

if args['a_width'] > 0:
    alphabet_width = args['a_width']
else:
    print 'Autotuning alphabet width'
    alphabet_width = int(math.ceil(math.sqrt(math.log(file_length, 2))))

if args['i_threshold'] != 0.:
    interestingness_threshold = args['i_threshold']
else:
    interestingness_threshold = 2.0

print 'Feature width: ', feature_width
print 'Alphabet width: ', alphabet_width
print 'Codon width: ', codon_width
print 'Granularity: ', granularity

# Generate random unit normal comparison vectors for random projection

PROJECTION_VECTORS = []
if args['hard']:
    for vector in range(0, alphabet_width):
        numpy.random.seed(vector)
        v = numpy.random.randint(2, size=feature_width)
        counter = 0
        for e in v:
            if e == 0:
                v[counter] = -0x1
                counter += 0x1
            else:
                v[counter] = 0x1
                counter += 0x1
        PROJECTION_VECTORS.append(v)
else:
    for vector in range(0, alphabet_width):
        numpy.random.seed(vector)
        PROJECTION_VECTORS.append(normalize(numpy.random.randn(0x1,
                                  feature_width)))

# CREATE GENOME
#
# Note irritating list method to guarantee unicode handling beyond
# UTF8 rather than concatenating a string which is nicer for STDOUT.
# To say I'm tempted to cap it at 8 is an understatement, it doesn't
# seem necessary for an effective tool

start = time.clock()
genome = []
for line in f:
    genome.append(locality_hash_vector(feature_hash_string(line,
                  vect_width, feature_width), alphabet_width))

genome_len = len(genome)
print 'File genome assembly time: ', round(time.clock() - start), 's'
print 'Genome length: ', genome_len

# SCORING
#

# Find codon counts over genome; create stream of windows for scoring display

counts = {}
codon_stream = []
for x in range(0, len(genome) - codon_width + 0x1):
    window = []
    for shift in range(0, codon_width):
        window.append(genome[x + shift])
    key = ''.join(window)
    if key in counts:
        counts[key] += 0x1
    else:
        counts[key] = 0x1
    codon_stream.append(key)

# Find coding cost for each codon and create scoring map.

scoring_map = {}
for key in counts.keys():
    scoring_map[key] = -math.log(float(counts[key]) / (len(genome)
                                 - codon_width + 0x1), 2)

# Normalize scores to 0..1
#
# for key in scoring_map:
#     scoring_map[key] = (scoring_map[key] - min) / max

# Floor scores; use if not normalizing to reduce visual noise

min = numpy.min(scoring_map.values())
for key in scoring_map:
    scoring_map[key] = scoring_map[key] - min

# Create stream of scores; Apply exponential decay to scores as symbols re-occur

codon_score_stream = []
decay = {}
score = 0

for e in codon_stream:
    if e in decay:
        decay[e] += 0x1
    else:
        decay[e] = 0x1
    if decay[e] < 16:
        score = 1.0 / 2 ** (decay[e] - 0x1) * scoring_map[e]
    else:
        score = 1.0 / 2 ** 16 * scoring_map[e]
    codon_score_stream.append(score)

max_anomaly = numpy.max(codon_score_stream)
print 'Max anomaly score: ', max_anomaly, 'bits'

# Find some robust statistics to pick out interesting items for printing
# Too expensive for files of any significant length
#
# start = time.clock()
# median_score = numpy.median(codon_score_stream)
# cr_s_estimator = 0
# med_dists = []
# for e in codon_score_stream:
#    e_dist = []
#    for all in codon_score_stream:
#        e_dist.append(e - all)
#    # print e_dist
#    med_dists.append(numpy.median(e_dist))
# cr_s_estimator = numpy.median(med_dists)
# print 'Median: ', median_score
# print 'CR Estimator: ', cr_s_estimator
# print 'Calc time: ', time.clock() - start

# Find some Gaussian statistics

average_score = numpy.average(codon_score_stream)
std_dev_score = numpy.std(codon_score_stream)
print 'Average score: ', average_score
print 'Std. deviation: ', std_dev_score

# PRINTING
#
# Should be replaced with some kind of auto scrolling line explorer
# where you highlight a region in the chart and it displays those lines + n extra
# from the file. It's more useful to view the lines in a viewer for scanning
#

# Scan over file and print out interesting lines. Does not use linecache to avoid memory caching file

tracker = 0
interesting = 0
for e in codon_score_stream:
    if e > average_score + interestingness_threshold * std_dev_score:
        f.seek(line_offset[tracker])
        print tracker, ',', f.readline().rstrip(), ',', e
        interesting += 0x1
    tracker += 0x1
    f.seek(0)
print '# of interesting things: ', interesting
f.close()

# CHARTING
#

if args['C'] == True:

    # Create an LSH bucket counter for charting

    buckets = {}
    for e in genome:
        if e in buckets:
            buckets[e] += 0x1
        else:
            buckets[e] = 0x1

    # Create moving average list for charting

    streamMA = moving_average(codon_score_stream, ma_window)

    # Chart stuff

    fig = plt.figure()

    score_plt = fig.add_subplot(211)
    score_plt.plot(streamMA)
    score_plt.set_ylabel('Anomaly Score')
    score_plt.set_ylim(0., numpy.max(streamMA) + 0x1)
    score_plt.set_xlabel('Line number')

    bucket_plt = fig.add_subplot(212)
    bucket_plt.bar(range(len(buckets.keys())), [math.log(float(y), 10)
                   for y in buckets.values()], 0x1)
    bucket_plt.set_ylim(0., numpy.max([math.log(float(y), 10) for y in
                        buckets.values()]) + 0x1)
    bucket_plt.set_xlim(0., 2 ** alphabet_width)
    bucket_plt.set_xlabel('LSH Bin number')
    bucket_plt.set_ylabel('Counts (log10)')

    plt.show()

