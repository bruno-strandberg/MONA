#!/usr/bin/env python
import numpy as np

def count_print(filename):
    with open(filename) as f:
        q = []
        total = []
        goodev = []
        goodtr = []
        goodsh = []
        badtr = []
        badsh = []
        badboth = []
        rejected = []
        shower = False
        for line in f:
            words = line.split(":")
            if "Quality" in line:
                q.append(float(words[1]))
            elif "Total events" in line:
                total.append(float(words[1]))
            elif "Good events" in line:
                goodev.append(float(words[1]))
            elif "Good tracks" in line:
                goodtr.append(float(words[1]))
            elif "Good showers" in line:
                goodsh.append(float(words[1]))
            elif "Bad tracks" in line:
                badtr.append(float(words[1]))
            elif "Bad showers" in line:
                badsh.append(float(words[1]))
            elif "Bad both" in line:
                badboth.append(float(words[1]))
            elif "Rejected" in line:
                rejected.append(float(words[1]))
            else:
                pass
    return [q, total, goodev, goodtr, goodsh, badtr, badsh, badboth, rejected]

def printlist(lst):
    for item in lst:
        #if item is not None:
        print item

numbers = count_print("./reconstruction_counts.txt")
q = numbers[0]
total = numbers[1]
goodev = numbers[2]
goodtr = numbers[3]
goodsh = numbers[4]
badtr = numbers[5]
badsh = numbers[6]
badboth = numbers[7]
rejected = numbers[8]

track_numbers = [l[:10] for l in numbers]
shower_numbers = [l[10:] for l in numbers]
#numbers = np.array(numbers)
print "Track numbers"
print "Total events"
print printlist(track_numbers[1])
print "Good events"
print printlist(track_numbers[2])
print "Good tracks"
print printlist(track_numbers[3])
print "Good showers"
print printlist(track_numbers[4])
print "Bad tracks"
print printlist(track_numbers[5])
print "Bad showers"
print printlist(track_numbers[6])
print "Bad both"
print printlist(track_numbers[7])
print "Rejected"
print printlist(track_numbers[8])

print "Shower numbers"
print "Total events"
print printlist(shower_numbers[1])
print "Good events"
print printlist(shower_numbers[2])
print "Good tracks"
print printlist(shower_numbers[3])
print "Good showers"
print printlist(shower_numbers[4])
print "Bad tracks"
print printlist(shower_numbers[5])
print "Bad showers"
print printlist(shower_numbers[6])
print "Bad both"
print printlist(shower_numbers[7])
print "Rejected"
print printlist(shower_numbers[8])
