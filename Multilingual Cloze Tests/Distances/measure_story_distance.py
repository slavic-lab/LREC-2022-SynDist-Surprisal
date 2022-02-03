#MEASURE SYNTACTIC DISTANCES BETWEEN SLAVIC LANGUAGES FROM MULTIPLE SYNTACTIC ALIGNMENT
#Written by Philip Georgis (2021-22) in cooperation with the SFB 1102 C4 Research Project

#Load required modules
import os
from pathlib import Path
from itertools import product
from statistics import mean
from collections import defaultdict

#Designate local and parent file paths
local_dir = Path(str(os.getcwd()))
parent_dir = local_dir.parent
grandparent_dir = parent_dir.parent

#Load syntactic distance tools
os.chdir(str(grandparent_dir) + "/Linguistic Distance/Syntactic_Distance")
from syntactic_distance import *
os.chdir(local_dir)

#%%
#Load North Wind and Sun alignment workbook
northwindsun = load_alignfile('Multilingual_Story_Alignments.xlsx')

#Calculate start lines of the fragment alignments
startlines = calculate_startlines(workbook=northwindsun,
                                  sheet='Multiple',
                                  index_row=1,
                                  n_alignments=9,
                                  n_spaces=1)    

#%%
#Calculate multiple syntactic distance for each aligned fragment by iterating through start lines of each
#Create nested dictionary of fragment numbers with multiple syntactic distances dictionaries as values
fragment_syndists = {}
for i in range(len(startlines)):
    fragment = i+1
    start_row = startlines[i]+1
    fragment_syndists[fragment] = multiple_syntactic_distance(workbook=northwindsun,
                                                              sheet='Multiple',
                                                              start_row=start_row,
                                                              start_column='B',
                                                              n_alignments=9,
                                                              
                                                              #Also specify a workbook/sheet with POS tags
                                                              #in order to account for the overlap in POS
                                                              POS_workbook=northwindsun,
                                                              POS_sheet='POS')

#Get list of language pairs
language_pairs = list(product(fragment_syndists[1].keys(), fragment_syndists[1].keys()))
language_pairs = [pair for pair in language_pairs if pair[0] != pair[1]]

#Designate abbreviations used for each language
lang_abbrev = {'RU':'Russian',
               'BE':'Belarusian',
               'UK':'Ukrainian',
               'PL':'Polish',
               'CS':'Czech',
               'SK':'Slovak',
               'SL':'Slovene',
               'HR':'Croatian',
               'BG':'Bulgarian'}

#%%
#Write CSV file with measured distances
#Simultaneously create a new dictionary to store the distances in a clearer format
distances = defaultdict(lambda:defaultdict(lambda:[]))


with open('Syntactic/NorthWindandSun_syntactic_distances.csv', 'w') as f:
    
    #Create all column names for various measures: InDel, Movement, Binary Movement
    #   - non-normalized and normalized version of each
    #   - measure for every fragment individually and the average over the whole text
    measure_labels = ['InDel', 'nInDel', 
                      'Movement', 'nMovement', 
                      'BinaryMovement', 'nBinaryMovement']
    fragment_labels = [f'Fragment{i}_{measure}' 
                       for i in range(1,8) 
                       for measure in measure_labels]
    mean_labels = [f'Mean_{measure}' for measure in measure_labels]
    columns = ['Language1', 'Language2'] + fragment_labels + mean_labels
    f.write(",".join(columns))
    f.write('\n')
    
    #Iterate through language pairs; extract and write distances for each measure
    for lang_pair in language_pairs:
        lang1, lang2 = lang_pair
        
        #Initialize empty lists containing the distances
        InDels, nInDels = [], []
        Movements, nMovements = [], []
        BinaryMovements, nBinaryMovements = [], []
        
        #Iterate through fragments of the text, extract non-normalized measures
        for fragment in fragment_syndists:
            align_length, indel, movement, moved_words = fragment_syndists[fragment][lang1][lang2]
            
            #Normalize measures by length of alignment
            nInDel = indel / align_length
            nMovement = movement / align_length
            nBinaryMovement = moved_words / align_length
            
            #Add measures to respective lists and to distance dictionary under associated label
            for measure, lis, label in zip([indel, nInDel, movement, nMovement, moved_words, nBinaryMovement], 
                                           [InDels, nInDels, Movements, nMovements, BinaryMovements, nBinaryMovements],
                                           measure_labels):
                lis.append(measure)
                distances[lang_pair][label].append(measure)
        
        #Write the measures for this language pair to a single row in the CSV
        lang1, lang2 = lang_abbrev[lang1], lang_abbrev[lang2]
        f.write(f'{lang1},{lang2},')
        
        #Write measures for each fragment
        fragment_measures = list(zip(InDels, nInDels, 
                                     Movements, nMovements, 
                                     BinaryMovements, nBinaryMovements))
        for i in fragment_measures:
            values = ",".join([str(v) for v in i])            
            f.write(f'{values},')
        
        #Add mean value for each measure
        measures =  [InDels, nInDels, 
                     Movements, nMovements, 
                     BinaryMovements, nBinaryMovements]
        for measure, label in zip(measures, measure_labels):
            mean_measure = mean(measure)
            f.write(f'{mean_measure},')            
        f.write('\n')