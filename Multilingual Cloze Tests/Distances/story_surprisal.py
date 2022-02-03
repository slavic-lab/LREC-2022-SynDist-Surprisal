#MEASURE ORTHOGRAPHIC AND PHONETIC SURPRISAL BETWEEN SLAVIC LANGUAGES
#BASED ON COGNATE ALIGNMENTS FROM MULTIPLE SYNTACTIC ALIGNMENT OF PARALLEL SENTENCES
#Written by Philip Georgis (2021-22) in cooperation with the SFB 1102 C4 Research Project

#Load required modules
import os
from pathlib import Path
from collections import defaultdict
from nltk import edit_distance
import pandas as pd
from statistics import mean
import seaborn as sns
sns.set(font_scale=1.0)

local_dir = Path(str(os.getcwd()))
parent_dir = local_dir.parent
grandparent_dir = parent_dir.parent

#Load syntactic distance tools
os.chdir(str(grandparent_dir) + "/Linguistic Distance/Syntactic_Distance")
from syntactic_distance import *
os.chdir(str(grandparent_dir))
from auxiliary_functions import normalize_dict, lidstone_smoothing, surprisal, draw_dendrogram, plot_distances
os.chdir(local_dir)

#%%
#Import phonetic distance tools from thesis directory
os.chdir('/Users/phgeorgis/Documents/School/MSc/Saarland_University/Courses/Thesis/Code')
from phonetic_distance import *
os.chdir(local_dir)

#%%
#ORTHOGRAPHIC ALIGNMENT

#Import orthographic equivalencies
from orthographic_equivalencies import *

#Load Slavic digraphs from CSV file
slavic_digraphs_csv = pd.read_csv('Slavic_digraphs.csv', sep='\t')
slavic_digraphs = {}
for i, row in slavic_digraphs_csv.iterrows():
    slavic_digraphs[row.Language] = row.Digraphs.split(', ')


def segment_orth(word, digraphs=[]):
    """
    Segments an orthographic string using optional language-specific digraph 
    conventions into a list of orthographic segments/units.

    Parameters
    ----------
    word : string
        Orthographic string, e.g. "cześć".
    digraphs : list, optional
        List of character sequences to be considered as digraphs/trigraphs, 
        e.g. ['ch', 'cz', 'dz', 'dź', 'dż', 'rz', 'sz'].
        By default there are no digraphs. 

    Returns
    -------
    segs : list
        Segmented list of orthraphic characters/units, e.g. ['cz', 'e', 'ś', 'ć'].

    """
    #If no digraphs are specified, simply segment each character into its own unit
    if len(digraphs) == 0:
        segs = list(word)
        return segs
    
    #Otherwise segment considering digraphs
    else:
        segs = []
        i = 0
        
        #Iterate through character indices of string
        #Check at each step whether the current plus following character(s) form a di/trigraph
        #Add the identified segmentable unit to list of segments
        while i < len(word):
            ch = word[i]
            try:
                nxt = word[i+1]
                if f'{ch.lower()}{nxt.lower()}' in digraphs:
                    
                    #Check if it forms part of a trigraph
                    try:
                        nxtnxt = word[i+2]
                        if f'{ch.lower()}{nxt.lower()}{nxtnxt.lower()}' in digraphs:
                            segs.append(f'{ch}{nxt}{nxtnxt}')
                            i += 3
                        else:
                            segs.append(f'{ch}{nxt}')
                            i += 2
                    except IndexError:
                        segs.append(f'{ch}{nxt}')
                        i += 2
                else:
                    segs.append(ch)
                    i += 1
                    
            except IndexError:
                segs.append(ch)
                i += 1
        
        return segs
        


def orth_align(word1, word2, 
               gop=-2, dist_func=orth_dist,
               lang1_digraphs=None, lang2_digraphs=None, 
               **kwargs):
    """
    Aligns two unsegmented orthographic strings.

    Parameters
    ----------
    word1 : string
        Orthographic string of first word to be aligned.
    word2 : string
        Orthographic string of second word to be aligned.
    gop : float, optional
        Value for the gap-opening penalty (GOP). The default is -2.
    dist_func : function, optional
        Function to use for pairwise alignment costs. The default is the orth_dist() function.
    lang1_digraphs : list, optional
        List of character sequences to be considered as digraphs/trigraphs
        for segmenting word1. By default there are no digraphs. 
    lang2_digraphs : list, optional
        List of character sequences to be considered as digraphs/trigraphs
        for segmenting word2. By default there are no digraphs. 
    **kwargs : dictionary, optional
        Dictionary of named keyword arguments to be passed to align_costs() 
        and the specified dist_func.

    Returns
    -------
    Best alignment using Needleman-Wunsch algorithm and pairwise alignment costs
    using specified orthographic distance function : list

    """
    
    #Segment using specific digraphs, if necessary (otherwise segmentation is unnecessary)
    if lang1_digraphs != None:
        word1 = segment_orth(word1, digraphs=lang1_digraphs)
    if lang2_digraphs != None:
        word2 = segment_orth(word2, digraphs=lang2_digraphs)

    #Calculate alignment cost for each pair of orthographic segments
    costs = align_costs(word1, word2, dist_func=dist_func, **kwargs)
    
    #Return least costly alignment according to Needleman-Wunsch alignment algorithm
    return best_alignment(SEQUENCE_1=list(word1), SEQUENCE_2=list(word2), 
                          SCORES_DICT=costs, GAP_SCORE=gop)


#%%
#EXTRACT CORRESPONDING WORDS/SYNTACTIC UNITS FROM SYNTACTIC ALIGNMENT

def find_matching_words(workbook, sheet, row, search_row, 
                        start_column, end_column=None, row_columns=None):
    """
    Identifies corresponding words/syntactic units in alignment sheet. 
    Units are considered to correspond if they are aligned in the same column
    or if their cells are highlighed in the same color.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook.
    sheet : string
        Label of sheet in workbook from which to extract matches.
    row : int
        Row number in worksheet to use as base for matches.
    search_row : int
        Row number in worksheet where matches with base row should be searched for.
    start_column : string
        Column string, e.g. 'B', specifying the first column of the alignment to search.
    end_column : string, optional
        Column string, e.g. 'Z', specifying the final column of the search range.
        By default all non-empty columns are searched.
    row_columns : list, optional
        Columns within row for which to search for matches, by default all columns.

    Returns
    -------
    matched : list
        List of tuples in format ('word1', 'word2') of matched units.

    """
    
    #Load workbook sheet
    sh = workbook[sheet]
    
    #Create list of all worksheet columns to search
    columns = [i[0].column_letter for i in list(sh.columns)]
    
    #Columns within row for which to search for matches, by default all columns
    if row_columns == None:
        row_columns = columns
    
    #By default, search all columns to the right of the start_column
    if end_column == None:
        columns = columns[columns.index(start_column):]
    
    #Otherwise search only the specified range of columns
    else:
        columns = columns[columns.index(start_column):columns.index(end_column)+1]
    
    #Create index of matched word pairs
    matched = []
    
    #Iterate through columns
    for c in row_columns:
        
        #Check if there is a word in this column in the first line
        entry1 = sh[f'{c}{row}'].value
        if entry1 != None:
            
            #Get the color of this cell
            color1 = get_color(workbook, sheet, f'{c}{row}')
            
            #Check whether the color is blank
            #If blank, simply check if there is an aligned word in the same column
            if color1 == blank:
                entry2 = sh[f'{c}{search_row}'].value
                
                #If the search row has an entry, add to index of matched words
                if entry2 != None:
                    matched.append((entry1.strip(), entry2.strip()))
                    
            
            #If not blank, find all words from search row with same color
            #And add each one to the matched word index
            else:
                matched_color = [sh[f'{cx}{search_row}'].value for cx in columns
                                 if get_color(workbook, sheet, f'{cx}{search_row}') == color1]
                for entry2 in matched_color:
                    matched.append((entry1.strip(), entry2.strip()))
    
    return matched

#Load North Wind and Sun alignment workbook
northwindsun = load_alignfile('Multilingual_Story_Alignments.xlsx')

#Calculate startlines of the fragment alignments
startlines = calculate_startlines(workbook=northwindsun,
                                  sheet='Multiple',
                                  index_row=1,
                                  n_alignments=9,
                                  n_spaces=1)

#Get list of all non-empty columns in the alignment sheet
columns = [i[0].column_letter for i in list(northwindsun['Multiple'].columns)]

#Create dictionary of orthographic and phonetic words per language
orth_words = defaultdict(lambda:defaultdict(lambda:0))
phonetic_words = defaultdict(lambda:defaultdict(lambda:0))

#Create dictionary of matched words per language pair
matched_words = defaultdict(lambda:[])

#Iterate through start lines
for line in startlines:
    #Get list of content rows within this fragment (given 9 languages)
    rows = list(range(line+1, line+10))
    
    #Iterate through unique combinations of rows i and j
    for i in range(len(rows)):
        row1 = rows[i]
        l1 = northwindsun['Multiple'][f'A{row1}'].value
        
        #Get all orthographic words in row
        l1_orth_words = [northwindsun['Multiple'][f'{col}{row1}'].value.lower()
                         for col in columns[1:] 
                         if northwindsun['Multiple'][f'{col}{row1}'].value != None]
        
        #Get all phonetic words in row
        l1_phon_words = [northwindsun['IPA'][f'{col}{row1}'].value.lower()
                         for col in columns[1:] 
                         if northwindsun['IPA'][f'{col}{row1}'].value != None]
        
        #Update counts of instances of each phonetic/orthographic word
        for words, count_dict in zip([l1_orth_words, l1_phon_words], [orth_words, phonetic_words]):
            for word in words:
                count_dict[l1][word] += 1
        
        #Iterate through second row for combination, j
        for j in range(i+1, len(rows)):
            row2 = rows[j]
            l2 = northwindsun['Multiple'][f'A{row2}'].value
            
            #Get the orthographic matches between row i and j
            orth_matches = find_matching_words(northwindsun, sheet='Multiple',
                                               row=row1, search_row=row2, 
                                               start_column='B')
            
            #Get the phonetic matches between row i and j
            ipa_matches = find_matching_words(northwindsun, sheet='IPA', 
                                              row=row1, search_row=row2, 
                                              start_column='B')
            
            #Save matches to dictionary of matched words/units
            matched_words[(l1, l2)].extend(list(zip(orth_matches, ipa_matches)))


def adjust_ipa(tr):
    """
    Removes spaces, extraneous unicode characters, and unneeded IPA diacritics
    for stress, pitch accents, and syllabic consonants from IPA strings.

    Parameters
    ----------
    tr : string
        IPA string.

    Returns
    -------
    adjusted : string
        IPA string with extraneous characters stripped away.

    """
    adjusted = strip_ch(tr, ["ˈ", "̂", "̌", "̩", " ", '\ufeff'])
    return adjusted


#Write a file containing all the matched word/unit pairs
with open('northwindsun_matched_words.csv', 'w') as f:
    f.write(','.join(['L1', 'Word1', 'IPA1', 'L2', 'Word2', 'IPA2', 'LD']))
    f.write('\n')
    
    for lang_pair in matched_words:
        l1, l2 = lang_pair
        
        for entry in matched_words[lang_pair]:
            orth, tr = entry
            word1, word2 = orth
            ipa1, ipa2 = tr
            
            #Lowercase the orthographic words
            word1 = word1.lower()
            word2 = word2.lower()
                        
            #Remove stress, tone, syllabic marking in IPA, also spaces
            ipa1, ipa2 = adjust_ipa(ipa1), adjust_ipa(ipa2)
            
            #Measure Levenshtein distance of IPA transcriptions
            LD = edit_distance(ipa1, ipa2)
            LD /= max(len(ipa1), len(ipa2))
            
            f.write(','.join([l1, word1, ipa1, l2, word2, ipa2, str(LD)]))
            f.write('\n')

#%%
#EXTRACT WORD TYPE COUNTS, TOKEN COUNTS, AND ORTHOGRAPHIC/ORTHOGRAPHIC CHARACTER FREQUENCIES

#Get statistics about number of unique tokens and word/character frequencies
text_stats = defaultdict(lambda:{})
ch_freqs = defaultdict(lambda:defaultdict(lambda:0))
phon_freqs = defaultdict(lambda:defaultdict(lambda:0))
for lang in orth_words:
    #Get numbers of unique words/units, and total numbers of word/unit tokens
    n_units = len(orth_words[lang])
    n_tokens = sum(orth_words[lang].values())
    n_words = sum([len(unit.split()) for unit in orth_words[lang]])
    n_word_tokens = sum([len(unit.split())*orth_words[lang][unit] for unit in orth_words[lang]])
    text_stats[lang]['N_unique_units'] = n_units
    text_stats[lang]['N_unit_tokens'] = n_tokens
    text_stats[lang]['N_unique_words'] = n_words
    text_stats[lang]['N_word_tokens'] = n_word_tokens
    
    #Get frequency of orthographic characters in each language
    lang_dg = slavic_digraphs[lang]
    for unit in orth_words[lang]:
        for word in unit.split():
            segments = segment_orth(word, digraphs=lang_dg)
            for ch in segments:
                ch_freqs[lang][ch] += 1
                
    #Get frequency of phonetic characters in each language
    for unit in phonetic_words[lang]:
        for word in unit.split():
            segments = segment_word(adjust_ipa(word))
            for seg in segments:
                phon_freqs[lang][seg] += 1

#Normalize character and phonetic segment frequencies
#for lang in orth_words:
#    ch_freqs[lang] = normalize_dict(ch_freqs[lang])
#    phon_freqs[lang] = normalize_dict(phon_freqs[lang])

#Write files with these statistics
for counts, file in zip([ch_freqs, phon_freqs], ['character_frequencies.csv', 
                                                 'phon_segment_frequences.csv']):
    header = ',' + ','.join(list(counts.keys()))    
    all_chs = sorted(list(set(ch for lang in counts for ch in counts[lang])))
    with open(file, 'w') as f:
        f.write(f'{header}\n')
        for ch in all_chs:
            vals = [counts[lang].get(ch, 0) for lang in counts]
            vals = ','.join([str(i) for i in vals])
            f.write(f'{ch},{vals}\n')

with open('text_statistics.csv', 'w') as f:
    cats = list(text_stats['RU'].keys())
    header = 'Language,' + ','.join(cats)
    f.write(f'{header}\n')
    for lang in text_stats:
        vals = ','.join([str(text_stats[lang][cat]) for cat in cats])
        f.write(f'{lang},{vals}\n')
            
#%%
#ALIGN ORTHOGRAPHIC AND PHONETIC SEQUENCES 
#AND EXTRACT CORRESPONDENCE PROBABILITIES FOR SURPRISAL CALCULATIONS
#NOTE: PRE-ANNOTATED FILE WITH COGNACY JUDGMENTS NEEDED

#Initialize dictionary containing phonetic correspondence counts
#e.g. phonetic_correspondences['RU']['UK']['l']['ʋ'] = 4
phonetic_correspondences = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:0))))
orthographic_correspondences = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:0))))
only_cognates = True #whether to use only cognates as basis for surprisal

try:
    #Load the file containing matched word/unit pairs with cognacy annotations
    all_aligned_matches = pd.read_csv('northwindsun_matched_words_annotated.csv')
    
    #Filter only cognate word pairs
    matched_cognates = all_aligned_matches[all_aligned_matches.Cognates == 'x']
    if only_cognates == True:
        matched_annotated = matched_cognates
    else:
        matched_annotated = all_aligned_matches
        
    #Iterate through rows of dataframe
    for index, row in matched_annotated.iterrows():
        l1, l2 = row.L1, row.L2
        
        #Extract the phonetic transcriptions, strip off extraneous characters
        ipa1, ipa2 = adjust_ipa(row.IPA1), adjust_ipa(row.IPA2)
        
        #Align the phonetic transcriptions and update counts of correspondences
        ph_alignment1 = phone_align(ipa1, ipa2)
        for pair in ph_alignment1:
            seg1, seg2 = pair
            phonetic_correspondences[l1][l2][seg1][seg2] += 1
        
        #Reverse the alignment to get the opposite direction and also update counts of correspondences
        ph_alignment2 = reverse_alignment(ph_alignment1)
        for pair in ph_alignment2:
            seg1, seg2 = pair
            phonetic_correspondences[l2][l1][seg1][seg2] += 1
        
        #Same process for orthographic forms
        #Extract orthograhic word forms
        orth1, orth2 = strip_ch(row.Word1, ['\ufeff']), strip_ch(row.Word2, ['\ufeff'])
        
        #Align orthographic word forms
        orth_alignment1 = orth_align(orth1, orth2, 
                                     lang1_digraphs=slavic_digraphs[l1], 
                                     lang2_digraphs=slavic_digraphs[l2])
        for pair in orth_alignment1:
            ch1, ch2 = pair
            orthographic_correspondences[l1][l2][ch1][ch2] += 1
        
        orth_alignment2 = reverse_alignment(orth_alignment1)
        for pair in orth_alignment2:
            ch1, ch2 = pair
            orthographic_correspondences[l2][l1][ch1][ch2] += 1
    
except FileNotFoundError:
    print('Error: CSV file with cognacy annotations not found!')


#Create dictionary of identified cognate pairs
cognate_dict = defaultdict(lambda:[])
for i, row in matched_cognates.iterrows():
    cognate_dict[(row.L1, row.L2)].append((row.Word1.strip(), row.Word2.strip()))
    cognate_dict[(row.L2, row.L1)].append((row.Word2.strip(), row.Word1.strip()))


#%%
#EXTRACT VALUES FOR d PARAMETER IN LIDSTONE SMOOTHING

#Extract number of attested phones in each language
d_values_phon = {lang1:len(set(seg for lang2 in phonetic_correspondences[lang1] 
                               for seg in phonetic_correspondences[lang1][lang2])) 
                 for lang1 in phonetic_correspondences}

#Extract number of attested characters/orthographic units in each language's orthography
d_values_orth = {lang1:len(set(ch for lang2 in orthographic_correspondences[lang1]
                               for ch in orthographic_correspondences[lang1][lang2]))
                 for lang1 in orthographic_correspondences}

#%%
#WRITE FILES WITH RAW COUNTS OF ORTHOGRAPHIC AND PHONETIC CORRESPONDENCES

for filename, d, d_values in zip(['orthographic_correspondences.csv', 'phonetic_correspondences.csv'],
                                 [orthographic_correspondences, phonetic_correspondences],
                                 [d_values_orth, d_values_phon]):
    with open(filename, 'w') as f:
        f.write(f'Language1,Character1,Language2,Character2,CorrCount,TotalCountCh1,TotalCountCh2,SmoothedSurprisal\n')
        for lang1 in d:
            for lang2 in d[lang1]:
                chs1 = sorted(list(d[lang1][lang2].keys()))
                for ch1 in chs1:
                    chs2 = sorted(list(d[lang1][lang2][ch1].keys()))
                    for ch2 in chs2:
                        corr_count = d[lang1][lang2][ch1][ch2]
                        if corr_count > 0:
                            total_ch1 = sum(d[lang1][lang2][ch1].values())
                            total_ch2 = sum(d[lang2][lang1][ch2].values())
                            smoothed_prob = lidstone_smoothing(x=corr_count, N=total_ch2,
                                                               d=d_values[lang1])
                            smoothed_surprisal = surprisal(smoothed_prob)
                            f.write(f'{lang1},{ch1},{lang2},{ch2},{corr_count},{total_ch1},{total_ch2},{smoothed_surprisal}\n')
    

#%%
#CALCULATE ORTHOGRAPHIC AND PHONETIC SURPRISAL BASED ON MATCHED/ALIGNED UNITS
#FROM SYNTACTIC ALIGNMENT SHEET AND EXTRACTED/SMOOTHED CORRESPONDENCE PROBABILITIES

#Initialize dictionaries with measures of phonetic and orthographic surprisal
#and counts of each type of aligned unit (cognate, non-cognate, gap)
phon_surprisal = defaultdict(lambda:[]) 
orth_surprisal = defaultdict(lambda:[])
cognate_type_counts = defaultdict(lambda:defaultdict(lambda:[]))

#Perform surprisal measurement for orthographic and phonetic units
#and in the same iteration loop also write a file with the results
for sheet, d_dict, corr_dict, align_func, surprisal_dict, label in zip(['Multiple', 'IPA'], 
                                                                       [d_values_orth, d_values_phon], 
                                                                       [orthographic_correspondences, phonetic_correspondences],
                                                                       [orth_align, phone_align],
                                                                       [orth_surprisal, phon_surprisal],
                                                                       ['Orthographic', 'Phonetic']):
    
    #Get list of columns in the alignment sheet
    northwindsun_columns = [i[0].column_letter for i in list(northwindsun[sheet].columns)]
    northwindsun_columns = northwindsun_columns[1:] #start from column B
    
    #Iterate through start lines
    for i in range(len(startlines)):
        print(f'Evaluating {label} fragment {i+1}...')
        line = startlines[i]
        
        #Get list of content rows within this fragment
        rows = list(range(line+1, line+10))
        
        #Iterate through all pairs of rows
        for row1 in rows:
            l1 = northwindsun[sheet][f'A{row1}'].value
            for row2 in rows:
                if row1 != row2:
                    l2 = northwindsun[sheet][f'A{row2}'].value
                    
                    alignments = []
                    
                    #Only count numbers of word alignment types once, either in orthographic or phonetic iteration
                    #Arbitrarily chosen to count these in orthographic iteration
                    if label == 'Orthographic':
                        wordtype_counts = defaultdict(lambda:[])
                    
                    #Iterate through columns
                    for c in northwindsun_columns:
                        
                        #Get the value of the cell in this row and column
                        word1 = northwindsun[sheet][f'{c}{row1}'].value
                        
                        #If not empty, find this cell's matches
                        if word1 != None:
                            word1 = word1.strip()
                            matches = find_matching_words(northwindsun, sheet, 
                                                          row=row1, search_row=row2,
                                                          start_column='B', end_column=None,
                                                          row_columns=[c])
                            
                            #If no matches are found, add empty string '' to list of matches
                            if len(matches) == 0:
                                matches.append((word1, ''))
                            
                            #Iterate through matches
                            for match in matches:
                                word1, word2 = match
                                
                                #Adjust IPA transcriptions and align phonetically
                                if sheet == 'IPA':
                                    word1, word2 = adjust_ipa(word1), adjust_ipa(word2).strip()
                                    alignment = align_func(word1, word2)
                                
                                #Lowercase and align orthographic words
                                else:
                                    word1, word2 = word1.lower().strip(), word2.lower().strip()
                                    alignment = align_func(word1, word2, 
                                                           lang1_digraphs=slavic_digraphs[l1], 
                                                           lang2_digraphs=slavic_digraphs[l2])
                                    
                                    #Check if the words are cognate, non-cognate, or a gap alignment
                                    #Save pair to appropriate list
                                    if word2 == '':
                                        wordtype_counts['gap'].append((word1, None))
                                    elif (word1, word2) in cognate_dict[(l1, l2)]:
                                        wordtype_counts['cognate'].append((word1, word2))
                                    else:
                                        wordtype_counts['non-cognate'].append((word1, word2))
                                
                                #Add to list of alignments
                                alignments.append(alignment)
                                
                                
                        #If there is no entry in the alignment slot for L1,
                        #check if there is an entry in this slot for L2
                        else:
                            word2 = northwindsun[sheet][f'{c}{row2}'].value
                            
                            #If there is a word in L2, align it with an empty
                            #string '' from L1, unless an equivalent in L1 can
                            #be found in a different columnn of the alignment
                            if word2 != None:
                                word2 = word2.strip()
                                
                                #Check whether there is no word in the L1
                                #position simply because it is elsewhere in the alignment
                                matches = find_matching_words(northwindsun, sheet, 
                                                          row=row2, search_row=row1,
                                                          start_column='B', end_column=None,
                                                          row_columns=[c])
                                
                                #If there is truly no L1 equivalent unit
                                if len(matches) == 0:
                                
                                    #Adjust IPA transcription 
                                    if sheet == 'IPA':
                                        word2 = adjust_ipa(word2)
                                    
                                    #Lowercase orthographic word
                                    else:
                                        word2 = word2.lower()
                                        
                                    #Align with gap and add to list of alignments
                                    alignment = align_func('', word2)
                                    alignments.append(alignment)
                                    
                                    #Add to list of gap alignments
                                    if label == 'Orthographic':
                                        wordtype_counts['gap'].append((word1, word2))
                                
                                #If matches are discovered elsewhere in the alignment,
                                #no further action is necessary; 
                                #this alignment is accounted for in previous block
                                else:
                                    pass
                            
                            #If there is no word in L2 either, skip this column
                            else:
                                pass
                    
                    #Iterate through list of alignments and score each
                    d_value = d_dict[l1]
                    alignment_scores = []
                    for alignment in alignments:
                        pair_scores = []
                        for pair in alignment:
                            seg1, seg2 = pair
                            
                            #Calculate Lidstone smoothed probability of correspondence
                            #from observed frequency of seg1-seg2 alignment,
                            #total occurrences of seg2 in L2, and 
                            #size of orthographic/phonetic alphabet in L1
                            corr_count = corr_dict[l1][l2][seg1].get(seg2, 0)
                            total_count = sum(corr_dict[l2][l1][seg2].values())
                            smoothed_prob = lidstone_smoothing(x=corr_count, 
                                                               N=total_count,
                                                               d=d_value)
                            
                            #Calculate surprisal of this correspondence from smoothed probability
                            pair_scores.append(surprisal(smoothed_prob))
                        
                        #Calculate the mean surprisal of the alignment
                        alignment_scores.append(mean(pair_scores))
                    
                    #Add all the alignment scores to the dictionary of surprisal
                    #scores for the L1-L2 pair
                    surprisal_dict[(l1, l2)].append(alignment_scores)
                    
                    #Add counts of each word type alignment to dictionary     
                    #during orthographic iteration
                    if label == 'Orthographic':
                        cognate_type_counts[(l1, l2)][i+1] = wordtype_counts
    
    
    #Write a file including all of the calculated scores
    with open(f'Surprisal/NorthWindandSun_Surprisal_{label}.csv', 'w') as f:
        
        #Prepare labels for each measurement in header of CSV
        fragment_labels = [f'Fragment{i}{x}' 
                           for i in range(1,8) 
                           for x in ['', '_normalized']]
        mean_labels = ['Mean', 'Mean_normalized']
        word_types = [f'{c}_F{i}' for i in range(1,8) 
                      for c in ['Cognate_Pairs', 'Non-Cognate_Pairs', 'Gap_Pairs', 'Stimulus_Cognate_Density']]
        columns = ['Language1_native', 'Language2_foreign'] + fragment_labels + mean_labels + word_types
        f.write(",".join(columns))
        f.write('\n')
        
        #Write data for each language pair according to each measurement
        for lang_pair in surprisal_dict:
            l1, l2 = lang_pair
            values = []
            totals, avgs = [], []
            wordtypes = []
            for i in range(len(surprisal_dict[lang_pair])):
                l = surprisal_dict[lang_pair][i]
                total, avg = sum(l), mean(l)
                totals.append(total)
                avgs.append(avg)
                values.extend([total, avg])
                cognate_count = len(cognate_type_counts[(l1, l2)][i+1]['cognate'])
                wordtypes.append(cognate_count)
                wordtypes.append(len(cognate_type_counts[(l1, l2)][i+1]['non-cognate']))
                wordtypes.append(len(cognate_type_counts[(l1, l2)][i+1]['gap']))
                
                #Calculate "cognate density" of stimulus sentence:
                #number of shared cognate words/units normalized by 
                #the number of words/units in the stimulus sentence
                l2_word_count = len([pair for word_type in cognate_type_counts[(l1, l2)][i+1]
                                     for entry in cognate_type_counts[(l1, l2)][i+1][word_type]
                                     if entry[1] != None])
                wordtypes.append(cognate_count/l2_word_count)
                
            values.append(mean(totals))
            values.append(mean(avgs))
            values.extend(wordtypes)
            values = [str(n) for n in values]
            
            f.write(','.join([l1, l2]+values))
            f.write('\n')

#%%
#CALCULATE PERCENTAGE OF COGNATE WORDS PER LANGUAGE PAIR
#ACCORDING TO NORTH WIND AND SUN PARALLEL TRANSLATIONS AND MANUAL COGNACY ANNOTATIONS

checked = []
cognate_counts = {}
with open('Cognate_counts.csv', 'w') as f:
    f.write('L1,L2,N_Cognates,N_UniqueCognates,TotalAligned,PercentCognate\n')
    for lang_pair in cognate_dict:
        l1, l2 = lang_pair
        if (l2, l1) not in checked:
            n_cognates = len(cognate_dict[lang_pair])
            n_unique_cognates = len(set(cognate_dict[lang_pair]))
            total_aligned_matches = all_aligned_matches[all_aligned_matches.L1==l1]
            total_aligned_matches = total_aligned_matches[total_aligned_matches.L2==l2]
            total_aligned_matches = len(total_aligned_matches)
            percent_cognates = round((n_cognates/total_aligned_matches)*100,1)
            f.write(','.join([l1, l2, str(n_cognates), 
                              str(n_unique_cognates), 
                              str(total_aligned_matches), 
                              str(percent_cognates)]))
            f.write('\n')
            checked.append(lang_pair)
            cognate_counts[(l1, l2)] = percent_cognates
            cognate_counts[(l2, l1)] = percent_cognates

#%%
#DRAW DENDROGRAMS BASED ON AVERAGE SURPRISAL BETWEEN LANGUAGE PAIRS

#Designate language names and abbreviations
langs = ['RU', 'BE', 'UK', 'PL', 'CS', 'SK', 'SL', 'HR', 'BG']
lang_names = ['Russian', 'Belarusian', 'Ukrainian', 
              'Polish', 'Czech', 'Slovak', 
              'Slovene', 'Croatian', 'Bulgarian']
lang_abbrev = {langs[i]:lang_names[i] for i in range(len(langs))}


#Load orthographic and phoentic surprisal measures
orth_surprisal = pd.read_csv('Surprisal/NorthWindandSun_Surprisal_Orthographic.csv')
phon_surprisal = pd.read_csv('Surprisal/NorthWindandSun_Surprisal_Phonetic.csv')


#Extract mean surprisal values for each pair in both directions
orth_surprisal_vals, phon_surprisal_vals = {}, {}
for source, dest in zip([orth_surprisal, phon_surprisal],
                        [orth_surprisal_vals, phon_surprisal_vals]):
    for i, row in source.iterrows():
        l1 = row.Language1_native
        l2 = row.Language2_foreign
        mean_nz_surprisal = row.Mean_normalized
        dest[(l1, l2)] = mean_nz_surprisal


def surprisal_dist(l1, l2, df):
    """
    Calculates the mean surprisal distance between the two languages
    given a dictionary of surprisal measurements.

    Parameters
    ----------
    l1 : string
        Name of language 1 as indexed in df.
    l2 : string
        Name of language 2 as indexed in df.
    df : dict
        Dictionary of surprisal measurements with keys in format ('L1', 'L2').

    Returns
    -------
    float
        Mean of L1-L2 and L2-L1 surprisal values.

    """
    return mean([df[(l1, l2)], df[(l2, l1)]])


#Draw dendrograms from mean surprisal measurements
for label, df in zip(['Phonetic', 'Orthographic'],
                     [phon_surprisal_vals, orth_surprisal_vals]):
    draw_dendrogram(group=langs, 
                    labels=lang_names, 
                    dist_func=surprisal_dist, 
                    sim=False,
                    title=f'Surprisal Distance ({label})',
                    save_directory='Surprisal/',
                    df=df)

#%%
#PLOT SCATTERPLOT OF PHONETIC AND ORTHOGRAPHIC SURPRISAL VALUES FOR LANGUAGE PAIRS

#Extract surprisal values
surprisal_data = {'Phonetic':{}, 'Orthographic':{}}
for data, label in zip([phon_surprisal, orth_surprisal], 
                       ['Phonetic', 'Orthographic']):
    for i, row in data.iterrows():
        l1, l2 = row.Language1_native, row.Language2_foreign
        s = row.Mean_normalized
        lang_pair = f'{l1}-{l2}'
        surprisal_data[label][lang_pair] = s

#Plot scatterplot
X = list(surprisal_data['Phonetic'].values())
Y = list(surprisal_data['Orthographic'].values())
labels = list(surprisal_data['Phonetic'].keys())
plt.figure(figsize=(15,12))
plt.scatter(X, Y)
plt.ylabel('Orthographic Adaptation Surprisal', fontsize=22)
plt.xlabel('Phonetic Adaptation Surprisal', fontsize=22)
for i, label in enumerate(labels):
    plt.annotate(label, (X[i], Y[i]), xytext = (X[i]+0.01, Y[i]-0.01))
plt.savefig('Orthographic and Phonetic Surprisal.png', dpi=200)
plt.show()
plt.close()
