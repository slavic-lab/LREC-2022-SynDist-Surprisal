#IMPLEMENTATION OF SYNTACTIC DISTANCE MEASURES
#Based on measured proposed in Heeringa et al. (2017)
#Code written by Philip Georgis (2020-22) in cooperation with SFB 1102 C4 Research Project

#Load required modules
from collections import defaultdict
from openpyxl import load_workbook
import math



def load_alignfile(excel_alignment):
    """
    Loads an  Excel file containing sentence/phrase alignments.

    Parameters
    ----------
    excel_alignment : string
        Path to Excel file.

    Returns
    -------
    openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file.

    """
    return load_workbook(excel_alignment, data_only = True)



def get_color(workbook, sheet, cell):
    """
    Returns the hexadecimal color value of a particular cell, e.g. cell 'K44'.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file.
    sheet : string
        Label of sheet within Excel workbook.
    cell : string
        Cell coordinates within Excel sheet, e.g. 'K44'.

    Returns
    -------
    string
        Hexadecimal color code of cell fill.

    """
    
    sh = workbook[sheet]
    return str(sh[cell].fill.start_color.index)



#Hexadecimal color code for white/blank cell
blank = '00000000'



def syntactic_dist(workbook, sheet, 
                   start_row, start_column, rows_down=1, 
                   print_movement=False, print_indel=False,
                   POS_workbook=None, POS_sheet=None):
    """
    Measures syntactic InDel and Movement distances between two rows in of a 
    syntactic alignment file. If a POS_workbook and POS_sheet are provided,
    measurements will take into account the overlap/match of associated POS tags.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file containing syntactic alignment.
    sheet : string
        Label of sheet within Excel workbook containing syntactic alignment.
    start_row : int
        Row number within Excel sheet to use as base row for comparison.
    start_column : string
        Column label (e.g. 'B') of first column in alignment.
    rows_down : int, optional
        Number of rows below to measure syntactic distance with respect to
        start_row. The default is 1 row down for simple pairwise alignents. 
        A negative number results in measuring a row above the start_row.
    print_movement : bool, optional
        Whether or not to print matched syntactic units which appear in 
        different parts of the alignment. The default is False.
    print_indel : bool, optional
        Whether or not to print instances of insertion/deletion. The default is False.
    POS_workbook : openpyxl.workbook.workbook.Workbook, optional
        Loaded Excel workbook containing syntactic alignment with POS tags. 
        The default is None.
    POS_sheet : TYPE, optional
        Label of sheet within Excel workbook containing syntactic alignment
        with POS tags. The default is None.

    Returns
    -------
    align_length : int
        Length of syntactic alignment.
    indel : int
        Number of insertions and deletions.
    movement : int
        Total number of positions by which corresponding words have been reordered.
    moved_words : int
        Number of corresponding words which have been reordered.

    """
    
    #Load syntactic alignment workbook sheet
    sh = workbook[sheet]
    
    #Load the POS annotation sheet if given
    eval_POS = False
    if POS_workbook != None:
        assert POS_sheet != None
        POS_sh = POS_workbook[POS_sheet]
        eval_POS = True
    
    #Calculate list of non-empty columns and rows in the alignment sheet
    columns = [i[0].column_letter for i in list(sh.columns)]
    rows = list(range(1,len(list(sh.rows))+1))
    
    #Calculate list of columns to the right of the specified starting column
    right_columns = [columns[i] for i in range(len(columns)) 
                     if i >= columns.index(start_column)]
    
    #Designate the row numbers of of L1 and L2 within alignment
    l1_row = start_row
    l2_row = start_row + rows_down
    
    #Iterate through columns: identify and count insertions/deletions, and moved units
    align_length, indel, movement, moved_words = 0, 0, 0, 0
    for c in range(len(right_columns)):
        column = right_columns[c]
        l1_cell = sh[f'{column}{l1_row}']
        l2_cell = sh[f'{column}{l2_row}']
        
        #If at least one of the two cells is not blank:
        if (l1_cell.value != None) or (l2_cell.value != None):
            align_length += 1
            
            #If the L2 cell is blank:
            if l2_cell.value == None:
                
                #Determine color of L1 cell and look for corresponding unit
                #in other columns, based on matching colors
                l1_color = get_color(workbook, sheet, f'{column}{l1_row}')
                
                #Count as deletion if the L1 cell has no color
                if l1_color == blank:
                    indel += 1
                    if print_indel == True:
                        print('InDel+=1: ', sh[f'{column}{l1_row}'].value, " | ", sh[f'{column}{l2_row}'].value)
                
                #Otherwise: count as movement, look for corresponding unit(s)
                else:
                    l2_matches = [i for i in range(len(right_columns))
                                  if get_color(workbook, sheet, f'{right_columns[i]}{l2_row}') == l1_color]
                    
                    #If there are no color matches, count it as InDel and 
                    #skip movement calculation (this situation would occur 
                    #if it is a multiple alignment and the cell is color coded 
                    #to show movement with respect to a different row)
                    if len(l2_matches) == 0:
                        indel += 1
                        if print_indel == True:
                            print('InDel+=1: ', sh[f'{column}{l1_row}'].value, " | ", sh[f'{column}{l2_row}'].value)
                        continue
                    
                    #Otherwise: calculate the amount of movement
                    distances = [abs(c - l2_matches[i]) for i in range(len(l2_matches))]
                    
                    #Iterate through cells between matched units and count 
                    #how many are empty in both languages
                    #Adjust distances by subtracting number of empty cells
                    #(relevant/necessary only for multiple alignment scheme)
                    for i in range(len(l2_matches)):
                        empty_cells = 0
                        l2_match = l2_matches[i]
                        for k in range(min(l2_match, c), max(l2_match, c)):
                            if (sh[f'{right_columns[k]}{l1_row}'].value == None) and (sh[f'{right_columns[k]}{l2_row}'].value == None):
                                empty_cells += 1
                        distances[i] -= empty_cells
                        
                    #Calculate final movement and InDel distances
                    if len(distances) > 0:
                        
                        #Calculate average distance among all matched units
                        distance = sum(distances) / len(distances) 
                        movement += distance
                        
                        #If POS annotations are specified, consider the POS of the 
                        #matched moved units for InDel calculation
                        if eval_POS == True:

                            for l2_c in l2_matches:
                                l2_c = right_columns[l2_c]
                                l1_POS = POS_sh[f'{column}{l1_row}'].value.strip()
                                l2_POS = POS_sh[f'{l2_c}{l2_row}'].value.strip()
                    
                                #Add InDel distance of 0.5 if the POS labels don't match
                                if l1_POS != l2_POS:
                                    indel += 0.5
                                    if print_indel == True:
                                        print('InDel+=0.5: ', sh[f'{column}{l1_row}'].value, " | ", sh[f'{column}{l2_row}'].value)
                                
                    
                    #Update calculation of binary movement distance
                    moved_words += len([l2_matches[i] for i in range(len(l2_matches))
                                        if l2_matches[i] != c])
                    
                    #Print matched units
                    if print_movement == True:
                        for i in range(len(l2_matches)):
                            l2_cell_value = sh[f'{right_columns[l2_matches[i]]}{l2_row}'].value
                            dist = abs(c - l2_matches[i])
                            print(f'Matching "{l1_cell.value}" with "{l2_cell_value}" (distance = {dist})')
            
            
            #If the L1 cell is blank:
            elif l1_cell.value == None:
                
                #Get the color of the L2 cell 
                l2_color = get_color(workbook, sheet, f'{column}{l2_row}')
                
                #Count as deletion if the L2 cell has no color
                if l2_color == blank:
                    indel += 1
                    if print_indel == True:
                        print('InDel+=1: ', sh[f'{column}{l1_row}'].value, " | ", sh[f'{column}{l2_row}'].value)
                
                #Otherwise: check for matches
                else: 
                    
                    #Check if this L2 cell color matches the color from one of the L1 cells
                    #If not, then count it as InDel
                    #This is only the case in multiple alignment scheme where 
                    #the L2 cell's color is meant to match with some other row
                    L1_colors = [get_color(workbook, sheet, f'{column}{l1_row}') for column in right_columns]
                    if l2_color not in L1_colors:
                        indel += 1
                        if print_indel == True:
                            print('InDel+=1: ', sh[f'{column}{l1_row}'].value, " | ", sh[f'{column}{l2_row}'].value)
                    
                    #Otherwise, a moved equivalent has been found elsewhere in the alignment
                    #But we only evaluate this case in one direction; pass
                    else:
                        pass
                    
            #If neither cell is blank, there is a matched equivalent syntactic unit
            else:
                
                #Consider the POS of the matched units if specified
                if eval_POS == True:
                    l1_POS = POS_sh[f'{column}{l1_row}'].value.strip()
                    l2_POS = POS_sh[f'{column}{l2_row}'].value.strip()
                    
                    #Add InDel distance of 0.5 if the POS labels don't match
                    if l1_POS != l2_POS:
                        indel += 0.5
                        if print_indel == True:
                            print('InDel+=0.5: ', sh[f'{column}{l1_row}'].value, " | ", sh[f'{column}{l2_row}'].value)
                
    return align_length, indel, movement, moved_words



def multiple_syntactic_distance(workbook, sheet, start_row, start_column, 
                                n_alignments,
                                print_movement=False,
                                **kwargs):
    """
    Calculate syntactic distance for all pairwise combinations of sentences
    in a multiple alignment scheme.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file containing syntactic alignment.
    sheet : string
        Label of sheet within Excel workbook containing syntactic alignment.
    start_row : int
        Row number within Excel sheet to use as base row for comparison.
    start_column : string
        Column label (e.g. 'B') of first column in alignment.
    n_alignments : int
        Number of rows in alignment (e.g. number of sentences/languages).
    print_movement : bool, optional
        Whether or not to print matched syntactic units which appear in 
        different parts of the alignment. The default is False.
    **kwargs : dict
        Dictionary of named keyword arguments to be passed to syntactic_dist().

    Returns
    -------
    syntactic_distances : collections.defaultdict
        Nested default dictionary of syntactic distances for each pair of sentences
        indexed by their IDs.

    """
    
    #Load workbook and sheet with multiple syntactic alignment
    sh = workbook[sheet]
    
    #Calculate list of non-empty columns in the alignment sheet
    columns = [i[0].column_letter for i in list(sh.columns)]
    
    #Designate column containing sentence/language IDs/names as the column
    #preceding the start column
    name_column = columns[columns.index(start_column)-1]
    names = [sh[f'{name_column}{str(start_row+i)}'].value for i in range(n_alignments)]
    
    #Iterate through pairs of rows and calculate pairwise syntactic distance between each
    syntactic_distances = defaultdict(lambda:{})
    for i in range(n_alignments):
        for j in range(i+1, n_alignments):
            dist = syntactic_dist(workbook, sheet, start_row+i, start_column, rows_down=j-i, print_movement=print_movement, **kwargs)
            syntactic_distances[names[i]][names[j]] = dist
            syntactic_distances[names[j]][names[i]] = dist  
    
    return syntactic_distances



def chunk_list(lis, n):
    """
    Splits a list into sublists of length n; if not evenly divisible by n,
    the final sublist contains the remainder.

    Parameters
    ----------
    lis : list
        List of items to be chunked.
    n : int
        Length of chunks.

    Returns
    -------
    list
        Chunked list.

    """

    return [lis[i * n:(i + 1) * n] for i in range((len(lis) + n - 1) // n)]



def calculate_startlines(workbook, sheet, index_row, n_alignments, n_spaces):
    """
    Returns a list of start lines (row numbers) for alignments in a workbook sheet,
    given the row number of the first alignment's index, the number of phrases
    per alignment, and the number of spaces between alignments.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file containing alignments.
    sheet : string
        Label of sheet within Excel workbook containing alignments.
    index_row : int
        Row number of the first alignment's index (line above alignment).
        If first alignment begins in first row of workbook, set to 0.
    n_alignments : int
        Number of sentences per alignment (e.g. 2 if pairwise).
    n_spaces : int
        Number of blank rows between alignments.

    Returns
    -------
    startlines : list
        List of alignment start rows.

    """
    
    #Load workbook sheet containing alignments
    sh = workbook[sheet]
    
    #Get list of rows in workbook sheet
    rows = list(sh.rows)
    
    #Disregard rows before the index row
    if index_row > 0:
        rows = rows[index_row-1:]
    
    #Split list of rows into chunks of specified size
    alignment_indices = chunk_list(list(range(len(rows))), n_alignments+n_spaces+1)
    
    #The start line for each alignment is the first row in each chunk +1 to account for the index row
    startlines = [i[0]+1 for i in alignment_indices]
    
    return startlines



def extract_phrases(workbook, sheet, index_row, label_column, n_alignments, 
                    remove_empty=True):
    """
    Extracts phrases from syntactic alignment format.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file containing alignments.
    sheet : string
        Label of sheet within Excel workbook containing alignments.
    index_row : int
        Row number of the first alignment's index (line above alignment).
        If first alignment begins in first row of workbook, set to 0.
    label_column : str
        Column (e.g. 'A') containing IDs/labels for each sentence.
    n_alignments : int
        Number of sentences per alignment (e.g. 2 if pairwise).
    remove_empty : bool, optional
        Whether to ignore/remove empty cells when extracting phrases. The default is True.

    Returns
    -------
    cell_values : dict
        Dictionary indexed by sentencee labels/IDs with the extracted
        text for each in string format.

    """
    #Load workbook sheet containing alignments
    sh = workbook[sheet]
    
    #Calculate list of non-empty columns in the alignment sheet
    columns = [i[0].column_letter for i in list(sh.columns)]
    
    #Designate columns for label and start of alignment
    start_column_index = columns.index(label_column)
    start_column = columns[start_column_index+1]
    
    #Extract labels/IDs for each sentence
    labels = [sh[f'{label_column}{i+1+index_row}'].value for i in range(n_alignments)]
    
    #Extract text from each row into dictionary
    cell_values = {labels[i]:[sh[f'{columns[start_column_index+j]}{index_row+1+i}'].value 
                              for j in range(start_column_index+1, len(columns))] 
                   for i in range(len(labels))}
    
    #Remove blank text from empty cells
    if remove_empty == True:
        cell_values = {label:' '.join([word for word in cell_values[label] 
                                       if word != None]) 
                       for label in cell_values}
    return cell_values
    


def write_measures(workbook, sheet, startlines, phrases, syndists,
                   label1, label2,
                   output_file):
    """
    Given an Excel alignment sheet, a list of start lines, extracted phrases, 
    measured syntactic distances, and labels for the languages, an output file
    is written displaying each sentence and its respective syntactic measurements.

    Parameters
    ----------
    workbook : openpyxl.workbook.workbook.Workbook
        Loaded Excel workbook file containing alignments.
    sheet : string
        Label of sheet within Excel workbook containing alignments.
    startlines : list
        List of row numbers marking the start of each alignment.
    phrases : dict
        Dictionary of extracted phrases indexed by sentence/language IDs.
    syndists : dict
        Dictionary of syntactic distances indexed by start lines and sentence/language IDs.
    label1 : str
        Label for the first sentence/language. 
        Must match label in phrase and syndists dictionaries.
    label2 : str
        Label for the second sentence/language.
        Must match label in phrase and syndists dictionaries.
    output_file : str
        File path to output file.

    Returns
    -------
    None.

    """
    
    with open(output_file, 'w') as f:
        f.write('\t'.join(['Phrase Index', label1, label2, 'Length', 
                           'INDEL', 'INDEL (Normalized)', 'Linear Movement', 
                           'Linear Movement (Normalized)', 'Log Movement', 
                           'Log Movement (Normalized)', 'Binary Movement', 
                           'Binary Movement (Normalized)']))
        f.write('\n')
        for i in range(len(startlines)):
            index = startlines[i]
            phrase1 = phrases[index][label1]
            phrase2 = phrases[index][label2]
            syndist = syndists[index][label1][label2]
            align_length, indel, movement, moved_words = syndist
            log_movement = math.log(movement) if movement > 0 else 0
            f.write('\t'.join([str(item) for item in [i+1, phrase1, phrase2, 
                                                      align_length, indel, 
                                                      indel/align_length, 
                                                      movement, movement/align_length, 
                                                      log_movement, log_movement/align_length, 
                                                      moved_words, moved_words/align_length]]))
            f.write(f'\n')
        
