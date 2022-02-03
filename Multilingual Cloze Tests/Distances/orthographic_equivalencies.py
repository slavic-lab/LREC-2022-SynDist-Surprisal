base_chs = {#Latin vowels
            'á':'a',
            'ä':'a',
            'ą':'a',
            'é':'e',
            'ě':'e',
            'ę':'e',
            'í':'i',
            'ó':'o',
            'ô':'o',
            'ú':'u',
            'ů':'u',
            'ý':'y',
            
            #Latin consonants
            'ć':'c',
            'č':'c',
            'ď':'d',
            'đ':'d',
            'ł':'l',
            'ĺ':'l',
            'ľ':'l',
            'ń':'n',
            'ň':'n',
            'ř':'r',
            'ŕ':'r',
            'š':'s',
            'ś':'s',
            'ť':'t',
            'ż':'z',
            'ź':'z',
            'ž':'z',
            
            #Cyrillic bases
            'ґ':'г',
            'ё':'е',
            'ї':'і',
            'щ':'ш',
            #'ў':'у' #don't associate <ў> with <у> because the former is a consonant and the latter is a vowel

            #Latin equivalents of Cyrillic characters
            'â':'a',
            'ë':'e',
            'ê':'e',
            'è':'e',
            'ï':'i',
            'ì':'i',
            'û':'u',
            'ŝ':'s'
            }

cyr_lat_eqs = {'а':'a',
               'б':'b',
               'в':'v',
               'г':'g',
               'ґ':'g',
               'д':'d',
               'е':'e',
               'ж':'ž',
               'з':'z',
               'и':'i',
               'й':'j',
               'к':'k',
               'л':'l',
               'м':'m',
               'н':'n',
               'о':'o',
               'п':'p',
               'р':'r',
               'с':'s',
               'т':'t',
               'у':'u',
               'ф':'f',
               'х':'h',
               'ц':'c',
               'ч':'č',
               'ш':'š',
               'щ':'ŝ',
               'ъ':'"',
               'ы':'y',
               'ь':"ʹ",
               'э':'è',
               'ю':'û',
               'я':'â',
               'ё':'ë',
               'є':'ê',
               'і':'ì',
               'ї':'ï',
               'ў':'ŭ',
               "'":'"',
               }

vowel_chs = {#Latin vowels 
             'a', 'á', 'ä', 'ą', 
             'ă', 'â', #equivalents of Bulgarian <ъ> and Cyrillic <я>
             'e', 'é', 'ě', 'ę',
             'ë', 'ê', 'è', #equivalents of Cyrillic <ё, є, э>
             'i', 'í', 
             'ï', 'ì', #equivalents of Cyrillic <ї, і>
             'o', 'ó', 'ô',
             'u', 'ú', 'ů',
             'û', #equivalent of Cyrillic <ю>
             'y', 'ý',
             
             #Cyrillic vowels
             'а', 'е', 'и', 'о', 'у', 'ы', 
             'э', 'я', 'ё', 'є', 'і', 'ї'
             }


consonant_chs = {#Latin consonants
                 'b', 
                 'c', 'ć', 'č', 
                 'd', 'ď', 'đ',
                 'f', 'g', 'h', 'j', 'k',
                 'l', 'ł', 'ĺ', 'ľ',
                 'm',
                 'n', 'ń', 'ň',
                 'p', 'q',
                 'r', 'ř', 'ŕ',
                 's', 'š', 'ś',
                 't', 'ť', 
                 'ŭ', #equivalent of Belarusian <ў>
                 'v', 'w', 'x', 
                 'ż', 'ź', 'ž',
                 
                 #Cyrillic consonants
                 'б', 'в', 'г', 'ґ', 'д'
                 'ж', 'з', 'й' ,'к', 'л',
                 'м', 'н', 'п', 'р', 'с', 
                 'т', 'ф', 'х', 'ц', 'ч', 
                 'ш', 'щ', 'ў'}

softhard_signs = {'ь', 'ъ', "'"}

#%%
def orth_dist(ch1, ch2, lang1=None, lang2=None):
    #First lower-case both characters
    ch1, ch2 = ch1.lower(), ch2.lower()
    
    #Check if both characters are from the same alphabet
    #Check using first element in case character is a digraph/trigraph
    ch1_cyr, ch2_cyr = False, False
    if ch1[0] in cyr_lat_eqs:
        ch1_cyr = True
    if ch2[0] in cyr_lat_eqs:
        ch2_cyr = True

    #Create accommodations for treating <ъ> as a vowel in Bulgarian only
    if lang1 != None:
        if lang1.lower() in ['bg', 'bulgarian']:
            if ch1 == 'ъ':
                ch1 = 'ă'
    
    if lang2 != None:
        if lang2.lower() in ['bg', 'bulgarian']:
            if ch2 == 'ъ':
                ch2 = 'ă'
    
    #Convert Cyrillic characters to Latin equivalents if alphabets don't match
    if (ch1_cyr, ch2_cyr).count(False) == 1:
        ch1 = ''.join([cyr_lat_eqs.get(ch1[i], ch1[i]) for i in range(len(ch1))])
        ch2 = ''.join([cyr_lat_eqs.get(ch2[i], ch2[i]) for i in range(len(ch2))])
        
    
    #Identical characters
    if ch1 == ch2:
        return 0
    
    #Characters with same base but different diacritics
    base1 = ''.join([base_chs.get(ch1[i], ch1[i]) for i in range(len(ch1))])
    base2 = ''.join([base_chs.get(ch2[i], ch2[i]) for i in range(len(ch2))])
    if base1 in base2:
        return -0.5
    elif base2 in base1:
        return -0.5    


    #Different characters: check using first element in case of digraphs/trigraphs
    if ch1[0] in vowel_chs:
        #Both vowels
        if ch2[0] in vowel_chs:
            return -1
        #Vowel and non-vowel
        else:
            return -4.5
    
    elif ch1[0] in consonant_chs:
        #Both consonants
        if ch2[0] in consonant_chs:
            return -1
        #Consonant and non-consonant
        else:
            return -4.5
    
    else: #ch1 is soft or hard sign
        if ch2 in softhard_signs:
            return -1
        #Soft/hard sign with vowel or consonant
        else:
            return -4.5
    