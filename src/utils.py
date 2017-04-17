'''
#####################################################################
#-------------------UTILITIES----------------------------------------
#####################################################################
'''

def lab2type(atomLabel):
    '''
    Convert atom label to atom type, e.g. 'C(2)' to 'C'. Return type.
    '''
    return atomLabel.split('(')[0].upper()

def spec2norm(atomLabel):
    '''
    Convert special atom label to normal atom label, e.g. 'C(2),asym' to 'C(2)'. Return label.
    '''
    return atomLabel.split(',')[0]

def rawInput2labels(rawInput):
    '''
    Convert raw input of unspecified format atom labels to list of correctly formatted labels. Return list.
    '''
    return formatLabels(labels2list(rawInput))

def labels2list(inputText):
    '''
    Convert raw input of atom labels to list of unformatted atom labels. Return list.
    '''
    inputText = inputText.upper()
    if ',' in inputText:
        inputText = inputText.replace(' ','').strip(',')
        inputAtomList = inputText.split(',')
    elif ' ' in inputText:
        inputAtomList = inputText.split()
    else:
        inputAtomList = [inputText]

return inputAtomList

def isfloat(x):
    '''
    Check if unknown value can be converted to a float. Return result as bool.
    '''
    floatBool = True

    try:
        float(x)
    except Exception:
        floatBool = False

    return floatBool
