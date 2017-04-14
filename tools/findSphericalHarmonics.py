#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 13:41:23 2017

@author: matt
"""
    
def getSpherHarms(mExp):
    '''
    Get key table line of spherical harmonics for a given site symmetry. Return line as string.
    '''
    spherHarms = []
    for l in range(0,5):            #l
        for mu in range(0,50):       #mu
        
            m = eval(mExp)
            
            if abs(m) > l:
                break
            
            else:
                spherHarms.append('{}{}'.format(l,m))
    
    return spherHarms
    
    
def makeKeyTableLine(spherHarms):
    '''
    Given spherical harmonics in the format '10' make the line in the key table. Return the line as a string.
    '''
    x = ['00', '00x', '11', '1-1', '10', '20', '21', '2-1', '22', '2-2', '30', '31', '3-1', '32', '3-2', '33', '3-3',
         '40', '41', '4-1', '42', '4-2', '43', '4-3', '44', '4-4']
    
    line = ''
    lastl = '0'
    
    for item in x:
        
        if lastl != item[0]:
            line+= ' '
        
        lastl = item[0]
            
        if item in spherHarms:
            line += '1'
        
        else:
            line += '0'
            
    return line
        
def getMultipoles(mExp):
    '''
    Get key table multipole line for a given m expression. Return the line as a string.
    '''
    return makeKeyTableLine(getSpherHarms(mExp))
    
def makeMultipoleBank():
    '''
    Make a dictionary of site symmetries and their correct multipole entries in the key table. Return this dictionary.
    '''
    bank = {'NO':'00 000 00000 0000000 000000000', '1':'10 111 11111 1111111 111111111',
            'cyl':'10 001 00000 0000000 000000000','cylX':'10 001 10000 1000000 000000000', 
            'cylXD':'10 001 10000 0000000 000000000'}
    
    symConds = {'2':'2*mu',
                'm':'l-2*mu',
                'mm2':'2*mu',
                '4':'4*mu',
                '4mm':'4*mu',
                '3':'3*mu',
                '6':'6*mu',
                '6mm':'6*mu'                
                }
    
    for sym, cond in symConds.items():
        bank[sym] = getMultipoles(cond)
        
    return bank
    
print(makeMultipoleBank())
