import os, sys
import re

## Class which reads a text file and parses the text into a 
## list of lists of strings

class ParseTxt:
    def __init__(self):
        self.inF = None
        self.parsed_data = []
        self.lines = []
    #

    # Parses text

    def Read(self, inFN, stride = 1):
        self.inF = None
        self.parsed_data = []
        self.inF = open(inFN, "r")
        all_lines = self.inF.readlines()
        all_lines = all_lines[::stride]
        
        j = -1
        for i in range(len(all_lines)):
            line = all_lines[i].rstrip('\n')
            if len(line) > 0:
                j = j + 1
                self.lines.append(line)
                self.parsed_data.append([])
                for word in re.findall(r'\S+', line.strip()):
                    self.parsed_data[j].append(word)
        
        self.inF.close()
    #

    #def Read(self, inFN):
    #    self.inF = None
    #    self.parsed_data = []
    #    self.inF = open(inFN, "r")
    #    all_lines = self.inF.readlines()
    #    lines = []
    #    
    #    j = -1
    #    for i in range(len(all_lines)):
    #        line = all_lines[i].rstrip('\n')
    #        if len(line) > 0:
    #            j = j + 1
    #            lines.append(line)
    #

    # Dump parsed data

    def Print(self):
        print (self.parsed_data)
#
