#!/usr/bin/env python3
# Name: Jonathan Lyda (jlyda)
# Group Members: None

'''
Reads a FASTQ file through std in, uses command line parameters to know what file type to convert the std in file to,
eliminates errors found in the file provided, and redirects the filtered and converted FASTQ file to std out.

input (command line)(illumina 1.8 to Phred+64):
python FastQ+translater.py <illumina1.8.fastq >newesttest1.8to1.3.fastq -P33in -P64out

output (std out) (redirected to newesttest1.8to1.3.fastq): output file is created in current directory
'''

import sys

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - allows user to specify input and output fastq formats.',
            epilog='Program epilog - will convert fastq file format based on user arguments.',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        # self.parser.add_argument('inFile', action='store', help='input file name')
        # self.parser.add_argument('outFile', action='store', help='output file name')
        self.parser.add_argument('-P33in', '--PHRED33input', action='store', nargs='?', const=True, default=False,
                                 help='input of Phred 33 (illumina 1.8)')

        self.parser.add_argument('-P64in', '--PHRED64input', action='store', nargs='?', const=True, default=False,
                                 help='input of Phred 64 (illumina 1.3)')

        self.parser.add_argument('-P64Bin', '--PHRED64withBoffsetinqualityvalues', action='store', nargs='?',
                                 const=True, default=False, help='input of Phred 64 with B offset in quality values')

        self.parser.add_argument('-P64SOLin', '--PHRED64withSOLEXAinterpretationsofQscore', action='store',
                                 nargs='?', const=True, default=False, help='input of Phred 64 Solexa (SOLEXA)')

        self.parser.add_argument('-P33out', '--PHRED33output', action='store', nargs='?', const=True, default=True,
                                 help='desired output of Phred 33')  # default set to True

        self.parser.add_argument('-P64out', '--PHRED64output', action='store', nargs='?', const=True, default=False,
                                 help='desired output of Phred 64')

        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class FastQTranslator:
    """
    Reads a Fastq file, parses it, and executes the specified format conversion using various conversion methods.
    """

    def __init__(self, filename):
        """
        initializes a dictionary of solexa values, the file object, the current Fastq entry information, and
        the current Q score value after conversion.
        """

        self.solValues = {-5: 1, -4: 1, -3: 2, -2: 2, -1: 3, 0: 3, 1: 4, 2: 4, 3: 5, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9,
                          9: 10, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18, 19: 19,
                          20: 20, 21: 21, 22: 22, 23: 23, 24: 24, 25: 25, 26: 26, 27: 27, 28: 28, 29: 29, 30: 30,
                          31: 31, 32: 32, 33: 33, 34: 34, 35: 35, 36: 36, 37: 37, 38: 38, 39: 39, 40: 40}
        self.filename = filename
        self.tempList = [[], [], [], []]   #stores each line in a FASTQ entry when reading one entry at a time
        self.QScore = ''

    def fastqReading(self):
        """
        Parses each FASTQ entry (4 lines), eliminates incorrect entries, and corrects invalid base calls (filters entry).
        """
        file1 = open(self.filename)
        counter = 0
        for line in file1:   #iterates through each line in the given file
            counter += 1
            if ((line.startswith('@')) and (counter == 1)):
                header = line.upper().rstrip()
                self.tempList[0].append(header)    #stores the header if it is correct
            elif counter == 2:
                seq = line.rstrip().replace('*', 'N').replace('.', 'N').replace('n', 'N')
                seq = seq.upper()
                self.tempList[1].append(seq)  #stores the sequence (after capitalizing and making all invalid bases "N"

            elif ((line.startswith('+')) and (counter == 3)):
                plus = line.rstrip()
                self.tempList[2].append(plus)  #stores the 3rd line of the entry (a "+" sign)

            elif ((counter == 4) and (self.tempList[0]) and (self.tempList[1]) and (self.tempList[2])):
                quality = line.rstrip()
                self.tempList[3].append(quality)  #stores quality line
                if (len(quality) == len(seq)):   #checks to make sure line 2 (sequence) and line 4 (quality) are equal
                    for item in self.tempList:
                        yield (item[0])  #gives current stored line data
                    self.tempList = [[], [], [], []]   #reset
                    header = ''
                    seq = ''
                    plus = ''
                    quality = ''
                    counter = 0
                else:
                    self.tempList = [[], [], [], []]  #reset
                    header = ''
                    seq = ''
                    plus = ''
                    quality = ''
                    counter = 0
            else:
                self.tempList = [[], [], [], []]  #reset
                header = ''
                seq = ''
                plus = ''
                quality = ''
                counter = 0

    def ill15ToIll(self):
        """
        Converts Fastq file in Illumina1.5 (PhredB+64) to Illumina 1.3 (Phred+64).
        """
        count = 0
        for line in self.fastqReading():  #iterates through each line of the parsed file
            if ((count == 0) or (count == 1) or (count == 2)):
                self.QScore = line
                count += 1
            else:
                for qual in line:  #iterates through each character in the quality line
                    if (qual == 'B'):
                        qual = '@'      #Phred+64 encoding
                        self.QScore += qual
                    else:
                        newQual = chr((ord(qual)-64)+64)
                        self.QScore += newQual
                count = 0
            yield (self.QScore)  #gives current line information
            self.QScore = ''  #reset

    def ill15ToSanger(self):
        """
        Converts Fastq file in Illumina1.5 (PhredB+64) to Sanger (Phred+33).
        """
        count = 0
        for line in self.fastqReading(): #iterates through each line of the parsed file
            if ((count == 0) or (count == 1) or (count == 2)):
                self.QScore = line
                count += 1
            else:
                for qual in line:   #iterates through each character in the quality line
                    if (qual == 'B'):
                        qual = '!'
                        self.QScore += qual
                    else:
                        newQual = chr((ord(qual)-64)+33)
                        self.QScore += newQual
                count = 0
            yield (self.QScore)
            self.QScore = ''

    def solToIll(self):
        """
        Converts Fastq file in Solexa (Phred+64 Solexa) to Illumina 1.3 (Phred+64).
        """
        count = 0
        for line in self.fastqReading(): #iterates through each line of the parsed file
            if ((count == 0) or (count == 1) or (count == 2)):
                self.QScore = line
                count += 1
            else:
                for qual in line:   #iterates through each character in the quality line
                    newQual = (ord(qual)-64)
                    newQual2 = self.solValues[newQual]    #uses solexa dictionary for conversion
                    qual = chr((newQual2)+64)
                    self.QScore += qual
                count = 0
            yield (self.QScore)
            self.QScore = ''

    def solToSanger(self):
        """
        Converts Fastq file in Solexa (Phred+64 Solexa) to Sanger (Phred+33).
        """
        count = 0
        for line in self.fastqReading(): #iterates through each line of the parsed file
            if ((count == 0) or (count == 1) or (count == 2)):
                self.QScore = line
                count += 1
            else:
                for qual in line:   #iterates through each character in the quality line
                    newQual = (ord(qual)-64)
                    newQual2 = self.solValues[newQual]   #uses solexa dictionary for conversion
                    qual = chr((newQual2)+33)
                    self.QScore += qual
                count = 0
            yield (self.QScore)
            self.QScore = ''

    def sangerToIll(self):
        """
        Converts Fastq file in Sanger (Phred+33) to Illumina 1.3 (Phred+64).
        """
        count = 0
        for line in self.fastqReading(): #iterates through each line of the parsed file
            if ((count == 0) or (count == 1) or (count == 2)):
                self.QScore = line
                count += 1
            else:
                for qual in line:    #iterates through each character in the quality line
                    newQual = chr((ord(qual)-33)+64)
                    self.QScore += newQual
                count = 0
            yield (self.QScore)
            self.QScore = ''

    def illToSanger(self):
        """
        Converts Fastq file in Illumina 1.3 (Phred+64) to Sanger (Phred+33).
        """
        count = 0
        for line in self.fastqReading(): #iterates through each line of the parsed file
            if ((count == 0) or (count == 1) or (count == 2)):
                self.QScore = line
                count += 1
            else:
                for qual in line:   #iterates through each character in the quality line
                    newQual = chr((ord(qual)-64)+33)
                    self.QScore += newQual
                count = 0
            yield (self.QScore)
            self.QScore = ''

    def sameType(self):
        """
        Returns the same FASTQ file format back to the user if specified to do so.
        """
        for line in self.fastqReading(): #iterates through each line of the parsed file
            self.QScore = line
            yield (self.QScore)
            self.QScore = ''

def main(inCL=None):
    '''
    Parse the given std in file, use commandline arguments to choose which conversion method to use, converts the Fastq
    file to the specified format, redirects the new translated information to a std out file.
    '''
    if inCL is None:  # if no parameters are entered on command line
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    filename = sys.stdin.fileno()   #uses std in
    myFastq = FastQTranslator(filename)   #creates a myFastq object using std in file

    if myCommandLine.args.PHRED33input is True:  # Phred33 to Phred64 or Phred33 (no format change)
        if myCommandLine.args.PHRED64output is True:
            for item in myFastq.sangerToIll():
                print(item)
        else:
            for item in myFastq.sameType():
                print(item)


    if myCommandLine.args.PHRED64input is True:  # Phred64 to Phred33 or Phred64 (no format change)
        if myCommandLine.args.PHRED33output is True:
            for item in myFastq.illToSanger():
                print(item)
        else:
            for item in myFastq.sameType():
                print(item)

    if myCommandLine.args.PHRED64withBoffsetinqualityvalues is True:  # Phred64B to Phred33 or Phred64
        if myCommandLine.args.PHRED64output is True:
            for item in myFastq.ill15ToIll():
                print(item)
        else:
            for item in myFastq.ill15ToSanger():
                print(item)

    if myCommandLine.args.PHRED64withSOLEXAinterpretationsofQscore is True:  # Phred64 (Solexa) to Phred33 or Phred64
        if myCommandLine.args.PHRED64output is True:
            for item in myFastq.solToIll():
                print(item)
        else:
            for item in myFastq.solToSanger():
                print(item)

if __name__ == "__main__":
    main()  # uses std in
