# CONTAINS THE FUNCTION FOR READING
# THE EYE AND TAIL MOVEMENTS
# FROM A CSV FILE
# LAST UPDATE 27 JUNE 2019
import os
import csv


def eye_reader(eyefile):
    '''
    :param eyefile: csv file containing the left and right eye angles
    :return: the left and right eye angles
    '''
    # READ EYE POSITION from the eye tracker

    LeftEye = [] # a list to store the angles
    RightEye = []

    with open(eyefile) as csvfile:

        filename, file_extension = os.path.splitext(eyefile) # split the filename
        readCSV = csv.reader(csvfile, delimiter=',') # declare the csv file separated by comma
        headers = readCSV.next() # store the headers of each column of interest
        left = headers.index('left') # first column with title left
        right = headers.index('right')

        for row in readCSV: # read the data of each rows of each column
            LeftEye.append(float(row[left]))
            RightEye.append(float(row[right]))

        eyes = [{'LeftEye': LeftEye, 'RightEye': RightEye, 'filename': filename}] # store the angles in dictionary

    return eyes



def tail_reader(tailfile):
    # READ TAIL POSITION from the 2p tail tracker

    tailangle = [] # list to store the tail angles

    with open(tailfile) as csvfile:

        filename, file_extension = os.path.splitext(tailfile) # split the filename
        readCSV = csv.reader(csvfile, delimiter=',') # declare the csv file separated by comma
        headers = readCSV.next() # store the headers of each column of interest
        t = headers.index('Tail') # first column with title Tail

        for row in readCSV: # read the data of each rows of each column
            tailangle.append(float(row[t]))

    return tailangle