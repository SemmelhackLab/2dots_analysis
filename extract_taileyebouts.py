import shelve
from bouts import *
from framemetrics import *
import peakdetector
from operator import add
import numpy as np
import operator


def extract_tail_eye_bout(eyes,vel, tailfits, Fs, bout_thresh, peakthres, delay):

    # --------------------- READ TAILFITS -----------------------------------
    shv = shelve.open(tailfits)
    print shv

    if shv:

        tailfit = shv[shv.keys()[0]].tailfit
        shv.close()

        tailfit = normalizetailfit(tailfit)
        tailangle = -(tail2angles(tailfit))
        boutedges, var = extractbouts(tailangle, bout_thresh)

        # --------------------- READ EYE ANGLES -----------------------------------
        LeftEye = eyes[0]['LeftEye']
        RightEye = eyes[0]['RightEye']
        LeftVel = vel[0]['LeftVel']
        RightVel = vel[0]['RightVel']
        filename = eyes[0]['filename']
        # --------------------- EXTRACT TAIL BOUTS -----------------------------------
        bouts = []
        for bout in boutedges:
            # if boutacceptable(tailfit[bout[0]:bout[1]]):  # tag
            bouts += [{'tail': tailfit[bout[0]:bout[1]],
                       'frames': [bout[0], bout[1]]}]
        '''
        bout_angles = []
        left_eyeangles = []
        right_eyeangles = []
        tailfreq = []
        peaks=[]
        '''
        sum_eyes = map(add, LeftEye, RightEye)
        taileye = []

        '''
        if not bouts:
            taileye += [{'tail': tailangle, 'left': LeftEye, 'right': RightEye, 'sum_eyes': sum_eyes,
                         'LeftVel': LeftVel, 'RightVel': RightVel, 'filename': filename, 'no_bouts': [1]}]
        '''
        # else:
        for i in range(len(bouts)):
            nFrames = len(bouts[i]['tail'])
            # frames = bouts[i]['tail']
            # print bouts[i]['frames'][0]

            left_eyeangle = LeftEye[bouts[i]['frames'][0]:bouts[i]['frames'][1]]
            right_eyeangle = RightEye[bouts[i]['frames'][0]:bouts[i]['frames'][1]]
            left_vel = LeftVel[bouts[i]['frames'][0]:bouts[i]['frames'][1]]
            right_vel = RightVel[bouts[i]['frames'][0]:bouts[i]['frames'][1]]

            if delay > bouts[i]['frames'][0]:
                left_eyeangle_delay = LeftEye[0:bouts[i]['frames'][1]]
                right_eyeangle_delay = RightEye[0:bouts[i]['frames'][1]]
                left_vel_delay = LeftVel[0:bouts[i]['frames'][1]]
                right_vel_delay = RightVel[0:bouts[i]['frames'][1]]
            else:
                left_eyeangle_delay = LeftEye[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]
                right_eyeangle_delay = RightEye[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]
                left_vel_delay = LeftVel[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]
                right_vel_delay = RightVel[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]

            sum_eyeangles = map(add, left_eyeangle, right_eyeangle)

            boutangle = tail2angles(bouts[i]['tail'])  # extract the tailfits of the bout frames
            boutangle = [-1*a for a in boutangle]
            peak = peakdetector.peakdetold(
                boutangle, peakthres)  # get the number of peaks which tells us about how many tail beats for the bout

            peak_new = []
            # print 'Length of peak: ', len(peak[0])
            # print 'Length of frames: ', nFrames
            # print 'Length of tail: ', len(bouts[i]['tail'])
            # print 'Length of tail: ', len(bouts[i]['tail'])
            for item in peak[0]:
                peak_new.append([boutedges[i][0] + item[0], item[1]])
            '''
            left_eyeangles.append(left_eyeangle)
            right_eyeangles.append(right_eyeangle)
            tailfreq.append(len(peak[0]) / float((Fs * nFrames)))
            bout_angles.append(boutangle)
            peaks.append(peak_new)
            '''
            # left_eyeangles = eye angles during bouts
            # left_eyeangles_delay = eye angles during bouts with delay
            # sum_eyeangles = sum of eye angles during bouts
            # sum_eyes = summation of eye angles during the whole trial
            # left_v = eye velocity during the bout
            # left_vel_delay = eye velocity during the bout with delay
            # LeftVel - eye velocity during the whole trial
            # frames = the initial and last frame of the bout
            # filename = the filename of the eye files

            taileye += [{'tail': tailangle, 'bout_angles': boutangle, 'tailfreq': len(peak[0])/(0.0033*nFrames),
                         'left_eyeangles': left_eyeangle, 'right_eyeangles': right_eyeangle,
                         'left_eyeangles_delay': left_eyeangle_delay, 'right_eyeangles_delay': right_eyeangle_delay,
                         'sum_eyeangles': sum_eyeangles, 'left': LeftEye, 'right': RightEye, 'sum_eyes': sum_eyes,
                         'left_v': left_vel, 'right_v': right_vel, 'LeftVel': LeftVel, 'RightVel': RightVel,
                         'frames': bouts[i]['frames'], 'filename': filename, 'left_vel_delay':
                             left_vel_delay, 'right_vel_delay': right_vel_delay, 'no_bouts': [0]}]

    else:
        print "Tailfit empty"
        taileye = []

    return taileye  # bout_angles,tailfreq, left_eyeangles, right_eyeangles


def extract_tail_eye_bout2(eyes,vel, tailangle, Fs, bout_thresh, peakthres, delay):

    if tailangle:

        boutedges, var = extractbouts(tailangle, bout_thresh)

        # --------------------- READ EYE ANGLES -----------------------------------
        LeftEye = eyes[0]['LeftEye']
        RightEye = eyes[0]['RightEye']
        LeftVel = vel[0]['LeftVel']
        RightVel = vel[0]['RightVel']
        filename = eyes[0]['filename']
        # --------------------- EXTRACT TAIL BOUTS -----------------------------------
        bouts = []
        for bout in boutedges:
            # if boutacceptable(tailfit[bout[0]:bout[1]]):  # tag
            bouts += [{'tail': tailangle[bout[0]:bout[1]],
                       'frames': [bout[0], bout[1]]}]

        sum_eyes = map(add, LeftEye, RightEye)
        taileye = []

        for i in range(len(bouts)):

            nFrames = len(bouts[i]['tail'])

            left_eyeangle = LeftEye[bouts[i]['frames'][0]:bouts[i]['frames'][1]]
            right_eyeangle = RightEye[bouts[i]['frames'][0]:bouts[i]['frames'][1]]
            left_vel = LeftVel[bouts[i]['frames'][0]:bouts[i]['frames'][1]]
            right_vel = RightVel[bouts[i]['frames'][0]:bouts[i]['frames'][1]]

            if delay > bouts[i]['frames'][0]:
                left_eyeangle_delay = LeftEye[0:bouts[i]['frames'][1]]
                right_eyeangle_delay = RightEye[0:bouts[i]['frames'][1]]
                left_vel_delay = LeftVel[0:bouts[i]['frames'][1]]
                right_vel_delay = RightVel[0:bouts[i]['frames'][1]]
            else:
                left_eyeangle_delay = LeftEye[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]
                right_eyeangle_delay = RightEye[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]
                left_vel_delay = LeftVel[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]
                right_vel_delay = RightVel[(bouts[i]['frames'][0] - delay):bouts[i]['frames'][1]]

            sum_eyeangles = map(add, left_eyeangle, right_eyeangle)

            boutangle = tailangle[bouts[i]['frames'][0]:bouts[i]['frames'][1]]  # extract the tail bouts
            boutangle = [-1 * a for a in boutangle]
            peak = peakdetector.peakdetold(
                boutangle, peakthres)  # get the number of peaks which tells us about how many tail beats for the bout

            peak_new = []

            for item in peak[0]:
                peak_new.append([boutedges[i][0] + item[0], item[1]])

            taileye += [{'tail': tailangle, 'bout_angles': boutangle, 'tailfreq': len(peak[0])/(0.0033*nFrames),
                         'left_eyeangles': left_eyeangle, 'right_eyeangles': right_eyeangle,
                         'left_eyeangles_delay': left_eyeangle_delay, 'right_eyeangles_delay': right_eyeangle_delay,
                         'sum_eyeangles': sum_eyeangles, 'left': LeftEye, 'right': RightEye, 'sum_eyes': sum_eyes,
                         'left_v': left_vel, 'right_v': right_vel, 'LeftVel': LeftVel, 'RightVel': RightVel,
                         'frames': bouts[i]['frames'], 'filename': filename, 'left_vel_delay':
                             left_vel_delay, 'right_vel_delay': right_vel_delay, 'no_bouts': [0]}]

    else:
        print "Tailfit empty"
        taileye = []

    return taileye  # bout_angles,tailfreq, left_eyeangles, right_eyeangles
