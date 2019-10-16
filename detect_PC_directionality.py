# THIS CODE WILL DETERMINE WHETHER THE TRIAL HAS A PREY CAPTURE LIKE KINEMATICS
# BASED ON THE ONSET OF THE CONTRA-LATERAL EYE CONVERGENCE
# The onset is defined as the point between the peak of eye velocity
# and the start of the tail bout.
# Last update: 27 JUNE 2019, Ivan

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.signal import savgol_filter

from csvfiles_reader import *
from extract_taileyebouts import *
from eye_movement_filters import *
from get_preylocation import *
import operator
from itertools import izip_longest
import seaborn as sns
from detect_preycapture import *

plt.style.use(['dark_background'])

# ============================SET PARAMETERS=======================================================

# ---------------------VIDEO RECORDING PARAMETERS------------------------------------------------

Fs = 300 # Sampling frequency used for the video
tfactor = 0.3  # convert frames to ms

# ---------------------PARAMETERS FOR EXTRACTING THE BOUTS------------------------------------------------
bout_thresh = 0.40  # 0.00 - 1.00 threshold value for extracting bouts, higher more bouts
peakthres = 4  # 0.00 - 20.00 lower vale more peaks, for calculating tail beat frequency

# ---------------------PARAMETERS FOR DETECTING THE PREY CAPTURE BOUTS, MINIMUM VALUES-----------------------------------
filter_avg_tail = 5.0 # threshold for average tail bout angle
filter_avg_binocular = 17.0 # threshold for eye angle binocular convergence, unit in degrees
filter_length_bout = 10 # length of bout based on frame number, unit in frame number
filter_eye_vel = -0.1 # diverging eye velocity in degrees/ms, unit in deg/ms
filter_eye_diverge = -1.0 # eye divergence by degrees, unit in deg
thresh_saccade_speed = 0.2 # the onset of eye movement should be greater than saccade threshold, unit in deg/ms
sensory_delay = 30 # delay of sensory transformation in number of frames, unit in number of frames

# ---------------------VISUAL PREY PARAMETERS------------------------------------------------
preypoints = [10, 70] # direction of prey in visual angle, e.g. [10,70] means going from 10deg to 70 deg
sensory_delay = 30 # how many frames to consider for sensory delay
delay = 50 # number of frames to be subtract or delay to the tail onset for finding the eye onset
speed = 40.0 # prey speed in degrees per second
padding_nonstimulus = [1800, 1650] # to add the frames as a waiting time to determine prey location, visual stimulus appears usually 3 s after trial
# ---------------------LOW PASS FILTER------------------------------------------------
order = 3
fs = 300.0  # sample rate, Hz
cutoff = 10  # desired cutoff frequency of the filter, Hz

# ============================ GENERATE ALL THE FILES =================================================================
# THIS PART IS ONLY FOR COMPILING ALL THE INPUT FILES. THIS IS SUBJECTIVE, IT DEPENDS ON HOW YOU ORGANIZE YOUR DATA
# I organize it this way home/experiment/fish/data (eye, tail), This part will combine all the experiments I have
# Directories
maindir = 'D:\\Semmelhack lab\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\'
dir = 'D:\\Semmelhack lab\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\'

maindir = 'E:\\new folder\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\'
dir = 'E:\\new folder\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\'

dir_output = dir

eye_files = [] # DECLARE A LIST TO STORE THE EYE ANGLE FILES
tail_files = []

fishIDs = list() # LIST FOR STORING EACH FISH NAME
fishIDs_dir = list() # LIST FOR SOME INFO FOR EACH FISH
fishIDs += [fish for fish in os.listdir(dir) if ".csv" not in fish]  # store the fish filenames
#fishIDs = fishIDs[0:-2]
# STORE [BIG PREY OR LEFT PREY FOR EACH FISH, RESPONSE TIME OF THE PREY CAPTURE BOUT]
fishIDs_dir += [{str(fish): [0, 0], str(fish) +'_bigorsmall': [], str(fish) +'_resptime': [] } for fish in fishIDs]

for fish in fishIDs: # READ EACH FISH EYE AND TAIL ANGLES

    eyefile = dir + fish + '\\results\\'
    tailfile = dir + fish + '\\tail\\tail results\\'

    for file in os.listdir(eyefile):
        if file.endswith(".csv"):
            # store the [file directory, directory output, prey direction, fish name, prey speed, filename of eye angles]
            eye_files.append([os.path.join(eyefile, file), dir_output, preypoints, fish, speed, file])
    for file in os.listdir(tailfile):
        if file.endswith(".csv"):
            tail_files.append([os.path.join(tailfile, file), dir_output, preypoints, fish, speed, file])

print len(eye_files), len(tail_files)
print len(fishIDs)

# print fishIDs
# ============================MACRO FOR PREY SELECTION ANALYSIS=======================================================

prey_positions1 = []
eye_diffs = []
tails = []
rtime = []
one_bout = 0
multi_bout = 0
nbouts = 0
eye_onsets = []
fish_movements = []
right_prey = 0
left_prey = 0
PC = 0
fish_PC = []
prey_side = []
response_time = []
duration = []
prey_size = []
sample_fish = []
correct = 0
wrong = 0
for j in range(0, len(eye_files)): # READ EACH TRIAL

    # Display which file is currently being processed
    print 'PROCESSING: '
    print eye_files[j][0]
    print tail_files[j][0]

    dir_output = str(eye_files[j][1]) # output for saving some info

    # 0s -----------------------------------DETERMINE PREY LOCATION ----------------------------

    prey_loc = prey_location(float(eye_files[j][4]), float(eye_files[j][2][0]), float(eye_files[j][2][1]), 0)
    first_padding = [preypoints[0]] * padding_nonstimulus[0] # visuall stim appeared after 3 s recording
    end_padding = [preypoints[0]] * padding_nonstimulus[1] # visual stim trial ends 3s before the recording
    prey_loc = first_padding + prey_loc + list(reversed(prey_loc)) + prey_loc + end_padding
    print 'PREY GOING FROM ', eye_files[j][2][0], ' TO', eye_files[j][2][1]

    # 0e ---------------------------------------------------------------------------------------

    print 'ONE DOT'

    smaller = None
    preysize = int(eye_files[j][5][eye_files[j][5].find("@") + -1])

    if 'L_' in eye_files[j][5]:
        preyside = 'l'
        prey_loc = list([-int(pr) for pr in prey_loc])
    elif 'R_' in eye_files[j][5]:
        preyside = 'r'
        prey_loc = prey_loc

    else:
        print "ONE DOT: ", eye_files[j][-1]
        continue
    '''
    if 'L_' in eye_files[j][5]: # FIND WHETHER THE PREY IS ON THE LEFT OR RIGHT
        preyside = 'l'
        prey_loc = list([-int(pr) for pr in prey_loc])
    elif 'R_' in eye_files[j][5]:
        preyside = 'r'
        prey_loc = prey_loc
    else:
        print "2 dots: ", eye_files[j][-1]
        continue
    '''
    # 1s -----------------read the eye and tail angles------------------------------------------

    eyes = eye_reader(str(eye_files[j][0]))
    eyes[0]['LeftEye'] = butter_freq_filter(eyes[0]['LeftEye'], cutoff, fs, order) # low pass filter
    eyes[0]['RightEye'] = butter_freq_filter(eyes[0]['RightEye'], cutoff, fs, order)
    tailangle = tail_reader(str(tail_files[j][0]))

    # 1e ----------------------------------------------------------------------------------------

    # 2s -----------------compute the velocity------------------------------------------

    LeftVel = savgol_filter(eyes[0]['LeftEye'], 3, 2, 1, mode = 'nearest')
    RightVel = savgol_filter(eyes[0]['RightEye'], 3, 2, 1, mode = 'nearest')
    eyes_vel = [{'LeftVel': LeftVel, 'RightVel': RightVel}]

    # 2e -------------------------------------------------------------------------------

    # 3s ----------------- DETECT THE TAIL BOUTS AND CORRESPONDING EYE MOVEMENTS -------------------------------

    TailEye = extract_tail_eye_bout2(eyes, eyes_vel, tailangle, Fs, bout_thresh, peakthres, delay) # extract the tail bouts

    if not TailEye: # if there's no bout detected
        print 'WARNING!! TailEye empty: ', eye_files[j]
        continue

    if len(TailEye) == 1: # determine whether it is a single bout trial or a multi-bout
        one_bout += 1
    else:
        multi_bout += 1

    # 3e --------------------------------------------------------------------------------------------------

    # 4s --------------------------------- SET TIME FOR EACH TRIAL TO seconds --------------------------------------------
    time = range(0, len(TailEye[0]['left']))
    time = [t / tfactor for t in time]
    # 4e ---------------------------------------------------------------------------------

    first_pc_bout = 1 # an ID for first prey capture bout
    pc_per_trial = []
    resp_time_per_trial = []

    for i in range(0, len(TailEye)): # READ THROUGH EACH BOUT
        print '============= BOUT #', i, 'of ', TailEye[i]['filename'] , '=================='
        t1 = TailEye[i]['frames'][0]
        t2 = TailEye[i]['frames'][1]

        if first_pc_bout == 0: # Check whether 1st prey capture bout has already been found
            print "First prey capture bout already found"
            continue

        # 5s ---------------------- DETECT PREY CAPTURE -------------------------------------------------------------
        bout_candid = {'frames': TailEye[i]['frames'], 'bout_angles': TailEye[i]['bout_angles'], 'sum_eyeangles': TailEye[i]['sum_eyeangles'],
                       'right_eyeangles': TailEye[i]['right_eyeangles'], 'left_eyeangles': TailEye[i]['left_eyeangles']}

        PC_threshold = {'sensory_delay': sensory_delay,  'bout_angles': filter_avg_tail, 'bout_length': filter_length_bout,
                        'sum_eyeangles': filter_avg_binocular, 'eye_vel': filter_eye_vel,
                        'eye_mag': filter_eye_diverge, 'tail_meanangle': 0.0}

        isPC = ispreycap(bout_candid, PC_threshold, tfactor)
        isPC = all(1 == pc for pc in isPC.values()) # check whether this bout satisfied all the filters set on the PC_threshold
        if not isPC:
            # if isPC is false, this bout is not prey capture, move to next one
            continue
        # 5e ---------------------------------------------------------------------------------------------------------

        print "DETECTED PREY CAPTURE BOUT"
        print eye_files[j][0]
        print tail_files[j][0]
        print 'Tail bout', TailEye[i]['frames']
        print 'Mean binocular angle', np.mean(TailEye[i]['sum_eyeangles'])
        print 'TAIL BOUT FREQ', TailEye[i]['tailfreq']

        if i == 0:
            rtime.append((TailEye[0]['frames'][0] / tfactor))

        tail = np.mean(TailEye[i]['bout_angles'])

        # get the index, and value of the saccade onset based on velocity
        # Get the velocity maximum peak
        # get the index, and value of the saccade onset based on velocity
        # Get the velocity maximum peak
        print TailEye[i]['frames'][0], TailEye[i]['frames'][1]
        r_maxpeak, r_max = max(enumerate(TailEye[i]['right_vel_delay']), key=operator.itemgetter(1))
        l_maxpeak, l_max = max(enumerate(TailEye[i]['left_vel_delay']), key=operator.itemgetter(1))

        # Get the minimum peak between the start of the bout to the peak velocity
        # You only want the minimum peak BEFORE the maximum peak
        # There are cases when the maximum peak is found within the first two points
        # one reason is the starting point of the detected tail bout is greatly delayed compared to
        # its corresponding eye movement

        if r_maxpeak == 0:
            r_min = TailEye[i]['right_vel_delay'][0]
            r_minpeak = 0
        else:
            r_minpeak, r_min = min(enumerate(TailEye[i]['right_vel_delay'][0:r_maxpeak]), key=operator.itemgetter(1))

        if l_maxpeak == 0:
            l_min = TailEye[i]['left_vel_delay'][0]
            l_minpeak = 0
        else:
            l_minpeak, l_min = min(enumerate(TailEye[i]['left_vel_delay'][0:l_maxpeak]), key=operator.itemgetter(1))

        r_sac_on = int(math.ceil((r_maxpeak + r_minpeak) / 2))
        l_sac_on = int(math.ceil((l_maxpeak + l_minpeak) / 2))

        r_sac = TailEye[i]['right_vel_delay'][r_sac_on]
        l_sac = TailEye[i]['left_vel_delay'][l_sac_on]

        right_con = TailEye[i]['right_eyeangles_delay'][r_sac_on]
        left_con = TailEye[i]['left_eyeangles_delay'][l_sac_on]

        right_dir = TailEye[i]['right_eyeangles_delay'][-1] - TailEye[i]['right_eyeangles_delay'][r_sac_on]
        left_dir = TailEye[i]['left_eyeangles_delay'][-1] - TailEye[i]['left_eyeangles_delay'][l_sac_on]

        right_eyebout = np.mean(TailEye[i]['right_eyeangles'])
        left_eyebout = np.mean(TailEye[i]['left_eyeangles'])

        if (TailEye[i]['frames'][0]-padding_nonstimulus[0]) < 0:
            print "NO STIMULUS YET"
            continue

        if (sensory_delay - int(TailEye[i]['frames'][0])) > 10:
            delayed_prey = 0
        else:
            delayed_prey = int(TailEye[i]['frames'][0]) - sensory_delay

        prey_pos1 = prey_loc[delayed_prey]  # delay the reaction time by 30 frames
        print prey_pos1
        # set the onset and contra eye to none, default option
        first_onset = 'none'
        contra_eye = 'none'
        eye_convergence = 0

        if r_max < thresh_saccade_speed and l_max < thresh_saccade_speed:
            print "======WARNING===== THE VELOCITY OF THE ONSET IS LESS THAN 150 deg/s"
            eye_convergence = 1
        # If there's a saccade, use it as the criteria else use eye convergence
        if eye_convergence == 0:
            if r_max < thresh_saccade_speed <= l_max:
                first_onset = 'left'
            elif l_max < thresh_saccade_speed <= r_max:
                first_onset = 'right'
            elif r_sac_on < l_sac_on:
                first_onset = 'right'
            elif l_sac_on < r_sac_on:
                first_onset = 'left'

            eye_onsets.append(first_onset)

        elif eye_convergence == 1:

            if right_eyebout > left_eyebout:
                contra_eye = 'right'
            elif left_eyebout > right_eyebout:
                contra_eye = 'left'

        prey_pos2 = -prey_pos1
        # Assign which prey was selected based on the detected contra-lateral eye
        if tail < -5.0:
            pc_onset = r_sac_on
            print "TAIL WINS", tail, first_onset, contra_eye, right_eyebout - left_eyebout
            if prey_pos1 < 0.0:
                prey_pos = prey_pos1
            elif prey_pos2 < 0.0:
                prey_pos = prey_pos2
        elif tail > 5.0:
            pc_onset = l_sac_on
            print "TAIL WINS", tail, first_onset, contra_eye, right_eyebout - left_eyebout
            if prey_pos1 > 0.0:
                prey_pos = prey_pos1
            elif prey_pos2 > 0.0:
                prey_pos = prey_pos2
        elif first_onset == 'right' or contra_eye == 'right':
            pc_onset = r_sac_on
            print "EYE WINS", tail, first_onset, contra_eye, right_eyebout - left_eyebout
            if prey_pos1 < 0.0:
                prey_pos = prey_pos1
            elif prey_pos2 < 0.0:
                prey_pos = prey_pos2
        elif first_onset == 'left' or contra_eye == 'left':
            pc_onset = l_sac_on
            print "EYE WINS", tail, first_onset, contra_eye, right_eyebout - left_eyebout
            if prey_pos1 > 0.0:
                prey_pos = prey_pos1
            elif prey_pos2 > 0.0:
                prey_pos = prey_pos2
        else:
            prey_pos = 'NONE'

        # Is prey on the right, left or center?
        if prey_pos < 0.0:
            prey_rl = 'left'
        elif prey_pos > 0.0:
            prey_rl = 'right'


        if first_pc_bout == 1:  # first bout
            PC += 1 # is there a prey capture in this trial?
            fish_PC.append(eye_files[j][3] + "_" + eye_files[j][0][-6:-4])
            prey_side.append(preyside)
            response_time.append((TailEye[i]['frames'][0]-padding_nonstimulus[0]) / 0.3)
            duration.append(TailEye[i]['frames'])
            prey_size.append(preysize)
            sample_fish.append(eye_files[j][3])
            first_pc_bout = 0 # set to 0 to notify that first PC bout is found

            if preyside == 'l' and prey_rl == 'left':
                correct += 1
            elif preyside == 'r' and prey_rl == 'right':
                correct += 1
            else:
                print 'WRONG--------------', preyside, prey_rl
                wrong += 1

    print TailEye[0]['filename']
    print 'done'

print set(sample_fish)
print len(set(sample_fish))
print PC
print np.mean(response_time), np.std(response_time)

print "D", correct, wrong
'''
results = izip_longest(fish_PC, prey_side, prey_size, response_time, duration)
header = izip_longest(['Fish'], ['Prey Side'], ['Prey Size'], ['Response time'], ['Duration'])
with open(maindir + "Summary of 1st PC_1dot_v4" + '.csv', 'wb') as myFile:
    # with open(dir_output + 'Velocity_Acceleration_' + filename + '.csv', 'wb') as myFile:
    wr = csv.writer(myFile, delimiter=',')
    for head in header:
        wr.writerow(head)
    for rows in results:
        wr.writerow(rows)
'''
'''
binwidth = 250
n1, bins1, patches1 = plt.hist(response_time, color=[0, 0, 1],edgecolor ='black',
                           bins=np.arange(min(response_time), max(response_time) + binwidth, binwidth))
response_time = [r/1000.0 for r in response_time]

sns.distplot(response_time,rug=True,hist=False)
plt.show()
'''

