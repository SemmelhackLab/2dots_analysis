# THIS CODE WILL DETERMINE THE SELECTED PREY AND ITS LOCATION
# BASED ON THE ONSET OF THE CONTRA-LATERAL EYE CONVERGENCE
# The onset is defined as the point between the peak of eye velocity
# and the start of the tail bout.
# Last update: 23 AUG 2018, Ivan

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
#plt.style.use(['dark_background'])
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
# ============================ MACRO FOR PREY SELECTION ANALYSIS =======================================================
# PARAMETERS FOR EXTRACTING THE BOUTS, AND PEAKS
Fs = 300
bout_thresh = 0.40  # 0.00 - 1.00 threshold value for extracting bouts, higher more bouts
tfactor = 0.3  # convert frames to ms
peakthres = 4  # 0.00 - 20.00 lower vale more peaks, for calculating tail beat frequency
speed = 40.0 # prey speed
filter_avg_binocular = 17.0 # threshold for eye angle binocular convergence
filter_length_bout = 10 # length of bout based on frame number
filter_eye_vel = -0.1 # diverging eye velocity in degrees/ms
filter_eye_diverge = -1.0 # eye divergence by degrees
thresh_saccade_speed = 0.2 # the onset of eye movement should be greater than saccade threshold
sensory_delay = 30 # delay of sensory transformation in number of frames
preypoints = [10, 70] # direction of prey in visual angle
delay = 50 # window delay for finding the peak
padding_nonstimulus = [1800,1650] # to add the frame for waiting time in prey location
prey_size = [2,3] # prey sizes to compare in ascending order

# =================== LOW PASS FILTER ========================
order = 3
fs = 300.0  # sample rate, Hz
cutoff = 10  # desired cutoff frequency of the filter, Hz

# ========================= GENERATE ALL THE FILES =================================================================
# Directories for 80 deg visual angle
maindir = 'E:\\new folder\\002 ANALYSIS\\2dots 2p\\2p_left_right_analysis\\'
dir = 'E:\\new folder\\002 ANALYSIS\\2dots 2p\\2p_left_right_analysis\\' #2dotsdiff_'+ str(prey_size[0]) + 'deg_' + str(prey_size[1]) + 'deg\\'
#maindir = 'D:\\Semmelhack lab\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\'
#dir = 'D:\\Semmelhack lab\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\' #2dotsdiff_'+ str(prey_size[0]) + 'deg_' + str(prey_size[1]) + 'deg\\'
dir_output = dir

eye_files = []
tail_files = []

fishIDs = list()
fishIDs_dir = list()

fishIDs += [fish for fish in os.listdir(dir) if ".csv" not in fish]  # store the fish filenames

fishIDs_dir += [{str(fish): [0, 0], str(fish) +'_bigorsmall': [], str(fish) +'_resptime': [] } for fish in fishIDs] # for storing the right or left information for each fish
neyes = [] # to check if there are missing files
ntail = [] # to check if there are missing files
for fish in fishIDs:

    eyefile = dir + fish + '\\results\\'
    tailfile = dir + fish + '\\tail\\tail results\\'

    for file in os.listdir(eyefile):
        if file.endswith(".csv"):
            eye_files.append([os.path.join(eyefile, file), dir_output, preypoints, fish, speed, file])
    for file in os.listdir(tailfile):
        if file.endswith(".csv"):
            tail_files.append([os.path.join(tailfile, file), dir_output, preypoints, fish, speed, file])

print len(eye_files), len(tail_files)
print len(fishIDs)

# print fishIDs

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
small_pc = 0
big_pc = 0
reyepos_rightprey = []
leyepos_rightprey = []
reyepos_leftprey = []
leyepos_leftprey = []
eyediff_leftpc = []
eyediff_rightpc = []
eye_b4pc = []
fish_PC = []
prey_side = []
response_time = []
duration = []
size_prey = []
eye_leftpc = []
eye_rightpc = []
curr_fish = -1
sample_fish = []
for j in range(0, len(eye_files)):
    curr_fish += 1
    # ============ PREY LOCATION ==================
    dir_output = str(eye_files[j][1])
    # prey = prey_location_osci(40.0,70.0,-50.0,-10.0,'ccw',5) # with oscillation (speed,p1,p2,p3,dir,times)
    prey_loc = prey_location(float(eye_files[j][4]), float(eye_files[j][2][0]), float(eye_files[j][2][1]), 0)  # without oscillation, (speed,p1,p2,times)
    first_padding = [preypoints[0]] * padding_nonstimulus[0]
    end_padding = [preypoints[0]] * padding_nonstimulus[1]
    prey_loc = first_padding + prey_loc + list(reversed(prey_loc)) + prey_loc + end_padding
    print 'PREY GOING FROM ', eye_files[j][2][0], ' TO', eye_files[j][2][1]

    if 'L-' in eye_files[j][5]:
        smaller = 'l'
        prey_size[0] = int(eye_files[j][5][eye_files[j][5].find("L") + 2])
        prey_size[1] = int(eye_files[j][5][eye_files[j][5].find("R") + 2])

    elif 'R-' + str(prey_size[0]) + 'deg' in eye_files[j][5]:
        smaller = 'r'
        prey_size[0] = int(eye_files[j][5][eye_files[j][5].find("R") + 2])
        prey_size[1] = int(eye_files[j][5][eye_files[j][5].find("L") + 2])
    else:
        print "ONE DOT: ", eye_files[j][-1]
        continue

    '''
    # WHERE IS THE SMALLER PREY?
    if 'L-' + str(prey_size[0]) +'deg' in eye_files[j][5]:
        smaller = 'l'
    elif 'R-' + str(prey_size[0]) + 'deg' in eye_files[j][5]:
        smaller = 'r'
    else:
        print "ONE DOT: ", eye_files[j][-1]
        continue
    '''
    print 'PROCESSING: '
    print eye_files[j][0]
    print tail_files[j][0]

    # read the left and right eye tracks
    eyes = eye_reader(str(eye_files[j][0]))
    tailangle = tail_reader(str(tail_files[j][0]))
    eyes[0]['LeftEye'] = butter_freq_filter(eyes[0]['LeftEye'], cutoff, fs, order)
    eyes[0]['RightEye'] = butter_freq_filter(eyes[0]['RightEye'], cutoff, fs, order)

    # compute the velocity
    LeftVel = savgol_filter(eyes[0]['LeftEye'], 3, 2, 1, mode = 'nearest')
    RightVel = savgol_filter(eyes[0]['RightEye'], 3, 2, 1, mode = 'nearest')
    eyes_vel = [{'LeftVel': LeftVel, 'RightVel': RightVel}]

    # DETECT THE TAIL BOUTS AND CORRESPONDING EYE MOVEMENTS
    TailEye = extract_tail_eye_bout2(eyes, eyes_vel, tailangle, Fs, bout_thresh, peakthres, delay)

    if not TailEye:
        print 'WARNING!! TailEye empty: ', eye_files[j]
        continue

    if len(TailEye) == 1:
        one_bout += 1
    else:
        multi_bout += 1

    time = range(0, len(TailEye[0]['left']))
    time = [t / tfactor for t in time]
    first_pc_bout = 1 # an ID for first prey capture bout
    pc_per_trial = []
    resp_time_per_trial = []

    for i in range(0, len(TailEye)):
        t1 = TailEye[i]['frames'][0]
        t2 = TailEye[i]['frames'][1]
        if t1 < sensory_delay:
            print "Too early bout"
            continue
        if first_pc_bout == 0:
            print "First prey capture bout already found"
            continue

        print '============= BOUT #', i, '=================='
        if i == 0:
            rtime.append((TailEye[0]['frames'][0] / tfactor))
        print eye_files[j][0]
        print tail_files[j][0]
        print 'Tail bout', TailEye[i]['frames']
        print 'Mean binocular angle', np.mean(TailEye[i]['sum_eyeangles'])
        print 'TAIL BOUT FREQ', TailEye[i]['tailfreq']
        # FILTER BASED ON BINOCULAR CONVERGENCE
        if np.mean(TailEye[i]['sum_eyeangles']) < filter_avg_binocular:  # [(mid_list  - tenth):(mid_list + tenth)]) < 30:
            print 'Not converge enough, bout #: ', i, 'Sample: ', TailEye[i]['filename']
            continue

        # FILTER BASED ON LENGTH OF BOUT
        if int(TailEye[i]['frames'][1] - TailEye[i]['frames'][0]) < filter_length_bout:
            print 'TOO SHORT BOUT'
            continue

        right_eyebout = np.mean(TailEye[i]['right_eyeangles'])
        left_eyebout = np.mean(TailEye[i]['left_eyeangles'])

        right_v = (TailEye[i]['right_eyeangles'][-1] - TailEye[i]['right_eyeangles'][0]) / (
                len(TailEye[i]['right_eyeangles']) / tfactor)
        left_v = (TailEye[i]['left_eyeangles'][-1] - TailEye[i]['left_eyeangles'][0]) / (
                len(TailEye[i]['left_eyeangles']) / tfactor)

        # FILTER BASED ON THE VELOCITY
        if right_v <= filter_eye_vel or left_v <= filter_eye_vel:
            print 'Diverging eyes'
            continue

        right_v = (TailEye[i]['right_eyeangles'][-1] - TailEye[i]['right_eyeangles'][0])# / (
                #len(TailEye[i]['right_eyeangles']) / tfactor)
        left_v = (TailEye[i]['left_eyeangles'][-1] - TailEye[i]['left_eyeangles'][0])# / (
                #len(TailEye[i]['left_eyeangles']) / tfactor)

        # FILTER BASED ON THE MAGNITUDE OF EYE DIVERGENCE
        if right_v < filter_eye_diverge or left_v < filter_eye_diverge:
            print 'Big divergence of eye'
            continue

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

        r_sac_on = int(math.ceil((r_maxpeak + r_minpeak)/2))
        l_sac_on = int(math.ceil((l_maxpeak + l_minpeak)/2))

        r_sac = TailEye[i]['right_vel_delay'][r_sac_on]
        l_sac = TailEye[i]['left_vel_delay'][l_sac_on]

        right_con = TailEye[i]['right_eyeangles_delay'][r_sac_on]
        left_con = TailEye[i]['left_eyeangles_delay'][l_sac_on]

        right_dir = TailEye[i]['right_eyeangles_delay'][-1] - TailEye[i]['right_eyeangles_delay'][r_sac_on]
        left_dir = TailEye[i]['left_eyeangles_delay'][-1] - TailEye[i]['left_eyeangles_delay'][l_sac_on]

        # Eye orientation before prey capture

        tail = np.mean(TailEye[i]['bout_angles'])

        EyeDiff = right_eyebout - left_eyebout
        VelDiff = right_v - left_v
        converge_thres = 2.0

        if EyeDiff > converge_thres:
            print 'Right eye converge more: ', 'R:', right_con, 'L:', left_con, \
                'Avg Velocities', right_v, left_v
        elif EyeDiff < -converge_thres:
            print 'Left eye converge more: ', 'L:', left_con, 'R:', right_con, \
                'Avg Velocities', left_v, right_v
        else:
            print 'Left and right eyes are equally converge: ', 'R:', right_con, 'L:', left_con, \
                'Avg Velocities', right_v, left_v

        if (sensory_delay - int(TailEye[i]['frames'][0])) > 10:
            delayed_prey = 0
        else:
            delayed_prey = int(TailEye[i]['frames'][0]) - sensory_delay

        if delayed_prey >= len(prey_loc):
            print '========= WARNING ========= No prey location info available for the given bout frame. ' \
                  'The frame location of the bout is greater than the length of frames where prey is available'
            continue

        # print 'delay', delayed_prey, len(prey_loc)
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

        prey_positions1.append(prey_pos)
        eye_diffs.append(EyeDiff)
        tails.append(tail)

        print 'Prey position', prey_pos
        print 'Eye Onset:', first_onset, 'R:', r_sac_on, 'L:', l_sac_on
        print 'Right saccade: ', r_maxpeak, r_max
        print 'Left saccade: ', l_maxpeak, l_max
        print 'right', right_con, right_dir
        print 'left', left_con, left_dir
        print 'tail angle', tail
        print 'tail beat frequency', TailEye[i]['tailfreq']

        stimsd = 1
        fi = fishIDs.index(eye_files[j][3])  # index of the current fish among all the fish samples in fishIDs
        if (TailEye[i]['frames'][0]-padding_nonstimulus[0]) < 0:
            print "NO STIMULUS YET"
            stimsd = 2
            continue

        pc_onset = pc_onset + t1 - delay
        lefteye_b4pc = np.mean(TailEye[i]['left'][pc_onset - 150: pc_onset])
        righteye_b4pc = np.mean(TailEye[i]['right'][pc_onset - 150: pc_onset])
        eye_b4pc.append([lefteye_b4pc, righteye_b4pc])

        '''
        if prey_rl == 'right':
            eye_rightpc.append(np.subtract(TailEye[i]['right'][pc_onset: pc_onset + 150],TailEye[i]['left'][pc_onset: pc_onset + 150]))
        elif prey_rl == 'left':
            eye_leftpc.append(np.subtract(TailEye[i]['right'][pc_onset: pc_onset + 150],TailEye[i]['left'][pc_onset: pc_onset + 150]))
        '''
        if prey_rl == 'right':
            eye_rightpc.append(np.subtract(TailEye[i]['left'][pc_onset-300: pc_onset+300],TailEye[i]['right'][pc_onset-300: pc_onset+300]))
        elif prey_rl == 'left':
            eye_leftpc.append(np.subtract(TailEye[i]['left'][pc_onset-300: pc_onset+300],TailEye[i]['right'][pc_onset-300: pc_onset+300]))

        #eye_b4pc.append(lefteye_b4pc)
        #eye_b4pc.append(np.mean(np.subtract(lefteye_b4pc, righteye_b4pc)))
        print "ONSET%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print TailEye[i]['frames'][0] - pc_onset

        #print stimsd
        if first_pc_bout == 1: # first bout
            if prey_pos > 0.0:
                if smaller == 'r':
                    small_pc += 1
                    pc_per_trial.append(1)
                    size_prey.append(prey_size[0])
                else:
                    big_pc += 1
                    pc_per_trial.append(0)
                    size_prey.append(prey_size[1])

                resp_time_per_trial.append(TailEye[i]['frames'][0]/0.3)
                right_prey += 1
                reyepos_rightprey.append(righteye_b4pc)
                leyepos_rightprey.append(lefteye_b4pc)
                eyediff_rightpc.append([(righteye_b4pc - lefteye_b4pc), prey_pos])
                fishIDs_dir[fi][fishIDs[fi]][1] += 1

                fish_PC.append(eye_files[j][3] + "_" + eye_files[j][0][-6:-4])
                prey_side.append(prey_pos)
                response_time.append((TailEye[i]['frames'][0] - padding_nonstimulus[0]) / 0.3)
                duration.append([pc_onset, TailEye[i]['frames'][1]])
                sample_fish.append(eye_files[j][3])

                first_pc_bout = 0 # turn off the first PC ID
            elif prey_pos < 0.0:
                if smaller == 'l':
                    small_pc += 1
                    pc_per_trial.append(1)
                    size_prey.append(prey_size[0])
                else:
                    big_pc += 1
                    pc_per_trial.append(0)
                    size_prey.append(prey_size[1])

                resp_time_per_trial.append(TailEye[i]['frames'][0]/0.3)
                reyepos_leftprey.append(righteye_b4pc)
                leyepos_leftprey.append(lefteye_b4pc)
                eyediff_leftpc.append([(righteye_b4pc - lefteye_b4pc), prey_pos])
                left_prey += 1
                fishIDs_dir[fi][fishIDs[fi]][0] += 1

                fish_PC.append(eye_files[j][3] + "_" + eye_files[j][0][-6:-4])
                #fish_PC.append(eye_files[j][0][-47:-38] + "_" + eye_files[j][0][-6:-4])
                prey_side.append(prey_pos)
                response_time.append((TailEye[i]['frames'][0] - padding_nonstimulus[0]) / 0.3)
                duration.append([pc_onset, TailEye[i]['frames'][1]])
                sample_fish.append(eye_files[j][3])


                first_pc_bout = 0 # turn off the first PC ID
        elif first_pc_bout == 0:
            if prey_pos > 0.0:
                resp_time_per_trial.append(TailEye[i]['frames'][0]/0.3)
                if smaller == 'r':
                    pc_per_trial.append(1)
                else:
                    pc_per_trial.append(0)

            elif prey_pos < 0.0:
                resp_time_per_trial.append(TailEye[i]['frames'][0]/0.3)
                if smaller == 'l':
                    pc_per_trial.append(1)
                else:
                    pc_per_trial.append(0)

        nbouts += 1
    if first_pc_bout == 0:
        fishIDs_dir[fi][fishIDs[fi] + '_resptime'].append(resp_time_per_trial)
    if pc_per_trial:
        fishIDs_dir[fi][fishIDs[fi] + '_bigorsmall'].append(pc_per_trial)
    print TailEye[0]['filename']
    print 'done'
'''
for i in eye_b4pc:
    plt.scatter(i[0], i[1], c = [1,0,1])
plt.xlim([-8,13])
plt.ylim([-8,13])

plt.show()
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(10, 8))

time = [(t/300.0)-1.0 for t in range(len(eye_rightpc[0]))]

for i in eye_rightpc:
    plt.plot(time, i, c = [1,0,1])

for i in eye_leftpc:
    plt.plot(time, i, c = [0,1,1])
#plt.fill_between([-1.0, 0], -15, 15, facecolor=[0.5,0.5,0.5], alpha=0.5, interpolate=True)
ax.axvline(0, c=[0.5,0.5,0.5], ls='--')
ax.set_ylabel("Eye difference ($^\circ$)", fontsize=30)
ax.set_xlabel("Time (s)", fontsize=30)
ax.set_ylim([-37,37])
ax.set_yticks([-20, -10, 10, 20])
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
fig.savefig('C:\\Users\\Semmelhack Lab\\Documents\\Presentations 2\\NOV 2019 ATTENTION GRANT\\eyeb4PC.png')
plt.show()
'''

print set(sample_fish)
print len(set(sample_fish))

print small_pc, big_pc
print np.mean(response_time), np.std(response_time)

results = izip_longest(fish_PC, prey_side, size_prey, response_time, duration)
header = izip_longest(['Fish'], ['Prey Side'], ['Prey Size'], ['Response time'], ['Duration'])
with open(maindir + "Summary of 1st PC_2dotsv4.csv", 'wb') as myFile:
     #with open(dir_output + 'Velocity_Acceleration_' + filename + '.csv', 'wb') as myFile:
    wr = csv.writer(myFile, delimiter=',')
    for head in header:
        wr.writerow(head)
    for rows in results:
        wr.writerow(rows)
