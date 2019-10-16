import numpy as np
import operator


def ispreycap(bout_candid, PC_threshold, tfactor):
    '''
    THIS FUNCTION CONTAINS SET OF FILTERS TO EXCLUDE NON PREY CAPTURE BOUTS
    :param bout_candid:
    :param PC_threshold:
    :param tfactor:
    :return:
    '''

    # Assign an ID for each filter for prey capture
    PC_filter = {'sensory_delay': 1,'tail_meanangle': 1, 'sum_eyeangles': 1,
                 'bout_length': 1, 'eye_vel': 1, 'eye_mag': 1}
    t1 = bout_candid['frames'][0]
    t2 = bout_candid['frames'][1]

    # 0s -------------------- FILTER BY SENSORY DELAY -------------------
    if t1 < PC_threshold['sensory_delay']:
        print "1. Too early bout"
        PC_filter['sensory_delay'] = 0
    # 0e -----------------------------------------------------------------

    # 1s ------------------- FILTER BY TAIL ANGLE ------------------------
    tail = np.mean(bout_candid['bout_angles']) # get the tail angle
    if np.abs(tail) < PC_threshold['tail_meanangle']:
        print '2. Not strong enough tail movement'
        PC_filter['tail_meanangle'] = 0

    # 1e ------------------------------------------------------------------

    # 2s ---------------------- FILTER BASED ON BINOCULAR CONVERGENCE -----------------------------------------------
    if np.mean(bout_candid['sum_eyeangles']) <  PC_threshold['sum_eyeangles']:
        print "3. Eyes are not converge enough"
        PC_filter['sum_eyeangles'] = 0
    # 2e -----------------------------------------------------------------------------------------------------------

    # 3s ---------------------------- FILTER BASED ON LENGTH OF BOUT ----------------------------------------------
    if int(bout_candid['frames'][1] - bout_candid['frames'][0]) < PC_threshold['bout_length']:
        print '4. TOO SHORT BOUT'
        PC_filter['bout_length'] = 0
    # 3e ---------------------------------------------------------------------------------------------------------

    # 4s ---------------------------- FILTER BASED ON EYE VELOCITY ----------------------------------------------

    right_v = (bout_candid['right_eyeangles'][-1] - bout_candid['right_eyeangles'][0]) / (
            len(bout_candid['right_eyeangles']) / tfactor)
    left_v = (bout_candid['left_eyeangles'][-1] - bout_candid['left_eyeangles'][0]) / (
            len(bout_candid['left_eyeangles']) / tfactor)

    if right_v <= PC_threshold['eye_vel'] or left_v <= PC_threshold['eye_vel']:
        print '5. Diverging eyes'
        PC_filter['eye_vel'] = 0
    # 4e ---------------------------------------------------------------------------------------------------------

    # 5s ---------------------------- FILTER BASED ON EYE VELOCITY ----------------------------------------------
    right_mag = (bout_candid['right_eyeangles'][-1] - bout_candid['right_eyeangles'][0])
    left_mag = (bout_candid['left_eyeangles'][-1] - bout_candid['left_eyeangles'][0])

    if right_mag <  PC_threshold['eye_mag'] or left_mag < PC_threshold['eye_mag']:
        print '6. Big divergence of eye'
        PC_filter['eye_mag'] = 0
    # 5e ---------------------------------------------------------------------------------------------------------

    return PC_filter

'''
def find_pc_onset(right_vel, left_vel):

    # Get the minimum peak between the start of the bout to the peak velocity
    # You only want the minimum peak BEFORE the maximum peak
    # There are cases when the maximum peak is found within the first two points
    # one reason is the starting point of the detected tail bout is greatly delayed compared to
    # its corresponding eye movement

    r_maxpeak, r_max = max(enumerate(right_vel), key=operator.itemgetter(1))
    l_maxpeak, l_max = max(enumerate(left_vel), key=operator.itemgetter(1))

    if r_maxpeak == 0:
        r_min = right_vel[0]
        r_minpeak = 0
    else:
        r_minpeak, r_min = min(enumerate(TailEye[i]['right_vel_delay'][0:r_maxpeak]), key=operator.itemgetter(1))

    if l_maxpeak == 0:
        l_min = TailEye[i]['left_vel_delay'][0]
        l_minpeak = 0
    else:
        l_minpeak, l_min = min(enumerate(TailEye[i]['left_vel_delay'][0:l_maxpeak]), key=operator.itemgetter(1))

    if (TailEye[i]['frames'][0] - padding_nonstimulus[0]) < 0:
        print "NO STIMULUS YET"
        continue

    return
'''