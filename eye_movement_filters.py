from scipy.signal import butter, filtfilt


def butter_filter(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_freq_filter(data, cutoff, fs, order=5):
    b, a = butter_filter(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def dynamic_filter(eye, sac, sac_cutoff, nonsac, nonsac_cutoff, order, fs):
    # eye = the eye angle data
    # sac = the index of saccade points
    # sac_cutoff = frequency cut off for saccade
    # nonsac = index of non saccade points
    # nonsac_cutoff = freq cut off for non saccade
    # order = order of filter
    # fs = sampling frequency

    sac_filtered = butter_freq_filter(eye, sac_cutoff, fs, order)  # high-pass filtered
    nonsac_filtered = butter_freq_filter(eye, nonsac_cutoff, fs, order)  # low-pass filtered

    eye_filtered = [None] * len(eye) # declare the eye angles array

    for ind, val in enumerate(sac):
        eye_filtered[val] = sac_filtered[val]

    for ind, val in enumerate(nonsac):
        eye_filtered[val] = nonsac_filtered[val]

    return eye_filtered
