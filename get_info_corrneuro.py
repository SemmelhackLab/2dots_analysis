import pickle


def get_info(pickle_file):

    pickle_in = open(pickle_file,"r")

    info = pickle.load(pickle_in)

    pickle_in.close()

    return info

