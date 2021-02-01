#! /bin/python3

import json
import math
import random
import copy
import matplotlib.pyplot as plt
import numpy as np
import argparse

from range import env, get_loss_for_env

### ARGS ###
ARGS={}

############## CHANGE THESE PARAMETERS TO CHANGE THE SIMULATION PARAMETERS ##################

AREA_LENGTH = AREA_BREADTH = 100
LAMBDA_A = 0.01                                      # AP density (unit: aps/meter2)
CHANNELS_AVAILABLE = [0]                             # list of all used for simulation for 2.4GHz channels

F_c = 2.44                                             # Center frequency  
P_t = {20:20,40:23}#100            # if dBm, should be 20      # transmit power at the ap in dBm with 20MHz bandwidth
P_t40 = 23                                          # transmit power at the ap in dBm with 40MHz bandwidth
NOISE_PSD = -174                                    # noise power spectal density in dBm/Hz
CHANNEL_BANDWIDTH = [20]                              # channel bandwidth
SPECTRAL_EFFICIENCY = {20:3.61,40:3.75}                          # maximum link spectral efficiency 20MHz
#SPECTRAL_EFFICIENCY = 3.61                          # maximum link spectral efficiency 20MHz
SPECTRAL_EFFICIENCY40 = 3.75                        # maximum link spectral efficiency 40MHz
BANDWIDTH_EFFICIENCY_1 = 1                          # bandwidth efficiency coefficient 1
BANDWIDTH_EFFICIENCY_2 = 1                          # bandwidth efficiency coefficient 2
# constants
MAX_RANGE = {
    'home': env(F_c, "home"),
    'mall': env(F_c, "mall"),
    'los': env(F_c, "los")
}
ENVIRONMENT = 'mall'

############################################################################################


# set matplotlib parameters 
plt.rcParams["figure.figsize"] = (10, 10)

CHANNEL_COLORS = [                                  # assign different colors for differnt channel
    '#FF0000',                                   # red 
    '#0000FF',                                   # blue 
    '#00FF00'                                   # green
    ]

# inline functions
DISTANCE = lambda x1, y1, x2, y2: math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
mW_TO_dBm = lambda p: 10 * math.log(p, 10)
dBm_TO_mW = lambda p: 10**(p/10)


def grid_arrange_aps(area_length=50, area_breath=50, lambda_a=1.85):
    """ 
    Arrange APs in a grid order
    This function is now limited to only arraging AP for a square environment with square ap service area.
    :TODO extend this function to work with rectangle environment aswell
    :TODO extend this function to work with rectangle AP service area.

    This function uses the AP density (lambda_a) to compute the number of APs required on each side
    """

    env_area = AREA_LENGTH * AREA_BREADTH
    num_of_aps = env_area * lambda_a                        # total number of AP required for the complete deployment
    aps_per_side = int(math.ceil(math.sqrt(num_of_aps)))    # since we are using a grid deployment calculate the perfect square of
    #dist = round(AREA_LENGTH / (2 * aps_per_side), 2)       # half the distance between APs
    dist = AREA_LENGTH / (2 * aps_per_side)       # half the distance between APs

    aps, ap_matrix = {}, []
    for i in range(aps_per_side):
        mat_tmp = list()                                                            # use this matix to get a APs matrix, this matrix could be usely later
        for j in range(aps_per_side):
            ap_id = aps_per_side * i + j                                            # AP id in sequencial order
            #y, x = round(dist + 2 * i * dist, 2), round(dist + 2 * j * dist, 2)     # x and y are in meters and the coordiante center is placed on top-left corner
            y, x = dist + 2 * i * dist, dist + 2 * j * dist     # x and y are in meters and the coordiante center is placed on top-left corner
            aps.update({ap_id: {'x': x, 'y': y }})
            mat_tmp.append(ap_id)
        ap_matrix.append(mat_tmp)

    return aps, ap_matrix


def randomly_arrange_aps(area_width=50, area_height=50, lambda_a=1.85, minimum_distance=1):
    """
    Area measured in m. minimum_distance of 1 meters to simulate sane worst case inter-apartment placement.
    Returns a list with points that are inside the env_area and follows minimum_distance.

    An idea for how to make this functions time complexity linear.
    - Divide the simulated area into squares that has the width of minimum_distance
    - Save this structure in a dictionary/hashmap
    - Every time a new point is to be tested
        - find which square the new point should be in
        - look in the 8 squares surounding the points square
        - look if any of the points in the surounding square are too close
        - if non are, add the point to it's square and to the list of APs
        - if there's a conflict try with a new point

    The challenge is to translate the coordiantes to the squares. The rest shouldn't be a problem.
    """

    env_area = area_height * area_width
    num_of_aps = env_area * lambda_a
    # Convert to cm from m.
    width = area_width*100
    height = area_height*100
    min = minimum_distance*100

    #print("Starting random AP assignment.")

    aps = []
    while num_of_aps > 0:
        ap = {'x': random.randint(0, width), 'y': random.randint(0, height)}
        too_close = False
        for a in aps:
            if math.sqrt(pow(ap['x'] - a['x'], 2) + pow(ap['y'] - a['y'], 2)) < min:
                too_close = True
                #print("Too close, retrying.")
                break
        if not too_close:
            aps.append(ap)
            num_of_aps -= 1
            #print("Assigned AP, %d of APs to go." % num_of_aps)

    # Convert back from cm to m.
    for a in aps:
        a['x'] = round(a['x'] / 100, 1)
        a['y'] = round(a['y'] / 100, 1)

    #print("Random AP assignment done.")

    return aps


def assign_channels_for_raa(aps):
    """
    Assign channels to output from randomly_assign_aps
    list of dictionaries -> dictionary with id as key
    """

    # Create list of channels to assign to aps
    channels = []
    channels_len = len(aps)
    segment_len = math.floor(channels_len/len(CHANNELS_AVAILABLE))
    for n in CHANNELS_AVAILABLE:
        for _ in range(segment_len):
            channels.append(n)

    diff = channels_len - len(channels)
    if diff != 0:
        for n in range(diff):
            channels.append(CHANNELS_AVAILABLE[n])      # diff should always be be less than len of CHANNELS_AVAILABLE

    random.shuffle(channels)

    # Create new dictionary and assign channels to aps
    aps_with_channel = {}
    ap_id = 0
    for ap in aps:
        ap['operating_channel'] = channels[ap_id]
        aps_with_channel[ap_id] = ap
        ap_id += 1

    return aps_with_channel


def assign_channels(aps):
    """
        Function to randomly assign channels to the APs. The strategy used here is to randomly assign channel
        to APs proportionally based on available channels
    """

    # shuffle aps to randomize channel assignment
    ap_ids = list(aps.keys())
    random.shuffle(ap_ids)
    aps_to_channel = math.ceil(len(ap_ids) / len(CHANNELS_AVAILABLE))
    # assign channels here
    channels = []
    for ch in CHANNELS_AVAILABLE:
        n = min(aps_to_channel, len(ap_ids)-len(channels))
        channels.extend([ch] * n)
    # update aps dictionary
    aps_tmp = copy.deepcopy(aps)
    for i, ap in enumerate(ap_ids):
        aps_tmp[ap].update({'operating_channel': channels[i]})

    return aps_tmp


def drop_users(aps, d_cs):
    """
    Randomly drop users on the given environment.
    Drop one user within the d_cs of all the aps 
    """

    users = dict()
    for i, ap in enumerate(aps):
        user_x = max(0,min(AREA_LENGTH,aps[ap]['x'] + random.uniform(-15,15)))                  # randomly sample user within an area of 30x30 m^2 originating from our ap, and inside our area
        user_y = max(0,min(AREA_LENGTH,aps[ap]['y'] + random.uniform(-15,15)))
        #user_x = random.uniform(aps[ap]['x'] - d_cs, aps[ap]['x'] + d_cs)/d_cs                 # randomly sample user coordinates
        #user_y = random.uniform(aps[ap]['y'] - d_cs, aps[ap]['y'] + d_cs)/d_cs
        users.update({i: {'x': user_x, 'y': user_y, 'associated_ap': ap} })

    return users


def get_transmitting_aps(aps, selected_channel, d_cs):
    """
    Get a list of non-contenting APs using the SSI process
    """

    realization = []
    # start with a random ap with the selected channel
    ap_ids = list(aps.keys())
    random.shuffle(ap_ids)
    start_ap = [ap for ap in ap_ids if aps[ap]['operating_channel'] == selected_channel][-1]
    realization.append(start_ap)

    # get a list of non-contenting aps
    for ap in aps:
        if ap != start_ap:
            is_contenting = False
            if aps[ap]['operating_channel'] == aps[start_ap]['operating_channel']:
                for rel_ap in realization:
                    d_ij = DISTANCE(aps[rel_ap]['x'], aps[rel_ap]['y'], aps[ap]['x'], aps[ap]['y'])
                    if d_ij <= d_cs:
                        is_contenting = True
                        break
                if not is_contenting:
                    realization.append(ap)

    return realization


def calculate_datarate(transmitter, all_transmitting_aps, environment, aps, users):
    receiver = [user for user in users if users[user]['associated_ap'] == transmitter][-1]
    
    bw = CHANNEL_BANDWIDTH[aps[transmitter]['operating_channel']]
    pt = P_t[bw]
    se = SPECTRAL_EFFICIENCY[bw]
    # signal power from selected transmitter and receiver pair
    dist_gii = DISTANCE(users[receiver]['x'], users[receiver]['y'], aps[transmitter]['x'], aps[transmitter]['y'])       # in meters
    loss_gii = get_loss_for_env(dist_gii, F_c, environment)                                                                  # in dBm 
    p_ii = dBm_TO_mW(pt - loss_gii)                                                                                    # in mW 

    # compute the signal power that reaches the received from all the other co transmitting APs
    p_ix = 0
    for ap_t in all_transmitting_aps:
        dist_gix = DISTANCE(users[receiver]['x'], users[receiver]['y'], aps[ap_t]['x'], aps[ap_t]['y'])                 # in meters
        loss_gix = get_loss_for_env(dist_gix, F_c, environment)                                                              # in dBm  
        p_ix = p_ix + dBm_TO_mW(pt - loss_gix)                                                                         # in mW 
    
    np.seterr(divide='ignore')
    # # sinr
    ###### check if this formula correct. 
    ###### here dBm_TO_mW(NOISE_PSD * CHANNEL_BANDWIDTH) is going to zero, check why
    SINR = np.divide(p_ii , (p_ix * dBm_TO_mW(NOISE_PSD) * bw))                                                     # in mW          

    # compute data rate
    data_rate = BANDWIDTH_EFFICIENCY_1 * math.log(1 + BANDWIDTH_EFFICIENCY_2 * SINR, 2) 
    data_rate = bw  * min(data_rate, se)                                                 # in Mbps

    return data_rate

def draw_transmitting(aps, transmitting_aps, users):
    """
    Plot transmitting ap
    """
    ap_xpoints,ap_ypoints,ap_color = [],[],[]
    i=0
    sending=0
    transmitting=()
    other=[]
    for ap in aps.values():
        ap_xpoints.append(ap['x'])
        ap_ypoints.append(ap['y'])
        if i==transmitting_aps[0]:
            ap_color.append('#00FF00')
            sending=i
            transmitting=(ap['x'],ap['y'])
        elif i in transmitting_aps:
            ap_color.append('#0000FF')
            other.append((ap['x'],ap['y']))
        else:
            ap_color.append('grey')
        i+=1

    # (x, y) for all the ordinates of users (x and y translate to distances in meters)
    user_xpoints, user_ypoints, user_color = [], [], []
    for user in users.values():
        if user['associated_ap'] == sending:
            user_xpoints.append(user['x'])
            user_ypoints.append(user['y'])
            user_color.append('#FF00FF')
    
    
    # plot
    cir=[]
    fig, ax = plt.subplots()
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.scatter(ap_xpoints, ap_ypoints, marker='o', s=30, c=ap_color)
    ax.scatter(user_xpoints, user_ypoints, marker='x', s=10, c=user_color)
    cir = plt.Circle(transmitting,MAX_RANGE[ENVIRONMENT],color='grey',fill=False,linestyle='--')
    #ax.set_aspect('equal', adjustable='datalim')
    ax.add_patch(cir)
    for tr in other:
        cir = plt.Circle(tr,MAX_RANGE[ENVIRONMENT],color='grey',fill=False,linestyle=':')
        ax.add_patch(cir)
    plt.show()


def draw(aps, users):
    """
    To plot all the given selection of parameters
    This function is used for visual verification and debugging the above functions 
    (Visual verification of above algorithms is easier in this case)
    """

    # (x, y) for all the AP (x and y translate to distances in meters)
    ap_xpoints, ap_ypoints, ap_channel_colors = [], [], []
    for ap in aps.values():
        ap_xpoints.append(ap['x'])
        ap_ypoints.append(ap['y'])
        ap_channel_colors.append(CHANNEL_COLORS[ap['operating_channel']])

    # (x, y) for all the ordinates of users (x and y translate to distances in meters)
    user_xpoints, user_ypoints, user_channel_colors = [], [], []
    for user in users.values():
        user_xpoints.append(user['x'])
        user_ypoints.append(user['y'])
        user_channel_colors.append(CHANNEL_COLORS[aps[user['associated_ap']]['operating_channel']])

    # plot
    plt.scatter(ap_xpoints, ap_ypoints, marker='o', s=30, c=ap_channel_colors)
    plt.scatter(user_xpoints, user_ypoints, marker='x', s=10, c=user_channel_colors)
    plt.show()


def simulate():

    # for access points
    if ARGS.random:
        aps = assign_channels_for_raa(randomly_arrange_aps(AREA_BREADTH, AREA_LENGTH, LAMBDA_A, 1))
    else:
        aps, ap_matrix = grid_arrange_aps(AREA_LENGTH, AREA_BREADTH, LAMBDA_A)
        aps = assign_channels(aps)

    # drop a random user
    users = drop_users(aps, d_cs=MAX_RANGE[ENVIRONMENT])

    average_throughput_for_all_channels = 0
    for channel in CHANNELS_AVAILABLE:
        # get transmitting aps
        transmitting_aps = get_transmitting_aps(aps, selected_channel=channel, d_cs=MAX_RANGE[ENVIRONMENT])
        # draw transmitting aps
        if ARGS.draw:
            draw_transmitting(aps, transmitting_aps, users)
        # calculate sinr and data rate
        average_data_rate = 0
        for transmitter in transmitting_aps:
            active_aps = list(set(transmitting_aps) - set([transmitter]))
            data_rate = calculate_datarate(transmitter, active_aps, ENVIRONMENT, aps, users)
            average_data_rate = average_data_rate + data_rate
        average_throughput_for_all_channels = average_throughput_for_all_channels + average_data_rate

    area_throughput = average_throughput_for_all_channels / float(AREA_LENGTH * AREA_BREADTH)
    # user_density = len(users) / float(AREA_LENGTH * AREA_BREADTH) # do not need to be calculated
    # user_throughput = area_throughput/ user_density # do not need to be calculated
    
    #print ("Mean Throughput per Area: %1.3f Mbps" % area_throughput)
    #print ("Total number of concurrent Transmitters: %f" % number_of_transmitters)
    #print ("Mean Throughput per User: %1.3f Mbps" % user_throughput)
    #print ()
    #plot
    #draw(aps, users)

    #return area_throughput, user_throughput
    return area_throughput

def get_args():
    parser = argparse.ArgumentParser(description='Simulate area throughput in 2.4GHz Wi-Fi depending on AP deployment density')
    parser.add_argument('-c', dest='bw', type=int, action='store', nargs='+', help='Add a channel and it\'s bandwidth | -c BW BW BW')
    parser.add_argument('-e', dest='environment', type=str, action='store', help='Set an environment | -e home | -e los | -e mall')
    parser.add_argument('-r', dest='random', action='store_true', help='Set use of random AP placement')
    parser.add_argument('-n', dest='loop', type=int, action='store', help='Set nr of loops')
    parser.add_argument('-d', dest='draw', action='store_true', help='Toggle drawing of transmitting APs positions')

    args = parser.parse_args()
    if args.bw != None:
        global CHANNELS_AVAILABLE
        global CHANNEL_BANDWIDTH
        if sum(args.bw) > 60:
            print("Too high bandwidth")
            exit(60)
        CHANNELS_AVAILABLE = np.arange(len(args.bw))
        CHANNEL_BANDWIDTH = args.bw
    if args.environment != None:
        global ENVIRONMENT
        ENVIRONMENT = args.environment
    if args.loop == None:
        args.loop = 1000
    global ARGS
    ARGS = args



if __name__ == "__main__":
    get_args()
    loop_nr=ARGS.loop
    l_a = []
    output = []
    area = AREA_LENGTH * AREA_BREADTH
    for i in range(4,50):
        l_a.append(math.pow(i,2)/area)
    for j in l_a:
        #LAMBDA_A = 0.0001+0.0003*j
        LAMBDA_A = j
        mean_throughput=0
        for i in range(loop_nr):
            area_throughput = simulate()
            mean_throughput += area_throughput
        mean_throughput = mean_throughput/loop_nr
        output.append([LAMBDA_A, mean_throughput])
        print(len(output))
        if len(output) > 2 and abs(mean_throughput - output[-2][1]) < 0.00001:
            break;
    filename = "output_data/"+ENVIRONMENT+str(CHANNEL_BANDWIDTH)
    if ARGS.random:
        filename = 'output_data/random_'+ENVIRONMENT+str(CHANNEL_BANDWIDTH)
    np.savetxt(filename, output, delimiter=',', header="AP density, Mean Throughput Mbps, Environment: "+ENVIRONMENT, fmt='%.4f %f')
