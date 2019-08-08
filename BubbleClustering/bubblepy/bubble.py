"""
This module is used to generate associator matrix based trees

note: all distances in the module are squared distance as only needed for comparisons


"""

import random
import numpy as np
from ete3 import Tree
import pandas as pd
import scipy
from scipy import spatial
import math


def dista(data_array, i, j):
    """
    distance between two samples in array

    Parameters
    ----------
    data_array : numpy array
                an array that contains the samples for the distance calculation
    i, j : int
            index values of the samples under consideration
    Returns
    -------
    result : float
        The Euclidean squared distance between two rows in the data matrix
    """
    result = sum((data_array[i, ] - data_array[j, ]) ** 2)
    return result


def find_min_max(distance_array):
    """
    find minimum and maximum distances between points in a data matrix

    Parameters
    ----------
    distance_array : numpy array
                an array that contains the distance matrix
    Returns
    -------
    minn, maxx : float
        The minimum distance between two data points and the maximum distance between points in the data matrix

    singleton_bubble: numpy array
        contains the minimum distance between point and next closest point
    """
    singleton_bubble = np.empty(distance_array.shape[0])
    minn = maxx = distance_array[0, 0]
    for i in range(0, (distance_array.shape[0] - 1)):
        singleton_bubble[i] = distance_array[i, i+1]
        for j in range((i + 1), distance_array.shape[0]):
            current_comparison = distance_array[i, j]
            if minn > current_comparison:
                minn = current_comparison  # check if new minimum?
            if maxx < current_comparison:
                maxx = current_comparison  # check if new maximum?
            if singleton_bubble[i] > current_comparison:
                singleton_bubble[i] = current_comparison  # if the
    return minn, maxx, singleton_bubble


def memberate(center, radius, sorted_distance):
    """
    find which points are within a radius of the selected center

    Parameters
    ----------
    center: int
            the index value of the point selected to be the center of the sphere

    radius: float
            the length of the radius selected for the sphere

    sorted_distance : list
                a list of lists of tuples that contains the sorted distance information
    Returns
    -------
    reward : float
        The association score to be added to the i,j pairs in the associator matrix

    membership : list
            list of which points are within the selected radial distance from the chosen center
    """
    membership = []
    '''
    #old code attempting to speed up
    for i in range(0, distance_array.shape[0]):
        current_distance = distance_array[center, i]
        # print(f"{current_distance} as compared to: {radius}")
        if current_distance < radius:
            membership.append(i)
            # print("Added!")
            # if a data point is within the radius distance of the center it is added to the membership list
    '''

    # beginning of log speed up but not sure if necessary
    ls = 0
    rs = len(sorted_distance)-1
    while rs-ls>1:
        new = math.floor((rs-ls)/2)+ls
        if sorted_distance[center][new][0] < radius:
            ls = new
        else:
            rs = new

    for i in range(0,ls+1):
        membership.append(sorted_distance[center][i][1])

    ''' # working not log speed up code
    for i in range(0,len(sorted_distance)):
        if sorted_distance[center][i][0] < radius:
            membership.append(sorted_distance[center][i][1])
        else:
            # print("not in bubble")
            break
    '''

    reward = 1.0 / len(membership)
    # Reward is inversely proportionate to the number of data items in the bubble
    # this allows the algorithm to compensate for over sampling of dense vs sparse regions
    return reward, membership


def associate(associator_matrix, membership, reward):
    """
    increase the association level between two points if they were found in the same sphere

    Parameters
    ----------
    associator_matrix : numpy array
                an array that contains the current association score matrix,
                higher scores indicate a tighter degree of association

    reward : float
        The association score to be added to the i,j pairs in the associator matrix

    membership : list
            list of which points are within the selected radial distance from the chosen center
    """
    if len(membership) == 1:
        k = membership[0]
        associator_matrix[k, k] += 1

    for i in range(0, len(membership) - 1):
        i_2 = membership[i]
        for j in range(i + 1, len(membership)):
            j_2 = membership[j]
            # grab pairs from the membership list and add the reward for the bubble to the already accumulated score
            associator_matrix[i_2, j_2] = associator_matrix[i_2, j_2] + reward


def snp_to_distance(data_name):
    """
        Takes a SNP profile file and generates the Hamming distance matrix

        Parameters
        ----------
        data_name : str
                    name of the file that contains the information
        Returns
        -------
        names: list
            Contains the ordered sample names

        data : numpy array
                the Hamming distance matrix calculated from the supplied allele profiles
        """
    # read the data into python, retaining row and column names, Date Entered is currently dropped from consideration
    sample_info = pd.read_csv(f'{data_name}', sep='\t', index_col=0, header=0)
    sample_info.drop(["Date Entered", "ST"], axis=1, inplace=True)

    data = scipy.spatial.distance.pdist(sample_info, metric='hamming', w=None, V=None, VI=None)
    data = scipy.spatial.distance.squareform(data, force='no', checks=True)
    names = sample_info.index.values

    return data, names


def cluster(data, level, is_distance):
    """
        Takes a SNP profile file and generates the Hamming distance matrix

        Parameters
        ----------
        data : numpy array
                the Hamming distance matrix calculated from the supplied allele profiles

        level: float
            the desired level of bubble coverage, number of iterations determined by n*10**level

        is_distance: boolean
            Indicates if the supplied data is already a distance matrix

        Returns
        -------
        associator_matrix: numpy array
            Contains the association data in an n by n array, higher numbers indicate stronger levels of association
        """
    d_size = data.shape[0]
    iteration = int(10 ** level)*d_size
    associator_matrix = np.full([d_size, d_size], 0.000001)
    distance_matrix = np.full([d_size, d_size], 0.00)
    if is_distance == False :
        # generate the euclidean distance matrix from the supplied data matrix
        for i in range(0, d_size):
            for j in range(0, d_size):
                distance_matrix[i, j] = dista(data, i, j)
    else:
        if is_distance == True :
        distance_matrix = data

        else:
            print("That is an incorrect value for if your data matrix is a distance matrix")

    # generate a data structure that contains the sorted distance information
    dist_from_center=[]
    for i in range(0, d_size):
        distance_by_point=[]
        for j in range(0, d_size):
            distance_by_point.append((distance_matrix[i, j], j))
        distance_by_point.sort()
        dist_from_center.append(distance_by_point)
    # find the global min/max and the local min distance between two points
    minn, maxx, singleton_distance = find_min_max(distance_matrix)
    maxx -= minn  # normalize the maximum

    perc = 0

    for i in range(0, iteration):
        # generate random center
        center = random.randint(0, d_size - 1)
        # generate random radius between the min distance between two points and the max distance between two points
        radius = random.uniform(0, 1) * maxx + minn

        # if the radius < min distance between center and the next closest point, bubble only has center
        if radius < singleton_distance[center]:
            reward = .00000000001
            membership = [center]
        else:
            # get the reward and the list of entries that get that association reward
            reward, membership = memberate(center, radius, dist_from_center)

        # update the associator matrix with the reward information
        # print(f"{i} center:{center} radius {radius} : {reward}, {membership}")
        associate(associator_matrix, membership, reward)
        # Give some indication that the code is doing something, * every 1% to completion and number every 5%

        if i % (iteration / 100) == 0:
            if i % (iteration / 20) == 0:
                print(perc, end='')
                perc += 5
            print("*", end='')

    # print(f"{minn}, {maxx}")
    return associator_matrix


def am_to_tree(associator_matrix, names, level):
    """
           Takes an associator matrix and generates an ete3 tree

           Parameters
           ----------
           associator_matrix : numpy array
                   the associator matrix as calculated by the bubble clustering algorithm
           names:
                the list of sample names to be used for the tree leaf identifiers

           level: float
               the level of bubble coverage used to generate the associator matrix
           Returns
           -------
           tree_structure: ete3 tree
               Contains the tree format of the tree constructed from the supplied association data
           """
    d_size = associator_matrix.shape[0]
    iteration = int(10 ** level)
    # change AM into list of tuples -> sum_reward_value, i,j
    suggested_joins = list()
    for i in range(0, d_size):
        for j in range(i + 1, d_size):
            suggested_joins.append((associator_matrix[i, j], i, j))

    suggested_joins.sort(reverse=True)  # sort the list to get a list of suggested join order

    dist = iteration/suggested_joins[0][0]
    # convert the reward into a distance like measure so that large numbers are now small distances
    joins = [(suggested_joins[0], f"({names[suggested_joins[0][1]]}:{dist},{names[suggested_joins[0][2]]}:{dist})")]
    # convert the first suggested join into a confirmed join and build the newick tree for that join
    last_subtree = list(range(-1, -1 * d_size - 1, -1))
    # update the last_subtree information
    last_subtree[suggested_joins[0][1]] = 0
    last_subtree[suggested_joins[0][2]] = 0

    # initialize iterator for the remainder of the joins
    k = 1

    while len(joins) < d_size-1:
        # extract the indexes of the next suggested join
        i = suggested_joins[k][1]
        j = suggested_joins[k][2]
        k += 1

        # check to see if the samples are already joined or not,
        # if they are already joined they will already be part of the same subtree

        if last_subtree[i] != last_subtree[j]:

            # extract the labels of the subtrees being joined
            before_i = last_subtree[i]
            before_j = last_subtree[j]
            # generate the identifier for the new subtree
            after = len(joins)
            # calculate the branch length based on the association score and the sample size used to generate that score
            dist = iteration/suggested_joins[k][0]

            # build the left branch of the Newick format tree for the new subtree
            if before_i < 0:  # indicates a leaf node
                identifier = (before_i*-1)-1  # converts the leaf node into the sample index
                left_branch = f"{names[identifier]}:{dist}"
            # if not a leaf then an internal node with a pre-existing Newick tree
            else:
                left_branch = f"{joins[before_i][1]}:{dist}"

            # build the right branch of the Newick format tree for the new subtree
            if before_j < 0:  # indicates a leaf node
                identifier = (before_j * -1)-1  # converts the leaf node into the sample index
                right_branch = f"{names[identifier]}:{dist}"

            # if not a leaf then an internal node with a pre-existing Newick tree
            else:
                right_branch = f"{joins[before_j][1]}:{dist}"
            # put the left and right trees together to generate the subtree
            newick_tree = f"({left_branch},{right_branch})"

            # update all instances of last_subtree from old subtree values to new ones
            for i in range(len(last_subtree)):
                if last_subtree[i] == before_j:
                    last_subtree[i] = after
            for i in range(len(last_subtree)):
                if last_subtree[i] == before_i:
                    last_subtree[i] = after
            # add the join to the list of joins along with the Newick tree for that subgroup
            joins.append((suggested_joins[k], newick_tree))
    # generate the final tree by pulling the last subtree and adding the ; used to denote the end of the Newick Tree
    final_tree = f"{joins[-1][1]};"

    # generate an ete3 tree from the Newick Tree
    tree_structure = Tree(f"{final_tree}")
    return tree_structure


def am_to_linkage(associator_matrix, level):
    """
           Takes an associator matrix and generates a scipy linkage matrix

           Parameters
           ----------
           associator_matrix : numpy array
                   the associator matrix as calculated by the bubble clustering algorithm

           level: float
               the level of bubble coverage used to generate the associator matrix
           Returns
           -------
           linkage_matrix: scipy linkage matrix
               Contains the linkage matrix format of the tree constructed from the supplied association data
           """
    d_size = associator_matrix.shape[0]
    iteration = int(10 ** level)
    # change AM into list of tuples -> sum_reward_value, i,j
    suggested_joins = list()
    for i in range(0, d_size):
        for j in range(i + 1, d_size):
            suggested_joins.append((associator_matrix[i, j], i, j))

    suggested_joins.sort(reverse=True)  # sort the list to get a list of suggested join order

    linkage_matrix = np.full([d_size-1, 4], 0.0)
    last_subtree = list(range(0, d_size))

    # convert the reward into a distance like measure so that large numbers are now small distances
    dist = iteration/suggested_joins[0][0]

    linkage_matrix[0, 0] = last_subtree[suggested_joins[0][1]]
    linkage_matrix[0, 1] = last_subtree[suggested_joins[0][2]]
    linkage_matrix[0, 2] = dist
    linkage_matrix[0, 3] = 2
    # convert the first suggested join into a confirmed join in the linkage format

    # update the last_subtree information
    last_subtree[suggested_joins[0][1]] = d_size
    last_subtree[suggested_joins[0][2]] = d_size

    # initialize iterator for the remainder of the joins
    k = 1
    joins = 1

    while linkage_matrix[-1, -1] == 0.0:
        i = suggested_joins[k][1]
        j = suggested_joins[k][2]

        # check to see if the samples are already joined or not,
        # if they are already joined they will already be part of the same subtree

        if last_subtree[i] != last_subtree[j]:

            # extract the labels of the subtrees being joined
            before_i = last_subtree[i]
            before_j = last_subtree[j]
            # generate the identifier for the new subtree
            after = joins + d_size
            # calculate the branch length based on the association score and the sample size used to generate that score
            dist = iteration / suggested_joins[k][0]

            # calculate the number of leaves in the sub tree
            if before_i < d_size:  # indicates a leaf node
                leaf_count = 1
            # if not a leaf then an internal node with a pre-existing subtree leaf count
            else:
                identifier = before_i - d_size
                leaf_count = linkage_matrix[identifier, 3]

            # build the right branch of the Newick format tree for the new subtree
            if before_j < d_size:  # indicates a leaf node
                leaf_count += 1

            # if not a leaf then an internal node with a pre-existing subtree leaf count
            else:
                identifier = before_j - d_size
                leaf_count += linkage_matrix[identifier, 3]

            # update the linkage_matrix
            linkage_matrix[joins, 0] = last_subtree[i]
            linkage_matrix[joins, 1] = last_subtree[j]
            linkage_matrix[joins, 2] = dist
            linkage_matrix[joins, 3] = leaf_count

            # update all instances of last_subtree from old subtree values to new ones
            for i in range(len(last_subtree)):
                if last_subtree[i] == before_j:
                    last_subtree[i] = after
            for i in range(len(last_subtree)):
                if last_subtree[i] == before_i:
                    last_subtree[i] = after

            # update iterators
            k += 1
            joins += 1
        else:
            # increase iterator if skipping line in suggested_joins because information already included
            k += 1
    return linkage_matrix
