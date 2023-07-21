#! /usr/bin/env python
# tools that cluster the contacts with the Ca as reference point.
from Bio.Cluster import *
from Vector_tools import *
import numpy as np

def Save_Distance_Matrix( Distance_Matrix, filename ):
    outfile = open( filename, 'w' )
    I,J = Distance_Matrix.shape
    for i in range(I):
        print(">> outfile", '\t'.join( map(str, Distance_Matrix[i]) ))
    outfile.close()

def Iterate_Cluster_vector(name, matrix, namelist, outputfile, cutoff_step, upperbound, lowerbound):
    '''distance cutoff increases from 0.3-0.8 if an increase of 0.1 does not merge more clusters, stop '''
    #    move the function of deciding the treshold to the treecutnumber function
    #    for distance_cutoff in np.arange(0.3,1.1,0.1):        
    #tree = treecluster(method = 'a',distancematrix = matrix)
    # https://biopython.org/docs/1.75/api/Bio.Cluster.html
    tree = treecluster(data=None, method = 'a',distancematrix = matrix)
    Nodelist = []
    #for eachnode in tree: 
    #    Nodelist.append(str(eachnode))
    for i in range(len(tree)):
        Nodelist.append(str(tree[i]))
    #    for distance_cutoff in np.arange(0.3,1,0.1):
    #    clusterNum_cutoff, distance_cutoff = Auto_TreeCutNumber( name, Nodelist, cutoff_step, upperbound, lowerbound )
    clusterNum_cutoff, distance_cutoff = Auto_TreeCutNumber_new( name, Nodelist, upperbound, lowerbound )
    clusterid = tree.cut( clusterNum_cutoff )
    print('distance_cut:%.2f, clusters:%d' %(distance_cutoff, clusterNum_cutoff ))
    out_file = open(outputfile,'w')
    rank_list = []
    for i in range(len(clusterid)):
        rank_list.append((clusterid[i],namelist[i]))
    for j in sorted(rank_list):
        out_file.write('%110s :%2d\n' %(j[1],j[0]))
    out_file.close()
    #  return a dictionary, to find each 
    return distance_cutoff, clusterNum_cutoff

def Vector_cosine_Distance_Matrix( Vector_set ):
    Distance_Matrix = np.zeros( len(Vector_set)*len(Vector_set) )
    Distance_Matrix.shape = len(Vector_set), len(Vector_set)
    scaled_VS = [ Scale_Vector(Vector) for Vector in Vector_set ]
    #print scaled_VS
    for i in range(len(scaled_VS)):
        for j in range(len(scaled_VS)):
            Distance_Matrix[i][j] = DistanceMeasure( scaled_VS[i], scaled_VS[j] )
    return Distance_Matrix

def Vector_position_Distance_Matrix( Vector_set ):
    '''set of tuple (x,y,z) stand for the positions'''
    Distance_Matrix = np.zeros( len(Vector_set)*len(Vector_set) )
    Distance_Matrix.shape = len(Vector_set), len(Vector_set)
    for i in range(len( Vector_set)):
        for j in range(len(Vector_set)):
            Distance_Matrix[i][j] = DistanceMeasure( Vector_set[i], Vector_set[j] )
    return Distance_Matrix

def Auto_TreeCutNumber(name, Nodelist, cutoff_step, upperbound, lowerbond ):
    # keep in a file for the record
    # cutoff step means that you know a merge larger than certain treshold is important enough to separate clusters
    number_nodelist = []
    tempfile = open('%s_tree_temp' %name,'w')
    for eachnode in Nodelist:
        number_nodelist.append(float(eachnode.split(':')[-1]))
        tempfile.write(eachnode+'\n')
    tempfile.close()
    # decide the cut
    I = len(number_nodelist) - 1
    for i in range( len(number_nodelist)-1 ):
        if number_nodelist[i+1] > upperbound :
            I = i
            break
        if ( number_nodelist[i+1] -number_nodelist[i]> cutoff_step and  number_nodelist[i+1] >= lowerbond ): 
            I = i
            break
        else: continue
    cutnumber = len(number_nodelist)-I
    #print 'for rmsd equals %f, we got %d clusters' %(rmsd_cutoff, cutnumber)
    return cutnumber, number_nodelist[I]


def Auto_TreeCutNumber_new(name, Nodelist, upperbound, lowerbond ):
    # keep in a file for the record
    # when you don't know which step make more sense, but have a idea about upper bond and lower bond
    number_nodelist = []
    tempfile = open('%s_tree_temp' %name,'w')
    for eachnode in Nodelist:
        number_nodelist.append(float(eachnode.split(':')[-1]))
        tempfile.write(eachnode+'\n')
    tempfile.close()
    # decide the cut
    I = len(number_nodelist) - 1
    max_step = 0
    for i in range( len(number_nodelist)-1 ):
        if number_nodelist[i+1] > upperbound :
            I = i
            break
        if number_nodelist[i+1] -number_nodelist[i]> max_step and number_nodelist[i+1] >= lowerbond: 
            max_step = number_nodelist[i+1] -number_nodelist[i]
            continue
        else: continue
    # now we know the max_step
    for i in range( len(number_nodelist)-1 ):
        if number_nodelist[i+1] < lowerbond: continue
        elif number_nodelist[i+1] -number_nodelist[i] == max_step:
           I = i
           break
        else: continue
    cutnumber = len(number_nodelist)-I
    # print 'for rmsd equals %f, we got %d clusters' %(utoff, cutnumber)
    return cutnumber, number_nodelist[I]
 
def treecutnumber(name, Nodelist, cutoff):
    # keep in a file for the record
    number_nodelist = []
    tempfile = open('%s_tree_temp' %name,'w')
    for eachnode in Nodelist:
        number_nodelist.append(float(eachnode.split(':')[-1]))
        tempfile.write(eachnode+'\n')
    tempfile.close()
    i  = 0
    while(cutoff > number_nodelist[i] and i < len(number_nodelist)-1):
        i += 1
    cutnumber = len(number_nodelist)+1-i
    #print 'for rmsd equals %f, we got %d clusters' %(rmsd_cutoff, cutnumber)
    return cutnumber

#if __name__ == '__main__':
#    main()
