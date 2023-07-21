#! /usr/bin/env python
# this is the statistic analysis program for the result in trimer comparison
#
import sys

def main():
    clusters = sys.argv[1]
    tag = '.'.join( clusters.split('.')[:-1] )
    Read_statistic(clusters,  0, 0, tag)

def Read_tree_temp( tag, max_num ):
    filein = open( tag+'_tree_temp' )
    cut_list = []
    for line in filein.readlines():
        cut_list.append( float( line.strip().split(':')[1].strip())  )
    cut_list.reverse()
    filein.close()
    return cut_list[max_num-1] 

def Read_statistic(clusters, rmsd_cut, max_num, tag ):
    clusterfile = open(clusters,'r')
    Stat_dict = {}
    for line in clusterfile.readlines():
        clusterno = int(line.strip().split(':')[-1])
        if '.pdb' in line:
            decoy = line.strip().split(':')[0].strip()[:-4]
        else: decoy = line.strip().split(':')[0].strip()
        Stat_dict.setdefault( clusterno, [])
        Stat_dict[ clusterno ].append( decoy )
    clusterfile.close()

    # output the '*' histogram
    for id in Stat_dict.keys():
        statics = '*'*len(Stat_dict[id])
        print('%3d\t%3d\t%s' %(id,len(Stat_dict[id]),statics))
        f = open( tag + '.cluster%d'%id,'w' ) 
        f.write( '\n'.join(Stat_dict[id]) )
        f.close()
    if rmsd_cut == 0 and rmsd_cut == 0:
    #############
        key_list = Stat_dict.keys()
        max_num = max( map( int, key_list ) ) + 1
        rmsd_cut = Read_tree_temp( tag, max_num )    
    #############

    # output the text file 
    f = open(tag+'.summary','w')
    for id in Stat_dict.keys():
        f.write('%3d\t%d\n' %(id, len(Stat_dict[id]) ) )
    f.write('distance_cut:%.2f, clusters:%d\n' %(rmsd_cut, max_num) )
    f.close()

    # output each cluster
    return max_num

if __name__ == '__main__':
    main()
