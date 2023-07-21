#! /usr/bin/env python
# cluster the ligand will 
# take into position and the conformation 
# or just take position into account
import sys, math, os
curr_dir = os.get_cwd()
sys.path.append("%s/clustering_scripts"%(curr_dir))
import numpy as np
from Vector_tools import *
from Structure_tools import *
import Island011_ClusterContact as ClusterContact
import Statistic_result as SR
import argparse

def main():
    parser = argparse.ArgumentParser( description =  '''
    Cluster_ligands_position.py -l docked_complex_structure_list -t output_tag    
       clustering step1 based on position of the ligand, average linkage hierarchy-clustering method is applied.
    The program will output three files, tree_temp, cluster, statistic_result
    ''')
    parser.add_argument('-l', nargs=1, help='required, add pdb file list')
    parser.add_argument('-t', nargs=1, help='required, add tag after this')
    parser.add_argument('-s', nargs=1, help='required, what is the protein type (e.g. active_des, inactive_wt)')
    args = parser.parse_args()
    try:
        listfilename = args.l[0]
        output_name = args.t[0]
        state = args.s[0].split("_")[0]
    except:
        parser.print_help()
        exit()
    structure_list = Readin_list( listfilename )
    ligand_list = Generate_ligandInfor_list( structure_list )
    ligand_Position_Matrix = ClusterContact.Vector_position_Distance_Matrix( [ ligand.Mass_Center for ligand in ligand_list] )

    #### HYERPARAMETERS ####
    if state == "active":
        cut_off_distance, cut_off_clusterno = ClusterContact.Iterate_Cluster_vector( output_name, ligand_Position_Matrix, structure_list, output_name+'.cluster', 0.5, 0.05, 0.04) # may need to play around with to get diversity
    elif state == "inactive":
        cut_off_distance, cut_off_clusterno = ClusterContact.Iterate_Cluster_vector( output_name, ligand_Position_Matrix, structure_list, output_name+'.cluster', 0.5, 1.00, 0.90)
    elif state == "crystal":
        cut_off_distance, cut_off_clusterno = ClusterContact.Iterate_Cluster_vector( output_name, ligand_Position_Matrix, structure_list, output_name+'.cluster', 0.5, 0.05, 0.04)
    else:
        cut_off_distance, cut_off_clusterno = ClusterContact.Iterate_Cluster_vector( output_name, ligand_Position_Matrix, structure_list, output_name+'.cluster', 0.5, 1.00, 0.90) # generic pdb file

    SR.Read_statistic( output_name+'.cluster', cut_off_distance, cut_off_clusterno, output_name )

def Generate_ligandInfor_list( structure_list ):
    ligand_list = []
    for structure_file in structure_list:
        Protein_infor = Protein_Structure( structure_file )
        ligand_list.append( Protein_infor.ligand )
    return ligand_list

def Readin_list( listfilename ):
    '''read in the list of structures'''
    file_in = open(listfilename,'r')
    structure_list = []
    for line in file_in.readlines():
        structure_list.append( line.strip() )
    file_in.close()
    return structure_list

if __name__ == '__main__':
    main()

