#! /usr/bin/env python
# 
import numpy as np
import re, sys, math
from Vector_tools import *

def main():
    if not len(sys.argv) == 4:  
        usage()
        exit(0)
    filename      = sys.argv[1]
    search_region = eval(sys.argv[2])
    search_list   = eval(sys.argv[3])
    #print search_list
    ProteinStructure = Protein_Structure( filename )
    # protein side
    atom_index, atom_matrix = Protein_AtomMatrix(ProteinStructure, search_region)
    # ligand side
    atominfor_list = Extract_Lig_Atom(ProteinStructure, search_list)     
    for atominfor in atominfor_list:
        subMatrix, indicator_range = Matrix_Select( atominfor.AtomLocation, atom_matrix, treshold )        
        indicator_which            = Distance_Search( atominfor.AtomLocation, subMatrix, treshold )
        #print indicator_range, indicator_which
        print(subMatrix[indicator_which,:])
        ProxAtom = np.array( atom_index )[ indicator_range ][ indicator_which]
        print(ProxAtom)

def PoDistance_Ligand(L1, L2):
#    CoordA = [ Atominfor.AtomLocation for Atominfor in L1.Atom_list ]
#    CoordB = [ Atominfor.AtomLocation for Atominfor in L2.Atom_list ]
    return DistanceMeasure( L1.Mass_Center,  L2.Mass_Center )

def RMSD_Ligand(L1, L2):
    assert len(L1.Atom_list) == len(L2.Atom_list), "%d %d" %(len(L1.Atom_list), len(L2.Atom_list))
    Squared_Distance_list = []
    for i in range( len(L1.Atom_list) ):
        Squared_Distance_list.append( Squared_DistanceMeasure( L1.Atom_list[i].AtomLocation, Find_Atom( L2.Atom_list, L1.Atom_list[i].AtomName ).AtomLocation ) )
    #print math.fsum( Squared_Distance_list )
    return math.sqrt( math.fsum( Squared_Distance_list )/ float( len(L1.Atom_list) )  )

def RMSD_Partial_Ligand( L1, L2, Select_Name ):
    #assert( len(L1.Atom_list) == len(L2.Atom_list))
    Squared_Distance_list = []
    for atom in Select_Name:
        Squared_Distance_list.append( Squared_DistanceMeasure( Find_Atom( L1.Atom_list, atom ).AtomLocation, Find_Atom( L2.Atom_list, atom ).AtomLocation ) )
    return math.sqrt( math.fsum( Squared_Distance_list)/ float( len(Squared_Distance_list) ) )

def Find_Atom(Ligand_atom_list, TAtomName):
    #print TAtomName
    for atom_infor in Ligand_atom_list:
        if atom_infor.AtomName == TAtomName: return atom_infor
    print('can not find %s' %TAtomName)
    exit(1)

def Generate_Ligand_RMSD_Matrix( Ligand_list ):
    Distance_Matrix = np.zeros( len(Ligand_list) * len(Ligand_list) )
    Distance_Matrix.shape = len(Ligand_list), len(Ligand_list)
    for i in range( len( Ligand_list ) ):
        for j in range( len( Ligand_list) ):
            Distance_Matrix[i][j] = RMSD_Ligand( Ligand_list[i], Ligand_list[j] )
#   Generate Partial Lrmsd matrix: coded as 5HT1. Need to automate the selection of partial Ligand
#            Distance_Matrix[i][j] = RMSD_Partial_Ligand( Ligand_list[i], Ligand_list[j], ['N5','O3','C19','C25','C24','C15','C16','C26','C32','N4','C20','C12','C10','C11','C5','C4','C8','C14','N1','C9'] )   
#            Distance_Matrix[i][j] = RMSD_Partial_Ligand( Ligand_list[i], Ligand_list[j], ['N4','O3','C19','C25','C24','C15','C16','C26','C32','N5','C20','C12','C10','C11','C5','C4','C8','C14','N1','C9'] )
#            Distance_Matrix[i][j] = RMSD_Partial_Ligand( Ligand_list[i], Ligand_list[j], ['C1','C7','C6','C10','O1','C11','N1','N2','N5','C12','N3','C13','N4','C14','N6','N7'] )
#   Generate Partial Lrmsd matrix: coded for 4PHU
#            Distance_Matrix[i][j] = RMSD_Partial_Ligand( Ligand_list[i], Ligand_list[j], ['C11','O4','C12','C22','C13','C28','C27','C23','O7','C26','C14','C24','C25','O5','O6']  )           
    return Distance_Matrix


def Extract_Lig_Atom( Protein, search_list ):
    ''' search the bonds to atoms you believe to be important '''
    Atominfor_list = []
    #print search_list
    for atom in search_list:
        #print atom
        Atominfor_list.append( Find_Atom( Protein.ligand.Atom_list, atom ) )
    return Atominfor_list

def Protein_AtomMatrix( Protein, ResiRegion=[] ):
    '''ResiRegion format [(resi_1,resi_2),(resi_3,resi_4)]'''
    atom_index = []
    atom_matrix = []
   # if ResiRegion == []:
   # else:
    Protein_res = Protein.Seq_coord.keys()
    Protein_res.sort()
    for resi_no, resi_name in Protein_res:
        # if we find resi in target region
        if Check_region(resi_no, ResiRegion ): 
            for atom in Protein.Seq_coord[(resi_no, resi_name)]: 
                atom_index.append( (resi_no, resi_name, atom.AtomName) ) 
                atom_matrix.append( atom.AtomLocation )
        else: continue
    return atom_index, np.array(atom_matrix)

def Check_region(resi, Region_list ):
    if Region_list == []: return True
    else:
        for region in Region_list:
            if resi > region[0] and resi < region[1]: return True
            else: continue
        return False

def Ligand_profile( ligand_profile ):
    polar_atom_list    = []
    Aromatic_ring_list = []
    rigid_part_list = []
    flexible_part_list = []
    with open(ligand_profile,'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.match('polar',line):
                polar_atom_list = line.strip().split()[1:]
            if re.match('ring', line):
                Aromatic_ring_list = [ eval(atom_list) for atom_list in line.strip().split()[1:] ]
            if re.match('rigid',line):
                #rigid_part_list = [eval(atom_list) for atom_list in line.strip().split()[1:] ]
                rigid_part_list = eval( line.strip().split()[1] )
            if re.match('flexible', line):
                #flexible_part_list = [eval(atom_list) for atom_list in line.strip().split()[1:]]  
                flexible_part_list = eval( line.strip().split()[1] )
    return polar_atom_list, Aromatic_ring_list, rigid_part_list, flexible_part_list

def Aromatic_interaction(RingA, RingB):
    '''define pi stacking between two rings'''
    angle = Angle_Between_Vector( RingA.direction, RingB.direction, 1,1 )
    # print 'angle = %.2f ' %(angle/math.pi)
    # print ( angle > math.pi/3 and angle < math.pi*2/3 )
    if angle <= math.pi/6: 
        # parallel mode, face to face
        T = DistanceMeasure( RingA.ringLocation, RingB.ringLocation) 
        #print T
        if T < 4.4:
            return 'FTF'
        else: return 'NoEffect'
    elif angle > math.pi/3 and angle < math.pi*2/3:
        # T-shape mode, edge to face
        T = DistanceMeasure( RingA.ringLocation, RingB.ringLocation)
        #print T
        if T < 5.5:
            return 'ETF'
        else: return 'NoEffect'
    else: return 'NoEffect'

def set_AA_ring( AA_name, Atom_infor_list ):
    '''AA_name could be one of Phe, Tyr, His, Trp'''
    if AA_name == 'F' or AA_name == 'Y' or AA_name == 'TYR' or AA_name == 'PHE':
        AA_ring = Aromatic_ring()
        AA_ring.Ring_setup( [ Find_Atom(Atom_infor_list, 'CE2').AtomLocation, Find_Atom(Atom_infor_list,'CG').AtomLocation, Find_Atom(Atom_infor_list, 'CE1').AtomLocation ] )
    #################
    elif AA_name == 'H' or AA_name == 'HIS':
        AA_ring = Aromatic_ring()
        AA_ring.Ring_setup( [ Find_Atom(Atom_infor_list, 'CG').AtomLocation, Find_Atom(Atom_infor_list, 'ND1').AtomLocation, Find_Atom(Atom_infor_list, 'CE1').AtomLocation, Find_Atom(Atom_infor_list, 'NE2').AtomLocation,Find_Atom(Atom_infor_list, 'CD2').AtomLocation ] )
    #################
    elif AA_name == 'TRP' or AA_name == 'W':
        AA_ring = Aromatic_ring()
        AA_ring.Ring_setup( [ Find_Atom(Atom_infor_list, 'CH2').AtomLocation, Find_Atom(Atom_infor_list, 'CE2').AtomLocation, Find_Atom(Atom_infor_list, 'CE3').AtomLocation ] )
    return AA_ring 

class Aromatic_ring():
    def __init__(self):
        # self.AtomName = []
        self.ringLocation = ()
        self.direction = ''
     
    def Ring_setup( self,Atom_list ): 
        # Atom_list = Extract_Lig_Atom( Protein_structure, search_list )
        normal_vector, center = PlaneFitTLSTSQ( np.vstack( Atom_list ) )    
        self.ringLocation = center
        self.direction = normal_vector

class Atom_Infor():
    def __init__(self):
        self.AtomName = ''
        self.AtomLocation = ()

    def Atom_readin(self,line):
        self.AtomName = line[11:16].strip()
        self.AtomLocation = np.array( [float(line[30:38]),float(line[38:46]),float(line[46:54])] )
    
    def __repr__(self):
        return '%s %.2f %.2f %.2f '%(self.AtomName, self.AtomLocation[0], self.AtomLocation[1], self.AtomLocation[2])

class Ligand_Structure():
    def __init__(self):
        self.resiNo    = 0
        self.Name      = ''
        self.Atom_list = []     
        self.Mass_Center = []
#        self.AromaticRing = ''

    def Coordinate_store(self, structure_file, HydrInclude ):
        filein = open( structure_file+'.pdb' ,'r')
        for line in filein.readlines():
            if not re.match('HETATM',line): continue
            else:
                if self.resiNo == 0 : self.resiNo = int(line[22:26])
                if self.Name == '' : self.Name = line[17:20].strip()
                AtomInfor = Atom_Infor()
                AtomInfor.Atom_readin(line)
                if HydrInclude == False and re.match('H', AtomInfor.AtomName): continue
                else:
                    self.Atom_list.append( AtomInfor )
        self.Mass_Center = Coordinate_Center( [ Atominfor.AtomLocation for Atominfor in self.Atom_list ] )
        filein.close()

    def Size_ligand( self ):
        ligandsize = []
        for axis in range(3):
            max = 0.0
            for i in range( len( self.Atom_list )-1 ):
                d = self.Atom_list[i].AtomLocation[axis] - self.Atom_list[i+1].AtomLocation[axis]
                if d > max: max = d
            ligandsize.append( max )
        assert len(ligandsize) == 3 , len( ligandsize )
        return ligandsize

class Protein_Structure():
    def __init__(self, filename, hydrogenInclude=False):
        if filename[-3:] == 'pdb':
            self.Id  = filename[:-4]
        else: self.Id = filename
        self.Seq_coord = {}
        self.Coordinate_store()
        self.ligand = Ligand_Structure()
        self.ligand.Coordinate_store( self.Id, hydrogenInclude )
        #print self.Seq_coord.keys()

    def Coordinate_store(self):
        ' storage coordinates from pdb file '
        filein = open( self.Id +'.pdb','r')
        for line in filein.readlines():
            #if not (re.match('ATOM',line) or re.match('HETATM',line)): continue
            if not (re.match('ATOM',line)): continue
            else:
                resiNo = int(line[22:26])
                resiName = line[17:20].strip()
                AtomInfor = Atom_Infor()
                AtomInfor.Atom_readin(line)
                if re.match('H', AtomInfor.AtomName): continue
                self.Seq_coord[(resiNo, resiName)] = self.Seq_coord.setdefault( (resiNo, resiName ),[] )
                self.Seq_coord[(resiNo, resiName)].append( AtomInfor )
        filein.close()

    def Find_ProteinAtom(self, AnyResiNo, AtomName ):
        ' find the atom location '
        for resi in self.Seq_coord.keys():
            if resi[0] == AnyResiNo: 
                for Atom in self.Seq_coord[resi]:
                    if Atom.AtomName == AtomName: 
                        return resi, Atom.AtomName, Atom.AtomLocation
        print(' %d %s not found ' %( AnyResiNo, AtomName))
        return AnyResiNo, AtomName, 0

#def usage():
#    print """ partial_ligand_comparison
#              %s x-ray_structure_file file_list
#              1. all structure files should only contain ligand on the HETATM part.
#""" %sys.argv[0]

if __name__ == '__main__':
    main()
