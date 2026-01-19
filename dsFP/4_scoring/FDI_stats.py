import argparse
import os.path

import json, time, os, sys, glob

import numpy as np
import itertools
import random
from collections import defaultdict

import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.conformation import form_disulfide

import sys
sys.path.insert(0, '/home/caolongxingLab/caolongxing/python_packages/xbin/')
sys.path.insert(0, '/home/caolongxingLab/caolongxing/python_packages/getpy/build/lib.linux-x86_64-3.8/')
sys.path.insert(0, '/home/caolongxingLab/caolongxing/python_packages/longxing_scripts/')

import xbin
import getpy
import hash_bt24_util
import basic_utils
import math_util
import neighbor_util

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet,fill_hbond_set

init( "-beta_nov16 -in:file:silent_struct_type binary -output_pose_energies_table false -output_virtual" +
#    " -holes:dalphaball /home/caolongxingLab/caolongxing/bin/DAlphaBall.gcc" +
#    " -extra_res_fa /home/caolx/Zinc_metalloenzyme/rosetta_scripts/ZNX.params" +
    " -extra_res_fa /home/caolongxingLab/caolongxing/IRFP/params/FDI.fa.params /home/caolongxingLab/caolongxing/IRFP/params/FDC.fa.params" +
    " -use_terminal_residues true -mute basic.io.database core.scoring" +
    " -dunbrack_prob_buried 0.8 -dunbrack_prob_nonburied 0.8" +
    " -dunbrack_prob_buried_semi 0.8 -dunbrack_prob_nonburied_semi 0.8" )

def parse_args( argv ):
    argv_tmp = sys.argv
    sys.argv = argv
    description = 'do protein sequence design using the MPNN model ...'
    parser = argparse.ArgumentParser( description = description )
    parser.add_argument('-pdbs', type=str, nargs='*', help='name of the input pdb file')
    parser.add_argument('-pdb_list', type=str, help='a list file of all pdb files')
    parser.add_argument('-probe_radius', type=float, default=1.4, help='the sasa probe radius')
    parser.add_argument('-bb_only', default=False, action='store_true', help='raw wefold bb')
    parser.add_argument('-hbond_energy', type=float, default=0.0, help='the cutoff value of hbond')
    parser.add_argument("-score_file", type=str, default='ABT_burial.sc', help='the output score file')
    args = parser.parse_args()
    sys.argv = argv_tmp

    return args

args = parse_args( sys.argv )

if args.pdbs == None:
    assert (args.pdb_list != None)
    with open(args.pdb_list) as f:
        all_pdbs = [line.strip() for line in f]
else:
    all_pdbs = args.pdbs

def to_numpy(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def to_rosetta(xyz):
    return pyrosetta.rosetta.numeric.xyzVector_double_t(xyz[0], xyz[1], xyz[2])

def transform_residue(rsd, m):
    for ii in range(1, rsd.natoms()+1):
        xyz = to_numpy(rsd.xyz(ii))
        xyz_1 = np.ones(4)
        xyz_1[:3] = xyz
        new_xyz = m @ xyz_1
        rsd.set_xyz(ii, to_rosetta(new_xyz[:3]))

def get_stub(rsd, dtype=np.float64):

    atoms = ['C3','C19','C1'] if rsd.name3() == 'FDC' else ['C8','C27','C6']

    N = to_numpy(rsd.xyz(atoms[0]))
    CA = to_numpy(rsd.xyz(atoms[1]))
    C = to_numpy(rsd.xyz(atoms[2]))

    m = np.zeros((3,3), dtype=dtype)

    m[0,:] = N
    m[1,:] = CA
    m[2,:] = C

    return hash_bt24_util.wefold_bb_stubs(m[None,:])[0]

# creat residues
rts = core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
rsd_gly = core.conformation.ResidueFactory.create_residue( rts.name_map('GLY') )
rsd_FDI = core.conformation.ResidueFactory.create_residue( rts.name_map('FDI') )
rsd_FDC = core.conformation.ResidueFactory.create_residue( rts.name_map('FDC') )

stub_FDC = get_stub(rsd_FDC)

# sfxn
sfxn = get_score_function()
opts = sfxn.energy_method_options()
hb_opts = opts.hbond_options()
hb_opts.decompose_bb_hb_into_pair_energies(True)
opts.hbond_options(hb_opts)
sfxn.set_energy_method_options(opts)

#xml filters
xml = "/home/caolongxingLab/caolongxing/IRFP/scripts/FDI_scoring.xml"
objs = protocols.rosetta_scripts.XmlObjects.create_from_file( xml )
# Load the movers we will need
ddg_norepack = objs.get_filter('ddg_norepack')
if ( isinstance(ddg_norepack, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    ddg_norepack = ddg_norepack.subfilter()
contact_FDC = objs.get_filter('contact_FDC')
if ( isinstance(contact_FDC, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    contact_FDC = contact_FDC.subfilter()
score_per_res = objs.get_filter('score_per_res')
if ( isinstance(score_per_res, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    score_per_res = score_per_res.subfilter()

scores = []
header_str = "description path sidechain_neighbors FDC_neighbors FDI_neighbors"
if not args.bb_only:
    header_str += " score_per_res FDC_sasa FDI_sasa FDI_hbonds ddg_norepack contact_FDC"
scores.append(header_str)

for pdb in all_pdbs:

    if pdb.endswith('.pdb'): tag = pdb.split('/')[-1][:-4]
    elif pdb.endswith('.pdb.gz'): tag = pdb.split('/')[-1][:-7]
    else: continue

    pose = pose_from_pdb(pdb)
    coords, seq = basic_utils.conformation_from_pose(pose, atoms=['N', 'CA', 'C', 'O'], residue_major=True)

    sidechain_neighbors_v = np.sum( neighbor_util.sidechain_neighbors(coords) ) / coords.shape[0]

    FDI_ires = seq.index('Z') + 1
    FDI_rsd = pose.residue(FDI_ires)

    FDC_coords = []
    for atom_n in ['N2', 'N3', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C27', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42']:
        xyz = FDI_rsd.xyz(atom_n)
        FDC_coords.append([xyz.x, xyz.y, xyz.z])
    FDC_coords = np.array(FDC_coords)
    FDC_neighbors = np.sum( neighbor_util.ligand_neighbors(coords, FDC_coords) )

    FDI_coords = FDC_coords.tolist()
    for atom_n in ['C24', 'C25', 'C26', 'C43', 'C44', 'C45']:
    #for atom_n in ['S1', 'O2', 'S2', 'O3', 'O4', 'O5', 'O6', 'O7', 'C24', 'C25', 'C26', 'C43', 'C44', 'C45']:
        xyz = FDI_rsd.xyz(atom_n)
        FDI_coords.append([xyz.x, xyz.y, xyz.z])
    FDI_coords = np.array(FDI_coords)
    FDI_neighbors = np.sum( neighbor_util.ligand_neighbors(coords, FDI_coords) )

    scores_str = f'{tag} {pdb} {sidechain_neighbors_v:.2f} {FDC_neighbors:.2f} {FDI_neighbors:.2f}'

    if not args.bb_only:
        
        score_per_res_v = score_per_res.compute(pose)

        atom_sasa = core.id.AtomID_Map_double_t()
        rsd_sasa  = utility.vector1_double()

        core.scoring.calc_per_atom_sasa(pose, atom_sasa, rsd_sasa, args.probe_radius)
        FDC_sasa = 0.0
        for atom_n in ['N2', 'N3', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C27', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42']:
            iatom = FDI_rsd.atom_index(atom_n)
            FDC_sasa += atom_sasa[core.id.AtomID(iatom, FDI_ires)]

        FDI_sasa = FDC_sasa
        for atom_n in ['C24', 'C25', 'C26', 'C43', 'C44', 'C45']:
        #for atom_n in ['S1', 'O2', 'S2', 'O3', 'O4', 'O5', 'O6', 'O7', 'C24', 'C25', 'C26', 'C43', 'C44', 'C45']
            iatom = FDI_rsd.atom_index(atom_n)
            FDI_sasa += atom_sasa[core.id.AtomID(iatom, FDI_ires)]

        # hbonds
        FDI_hbonds = 0
        sfxn.score(pose)
        hbset = HBondSet()
        derivatives=False
        exclude_bb=True
        exclude_bsc=False
        exclude_scb=False
        exclude_sc=False
        fill_hbond_set(
            pose, derivatives, hbset, exclude_bb, exclude_bsc, exclude_scb,
            exclude_sc
        )
        for bond in hbset.residue_hbonds(FDI_ires):
            if bond.acc_res() == FDI_ires and (not bond.acc_atm_is_backbone()):
                if bond.energy() < args.hbond_energy:
                    print(f"Don Res: {bond.don_res()}; Hbond energy: {bond.energy()}")
                    FDI_hbonds += 1

        # mutate the FDI to GLY
        stub_FDI = get_stub(FDI_rsd)
        m = stub_FDI @ np.linalg.inv(stub_FDC)
        rsd_FDC_work = rsd_FDC.clone()
        transform_residue(rsd_FDC_work, m)

        pose.replace_residue(FDI_ires, rsd_gly, True)
        pose.append_residue_by_jump(rsd_FDC_work, pose.size())
        sfxn.score(pose)
        ddg_v = ddg_norepack.compute(pose)
        contact_v = contact_FDC.compute(pose)

        scores_str += f' {score_per_res_v:.2f} {FDC_sasa:.2f} {FDI_sasa:.2f} {FDI_hbonds} {ddg_v:.2f} {contact_v:.2f}'

    scores.append(scores_str)


with open(args.score_file, 'w') as f:
    for line in scores:
        f.write(line+'\n')
