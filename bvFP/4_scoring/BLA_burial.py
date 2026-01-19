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
    " -extra_res_fa /home/caolongxingLab/caolongxing/wuRFP/scripts/BLA.fa.params /home/caolongxingLab/caolongxing/wuRFP/scripts/FLUOROPHORE.params" +
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
    parser.add_argument('-max_hbond_energy', type=float, default=0.0, help='the sasa probe radius')
    parser.add_argument('-bb_only', default=False, action='store_true', help='raw wefold bb')
    parser.add_argument('-Pfr', default=False, action='store_true', help='raw wefold bb')
    parser.add_argument("-score_file", type=str, default='BLA_burial.sc', help='the output score file')
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
    N = to_numpy(rsd.xyz('C14'))
    CA = to_numpy(rsd.xyz('C8'))
    C = to_numpy(rsd.xyz('C17'))

    m = np.zeros((3,3), dtype=dtype)

    m[0,:] = N
    m[1,:] = CA
    m[2,:] = C

    return hash_bt24_util.wefold_bb_stubs(m[None,:])[0]

# creat residues
rts = core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
rsd_gly = core.conformation.ResidueFactory.create_residue( rts.name_map('GLY') )
rsd_BLA = core.conformation.ResidueFactory.create_residue( rts.name_map('BLA') )
rsd_BLB = core.conformation.ResidueFactory.create_residue( rts.name_map('BLB') )

stub_BLB = get_stub(rsd_BLB)

# scoring func
score_func = get_score_function()

# sfxn
sfxn = get_score_function()
opts = sfxn.energy_method_options()
hb_opts = opts.hbond_options()
hb_opts.decompose_bb_hb_into_pair_energies(True)
opts.hbond_options(hb_opts)
sfxn.set_energy_method_options(opts)

#xml filters
xml = "/home/caolongxingLab/caolongxing/wuRFP/scripts/BLA_scoring.xml"
objs = protocols.rosetta_scripts.XmlObjects.create_from_file( xml )
# Load the movers we will need
ddg_norepack = objs.get_filter('ddg_norepack')
if ( isinstance(ddg_norepack, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    ddg_norepack = ddg_norepack.subfilter()
contact_BLA = objs.get_filter('contact_BLA')
if ( isinstance(contact_BLA, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    contact_BLA = contact_BLA.subfilter()

scores = []
header_str = "description path sidechain_neighbors bla_neighbors D_ring_neighbors"
if not args.bb_only:
    header_str += " score_per_res_with_BLA score_per_res_without_BLA bla_sasa ddg_norepack contact_BLA BLA_hbonds DistalRing_O_hbonds DistalRing_NH_hbonds"
scores.append(header_str)

for pdb in all_pdbs:

    if pdb.endswith('.pdb'): tag = pdb.split('/')[-1][:-4]
    elif pdb.endswith('.pdb.gz'): tag = pdb.split('/')[-1][:-7]
    else: continue

    pose = pose_from_pdb(pdb)
    coords, seq = basic_utils.conformation_from_pose(pose, atoms=['N', 'CA', 'C', 'O'], residue_major=True)

    sidechain_neighbors_v = np.sum( neighbor_util.sidechain_neighbors(coords) ) / coords.shape[0]

    BLA_ires = seq.index('Z') + 1
    BLA_rsd = pose.residue(BLA_ires)

    BLA_coords = []
    for iatom in range(6, BLA_rsd.nheavyatoms()+1):
        xyz = BLA_rsd.xyz(iatom)
        BLA_coords.append([xyz.x, xyz.y, xyz.z])
    BLA_coords = np.array(BLA_coords)
    BLA_neighbors = np.sum( neighbor_util.ligand_neighbors(coords, BLA_coords) )

    DRing_coords = []
    for atm in ['O5','C25','N3','C18','C20','C27','C32','C21','C28']:
        xyz = BLA_rsd.xyz(atm)
        DRing_coords.append([xyz.x, xyz.y, xyz.z])
    DRing_coords = np.array(DRing_coords)
    DRing_neighbors = np.sum( neighbor_util.ligand_neighbors(coords, DRing_coords) )

    scores_str = f'{tag} {pdb} {sidechain_neighbors_v:.2f} {BLA_neighbors:.2f} {DRing_neighbors:.2f}'

    if not args.bb_only:

        score_per_res_with_BLA = score_func.score(pose)/pose.size()

        atom_sasa = core.id.AtomID_Map_double_t()
        rsd_sasa  = utility.vector1_double()

        core.scoring.calc_per_atom_sasa(pose, atom_sasa, rsd_sasa, args.probe_radius)

        v = 0.0
        for iatom in range(6, BLA_rsd.nheavyatoms()+1):
            v += atom_sasa[core.id.AtomID(iatom, BLA_ires)]

        # hbonds
        BLA_hbonds = 0
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
        for bond in hbset.residue_hbonds(BLA_ires):
            if bond.acc_res() == BLA_ires and (not bond.acc_atm_is_backbone()):
                BLA_hbonds += 1
            if bond.don_res() == BLA_ires and (not bond.don_hatm_is_backbone()):
                BLA_hbonds += 1


        O5_hbonds = 0
        H9_hbonds = 0
        for hbond in hbset.residue_hbonds(BLA_ires):
            if hbond.energy() > args.max_hbond_energy: continue
            if hbond.don_res() == BLA_ires and pose.residue(BLA_ires).atom_name(hbond.don_hatm()).strip() == "H9":
                H9_hbonds += 1
            elif hbond.acc_res() == BLA_ires and pose.residue(BLA_ires).atom_name(hbond.acc_atm()).strip() == "O5":
                O5_hbonds += 1


        # mutate the AFT to GLY
        stub_BLA = get_stub(BLA_rsd)
        m = stub_BLA @ np.linalg.inv(stub_BLB)
        rsd_BLB_work = rsd_BLB.clone()
        transform_residue(rsd_BLB_work, m)

        if args.Pfr:
            for atm in ['C18', 'N3', 'H20', 'H9', 'C25', 'O5', 'C21', 'C28', 'H28', 'H29', 'H30', 'C20', 'C27', 'C32', 'H35', 'H36', 'H27']:
                rsd_BLB_work.set_xyz(atm, BLA_rsd.xyz(atm))


        pose.replace_residue(BLA_ires, rsd_gly, True)
        score_per_res_without_BLA = score_func.score(pose)/pose.size()
        pose.append_residue_by_jump(rsd_BLB_work, pose.size())
        sfxn.score(pose)
        ddg_v = ddg_norepack.compute(pose)
        contact_v = contact_BLA.compute(pose)

        scores_str += f' {score_per_res_with_BLA:.2f} {score_per_res_without_BLA:.2f} {v:.2f} {ddg_v:.2f} {contact_v:.2f} {BLA_hbonds} {O5_hbonds} {H9_hbonds}'

    scores.append(scores_str)


with open(args.score_file, 'w') as f:
    for line in scores:
        f.write(line+'\n')
