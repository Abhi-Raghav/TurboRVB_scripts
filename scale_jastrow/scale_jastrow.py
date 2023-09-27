#!/usr/bin/env python3

from turbogenius.pyturbo.io_fort10 import IO_fort10
import numpy as np
import matplotlib.pyplot as plt
import shutil
import argparse
import logging
import scienceplots
from matplotlib import cm
plt.style.use('science')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument("-wfu", type=str, default="fort.10_unit", help="WF file for unit cell")
parser.add_argument("-wfs", type=str, default="fort.10_sp", help="WF file for the supercell")
parser.add_argument("-j", "--jobtype", type=str, choices=["scale", "plot"], help="Choose scale to scale jastrow and plot for plotting jastrow matrices")
args = parser.parse_args()
logging.info(args)

fort10_unit = IO_fort10(args.wfu)
fort10_super = IO_fort10(args.wfs, in_place=True)
jobtype = args.jobtype

def scale_jas_matrix(fort10_unit, fort10_super):
    nel_tot_unit = fort10_unit.f10header.nel
    logging.info(f'Electrons in the unit cell = {nel_tot_unit}')
    nel_tot_super = fort10_super.f10header.nel
    logging.info(f'Electrons in the supercell = {nel_tot_super}')
    unit_cell_repetition = int(nel_tot_super / nel_tot_unit)
    logging.info(f'Number of unit cell repetitions = {unit_cell_repetition}')
    scaling_inh_one_body_j = (nel_tot_unit - 1) / (nel_tot_super - 1)
    logging.info(f'Scaling factor for inh. one body J = {scaling_inh_one_body_j}')
    jas_occ_unit = len(fort10_unit.f10jasocc.occ)
    jas_occ_super = len(fort10_super.f10jasocc.occ)
    jas_mat_nonzero_unit = fort10_unit.f10header.jas_mat_nonzero
    jas_mat_nonzero_coeff_super = fort10_super.f10jasmatrix.coeff

    # main part of algorithm starts here
    logging.info(f"Jastrow scaling starts ...")
    for index, colval in enumerate(fort10_unit.f10jasmatrix.col):
        if colval == jas_occ_unit:   # inhomogeneous one body jastrow
            jas1body_unit = fort10_unit.f10jasmatrix.coeff[index]
            jas1body_scaled = scaling_inh_one_body_j * jas1body_unit
            if index == (jas_mat_nonzero_unit - 1): # The last element is repeated only once
                jas_mat_nonzero_coeff_super[index + (unit_cell_repetition - 1)*(jas_mat_nonzero_unit - 1)] = jas1body_scaled
            else:
                for i in range(0, unit_cell_repetition):
                    jas_mat_nonzero_coeff_super[index + i*(jas_mat_nonzero_unit - 1)] = jas1body_scaled

        else:   # 3body jastrow
            jas3body_unit = fort10_unit.f10jasmatrix.coeff[index]
            for i in range(0, unit_cell_repetition):
                jas_mat_nonzero_coeff_super[index + i*(jas_mat_nonzero_unit - 1)] = jas3body_unit

    fort10_super.f10jasmatrix.coeff = jas_mat_nonzero_coeff_super
    logging.info(f"Jastrow scaling ends ...")


def plot_jas_matrix(fort10_unit, fort10_super):
    logging.info(f'Plotting Unit cell jasmatrix ...')
    fig1, axes1 = plt.subplots()
    axes1.set_title('Unit cell jasmatrix')
    axes1.grid()
    axes1.set_xlabel('Col')
    axes1.set_ylabel('Row')

    max_coeff_unit = np.max(fort10_unit.f10jasmatrix.coeff)
    max_coeff_super = np.max(fort10_super.f10jasmatrix.coeff)
    min_coeff_unit = np.min(fort10_unit.f10jasmatrix.coeff)
    min_coeff_super = np.min(fort10_super.f10jasmatrix.coeff)
    vmax = np.max([max_coeff_unit, max_coeff_super])
    vmin = np.min([min_coeff_unit, min_coeff_super])

    axes1.scatter(x=fort10_unit.f10jasmatrix.col, y=fort10_unit.f10jasmatrix.row, c=fort10_unit.f10jasmatrix.coeff, vmin=vmin, vmax=vmax, cmap=cm.coolwarm)
    axes1.plot()
    fig1.set_size_inches(5, 5)
    fig1.set_dpi(100)
    #fig1.colorbar(cm.ScalarMappable(cmap=cm.coolwarm), ax=axes1)
    plt.savefig('unit_jasmatrix.pdf')
    logging.info(f'Plotting Unit cell jasmatrix ... Done !')

    logging.info(f'Plotting supercell jasmatrix ...')
    fig2,axes2 = plt.subplots()
    axes2.set_title('Supercell jasmatrix')
    axes2.grid()
    axes2.set_xlabel('Col')
    axes2.set_ylabel('Row')
    axes2.scatter(x=fort10_super.f10jasmatrix.col, y=fort10_super.f10jasmatrix.row, c=fort10_super.f10jasmatrix.coeff, vmin=vmin, vmax=vmax, cmap=cm.coolwarm)
    axes2.plot()
    fig2.set_size_inches(10, 10)
    fig2.set_dpi(100)
    plt.savefig('super_jasmatrix.pdf')
    logging.info(f'Plotting supercell jasmatrix ... Done !')
    plt.show()

if jobtype == "scale":
    # First copy the fort10_sp for the supercell to fort.10_sp0. Keep the original files unchanged
    logging.info(f'Jastrow matrix from {args.wfu} will be scaled to {args.wfs} !')
    scale_jas_matrix(fort10_unit, fort10_super)
    logging.info('Job done !')
else:
    plot_jas_matrix(fort10_unit, fort10_super)


