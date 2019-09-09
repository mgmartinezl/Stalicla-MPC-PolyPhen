# -*- coding: utf-8 -*-

"""
STALICLA
MPC and pph2 annotations module (2)
This module contains the functions to annotate MPC and pph2 predictions to a dataset of patients and pathways.
It takes three main input files: patients and mutations, pathways and mpc official values document.
Author: Gabriela Martinez - airamgabriela17@gmail.com
Last update: September 2019
"""

from src.getMPC import *


def main():

    """ Main function to execute the getMPC protocol """

    # Arguments parser
    args = create_arg_parser_mpc()

    directory = os.path.dirname(os.path.abspath(__file__))

    chunks = 'data/processed/chunks'
    dir_chunks = os.path.join(directory, chunks)
    chunk_size = 10

    default_folder_output = 'analysis'
    dir_default_folder_output = os.path.join(directory, default_folder_output)
    dir_default_logs = os.path.join(directory, default_folder_output)

    if args['path'] is not None:
        dir_output = args['path']
        dir_logs = args['path']
    else:
        dir_output = dir_default_folder_output
        dir_logs = dir_default_logs

    # ------ Read patients and mutations ------ #

    # Input files reading, filtering and cleaning
    patients = drop_na(read_mutations(args['inputFile']))

    patients = create_key(patients)

    # ------ Filtering data ------ #

    # Restrict analysis to specific patients
    patients = filter_patients(patients, args['patient'])

    # Restrict analysis to specific genes
    patients = filter_genes(patients, args['gene'])

    # Restrict analysis to specific consequences
    patients = filter_consequences(patients, args['csq'])

    # ------ MPC annotations ------ #

    # Partition MPC file for easier processing
    # Join MPC, pph2 and adjusted consequence to patients
    # Also, a unique ID is created for every patient and mutation

    patients_mpc = blocks_mpc(args['inputMPC'], dir_chunks, chunk_size, patients)
    patients_mpc_pph2 = annotate_pph2_pred(patients_mpc)
    patients_mpc_pph2_adj = format_df(create_id(annotate_adj_cons(patients_mpc_pph2)))

    # ------ Pos-annotations filtering ------ #

    # Restrict analysis to specific values of pph2 predictions
    data = filter_PolyPhen(patients_mpc_pph2_adj, args['pph2'])

    # Restrict analysis to specific values of MPC
    data = filter_mpc(data, args['mpc_gt'])

    # Restrict analysis to specific values of adjusted consequence
    data = filter_adj_consequence(data, args['adj_csq'])

    # Export annotated patients
    data.to_csv(os.path.join(dir_output,
                'MPC-pph2-annotations-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime()))),
                index=False)

    # ------ Logging ------ #

    logging_info_mpc(args, dir_logs, dir_output)


if __name__ == "__main__":
    main()
