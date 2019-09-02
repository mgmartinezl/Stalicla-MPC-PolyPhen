# -*- coding: utf-8 -*-

"""

STALICLA

MPC and pph2 annotations module (2)

This module contains the functions to annotate MPC and pph2 predictions to a dataset of patients and pathways.

It takes three main input files: patients and mutations, pathways and mpc official values document.

Author: Gabriela Martinez - airamgabriela17@gmail.com
Last update: September 2019

"""

from getMPC import *
import getpass
getpass.getuser()


def main():

    """ Main function to execute the getMPC protocol """

    directory = os.path.dirname(os.path.abspath(__file__))
    chunks = 'InputFiles/Chunks'
    logs = 'Logs'
    output = 'Output'
    dir_chunks = os.path.join(directory, chunks)
    dir_logs = os.path.join(directory, logs)
    dir_output = os.path.join(directory, output)
    chunksize = 10

    create_directory(dir_logs)
    create_directory(dir_output)
    create_directory(dir_chunks)

    # Arguments parser
    args = create_arg_parser_mpc()

    # ------ Read patients and mutations ------ #

    # Input files reading, filtering and cleaning
    patients = drop_na(read_mutations(args['inputFile']))
    pathways = drop_na(read_pathways(args['pathwaysDirectory'], args['pathway']))

    # Joining data: patient mutations and pathways
    joint_data = join_input_data(patients, pathways)

    # Create raw base matrix
    base_matrix = create_base_matrix(joint_data)

    # Process patients data
    patients = process_patients(base_matrix)

    # ------ Filtering data ------ #

    # Restrict analysis to specific patients
    patients = analyze_patients(patients, args['patient'])

    # Restrict analysis to specific genes
    patients = analyze_genes(patients, args['gene'])

    # Restrict analysis to specific mutations
    patients = analyze_mutations(patients, args['mutation'])

    # ------ MPC annotations ------ #

    # Partition MPC file for easier processing
    # Join MPC, pph2 and adjusted consequence to patients
    # Also, a unique ID is created for every patient and mutation

    patients_mpc = blocks_mpc(args['inputMPC'], dir_chunks, chunksize, patients)
    patients_mpc_pph2 = annotate_pph2_pred(patients_mpc)
    patients_mpc_pph2_adj = create_id(annotate_adj_cons(patients_mpc_pph2))

    # Export annotated patients
    patients_mpc_pph2_adj.to_csv(
        os.path.join(dir_output, 'MPC-pph2-annotations-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime()))),
        index=False)

    # ------ MPC + binary matrix for clustering part I ------#

    # Compute binary variables from intermediate matrix from the PBPM protocol
    binary_features = compute_binary(joint_data)
    joint_binary = join_dfs(patients_mpc_pph2_adj, binary_features)

    # Export annotated patients with binary variables
    joint_binary.to_csv(os.path.join(dir_output, 'MPC-pph2-pathways-annotations-{}.csv'.format(
        strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime()))), index=False)

    # ------ Logging ------ #

    logging_info_mpc(args, dir_logs, dir_output)


if __name__ == "__main__":
    main()