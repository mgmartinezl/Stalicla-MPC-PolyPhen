# -*- coding: utf-8 -*-

"""
STALICLA
MPC and pph2 annotations module (1)
This module contains the functions to annotate MPC and pph2 predictions to a dataset of patients and pathways.
It takes three main input files: patients and mutations, pathways and mpc official values document.
Author: Gabriela Martinez - airamgabriela17@gmail.com
Last update: September 2019
"""

import pandas as pd
from time import gmtime, strftime
import argparse
import os
from pathlib import Path
import logging
import datetime
import getpass

local_time = datetime.datetime.now()


def drop_na(df):

    """ Drops records with null or 'NA' genes given by the HGNC_symbol.

    Parameters:
        df (pandas dataframe):Pandas dataframe with the column HGNC_symbol.

    Returns:
        drop_na(df):Pandas dataframe without null or invalid genes

    """

    df = df.dropna(subset=['HGNC_symbol'])
    df = df[pd.notnull(df['HGNC_symbol'])]  # Make sure that null columns are deleted
    df = df[~df['HGNC_symbol'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
    return df


def str_to_int(s):

    """ Function to convert a String object to integer 0.

    Parameters:
        s (None): a None Python object.

    Returns:
        none_to_str(s): an integer Python object.
    """

    import numbers

    convert = lambda i: 0
    if isinstance(s, str):
        s = convert(s)
        return s
    if isinstance(s, numbers.Numbers):
        return s


def none_to_str(s):

    """ Function to convert a None object to string.

    Parameters:
        s (None): a None Python object.

    Returns:
        none_to_str(s): a string Python object.

    """

    convert = lambda i: i or ''
    if s is None:
        s = convert(s)
        return s
    if s is not None:
        return s


def create_arg_parser_mpc():

    """ Creates a parser to read arguments via the command line.

    Parameters:
        parser object from command line.

    Returns:
        create_arg_parser():Args dictionary mapping user inputs.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', help='(Mandatory) path to the input file containing patient mutations file.')
    parser.add_argument('inputMPC', help='(Mandatory) path to file that contains the MPC official files.'
                                         'It can also be a directory with chunks of MPC values.')
    parser.add_argument('-g',
                        '--gene',
                        '--genes',
                        dest='gene',
                        help='Gene(s) to be extracted (default: all). A single gene, a subset of genes separated '
                             'by comma without spaces or a filepath can be provided.',
                        default=None)
    parser.add_argument('-id',
                        '--patient',
                        '--patients',
                        dest='patient',
                        help='Patient(s) id(s) to be extracted (default: all). A single patient id, a subset of '
                             'patients separated by comma without spaces or a filepath can be provided.',
                        default=None)
    parser.add_argument('-csq',
                        '--consequence',
                        '--consequences',
                        dest='csq',
                        help='Consequences(s) to be extracted (default: all). A single consequence, a subset '
                             'of consequences separated by comma without spaces or a filepath can be provided.',
                        default=None)
    parser.add_argument('-mpc',
                        '--mpc_gt',
                        dest='mpc_gt',
                        help='Filters records with values greater than or equal to an MPC threshold',
                        default=None)
    parser.add_argument('-pph2',
                        '--polyphen',
                        '--polyphen2',
                        dest='pph2',
                        help='Filters records for qualifiers of pph2 predictions. Available options: "benign", '
                             '"possibly damaging", "probably damaging".',
                        default=None)
    parser.add_argument('-adj_csq',
                        '--adjusted_consequence',
                        '--adjusted_consequences',
                        dest='adj_csq',
                        help='Filters records by specific value(s) of adjusted consequence. Available '
                             'options: PTV, Missense3, Missense, etc.',
                        default=None)
    parser.add_argument('--path',
                        help='Specify the path where the final annotations will be saved',
                        default=None)
    args = parser.parse_args()
    args_dict = vars(args)
    return args_dict


def read_mutations(filepath):

    """ Reads the mutations file with gene annotations.
        Column names must be respected.

    Parameters:
        filepath (str):The path to the mutations file.

    Returns:
        read_mutations(filepath):Pandas dataframe with mutations file information.

    """

    if filepath:
        input_file = pd.read_csv(filepath, header=0, sep="\t")
        return input_file


def filter_patients(df, patient_id):

    """ Filters a matrix for a set of patients.
    The patient parameter should be pass like: -id Patient_XXXX or --patient Patient_XXXX or --patients Patient_XXXX
    In case more than one patient wants to be included, type: -id Patient_XXXX, Patient_YYYY, Patient_ZZZZ
    This function also accepts a path to a file of patients (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        patient_id (str): single id or subset of ids to be filtered.

    Returns:
        filter_patients(df, patient_id): pandas dataframe object with filter applied.

    """

    df = df[pd.notnull(df['child_id'])]  # Make sure that null columns are deleted
    df = df[~df['child_id'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python

    if patient_id:
        if os.path.isfile(patient_id):
            file_name = Path(patient_id)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['child_id'].isin(file_data)]
            return df
        else:
            if ',' in patient_id:
                patients = patient_id.split(',')
                df = df[df['child_id'].isin(patients)]
                return df
            elif ',' not in patient_id:
                df = df[df.child_id == patient_id]
                return df
    elif patient_id is None:
        return df


def filter_genes(df, gene_id):

    """ Filters a matrix for a set of genes.
    The gene parameter should be pass like: -g XXXXX or --gene XXXXX or --genes XXXXX
    In case more than one gene wants to be included, type: -g XXXXX,YYYYY,ZZZZZ
    This function also accepts a path to a file of genes (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        gene_id (str): single id or subset of ids to be filtered.

    Returns:
        filter_genes(df, gene_id): pandas dataframe object with filter applied.

    """

    if gene_id:
        if os.path.isfile(gene_id):
            file_name = Path(gene_id)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['HGNC_symbol'].isin(file_data)]
            return df
        else:
            if ',' in gene_id:
                patients = gene_id.split(',')
                df = df[df['HGNC_symbol'].isin(patients)]
                return df
            elif ',' not in gene_id:
                df = df[df.HGNC_symbol == gene_id]
                return df
    elif gene_id is None:
        return df


def filter_consequences(df, csq):

    """ Filters a matrix for a set of consequences.
    The consequence parameter should be pass like: -csq desiredcsq or --consequence desiredcsq or --consequences desiredcsq
    In case more than one consequence wants to be included, type: -csq desiredmutation1,desiredmutation2
    This function also accepts a path to a file of consequences (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        csq (str): subset of consequences to be filtered.

    Returns:
        filter_consequences(df, csq): pandas dataframe object with filter applied.

    """

    if csq:
        if os.path.isfile(csq):
            file_name = Path(csq)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['consequence'].isin(file_data)]
            return df
        else:
            if ',' in csq:
                patients = csq.split(',')
                df = df[df['consequence'].isin(patients)]
                return df
            elif ',' not in csq:
                df = df[df.consequence == csq]
                return df
    elif csq is None:
        return df


def filter_PolyPhen(df, pph2):

    """ Filters a matrix by a subset of pph2 labels.
    The parameter should be pass like: -pph2 qualifier or --polyphen qualifier or --polyphen2 qualifier
    Example: -pph2 possibly damaging
             -pph2 "possibly damaging, probably damaging"

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        pph2 (str): subset of labels to be filtered.

    Returns:
        filter_polyphen(df, pph2): pandas dataframe object with filter applied.

    """

    if pph2:
        df = df[pd.notnull(df['pph2_prediction'])]  # Make sure that null columns are deleted
        df = df[~df['pph2_prediction'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        if ',' in pph2:
            pph2_labels = pph2.split(',')
            df = df[df['pph2_prediction'].isin(pph2_labels)]
            return df
        elif ',' not in pph2:
            df = df[df.pph2_prediction == pph2]
            return df
    elif pph2 is None:
        return df


def filter_mpc(df, mpc):

    """ Filters a matrix by an mpc value.
    The parameter should be pass like: -mpc_gt desiredmpc or --mpc desiredmpc
    Results will be shown for mutations with mpc greater or equal to the specified threshold.

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        mpc (float): mpc threshold.

    Returns:
        filter_mpc(df, mpc): pandas dataframe object with filter applied.

    """

    if mpc:
        mpc = float(mpc)
        df = df[pd.notnull(df['MPC'])]  # Make sure that null columns are deleted
        df = df[~df['MPC'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        df = df[df.MPC >= mpc]
        return df
    elif mpc is None:
        return df


def filter_adj_consequence(df, adj_cons):

    """ Filters a matrix by a subset of adjusted consequence labels.
    The parameter should be pass like: -adj_csq desiredcsqs or --adjusted_consequence desiredcsqs or
    --adjusted_consequences desiredcsqs
    Example: -adj_csq Missense3
             -adj_csq Missense3, PTV

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        adj_cons (str): subset of labels to be filtered.

    Returns:
        filter_adj_consequence(df, adj_cons): pandas dataframe object with filter applied.

    """

    if adj_cons:
        df = df[pd.notnull(df['adjusted_consequence'])]  # Make sure that null columns are deleted
        df = df[~df['adjusted_consequence'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        if ',' in adj_cons:
            cons_labels = adj_cons.split(',')
            df = df[df['adjusted_consequence'].isin(cons_labels)]
            return df
        elif ',' not in adj_cons:
            df = df[df.adjusted_consequence == adj_cons]
            return df
    elif adj_cons is None:
        return df


def create_key(df):

    """ Creates a single key field for every mutation found in a file.

    Parameters:
        df (pandas dataframe): matrix containing the mutations defined by Chr | Pos | Ref | Alt

    Returns:
        create_key(df): matrix with a new key column for every mutation.

    """

    df['key'] = df['Chr'].map(str) + "|" + df['Position'].map(str) + "|" + df['Ref'].map(str) + "|" + df['Alt'].map(str)
    return df


def partition_mpc(dir_mpc, dir_chunks, chunk_size):

    """ Function to read big file with MPC official values.
    Partitions MPC 16GB file into X lines chunks.

    Parameters:
        dir_mpc(str, path): directory where the original MPC large file is located.
        dir_chunks (str, path): directory where chunks of the big file will be created and stored.
        chunk_size (integer): number of lines to include in each chunk.

    Returns:
        partition_mpc(dir_mpc, dir_chunks, chunk_size): folder with chunks of MPC values.

    """

    i = 1
    for chunk in pd.read_csv(dir_mpc, header=0, chunksize=chunk_size,
                             usecols=["chrom", "pos", "ref", "alt", "PolyPhen", "MPC"], sep='\t'):
        chunk.to_csv(os.path.join(dir_chunks, "chunk_{}.csv".format(i)), index=False)
        i = i + 1


def annotate_mpc(dir_chunks, patients):

    """ Function to join MPC chunks with patient records.

    Parameters:
        dir_chunks (str, path): folder with chunks of MPC values.
        patients (pandas dataframe): patients dataframe already read.

    Returns:
        annotate_mpc(): patient records with annotated MPC official values.

    """

    chunks_to_read = os.listdir(dir_chunks)
    df = []

    for chunk in chunks_to_read:
        file = chunk
        path = os.path.join(dir_chunks, file)
        small_chunk = pd.read_csv(path, header=0)
        small_chunk['key'] = small_chunk['chrom'].map(str) + "|" + small_chunk['pos'].map(str) + "|" + small_chunk['ref'].map(str) + "|" + small_chunk['alt'].map(str)
        small_chunk = small_chunk[['key', 'PolyPhen', 'MPC']]
        joint_data = pd.merge(patients, small_chunk, on=['key'])
        if not joint_data.empty:
            df.append(joint_data)
        elif joint_data.empty:
            patients['MPC'] = 'NA'
            patients['PolyPhen'] = 'NA(NA)'
            df.append(patients)
    final_df = pd.concat(df)
    return final_df


def blocks_mpc(path, dir_chunks, chunk_size, df):

    """ Function to process the MPC official values entry.
    If it is a file, it will partition it (input file must be a txt) to then annotate values to patients.
    If it is a folder containing chunks of values, the function will iterate them while annotating information.
    This function will call the function annotate_mpc(dir_chunks, patients) on the input provided.

    Parameters:
        path (str, path):
        dir_chunks (str, path):
        chunk_size (integer):
        df (pandas dataframe):

    Returns:
        blocks_mpc(path, dir_chunks, chunk_size, df): patient records with annotated MPC official values.

    """

    if os.path.isfile(path):
        partition_mpc(path, dir_chunks, chunk_size)
        patients_mpc = annotate_mpc(dir_chunks, df)
        return patients_mpc
    elif os.path.isdir(path):
        patients_mpc = annotate_mpc(path, df)
        return patients_mpc


def annotate_pph2_pred(df):

    """ Function to annotate PolyPhen prediction label to a dataframe of patients.

    Parameters:
        df (pandas dataframe): patients dataframe with annotated MPC official values and PolyPhen information.

    Returns:
        annotate_pph2_pred(df): patients dataframe with annotated PolyPhen prediction label

    """

    df['pph2_prediction'] = df.apply(lambda row: row['PolyPhen'].replace('(', ' ').replace(')', '').split(' ')[0], axis=1)
    df['pph2_value'] = df.apply(lambda row: row['PolyPhen'].replace('(', ' ').replace(')', '').split(' ')[1], axis=1)
    df = df.drop(['PolyPhen'], axis=1)
    return df


def adj_consequence(row):

    """ Function to return an adjusted consequence for a patient according to a set of conditions.

    Parameters:
        row (pandas dataframe row): every row of the patients dataframe with annotated MPC official values.

    Returns:
        adj_consequence(row): patients with annotated adjusted consequence.

    """

    if row['consequence'] == 'missense_variant' and (row['pph2_prediction'] == 'probably damaging' or str_to_int(row['MPC']) >= 2):
        return 'Missense3'
    if row['consequence'] == 'missense_variant':
        return 'Missense'
    elif row['consequence'] == 'frameshift_variant' or row['consequence'] == 'splice_acceptor_variant' or \
            row['consequence'] == 'splice_donor_variant' or row['consequence'] =='stop_gained':
        return 'PTV'
    else:
        return row['consequence']


def annotate_adj_cons(df):

    """ Function to annotate adjusted consequence based on the adj_consequence(row) function.

    Parameters:
        df (pandas dataframe): patients dataframe with annotated MPC official values.

    Returns:
        adj_consequence(row): patients with annotated adjusted consequence.

    """

    df['adj_cons'] = df.apply(adj_consequence, axis=1)
    return df


def clean_missense(df):

    """ Function to drop missense records that do not hold MPC or PolyPhen prediction.

    Parameters:
        df (pandas dataframe): patients dataframe with consequences and PolyPhen prediction labels.

    Returns:
        clean_missense(df): patients dataframe witn valid missense records.

    """

    df = df.drop(df[(df.consequence == 'missense_variant') & (df.pph2_prediction == 'unknown')
                    & (df.pph2_value == 0)].index)
    df = df.drop(df[(df.consequence == 'missense_variant') & (df.MPC.isnull())].index)
    return df


def create_id(df):

    """ Function to create a unique id for every patient and mutation.

    Parameters:
        df (pandas dataframe): patients dataframe with mutations.

    Returns:
        create_id(df): patients dataframe with mutations and unique id per register patient | mutation.

    """

    df['_id'] = df.index + 1
    df = df.drop_duplicates(subset='key')
    return df


def format_df(df):

    """ Function to change column names and order.

    Parameters:
        df (pandas dataframe): patients dataframe with mutations and annotated data.

    Returns:
        format_df(df): pandas dataframe with columns names changed.

    """

    df = df.rename(columns={"child_id": "Child_id", "Chr": "Chr", "Position": "Pos", "Ref": "Ref", "Alt": "Alt",
                            "consequence": "Consequence", "HGNC_symbol": "HGNC_Symbol", "key": "Key", "MPC": "MPC",
                            "pph2_prediction": "PolyPhen_pred", "pph2_value": "PolyPhen_Value",
                            "adj_cons": "Adj_Consequence", "_id": "id"})
    df = df[['id', 'Child_id', 'Key', 'Chr', 'Pos', 'Ref', 'Alt', 'Consequence', 'HGNC_Symbol', 'MPC',
             'PolyPhen_pred', 'PolyPhen_Value', 'Adj_Consequence']]
    return df


def export_annotations(df, folderpath):

    """ Exports the MPC annotated dataframe of patients.

    Parameters:
        df (pandas dataframe): dataframe to be exported as a csv file.
        folderpath (str, path): path where the final csv will be stored.

    Returns:
        export_annotations(df, folderpath): annotated dataset of patients.

    """

    df.to_csv(os.path.join(folderpath,
                           'MPC-pph2-annotations-{}.csv'.format(local_time.strftime("%Y-%m-%d_%H꞉%M꞉%S"))),
                           index=False)
    return df


def logging_info_mpc(args, folderpath, folderpath2):

    """ Records an INFO type log.

    Parameters:
        args (dictionary): parsed arguments from command line.
        folderpath (str): destination of the info log.
        folderpath2 (str): folder where the annotated patients file is downloaded.

    Returns:
        logging_info(args, folderpath): info log.

    """

    is_base = "Annotated MPC and PolyPhen2 data downloaded at: {}".format(folderpath2)
    log_filename = os.path.join(folderpath, 'logINFO_{}.log'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    logging.basicConfig(filename=log_filename, level=logging.INFO)

    logging.info('MPC-PolyPhen protocol started generation on: ' + local_time.strftime("%Y-%m-%d_%H꞉%M꞉%S"))
    logging.info('Process initiated by the user ' + getpass.getuser())
    logging.info('\n' + '\n' +
                 "Protocol successfully generated! " + '\n' + '\n' +
                 " ------ Parameters ------ " + '\n' +
                 "Mutations input file: " + args['inputFile'] + '\n' +
                 "MPC input file/folder: " + args['inputMPC'] + '\n' +
                 "Gene(s): " + none_to_str(args['gene']) + '\n' +
                 "Patient(s): " + none_to_str(args['patient']) + '\n' +
                 "Consequence(s): " + none_to_str(args['csq']) + '\n' +
                 "MPC (greater/equal than): " + none_to_str(args['mpc_gt']) + '\n' +
                 "PolyPhen: " + none_to_str(args['pph2']) + '\n' +
                 "Adjusted consequence(s): " + none_to_str(args['adj_csq']) + '\n' + is_base)

    print("---")
    print("Protocol successfully executed!")
