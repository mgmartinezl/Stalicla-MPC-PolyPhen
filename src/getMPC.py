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
from datetime import datetime
import getpass


def create_directory(folderpath):

    """ Creates a folder.

    Parameters:
        folderpath (str):The path to the folder that will be create.

    Returns:
        create_directory(folderpath):The created folder in the specified path.

    """

    if not os.path.exists(folderpath):
        os.makedirs(folderpath)


def create_arg_parser_mpc():

    """ Creates a parser to read arguments via the command line.

    Parameters:
        parser object from command line.

    Returns:
        create_arg_parser():Args dictionary mapping user inputs.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', help='(Mandatory): Path to the input file containing patient mutations file')
    parser.add_argument('pathwaysDirectory', help='(Mandatory ):Path to the directory that contains the pathway files '
                                                  'with annotated genes')
    parser.add_argument('inputMPC', help='Path of the file containing the official MPC values ')
    parser.add_argument('-pathway', help='Pathway(s) to be extracted. If not given, all by default. '
                                         'A single pathway, a subset of pathways separated '
                                         'by comma without spaces or a filepath (tab text with no headers '
                                         'and a single column) can be provided. '
                                         'Examples: -pathway R-HSA-69620 | -pathway R-HSA-69620,0051705 | '
                                         '-pathway path-to-my-file.txt', default=None)
    parser.add_argument('-gene', help='Gene(s) to be extracted. '
                                      'If not given, all by default. '
                                      'A single gene, a subset of genes separated by comma without spaces or a '
                                      'filepath (tab text with no headers and a single column) can be '
                                      'provided. Examples: -gene CTR9 | -gene CTR9,NOCL2 | '
                                      '-gene path-to-my-file.txt', default=None)
    parser.add_argument('-patient', help='Patient(s) id(s) to be extracted. '
                                         'If not given, all by default. '
                                         'A single patient id, a subset of patients separated '
                                         'by comma without spaces or a filepath (tab text with no headers and a '
                                         'single column) can be provided. '
                                         'Examples: -patient Patient_X | -patient Patient_X,Patient_Y | '
                                         '-patient path-to-my-file.txt', default=None)
    parser.add_argument('-mutation', help='Mutation(s) to be extracted. '
                                          'If not given, all by default. '
                                          'A single mutation, a subset of mutations separated '
                                          'by comma without spaces or a filepath (tab text with no headers and a '
                                          'single column) can be provided. '
                                          'Examples: -mutation missense_variant | -mutation missense_variant,Intron | '
                                          '-mutation path-to-my-file.txt', default=None)

    args = parser.parse_args()
    args_dict = vars(args)
    return args_dict


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


def str_to_int(s):

    """ Function to convert a String object to integer 0

    Parameters:
        s (None): a None Python object.

    Returns:
        none_to_str(s): a int Python object.

    """

    import numbers

    convert = lambda i: 0
    if isinstance(s, str):
        s = convert(s)
        return s
    if isinstance(s, numbers.Numbers):
        return s


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
        # input_file['MPC'] = np.random.randint(0, 5, size=(len(input_file), 1)) # VALIDAR CON LAURA SI ES CALCULADO
        # input_file['adjusted_consequence'] = input_file.apply(lambda row: adj_consequence(row), axis=1) # VALIDAR
        # input_file.to_csv('mutations-last.csv', index=False)
        return input_file


def extract_pathways(folderpath, pathway):

    """ Extracts directory locations for a subset of pathways to be analyzed.
        The pathway parameter should be pass like: -pathway yourdesiredpathways.
        This function process the pathway argument as follows:

        * if you type '-pathway PathwayName' only information for this pathway will be computed
        * if you type '-pathway PathwayName1,PathwayName2' pathways with those two names will be computed
        * if you type nothing, by default, all the pathways will be computed

    Parameters:
        folderpath (str):The path to the folder containing the pathways files.
        pathway (str): single string or collection of paths to be analyzed.

    Returns:
        read_mutations(folderpath, pathway):List or string containing the pathways files' paths to process

    """

    if folderpath and pathway:
        if os.path.isfile(pathway):
            file_name = Path(pathway)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            files = []
            for i in file_data:
                file_name = '{}_with_gene_annotations'.format(i)
                path_to_find = os.path.join(folderpath, file_name + '.' + 'txt')
                files.append(path_to_find)
            return files
        elif ',' not in pathway:
            files = []
            file_name = '{}_with_gene_annotations'.format(pathway)
            path_to_find = os.path.join(folderpath, file_name + '.' + 'txt')
            files.append(path_to_find)
            return files
        elif ',' in pathway:
            paths = pathway.split(',')
            files = []
            for i in paths:
                file_name = '{}_with_gene_annotations'.format(i)
                path_to_find = os.path.join(folderpath, file_name + '.' + 'txt')
                files.append(path_to_find)
            return files
    elif folderpath and pathway is None:
        files = []
        for file in os.listdir(folderpath):
            if file.endswith('_with_gene_annotations.txt'):
                files.append(os.path.join(folderpath, file))
        return files


def read_pathways(folderpath, pathway):

    """ Reads a subset of pathways files.

    Parameters:
        folderpath (str):The path to the folder containing the pathways files.
        pathway (str): single string or collection of paths to be analyzed.

    Returns:
        read_pathways(folderpath, pathway):Pandas dataframe containing the pathways files' joint data

    """

    frames = []
    paths_to_find = extract_pathways(folderpath, pathway)
    for path in paths_to_find:
        g = pd.read_csv(path, header=0, sep='\t', engine='python')
        g['Pathway'] = '{}'.format(path.split('_with_gene_annotations')[0].split('/')[-1])
        g.columns = ['Gene', 'HGNC_symbol', 'HGNC_mapping', '%RVIS_ESP_0.1%', '%RVIS_ExAC_0.01%',
                     '%RVIS_ExAC_0.05%popn', '%RVIS_ExAC_0.1%popn', 'constraint_score', 'pLI', 'pRecessive', 'pNull',
                     'pHI', 'chrN:N-N', 'Pathway']
        g = g[['Pathway', 'Gene', 'HGNC_symbol', 'HGNC_mapping', '%RVIS_ESP_0.1%', '%RVIS_ExAC_0.01%',
               '%RVIS_ExAC_0.05%popn', '%RVIS_ExAC_0.1%popn', 'constraint_score', 'pLI', 'pRecessive', 'pNull', 'pHI',
               'chrN:N-N']]
        frames.append(g)
    data = pd.concat(frames, ignore_index=True)
    data = data.sort_values(by='Pathway', ascending=True)
    return data


def join_input_data(df1, df2):

    """ Joins data coming from mutations file and pathways, on HGNC_symbol column.

    Parameters:
        df1 (pandas dataframe):Pandas dataframe with mutations and patients data.
        df2 (pandas dataframe):Pandas dataframe with the pathways data.

    Returns:
        join_input_data(df1, df2):joint pandas dataframe.

    """

    pathways = df2[['Pathway', 'HGNC_symbol']]
    joint_data = pd.merge(df1, pathways, on=['HGNC_symbol'])
    return joint_data


def create_base_matrix(df):

    """ Creates base data matrix containing summarized information per patient and mutation (given by
    chromosome, position, reference and alteration).

    Parameters:
        df (pandas dataframe):Pandas dataframe with mutations, patients and pathways data.

    Returns:
        create_base_matrix(df):joint pandas dataframe.

    """

    pivot_table = pd.pivot_table(df,
                                 index=['child_id', 'Chr', 'Position', 'Ref', 'Alt', 'consequence'],
                                 values=['HGNC_symbol', 'Pathway'],
                                 aggfunc=lambda x: ', '.join(set(x)))
    flattened = pd.DataFrame(pivot_table.to_records())
    return flattened


def analyze_patients(df, patient_id):

    """ Filters a matrix for a set of patients.
    The patient parameter should be pass like: -patient Patient_XXXX
    In case more than one patient wants to be included, type: -patient Patient_XXXX, Patient_YYYY, Patient_ZZZZ
    This function also accepts a path to a file of patients (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        patient_id (str): single id or subset of ids to be filtered.

    Returns:
        analyze_patients(df, patient_id): pandas dataframe object with filter applied.

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


def analyze_genes(df, gene_id):

    """ Filters a matrix for a set of genes.
    The gene parameter should be pass like: -gene XXXXX
    In case more than one gene wants to be included, type: -gene XXXXX,YYYYY,ZZZZZ
    This function also accepts a path to a file of genes (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        gene_id (str): single id or subset of ids to be filtered.

    Returns:
        analyze_genes(df, gene_id): pandas dataframe object with filter applied.

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


def analyze_mutations(df, mutation):

    """ Filters a matrix for a set of mutations.
    The mutation parameter should be pass like: -mutation desiredmutation
    In case more than one mutation wants to be included, type: -mutation desiredmutation1,desiredmutation2
    This function also accepts a path to a file of mutations (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        mutation (str): subset of mutations to be filtered.

    Returns:
        analyze_mutations(df, mutation): pandas dataframe object with filter applied.

    """

    if mutation:
        if os.path.isfile(mutation):
            file_name = Path(mutation)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['consequence'].isin(file_data)]
            return df
        else:
            if ',' in mutation:
                patients = mutation.split(',')
                df = df[df['consequence'].isin(patients)]
                return df
            elif ',' not in mutation:
                df = df[df.consequence == mutation]
                return df
    elif mutation is None:
        return df


def process_patients(df):

    # Create a unique key for every row
    df['key'] = df['Chr'].map(str) + "|" + df['Position'].map(str) + "|" + df['Ref'].map(str) + "|" + df['Alt'].map(str)

    # Count number of affected pathways
    df['affected_pathways'] = df.Pathway.map(lambda x: [i.strip() for i in x.split(",")])
    df['#_affected_pathways'] = df.affected_pathways.apply(len)
    df = df.drop(['affected_pathways'], axis=1)
    return df


def partition_mpc(dir_mpc, dir_chunks, chunksize):

    """

    Function to read big file with MPC annotations.
    Partition MPC 16GB file into X lines chunks.

    """

    i = 1
    for chunk in pd.read_csv(dir_mpc, header=0, chunksize=chunksize,
                             usecols=["chrom", "pos", "ref", "alt", "PolyPhen", "MPC"], sep='\t'):
        chunk.to_csv(os.path.join(dir_chunks, "chunk_{}.csv".format(i)), index=False)
        i = i + 1


def annotate_mpc(dir_chunks, patients):

    """

    Function to join mpc chunks with patients records.

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


def annotate_pph2_pred(df):

    """ Function to annotate adjusted consequence """

    df['pph2_prediction'] = df.apply(lambda row: row['PolyPhen'].replace('(', ' ').replace(')', '').split(' ')[0], axis=1)
    df['pph2_value'] = df.apply(lambda row: row['PolyPhen'].replace('(', ' ').replace(')', '').split(' ')[1], axis=1)
    df = df.drop(['PolyPhen'], axis=1)
    return df


def adj_consequence(row):

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

    """ Function to annotate adjusted consequence """

    df['adj_cons'] = df.apply(adj_consequence, axis=1)
    return df


def create_id(df):

    # Drop missense records that do not hold MPC or PolyPhen prediction
    df = df.drop(df[(df.consequence == 'missense_variant') & (df.pph2_prediction == 'unknown') & (df.pph2_value == 0)].index)
    df = df.drop(df[(df.consequence == 'missense_variant') & (df.MPC.isnull())].index)

    # Create an id (we have duplicated patients)
    df['_id'] = df.index + 1
    df = df.drop_duplicates(subset='key')

    # Reorder data and change some names
    df = df.rename(columns={"child_id": "Child_id", "Chr": "Chr", "Position": "Pos", "Ref": "Ref", "Aly": "Alt",
                       "consequence": "Consequence", "HGNC_symbol": "HGNC_Symbol", "Pathway": "Pathway",
                       "key": "Key", "#_affected_pathways": "Num_Affected_Pathways", "MPC": "MPC",
                       "pph2_prediction": "pph2_Prediction", "pph2_value": "pph2_Value",
                       "adj_cons": "Adj_Consequence", "_id": "id"})
    df = df[['id', 'Child_id', 'Key', 'Chr', 'Pos', 'Ref', 'Alt', 'Consequence', 'HGNC_Symbol', 'Pathway',
                 'Num_Affected_Pathways', 'MPC', 'pph2_Prediction', 'pph2_Value', 'Adj_Consequence']]
    return df


# ---------- Extract dataset for clustering Part I ----------- #


def compute_binary(df):

    df['key'] = df['Chr'].map(str) + "|" + df['Position'].map(str) + "|" + df['Ref'].map(str) + "|" + df['Alt'].map(str)

    # Compute again the binary variables
    pivot = pd.concat([df.drop('Pathway', 1), pd.get_dummies(df.Pathway).mul(1)], axis=1)
    df = pd.DataFrame(pivot.to_records())
    cols_to_filter = [col for col in df if (col.startswith('path') or col.startswith('child') or col.startswith('key'))]
    df = df[cols_to_filter]
    df = df.groupby(['child_id', 'key']).max()
    df = df.sort_values(by='child_id', ascending=True)
    df = pd.DataFrame(df.to_records())
    df = df.rename(columns={"child_id": "Child_id", "key": "Key"})
    return df


def join_dfs(df1, df2):

    df = pd.merge(df1, df2, on=['Child_id', 'Key'])
    return df


def blocks_mpc(path, dir_chunks, chunksize, df):

    if os.path.isfile(path):
        partition_mpc(path, dir_chunks, chunksize)
        patients_mpc = annotate_mpc(dir_chunks, df)
        return patients_mpc
    elif os.path.isdir(path):
        patients_mpc = annotate_mpc(path, df)
        return patients_mpc


def logging_info_mpc(args, folderpath, folderpath2):

    """ Records an INFO type log.

    Parameters:
        args (dictionary): parsed arguments from command line.
        folderpath (str): destination of the info log.
        folderpath2 (str): folder where the annotated patients matrix is downloaded.

    Returns:
        logging_info(args, folderpath): info log.

    """

    is_base = "Annotated MPC data downloaded at: {}".format(folderpath2)

    log_filename = os.path.join(folderpath, 'MPC-pph2_logINFO_{}.log'.format(strftime("%Y-%m-%d_%H:%m:%S", gmtime())))
    logging.basicConfig(filename=log_filename, level=logging.INFO)

    logging.info('MPC annotations protocol started generation on: ' + datetime.now().strftime('%d-%m-%Y %H:%M:%S'))
    logging.info('Process initiated by the user ' + getpass.getuser())
    logging.info('\n' + '\n' +
                 "Script generated with the following parameters --- " + '\n' +
                 " mutations-File: " + args['inputFile'] + '\n' +
                 " pathways-Directory: " + args['pathwaysDirectory'] + '\n' +
                 " MPC-File: " + args['inputMPC'] + '\n' +
                 "-pathway: " + none_to_str(args['pathway']) + '\n' +
                 "-gene: " + none_to_str(args['gene']) + '\n' +
                 "-patient: " + none_to_str(args['patient']) + '\n' +
                 "-mutation: " + none_to_str(args['mutation']) + '\n' +
                 is_base)


"""    
# Compute again normalized variables

df = pd.read_csv('/home/gabi/Documents/GitHub-PBPM/PBPM/Clustering/int-matrix-nn.csv', header=0)
df['key'] = df['Chr'].map(str) + "|" + df['Position'].map(str) + "|" + df['Ref'].map(str) + "|" + df['Alt'].map(str)
pivot = df.pivot_table(values='consequence', index=['child_id', 'key'], columns=['Pathway'], aggfunc='count', fill_value=0)
flattened = pd.DataFrame(pivot.to_records())

# 2753

p = list(flattened.columns.values)[1:]
g = pd.read_csv('/home/gabi/Documents/GitHub-PBPM/PBPM/Clustering/all_genes_to_norm.csv', header=0)
df_to_norm = g[g['Pathway'].isin(p)]

df2 = flattened.set_index(['child_id', 'key']).div(df_to_norm.set_index('Pathway')['HGNC_symbol']).reset_index()

# Join both matrices
final = pd.merge(data, df2, on='key')
"""