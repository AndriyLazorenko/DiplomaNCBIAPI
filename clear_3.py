import urllib

import requests
import time
from Bio import Entrez
from functools import wraps
import errno
import os
import signal  # would NOT work properly on WINDOWS. UNIX only


class TimeoutError(Exception):
    pass


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    '''
    Function creates a decorator, which applied to other function would break it if that function
    takes too long to execute.

    Note that this will only work on UNIX.
    The basic idea is to use signal handlers to set an alarm for some time interval and raise
    an exception once that timer expires.
    The process for timing out an operations is described in the documentation for signal module.
    :param seconds: specifies "too long" in seconds
    :param error_message: custom error message
    :return:
    '''

    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator


def get_search_query_options(database, term):
    '''
    Returns result of eSearch query on NCBI server: № of results, QueryKey and WebEnv.

    See details in docs: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    :param database: database name on NCBI server
    :param term: desired search query on NCBI, one can obtain it from 'Search details' window
     on NCBI cite to the right (ex.: http://www.ncbi.nlm.nih.gov/snp )
    :return:
    '''
    handle = Entrez.esearch(db=database, retmode="xml", term=term, usehistory='y')
    records = Entrez.read(handle)
    count = int(records['Count'])
    query_key = records['QueryKey']
    web_env = records['WebEnv']
    return count, query_key, web_env


@timeout(20, os.strerror(errno.ETIMEDOUT))
def write_results(datab, retmax, fh, query_key, web_env, retStart=10):
    '''
    Writes results obtained from NCBI (use QueryKey and WebEnv from eSearch) to file specified.

    See details in docs: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    :param datab: database name on NCBI server
    :param retmax: Total number of records from the input set to be retrieved, up to
    a maximum of 10,000. Optionally, for a large set the value of retstart can be
    iterated while holding retmax constant, thereby downloading the entire set in
    batches of size retmax.
    :param fh: file handler to write results into
    :param query_key: QueryKey from eSearch
    :param web_env: WebEnv from eSearch
    :param retStart: Sequential index of the first record to be retrieved (default=0,
    corresponding to the first record of the entire set). This parameter can be used
    in conjunction with retmax to download an arbitrary subset of records from the input set.
    :return:
    '''
    handle_fetch = Entrez.efetch(db=datab, rettype='txt', retmode='text', retmax=retmax, query_key=query_key,
                                 WebEnv=web_env, retstart=retStart)
    req = requests.post(handle_fetch.url)
    print(handle_fetch.url)
    fh.write(req.text)


def specify_parameters():
    print("Please specify some parameters. Hit 'Enter' for default values.")
    database_def = 'SNP'
    term_def = '("Homo sapiens"[Organism] OR homo sapiens[All Fields]) AND (snp[Snp_Class] AND (00000.0500[GLOBAL_MAF] : 00000.5000[GLOBAL_MAF]) AND by 1000g[Validation] NOT merged rs[Filter] AND Y[ALLELE])'
    retMax_def = 1300
    database = input("Specify NCBI Database [def: '{}']: ".format(database_def))
    if not database: database = database_def
    print('Database: {}'.format(database))
    term = input("Specify desired NCBI eSearch query [default]: ".format(term_def))
    if not term: term = term_def
    print('Term: {}'.format(term))
    retMax = input("Specify retmax parameter [def: {}]: ".format(retMax_def))
    if not retMax: retMax = retMax_def
    print('retMax: {}'.format(retMax))
    while True:
        email = input("Enter your email (NCBI will use it to spy and spam you down): ".format(
            retMax_def))  # пошта необхідна для API-шки
        if email:
            if len(email.split('@')) != 2: continue  # just sanity check
            break
    return database, term, retMax, email


if __name__ == "__main__":
    base_url_esearch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    base_url_efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

    database, term, retMax, Entrez.email = specify_parameters()

    dataset = 'dataset_raw_{}.txt'.format(database)  # destination filename
    print("Results will be downloaded in '{}' file. Be patient.".format(dataset))

    count, query_key, web_env = get_search_query_options(database, term)
    print('QueryKey: {}, WebEnv: {}'.format(query_key, web_env))
    with open(dataset, 'w') as fh:
        # for i in range(count // retMax + 1):
        i = 0
        while True:
            try:
                t0 = time.time()
                write_results(database, retStart=i * retMax, retmax=retMax, fh=fh,
                              query_key=query_key,
                              web_env=web_env)
                t1 = time.time() - t0
                print('batch #{} ({} of {} lines): {} sec.'.format(i + 1, (i + 1) * retMax, count, t1))
                i += 1
            except TimeoutError as err:
                # retMax -=100
                print(err, 'Trying again: ret_max = {}'.format(retMax))
                pass
            except urllib.error.HTTPError as err2:
                print(err2)
                print("Database may be downloaded. Check '{}' file. Have fun.".format(dataset))
                break
    with open(dataset) as fh:
        print(len(list(fh.readlines())))
