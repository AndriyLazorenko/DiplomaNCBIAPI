import requests
import time
from Bio import Entrez


def read_temp_txt(temp, data):
    with open(temp) as temp_fh:
        with open(data, 'a+') as data_fh:
            for line in temp_fh:
                line = line.strip()
                if not line.startswith('rs'):
                    data_fh.write(line[line.index('rs'):] + '\n')
                else:
                    data_fh.write(line + '\n')


def get_IDs(database, term, retmax, retstart):
    handle = Entrez.esearch(db=database, retmode="xml", term=term, usehistory='y', retmax=retmax, retstart=retstart)
    records = Entrez.read(handle)
    id_list = records['IdList']
    count = records['Count']
    print('len', len(id_list))
    return id_list, ','.join(id_list), count


def get_IDs_to_file(database, term, retMax):
    count = int(get_IDs(database, term, 5, retstart=5)[2])
    rng = count // retMax + 1
    # print(rng)
    with open(ids_file, 'a+') as fh:
        for i in range(rng):
            ids_str = get_IDs(database, term, retMax, i * retMax)[0]
            for it in ids_str:
                fh.write(it + '\n')


def get_strings_from_IDs_to_file(datab, ids, retmax, fh):
    handle_fetch = Entrez.efetch(db=datab, id=ids, rettype='txt', retmode='text', retmax=retmax)
    req = requests.get(handle_fetch.url)
    # with open(file, 'a+') as fh:
    print(handle_fetch.url)
    fh.write(req.text)


def get_downloaded_dataset_length(data):
    try:
        with open(data) as data_fh:
            # print('list',len(list(data_fh)))
            return len(list(data_fh))
    except FileNotFoundError as err:
        return 0

        # if data_fh: return (len(data_fh.readlines()))
        # else: return 0


if __name__ == "__main__":
    # параметри API:
    Entrez.email = 'ivan.shpotenko@gmail.com'  # пошта необхідна для API-шки
    base_url_esearch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    base_url_efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    # параметри запиту:
    database = 'SNP'
    term = '("Homo sapiens"[Organism] OR homo sapiens[All Fields]) AND (snp[Snp_Class] AND (00000.0500[GLOBAL_MAF] : 00000.5000[GLOBAL_MAF]) AND by 1000g[Validation] NOT merged rs[Filter] AND Y[ALLELE])'
    # дії:
    temp = 'temp.txt'
    ids_file = 'ids_list_{}.txt'.format(database)
    dataset = 'dataset_raw_{}.txt'.format(database)

    retMax = 100  # TODO: trying with bigger then 200 causes errors. Why?

    #     скачування ID-шок:
    # get_IDs_to_file(database,term,retMax)

    count = int(get_IDs(database, term, 5, retstart=5)[2])
    rng = count // retMax + 1

    with open(dataset,'a+') as fh:
        with open(ids_file) as idsh:
            id_list = idsh.readlines()
        for i in range(rng):
            get_strings_from_IDs_to_file(database,id_list[i*retMax:i*retMax+retMax],retMax,fh)

