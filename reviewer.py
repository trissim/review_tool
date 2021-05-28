from manager import dataset, merge_datasets
import pandas as pd
import argparse
import math
import os.path as osp
import mygene
import sys
from pyorthomap import FindOrthologs

mg = mygene.MyGeneInfo()
species_id = {"human":9606,"rat":10116, "mouse":10090}
#species_ensembl_name = { "human": {
#                                   'from_dataset' : 'hsapiens_gene_ensembl',
#                                   'to_dataset' : 'mmusculus_gene_ensembl',
#                                   'to_homolog_attribute' : 'mmusculus_homolog_ensembl_gene',
#                                   'from_gene_id_name' : 'human_ensembl_gene_id',
#                                   'to_gene_id_name' : 'mouse_ensembl_gene_id'
#                                   },
#                       "rat": {
#                                   'from_dataset' : 'rnorvegicus_gene_ensembl',
#                                   'to_dataset' : 'rnorvegicus_gene_ensembl',
#                                   'to_homolog_attribute' : 'rnorvegicus_homolog_ensembl_gene',
#                                   'from_gene_id_name' : 'rnorvegicus_ensembl_gene_id',
#                                   'to_gene_id_name' : 'rnorvegicus_ensembl_gene_id'
#                                   },
#                       "mouse": {
#                                   'from_dataset' : 'mmusculus_gene_ensembl',
#                                   'to_dataset' : 'mmusculus_gene_ensembl',
#                                   'to_homolog_attribute' : 'mmusculus_homolog_ensembl_gene',
#                                   'from_gene_id_name' : 'hun_ensembl_gene_id',
#                                   'to_gene_id_name' : 'mouse_ensembl_gene_id'
#                                   }
#                       }
#"mouse":10090}
converted_genes_id = {"human":{},"rat":{},"mouse":{}}

def parse_args():
    parser = argparse.ArgumentParser(description="compare RNA-Seq results")
    parser.add_argument("-m","--master_sheet",required=True,type=str,help="path to master sheet")
    parser.add_argument("-e","--excels_path", default="./Excel Files/astrocytes/",type=str,help="path to folder with excel sheets")
    parser.add_argument("-d","--disease",nargs="*",help="filter by disease/model")
    parser.add_argument("-s","--stage",nargs="*",help="filter by stage of disease/model")
    parser.add_argument("-c","--cell",nargs="*",help="filter by cell type")
    parser.add_argument("-C","--convert",nargs="*",help="organism name to convert gene symbols to")
    parser.add_argument("-t","--tissue",nargs="*",help="filter by tissue")
    parser.add_argument("-p","--species",nargs="*",help="filter by species")
    parser.add_argument("-x","--excel", action="store_true", help="use this flag if the tables for each paper are excel files")
    parser.add_argument("--percent", type=int, help="fiter by minimum percent of papers that have the same gene")
    parser.add_argument("--absolute",type=int, help="filter by minimum number of papers that have the same gene")
    parser.add_argument("-u","--unique_papers", action="store_true", help="count multiple sheets from the same paper as one sheet")
    parser.add_argument("-o","--output",required=True,help="output file name")
    return parser.parse_args()

def filter_master_sheet(master_sheet,args):
    if args.disease:
        master_sheet = filter_sheets(master_sheet,"Disease/model", args.disease)
    if args.stage:
        master_sheet = filter_sheets(master_sheet,"Disease stage", args.stage)
    if args.cell:
        master_sheet = filter_sheets(master_sheet,"Cell type", args.cell)
    if args.tissue:
        master_sheet = filter_sheets(master_sheet,"Tissue", args.tissue)
    if args.species:
        master_sheet = filter_sheets(master_sheet,"Species", args.species)
    return master_sheet

def filter_sheets(master_sheet,column, values, operator=any):
    def filter_row(row):
        col_values = str(row[column]).split(",")
        return any(value in values for value in col_values)
    master_sheet = master_sheet[master_sheet.apply(filter_row, axis=1)]
    return master_sheet

def add_sheet(index,prefix,use_excel):
    if use_excel:
        extension = ".xlsx"
    else:
        extension = ".csv"
    path = osp.join(prefix,str(index)+extension)
    if osp.exists(path):
        try:
            return dataset(path,"Gene","logFC",title=osp.basename(path))
        except Exception:
            return None
    else:
        return None

def add_sheets(master_sheet,prefix,use_excel):
    all_sheets = {}
    for i ,row in master_sheet.iterrows():
        specie = class_sheet_by_species(row)
        if not specie is None:
            if not specie in all_sheets.keys():
                all_sheets[specie] = []
            sheet = add_sheet(int(row['Index']),prefix,use_excel)
            if not sheet is None:
                all_sheets[specie].append(sheet)
    return all_sheets

def class_sheet_by_species(master_sheet_row):
    species = master_sheet_row['Species']
    if not type(species) is str: return None
    if "human" in species: return "human"
    if "mouse" in species: return "mouse"
    if "rat" in species: return "rat"
    return None

def convert_all_sheets_to_species(all_sheets,target_species):
    conversion_dict = {}
    final_sheets=[]
    for origin_species,sheet_list in all_sheets.items():
        if not target_species == origin_species:
            all_genes = []
            for sheet in sheet_list:
                all_genes+=(sheet.df['Gene ID'].tolist())
            all_genes = [str(gene) for gene in all_genes]
            all_genes = set(all_genes)
    #       converted = {**conversion_dict, **convert_gene_name(origin_species, target_species, all_genes)}
            ortholog_dict = convert_gene_name(origin_species, target_species, all_genes)
            final_sheets += convert_sheet_list_to_species(sheet_list,ortholog_dict)

#            converted_genes_entrezid = convert_genes_to_entrezids(all_genes,origin_species,target_species)
#            entrezids_to_geneids = convert_entrezids_to_gene_ids(list(converted_genes_entrezid.values()))
#            converted += convert_sheet_list_to_species(sheet_list,converted_genes_entrezid,entrezids_to_geneids)
        else:
            final_sheets += sheet_list
        #print(len(converted))
    return final_sheets

def convert_entrezids_to_gene_ids(entrezids):
    #results = mg.getgenes(entrezids,returnall=True)
    results = mg.getgenes(entrezids)
    entrezid_to_geneid = {}
    for result in results:
        if 'notfound' in result.keys():
            entrezid_to_geneid[result["query"]] = None
        else:
            entrezid_to_geneid[result["query"]] = result['symbol']
    return entrezid_to_geneid

def convert_genes_to_entrezids(gene_names,origin_species,target_species):
    converted_genes_id = {}
    #results = mg.querymany(gene_names,species=origin_species,fields="homologene",scopes="symbol",returnall=True)
    results = mg.querymany(gene_names,species=origin_species,fields="homologene",scopes="symbol")
    for result, gene_name in zip (results,gene_names):
        converted_id = None
        homologues = None
        if not "notfound" in result.keys():
            if 'homologene' in result.keys():
                homologues = result['homologene']['genes']
            if not homologues is None:
                homologues = [tup for tup in homologues if tup[0] == species_id[target_species]]
                if len(homologues) > 0:
                    converted_id = homologues[0][1]
        converted_genes_id[gene_name] = converted_id
    return converted_genes_id


def convert_sheet_list_to_species(sheets,ortholog_dict):
    converted_sheets = []
    for sheet in sheets:
        to_remove = []
        #sheet.df.insert(0, 'Original Gene ID', sheet.df['Gene ID'])
        for i, row in sheet.df.iterrows():
            if row['Gene ID'] in ortholog_dict.keys():
                sheet.df.loc[i, 'Gene ID'] = ortholog_dict[row['Gene ID']]
            else:
                to_remove.append(i)
        sheet.df = sheet.df.drop(to_remove)
        converted_sheets.append(sheet)
    return converted_sheets

#def convert_sheet_list_to_species(sheets,converted_genes_to_entrezid,entrezids_to_geneids):
#   converted_sheets = []
#   for sheet in sheets:
#       for i, row in sheet.df.iterrows():
#           to_remove = []
#           entrezid = str(converted_genes_to_entrezid[row["Gene ID"]])
##           if entrezid is None or not entrezid in entrezids_to_geneids.keys():
#           try:
#               if entrezid == "None":
#                   raise
#               sheet.df.loc[i, 'Gene'] = entrezids_to_geneids[entrezid]
#           except Exception:
#               #print("missing",entrezid)
#               to_remove.append(i)
##           else:
##               converted_symbol = entrezids_to_geneids[entrezid]
##               sheet.df.loc[i, 'Gene'] = gene_name
#       sheet.df = sheet.df.drop(to_remove)
#       converted_sheets.append(sheet)
#   return converted_sheets

#def convert_gene_name(gene_name,origin_species,target_species):
#        hits = results['hits']
#        if len(hits) > 0:
#            if 'homologene' in hits:
#                homologues = hits[0]['homologene']['genes']
#            else: return None
#        else: return None
#        homologues = [tup for tup in homologues if tup[0] == species_id[target_species]]
#        if len(homologues) > 0:
#            return mg.getgene(homologues[0][1])['symbol']
#        else: return None

#orthologsbiomart way
def convert_gene_name(origin_species,target_species, *gene_names):

    def hs2mm(genes):
        searcher = FindOrthologs(
        host = 'http://www.ensembl.org',
        mart = 'ENSEMBL_MART_ENSEMBL',
        from_dataset = 'hsapiens_gene_ensembl',
        to_dataset = 'mmusculus_gene_ensembl',
        from_filters = 'hgnc_symbol',
        from_values = genes,
        to_attributes = ['external_gene_name'],
        to_homolog_attribute = 'mmusculus_homolog_ensembl_gene',
        from_gene_id_name = 'human_ensembl_gene_id',
        to_gene_id_name = 'mouse_ensembl_gene_id'
        )
        return searcher.map()

    def mm2hs(genes):
        searcher = FindOrthologs(
        host = 'http://www.ensembl.org',
        mart = 'ENSEMBL_MART_ENSEMBL',
        from_dataset = 'mmusculus_gene_ensembl',
        to_dataset = 'hsapiens_gene_ensembl',
        from_filters = 'hgnc_symbol',
        from_values = genes,
        to_attributes = ['hgnc_symbol'],
        to_homolog_attribute = 'hsapiens_homolog_ensembl_gene',
        from_gene_id_name = 'mouse_ensembl_gene_id',
        to_gene_id_name = 'human_ensembl_gene_id'
        )
        return searcher.map()

    gene_names = list(gene_names[0])
    if origin_species  == "human" and  target_species  == "mouse":
        orthologs = hs2mm(gene_names)


    if origin_species  == "mouse" and  target_species  == "human":
        orthologs = mm2hs(gene_names)
        print(orthologs)

    mapped = {}
    orthologs = orthologs[orthologs['external_gene_name'].notna()]
    for i,r in orthologs.iterrows():
        mapped[str(r['external_gene_name']).upper()]=r['hgnc_symbol']

    return mapped

def num_percent_common(series):
    total = len(series)-2
    remain = len(series.dropna())-2 #original gene ID and gene ID excluded
    return remain, float(remain)/float(total)

def get_unique_papers(master_sheet):
    num_papers=len(set(master_sheet["DOI"]))
    index_to_doi = {}
    for i, row in master_sheet.iterrows():
        if not pd.isna(row['Index']):
            index_to_doi[str(int(row['Index']))] = row["DOI"]
    return num_papers, index_to_doi

def num_percent_paper(series,num_papers, index_to_doi):
    series = series.dropna()
    indexes = list(series.keys())
    def contains_digit(string):
        for c in string:
            if c.isdigit():
                return True
        return False
    indexes = ["".join(c for c in index if c.isdigit()) for index in indexes if contains_digit(index)]
    papers = len({index_to_doi[index] for index in indexes})
    #remain = len(series.dropna())-2 #original gene ID and gene ID excluded

    return papers, float(papers)/float(num_papers)

def filter_min_row(final_csv,min_val,master_sheet,absolute=False,unique_papers=False):
    to_remove = []
    if unique_papers:
        num_papers, index_to_doi = get_unique_papers(master_sheet)
    for i, row in final_csv.iterrows():
        if unique_papers:
            remain, percent = num_percent_paper(row, num_papers, index_to_doi)
        else:
            remain, percent = num_percent_common(row)
        if absolute:
            if not min_val <= remain: to_remove.append(i)
        else:
            if not min_val <= percent: to_remove.append(i)
    final_csv = final_csv.drop(to_remove)
    return final_csv

def backup_original_gene_id(all_sheets):
    for k,sheets in all_sheets.items():
        for sheet in sheets:
            sheet.df.insert(0, 'Original Gene ID', sheet.df['Gene ID'])
    return all_sheets

def main():
    args = parse_args()
    print("reading master sheet")
    master_sheet = pd.read_excel(args.master_sheet)
    print("filtering master sheet")
    master_sheet = filter_master_sheet(master_sheet,args)
    all_sheets = add_sheets(master_sheet,args.excels_path,args.excel)
    if len(all_sheets.keys()) > 1:
        if not args.convert:
            sys.exit("More than one species detected.\nSpecify species to convert to using --convert")
        print("converting gene symbols")
        all_sheets = backup_original_gene_id(all_sheets)
        all_sheets = convert_all_sheets_to_species(all_sheets,args.convert[0])
    else:
        all_sheets = list(all_sheets.values())[0]
    print("merging datasets")
    merged = merge_datasets(all_sheets)
    merged.drop(merged.loc[merged['Gene ID']=='NAN'].index, inplace=True)
    print("filtering merged dataset")
    if args.percent:
        if not args.percent == 100:
            merged = filter_min_row(merged,float(args.percent)/100.0, master_sheet, unique_papers=args.unique_papers)
    elif args.absolute:
        merged = filter_min_row(merged, args.absolute, master_sheet, absolute=True, unique_papers=args.unique_papers)
    print("writing csv")
    merged.to_csv(args.output)

if __name__ == "__main__":
    main()
