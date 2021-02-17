from manager import dataset, merge_datasets
import pandas as pd
import argparse
import math
import os.path as osp
import mygene

def parse_args():
    parser = argparse.ArgumentParser(description="compare RNA-Seq results")
    parser.add_argument("-m","--master_sheet",required=True,type=str,help="path to master sheet")
    parser.add_argument("-e","--excels_path", default="./Excel Files/astrocytes/",type=str,help="path to folder with excel sheets")
    parser.add_argument("-d","--disease",nargs="*",help="filter by disease/model")
    parser.add_argument("-s","--stage",nargs="*",help="filter by stage of disease/model")
    parser.add_argument("-c","--cell",nargs="*",help="filter by cell type")
    parser.add_argument("-t","--tissue",nargs="*",help="filter by tissue")
    parser.add_argument("-p","--species",nargs="*",help="filter by species")
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

def add_sheet(index,prefix):
    path = osp.join(prefix,str(index)+".xlsx")
    if osp.exists(path):
        try:
            return dataset(path,"Gene","logFC",title=osp.basename(path))
        except Exception:
            return None
    else:
        return None

def add_sheets(master_sheet,prefix):
    all_sheets = []
    for i ,row in master_sheet.iterrows():
        all_sheets.append(add_sheet(int(row['Index']),prefix))
    all_sheets = list(filter(None,all_sheets))
    return all_sheets

def count_species(master_sheet):
    unique_species = set(master_sheet['species'])
    if unique_species > 1:
        most_species = ""
        max_count = 0
        for specie in unique_species:
             count = list(master_sheet['species']).count(v)
             if count > max_count:
                 max_count = count
                 most_species = specie

def class_sheet_by_species(master_sheet):
    sheets_by_species = {}
#    for i,row in master_sheet.iterrows():
#        if 
#        sheets_by_species['



def convert_gene_name(gene_name,origin_species,target_species):
        results = mg.query(gene_name,species=origin_species,fields="homologene",size=1)
        homologues = results['hits'][0]['homologene']['genes']
        homologues = [tup for tup in homologues if tup[0] == species_id[target_species]]
        if homologues > 0:
            return mg.get_gene(homologene[0][1])['symbol']
        else: return None


def main():
    args = parse_args()
    master_sheet = pd.read_excel(args.master_sheet)
    master_sheet = filter_master_sheet(master_sheet,args)
    all_sheets = add_sheets(master_sheet,args.excels_path)
    merged = merge_datasets(all_sheets)
    merged.to_csv(args.output)

if __name__ == "__main__":
    main()
